
import os
import sys
import json
import argparse
from functools import partial
from concurrent.futures import ProcessPoolExecutor

import uuid
import numpy as np
from tqdm import tqdm

import matplotlib
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from PIL import Image, ImageOps

import ase.db
from ase.io import read, write
from xtb.ase.calculator import XTB
from ase.optimize.bfgs import BFGS
from ase.utils.forcecurve import fit_images
from ase.neb import NEB, NEBOptimizer, NEBTools

class SuppressStderr:
    def __enter__(self):
        self._original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stderr.close()
        sys.stderr = self._original_stderr

def plot_mep(fit_list):
    fit_list[:,1,:] *= 23.0609 #to kcal/mol
    gray_scale = matplotlib.colormaps.get('binary', len(fit_list))
    fig, ax = plt.subplots()
    for i in range(len(fit_list)):
    
        if i+1 == len(fit_list):
            ax.plot(fit_list[i,0,:], fit_list[i,1,:], color='red', linewidth=3)
            break
    
        color = gray_scale(max(i / len(fit_list), 0.1))
        ax.plot(fit_list[i,0,:], fit_list[i,1,:], color=color)
    
    ax.set_title(f'Iter {len(fit_list)}')    
    ax.set_axisbelow(True)
    ax.set_ylabel("Energy [kcal/mol]")
    ax.set_xlabel("Reaction Coordinate [AA]")
    return fig

def get_fit(neb_tools):
    fit = fit_images(neb_tools.images)
    return fit.fit_path, fit.fit_energies

class CalculationChecker:
    def __init__(self, neb):
        self.neb = neb

    def check_calculations(self):
        missing_calculations = []
        for i, image in enumerate(self.neb.images[1:-1]):
            if {"forces", "energy"} - image.calc.results.keys():
                missing_calculations.append(i)

        if missing_calculations:
            raise ValueError(f"missing calculation for image(s) {missing_calculations}")


class DBWriter:
    def __init__(self, db_path, atomss):
        self.atomss = atomss
        self.db_path = db_path

    def write(self):
        with ase.db.connect(self.db_path) as db:
            for atoms in self.atomss:
                if atoms.calc.results:
                    db.write(atoms, data=atoms.calc.results)


def interpolate_band(atom_configs, transition_state=None):
    if transition_state:
        transition_state = read(transition_state)
        ts_positions = transition_state.get_positions()
        middle_idx = len(atom_configs) // 2
        atom_configs[middle_idx].set_positions(ts_positions)
        first_band = NEB(atom_configs[: middle_idx + 1])
        second_band = NEB(atom_configs[middle_idx:])
        first_band.interpolate("idpp")
        second_band.interpolate("idpp")
    else:
        band = NEB(atom_configs)
        band.interpolate("idpp")
    return atom_configs


def max_dimensions(frames):
    """Get the maximum width and height among a list of images."""
    max_width = max_height = 0
    for frame in frames:
        with Image.open(frame) as img:
            width, height = img.size
            max_width = max(max_width, width)
            max_height = max(max_height, height)
    return max_width, max_height

def pad_image(image_path, target_size):
    """Pad an image to the target size."""
    with Image.open(image_path) as img:
        img = ImageOps.expand(img, border=((target_size[0]-img.size[0])//2,
                                            (target_size[1]-img.size[1])//2,
                                            (target_size[0]-img.size[0]+1)//2,
                                            (target_size[1]-img.size[1]+1)//2),
                                fill='white')  # or another suitable color for your images
        return img

def frames_to_gif(frames, output_gif):
    # First, render each Atoms frame to an image
    image_paths = []
    for i, frame in enumerate(frames):
        img_path = f"tmp_frame_{i}_{uuid.uuid4()}.png"
        write(img_path, frame) 
        image_paths.append(img_path)

    # Determine the max dimensions
    max_width, max_height = max_dimensions(image_paths)

    # Create a list to store processed frames
    processed_frames = []

    # Pad each frame, ensuring a non-transparent background
    for img_path in image_paths:
        with Image.open(img_path) as opened_img:
            padded_frame = pad_image(img_path, (max_width, max_height))
            
            # Create a white background and paste the frame onto it to ensure non-transparency
            background = Image.new('RGB', padded_frame.size, (255, 255, 255))
            background.paste(padded_frame, mask=(padded_frame.split()[3] if len(padded_frame.split()) == 4 else None))
            processed_frames.append(np.array(background))

    # Extend the list of processed frames with a reversed copy (excluding the last frame)
    extended_frames = processed_frames + processed_frames[-2::-1]

    # Save the gif using imageio
    with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
        for processed_frame in extended_frames:
            writer.append_data(processed_frame)

    # Cleanup the temporary image files
    for img_path in image_paths:
        os.remove(img_path)


def process_seed(seed, n_images, neb_fmax, cineb_fmax, steps, output_path):
    
    with SuppressStderr(): # xTB is so noisy when not converged
        try:
            #print(f"Starting from seed : {seed}")
            reactant         = os.path.join(seed, 'reactant.xyz')
            product          = os.path.join(seed, 'product.xyz')
            transition_state = os.path.join(seed, 'ts.xyz')
            product = read(product)
            reactant = read(reactant)
            
            output = os.path.join(output_path, seed.split('/')[-2]+'-'+seed.split('/')[-1])
            os.makedirs(output, exist_ok=True)
            atom_configs = [reactant.copy() for i in range(n_images - 1)] + [product]
            
            for i, atom_config in enumerate(atom_configs):
                atom_config.calc = XTB(method='GFN2-xTB')

            #print("Relaxing endpoints ... ")
            BFGS(atom_configs[0], logfile=None).run()
            BFGS(atom_configs[-1], logfile=None).run()

            #print("Interpolating band ... ")
            interpolate_band(atom_configs, transition_state)

            #print("Running NEB ... ")
            neb = NEB(atom_configs, climb=True, parallel=False)
            calculation_checker = CalculationChecker(neb)
            neb_tools = NEBTools(neb.images)

            relax_neb = NEBOptimizer(neb, logfile=None)
            db_writer = DBWriter(os.path.join(output, "neb.db"), atom_configs)
            fmaxs = []
            fit_list = []
            relax_neb.attach(calculation_checker.check_calculations)
            relax_neb.attach(db_writer.write)
            relax_neb.attach(lambda: fmaxs.append(neb_tools.get_fmax()))
            relax_neb.attach(lambda: fit_list.append(get_fit(neb_tools)))
        
            converged = relax_neb.run(fmax=neb_fmax, steps=steps)

            if not converged:
                raise 
            
            #print("NEB has converged, turn on CI-NEB ...")
            neb.climb = True
            ci_converged = relax_neb.run(fmax=cineb_fmax, steps=steps)
                
            if ci_converged:
                open(os.path.join(output, "converged"), "w")
                #print("Reaction converged ... ")
            fit_list = np.array(fit_list)
            fig = plot_mep(fit_list)
            if ci_converged:
                np.save(os.path.join(output, "fitlist.npy"), fit_list)

            fig.savefig(os.path.join(output, "mep.png"))
            json.dump(fmaxs, open(os.path.join(output, "fmaxs.json"), "w"), indent=4)
            transition_state = max(atom_configs, key=lambda x: x.get_potential_energy())
            write(os.path.join(output, "transition_state.xyz"), transition_state)
            write(os.path.join(output, "transition_state.png"), transition_state)
            write(os.path.join(output, "reactant.xyz"), atom_configs[0])
            write(os.path.join(output, "reactant.png"), atom_configs[0])
            write(os.path.join(output, "product.xyz"), atom_configs[-1])
            write(os.path.join(output, "product.png"), atom_configs[-1])
            write(os.path.join(output, "mep.xyz"), atom_configs)        
            frames_to_gif(atom_configs, os.path.join(output, "mep.gif"))

            return seed
        
        except Exception as e:
            #print(f"Error processing seed {seed}: {e}")
            return None

def main(args):
    
    print_args(args)
    
    input_path  = args.input_path
    max_workers = args.max_workers
    n_images    = args.n_images
    neb_fmax    = args.neb_fmax
    cineb_fmax  = args.cineb_fmax
    steps       = args.steps
    output_path = args.output_path

    
    seeds = [dirpath for dirpath, _, filenames in os.walk(input_path) if "ts.png" in filenames]
    
    bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
    # Use a partial function to pass the extra arguments to process_seed
    process_with_args = partial(process_seed, n_images=n_images, neb_fmax=neb_fmax, 
                                cineb_fmax=cineb_fmax, steps=steps, output_path=output_path)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(tqdm(executor.map(process_with_args, seeds), 
                            desc='Seeds', total=len(seeds), smoothing=0, bar_format=bar_format, ncols=70))

    print('xTB-NEB completed!')

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()
    

def get_parser():
    parser = argparse.ArgumentParser(description="Run NEB calculations on filtered gsm jobs")

    parser.add_argument('-i', '--input_path', type=str, required=True,
                        help='Path of input directory containing filtered gsm jobs.')
    parser.add_argument('-o', '--output_path', type=str, required=True,
                        help='Path of output directory to store results.')
    parser.add_argument('-n', '--max_workers', type=int, default=1,
                        help='Number of processes to use for parallel execution.')
    parser.add_argument('--n_images', type=int, default=10,
                        help='Number of images for NEB.')
    parser.add_argument('--neb_fmax', type=float, default=0.5,
                        help='Fmax threshold for NEB.')
    parser.add_argument('--cineb_fmax', type=float, default=0.05,
                        help='Fmax threshold for CI-NEB.')
    parser.add_argument('--steps', type=int, default=500,
                        help='Maximum number of optimization steps.')

    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
