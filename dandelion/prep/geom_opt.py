import os
import shutil
import argparse
import warnings
from ase.io import read
from ase.optimize import BFGS
from xtb.ase.calculator import XTB
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def write_xyz(filename, atoms):
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n\n")
        for atom in atoms:
            f.write(f"{atom.symbol:<2}  {atom.position[0]:15.8f}  {atom.position[1]:15.8f}  {atom.position[2]:15.8f}\n")

def generate_eq_struc(atoms):
    atoms.calc = XTB(method="GFN2-xTB")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        opt = BFGS(atoms, logfile=None)
        opt.run(fmax=1e-4)
    return atoms

def process_file(input_file, output_dir):
    filename = os.path.basename(input_file)
    mol_dir = os.path.join(output_dir, os.path.splitext(filename)[0])
    os.makedirs(mol_dir, exist_ok=True)
    
    # Copy original file
    shutil.copy(input_file, mol_dir)
    
    # Generate and save optimized structure
    atoms = read(input_file)
    optimized_atoms = generate_eq_struc(atoms)
    write_xyz(os.path.join(mol_dir, 'struc.xyz'), optimized_atoms)
    
    # Remove the original copied file
    os.remove(os.path.join(mol_dir, filename))

def main(args):
    print_args(args)
    
    input_path = os.path.abspath(args.input_path)
    output_path = os.path.abspath(args.output_path)
    max_workers = args.max_workers
    
    # Get list of all .xyz files
    xyz_files = []
    for root, _, files in os.walk(input_path):
        xyz_files.extend([os.path.join(root, f) for f in files if f.endswith('.xyz')])
    
    # Process files in parallel with progress bar
    bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        list(tqdm(executor.map(process_file, xyz_files, [output_path]*len(xyz_files)), 
                  total=len(xyz_files), desc="Optimizing structures", smoothing=0, bar_format=bar_format, ncols=70))

def print_args(args):
    print("\nArguments provided:")
    for key, value in vars(args).items():
        print(f"  {key}: {value}")
    print()

def get_parser():
    parser = argparse.ArgumentParser(description="Optimize geometries using xTB")
    parser.add_argument('-i', '--input_path', required=True, 
                        help="Path of the input reactants directory")
    parser.add_argument('-o', '--output_path', required=True,
                        help='Path of output directory to store optimized geometries')
    parser.add_argument('-n', '--max_workers', type=int, default=1,
                        help='Number of processes to use for parallel execution.')
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)