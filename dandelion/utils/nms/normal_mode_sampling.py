import os
import shutil
import numpy as np

from ase import Atoms
from ase.io import read
from ase.db import connect
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.build import minimize_rotation_and_translation
from xtb.ase.calculator import XTB

# Constants
kB = 8.617333262145e-5  # Boltzmann constant in eV/K
cm_to_eV = 1.23981e-4

def align_non_eq_strucs(eq_struc, non_eq_strucs):
    reference = eq_struc  # Use the first structure as the reference  
    for i in range(0, len(non_eq_strucs)):
        minimize_rotation_and_translation(reference, non_eq_strucs[i])  # Align each structure to the reference

def write_xyz(filename, atoms):
    with open(filename, 'w') as f:
        f.write(str(len(atoms)) + '\n\n')
        for atom in atoms:
            f.write(f'{atom.symbol:<8}   {atom.x:12.8f}   {atom.y:12.8f}   {atom.z:12.8f}\n')

def generate_eq_and_non_eq_strucs(eq_struc, temp=1500):
    non_eq_strucs = []
    non_eq_energy = []
    energy_histogram = [400, 240, 150, 90, 50, 30, 20, 10, 10]  # Maximum room for each energy bin
    
    ############# Optimize the eq_struc and get normal modes ############
    calc = XTB(method="GFN2-xTB", atoms=eq_struc)
    eq_struc.calc = calc
    opt = BFGS(eq_struc)
    opt.run(fmax=1e-4)
    
    eq_energy = eq_struc.get_potential_energy()
    eq_positions = eq_struc.get_positions()

    N_atoms = len(eq_struc.get_atomic_numbers())
    NkBT = N_atoms*kB*temp  # Thermal energy

    vib = Vibrations(eq_struc)
    vib.run()
    modes = vib.get_vibrations().get_modes()
    frequencies = vib.get_vibrations().get_frequencies() * cm_to_eV
    vib.clean() # clear the json file
    N_modes = len(frequencies)
    #####################################################################


    # loop over non-eq candidates
    for iter in range(100000):
        displaced_positions = eq_positions.copy()


        ########## Linear combination of 3N-6 normal modes ##########
        for i in range(6, N_modes):
            
            frequency = frequencies[i]
            
            # gaussian_widths for different iteration ranges
            gaussian_widths = [0.0001, 0.0002, 0.0005, 0.001,0.003, 0.00001, 0.00005, 0.000005]
            index = iter % len(gaussian_widths)
            width = gaussian_widths[index]
            
            delta_e = np.random.normal(0, width, 1) #controll middle number
            amplitude = np.sqrt(2 * abs(delta_e)) / abs(frequency)
            mode = modes[i]
            random_sign = np.random.choice([-1, 1])
            
            displaced_positions += amplitude * mode * random_sign
            
        
        
        ###### Calculate the single point energy of the non_eq struc #######
        non_eq_struc = Atoms(symbols=eq_struc.get_chemical_symbols(),
                             positions=displaced_positions)
        try:
            calc = XTB(method="GFN2-xTB", atoms=eq_struc)
            non_eq_struc.set_calculator(calc)
            energy = non_eq_struc.get_potential_energy()
        except: #If not converged as perturbation was too big
            print(f'{iter}th iteration', energy_histogram, f'Current : error!')
            continue
        
        
        
        ###### Determine the energy bin and update the energy histogram ##########
        
        energy_bin = int( (energy - eq_energy) / (0.75 * NkBT) )
        print(f'{iter}th iteration', energy_histogram, f'Current : {energy_bin}')
        
        if energy_bin >= 9: #Too high energy structure
            continue
            
        if energy_histogram[energy_bin] > 0:
            energy_histogram[energy_bin] -= 1 #fill the bin
            non_eq_strucs.append(non_eq_struc)
            non_eq_energy.append(energy - eq_energy)
            
        
        
        ############### Check whether energy bin is full ####################
        if sum(energy_histogram) == 0:
            
            # sort non_eq_strucs based on the relative energy
            non_eq_energy_idx = np.argsort(non_eq_energy)
            non_eq_strucs = [non_eq_strucs[i] for i in non_eq_energy_idx]
            
            print('''     
          ******************************HURRAY****************************
          *** 1000 Non-eq geometries in Boltzmann Dist. were generated ***
          ****************************************************************     
                  ''')
            break
        
    return eq_struc, non_eq_strucs


def main(args):
    
    print_args(args)
    
    input_path  = args.input_path
    if not os.path.isdir(input_path):
        sys.exit(f"Error: '{input_path}' is not a directory.")
    output_path  = args.output_path

    with connect(output_path) as db:
        for root, dirs, files in os.walk(input_path):
            for filename in files:
                if filename.endswith('.xyz'):
                    
                    #Get the premature xyz
                    xyz_path = os.path.join(root, filename)
                    premature = read(xyz_path)

                    eq_struc, non_eq_strucs = generate_eq_and_non_eq_strucs(premature)
                    align_non_eq_strucs(eq_struc, non_eq_strucs)
                    db.write(eq_struc)

                    ####### Write the non-equilibrium structure to the subdirectory #######
                    for i, non_eq_struc in enumerate(non_eq_strucs):
                        db.write(non_eq_struc)


def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()
    

def get_parser():
    parser = argparse.ArgumentParser(description="Run normal mode sampling")

    parser.add_argument('-i', '--input_path', type=str, required=True,
                        help='Path of input directory containing reactants')
    parser.add_argument('-o', '--output_path', type=str, required=True,
                        help='Path of output ase db file.')

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
