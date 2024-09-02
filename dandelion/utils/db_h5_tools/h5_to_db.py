import os
import sys
import h5py
import argparse

from tqdm import tqdm
from ase import Atoms
from ase.db import connect
from ase.calculators.singlepoint import SinglePointCalculator


def main(args):
    
    print_args(args)
    
    input_path  = args.input_path
    if not os.path.isfile(input_path):
        sys.exit(f"Error: '{input_path}' is not a file.")
    output_path = args.output_path


    with h5py.File(input_path, 'r') as h5_file:
        data_group = h5_file['data']
        
        # Count total number of configurations 
        total_configs = sum(
            rxn_group['wB97x_6-31G(d).energy'].shape[0]
            for chem_group in data_group.values()
            for rxn_group in chem_group.values()
        )
        
        with connect(output_path) as db:
            with tqdm(total=total_configs, desc="Converting", unit="config") as pbar:
                for chem_group_name, chem_group in data_group.items():
                    for rxn_group_name, rxn_group in chem_group.items():
                        atomic_numbers = rxn_group['atomic_numbers'][:]
                        positions = rxn_group['positions'][:]
                        energies = rxn_group['wB97x_6-31G(d).energy'][:]
                        forces = rxn_group['wB97x_6-31G(d).forces'][:]
                        
                        for i in range(len(energies)):
                            atoms = Atoms(
                                numbers=atomic_numbers,
                                positions=positions[i],
                            )
                            atoms.set_calculator(SinglePointCalculator(
                                atoms,
                                energy=energies[i],
                                forces=forces[i]
                            ))
                            
                            unique_id = f"{chem_group_name}_{rxn_group_name}_{i}"
                            db.write(atoms, data={'unique_id': unique_id})
                            
                            pbar.update(1)

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()

def get_parser():
    parser = argparse.ArgumentParser(description="Translate hdf5 file into ase db file.")
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help="Path of the input wB97X hdf5 file")
    parser.add_argument('-o', '--output_path', required=True, 
                        help="Path of the output wB97X db file")

    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
