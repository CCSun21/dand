import argparse

import h5py
from ase.db import connect


def main(args):
    
    print_args(args)
    
    input_path  = args.input_path
    output_path = args.output_path
    
    # Data structure to hold the computed results
    rxn_data = {}

    # Extract data from ASE database
    with connect(input_path) as db:
        for row in db.select():
            # Extract unique_id and other data
            unique_id = row.data['unique_id']
            chem_group_name, rxn_group_name, index = unique_id.split('_')
            
            if chem_group_name not in rxn_data:
                rxn_data[chem_group_name] = {}
            
            if rxn_group_name not in rxn_data[chem_group_name]:
                rxn_data[chem_group_name][rxn_group_name] = {
                    'atomic_numbers': row.toatoms().numbers,
                    'energies': [],
                    'forces': [],
                    'positions': []
                }
            rxn_data[chem_group_name][rxn_group_name]['energies'].append(row.energy)
            rxn_data[chem_group_name][rxn_group_name]['forces'].append(row.forces)
            rxn_data[chem_group_name][rxn_group_name]['positions'].append(row.toatoms().positions)

    # Save the data to an h5 file
    with h5py.File(output_path, 'w') as h5file:
        # Ensure the 'data' group exists
        if 'data' not in h5file:
            data_group = h5file.create_group('data')
        else:
            data_group = h5file['data']
        
        # Iterate through the rxn_data dictionary to save datasets
        for chem_group_name in rxn_data:
            if chem_group_name not in data_group:
                chem_group = data_group.create_group(chem_group_name)
            else:
                chem_group = data_group[chem_group_name]
            
            for rxn_group_name, rxn_entry in rxn_data[chem_group_name].items():
                if rxn_group_name not in chem_group:
                    rxn_group = chem_group.create_group(rxn_group_name)
                else:
                    rxn_group = chem_group[rxn_group_name]
                
                # Add datasets to the reaction group
                rxn_group.create_dataset('atomic_numbers', data=rxn_entry['atomic_numbers'])
                rxn_group.create_dataset('wB97x_6-31G(d).energy', data=rxn_entry['energies'])
                rxn_group.create_dataset('wB97x_6-31G(d).forces', data=rxn_entry['forces'])
                rxn_group.create_dataset('positions', data=rxn_entry['positions'])

    print('Compiled successfully!')

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()

def get_parser():
    parser = argparse.ArgumentParser(description="Translate ase db file into hdf5 file.")
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help="Path of the input wB97X ASE db file")
    parser.add_argument('-o', '--output_path', required=True, 
                        help="Path of the output wB97X hdf5 file")

    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)


