import os
import sys
import argparse
from itertools import repeat

import h5py
from tqdm import tqdm
import glob 


def main(args):
    
    print_args(args)
    
    input_path           = args.input_path
    if not os.path.isdir(input_path):
        sys.exit(f"Error: '{input_path}' is not a directory.")
    output_path          = args.output_path

    # Open the output file
    with h5py.File(output_path, 'w') as h5file_out:
        # Ensure the 'data' group exists in the output file
        if 'data' not in h5file_out:
            data_group_out = h5file_out.create_group('data')
        else:
            data_group_out = h5file_out['data']
        
        # Iterate through each input file
        for input_path in glob.glob(os.path.join(input_path, '**/wb97x.h5'), recursive=True):
            print(input_path)
            # Determine the prefix ('a' or 'b') based on the input file name
            prefix = os.path.basename(os.path.dirname(input_path))  # Assumes file name is 'a.h5' or 'b.h5'
            
            # Open the input file
            with h5py.File(input_path, 'r') as h5file_in:
                # Iterate through chemical groups in the input file
                for chem_group_name, chem_group in tqdm(h5file_in['data'].items(), desc="Formulas"):
                    # Ensure the chemical group exists in the output file
                    if chem_group_name not in data_group_out:
                        chem_group_out = data_group_out.create_group(chem_group_name)
                    else:
                        chem_group_out = data_group_out[chem_group_name]
                    
                    # Iterate through reaction groups in the chemical group
                    for rxn_group_name, rxn_group in tqdm(chem_group.items(), desc=f"Rxns in {chem_group_name}", leave=False):
                        # Prefix the reaction group name with 'a' or 'b'
                        rxn_group_name_prefixed = f"{prefix}_{rxn_group_name}"
                        
                        # Ensure the reaction group exists in the output file
                        if rxn_group_name_prefixed not in chem_group_out:
                            rxn_group_out = chem_group_out.create_group(rxn_group_name_prefixed)
                        else:
                            rxn_group_out = chem_group_out[rxn_group_name_prefixed]
                        
                        # Copy datasets from input to output, creating new datasets
                        for dset_name, dset in rxn_group.items():
                            data = dset[:]
                            rxn_group_out.create_dataset(dset_name, data=data)

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()
    
def get_parser():
    parser = argparse.ArgumentParser(description='Merge h5 files in input directory')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Input path of directory containing h5 files to merge')    
    parser.add_argument('-o', '--output_path', required=True, 
                        help='Output path of the merged h5 file.')
    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)

