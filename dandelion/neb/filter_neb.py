import os
import json
import shutil
import argparse
from collections import defaultdict

from tqdm import tqdm
import numpy as np
from ase.io import read
from ase.vibrations import Vibrations
from xtb.ase.calculator import XTB


def is_transition_state(struc, threshold=50): #cm-1
    struc.calc = XTB(method="GFN2-xTB")
    
    try:
        vib = Vibrations(struc)
        vib.run()
        frequencies = vib.get_frequencies()
        vib.clean()
        
        # Filter out imaginary frequencies below the threshold
        significant_imaginary_freqs = np.count_nonzero(np.abs(np.imag(frequencies)) > threshold)

        return significant_imaginary_freqs == 1
    except:
        return False
    
def main(args):

    print_args(args)

    input_path  = args.input_path
    output_path = args.output_path
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    grown_seeds = [dirpath for dirpath, _, filenames in os.walk(input_path) if "converged" in filenames]
    grown_seeds_copy = grown_seeds
    # Group by mother string
    grouped_seeds = defaultdict(list)
    for seed in grown_seeds:
        mother_string = os.path.basename(seed)[:-8] # gsmGeom-m1-i1-c1-opt-gsm0044 -> gsmGeom-m1-i1-c1-opt
        grouped_seeds[mother_string].append(seed)
    rxn_list = []

    bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
    for mother_string, seeds in tqdm(grouped_seeds.items(), desc="Mothers", position=0, bar_format=bar_format, ncols=70):
        idx = 0
        for f in tqdm(seeds, desc=f"Rxns in {mother_string}", position=1, bar_format=bar_format, ncols=70, leave=False):

            ts_file_path = os.path.join(f, 'transition_state.xyz')
            if not os.path.exists(ts_file_path):
                continue
            
            struc = read(ts_file_path)
            if not is_transition_state(struc):
                #print(f"Directory {f} is not a genuine transition state. Skipping...")
                continue

            # If True, copy the directory
            new_name = os.path.join(output_path, f'{mother_string}-rxn{idx:04}')
            shutil.copytree(f, new_name)
            rxn_list.append(new_name)
            idx += 1

    with open(os.path.join(output_path, 'reactions.json'), 'w') as f:
        json.dump(rxn_list, f, indent=4)

    print(f'\n{len(rxn_list)}/{len(grown_seeds_copy)} rxns were saved to {output_path}/reactions.json')
    print('Filtering NEB finished!')


def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()

def get_parser():
    parser = argparse.ArgumentParser(description='Filter neb jobs and make reactions.json')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Input path of finished neb jobs')    
    parser.add_argument('-o', '--output_path', required=True, 
                        help='Output path of filtered neb jobs')

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)




