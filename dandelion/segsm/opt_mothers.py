import os
import sys
import shutil
import argparse

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from xtb.ase.calculator import XTB


def generate_eq_strucs(eq_struc):
    calc = XTB(method="GFN2-xTB", atoms=eq_struc)
    eq_struc.calc = calc
    opt = BFGS(eq_struc, logfile=None)
    opt.run(fmax=1e-4)
    
    return eq_struc
  

def main(args):
    
    print_args(args)
    
    input_path  = args.input_path
    output_path = args.output_path

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    

    for root, dirs, files in os.walk(input_path):
        for filename in files:
            if filename.endswith('.xyz'):
                
                #Get the premature xyz
                premature = read(os.path.join(root, filename))
                print(f'Handling {os.path.join(root, filename)}')
                            
                try:
                    eq_struc = generate_eq_strucs(premature)
                except:
                    print(f'{os.path.join(root, filename)} not converged')
                    continue
                
                print(f'{os.path.join(root, filename)} converged.')
                print(f'saved to {os.path.join(output_path, os.path.dirname(filename), "struc.xyz")}\n')
                write(os.path.join(output_path, os.path.dirname(filename), 'struc.xyz'), eq_struc)
                
def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()
    
def get_parser():
    parser = argparse.ArgumentParser(description='Optimize crude structures to generate mother structures')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Input path of crude sturctures')    
    parser.add_argument('-o', '--output_path', required=True, 
                        help='Output path of mother structures')

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
