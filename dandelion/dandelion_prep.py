#!/usr/bin/env python3
import os
import sys
import time
import argparse

from dandelion import __version__, print_separator, merge_args_with_defaults
from dandelion.prep.smiles_to_isoconfs import main as smiles_to_isoconfs, get_parser as smiles_to_isoconfs_parser
from dandelion.prep.geom_opt import main as geom_opt, get_parser as geom_opt_parser


def print_header(width=70):
    
    print(f'''

              H                 H                          
               \\\\               -                        
                 \\\\            -                         
                   C──────────C\                   H       
                  -              \\\\               /      
                 -                 \\\\            /       
       H────────C     O=Cc1ccccc1    C──────────C          
                 \\\\                 -            \\\\    
                   \\\\              -              \\\\   
                     \\C─────────C-                  O     
                      -           \\\\                     
                     -              \\\\                   
                    H                 H                    

{"Prepare Iso/Conformers from SMILES stirngs".center(width)}    
{("Dandelion " + __version__  + " by mlee").center(width)}
''')


def main():
    args = parse_arguments()
    
    input_path = args.input_path
    if not os.path.isfile(input_path):
        print(f"Error: The specified input path '{input_path}' is not a file.", file=sys.stderr)
        sys.exit(1)
    working_path = os.path.dirname(input_path)
    max_workers = args.max_workers

    phases = [
        ("1. Sample iso/conformers from SIMLES strings", smiles_to_isoconfs, smiles_to_isoconfs_parser, {
            'input_path': input_path,
            'output_path': os.path.join(working_path, '-1_isoconfs'),
        }),
        ("2. Optimize geometries", geom_opt, geom_opt_parser, {
            'input_path': os.path.join(working_path, '-1_isoconfs'),
            'output_path': os.path.join(working_path, '0_reactants'),
            'max_workers': max_workers
        }),                  
    ]

    print_header()
    
    for title, function, parser, custom_args in phases:
        time.sleep(3)
        print_separator(title)
        merged_args = merge_args_with_defaults(parser, custom_args)
        function(merged_args)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Prepare optimized iso/conformers from SMILES,\
                                     Other parameters can be set in each modules')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Input path of a file containing SMILES')    
    parser.add_argument('-n', '--max_workers', type=int, default=1,
                        help='Number of processes to use for parallel execution.')
    
    return parser.parse_args()


if __name__ == "__main__":
    main()
