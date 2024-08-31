import os
import sys
import glob
import argparse

from ase.io import read
from .ard_gsm.mol import MolGraph
from .ard_gsm.limits import connection_limits
from .ard_gsm.driving_coords import generate_driving_coords


def main(args):
    
    print_args(args)
    
    input_path           = args.input_path
    if not os.path.isdir(input_path):
        sys.exit(f"Error: '{input_path}' is not a directory.")
    output_path          = args.output_path
    maxbreak             = args.maxbreak
    maxform              = args.maxform
    maxchange            = args.maxchange
    minbreak             = args.minbreak
    minform              = args.minform
    minchange            = args.minchange
    ignore_single_change = args.ignore_single_change
    equiv_Hs             = args.equiv_Hs

    pdir = output_path
    if not os.path.exists(pdir):
        os.makedirs(pdir)

    with open(os.path.join(pdir, 'params.log'), 'w') as f:
        f.write('Connection limits:\n')
        for symbol in connection_limits:
            ll = connection_limits[symbol][0]
            ul = connection_limits[symbol][1]
            f.write('  {}: {}, {}\n'.format(symbol, ll, ul))
        f.write('maxbreak = {}\n'.format(maxbreak))
        f.write('maxform = {}\n'.format(maxform))
        f.write('maxchange = {}\n'.format(maxchange))
        f.write('single_change = {}\n'.format(not ignore_single_change))
        f.write('equiv_Hs = {}\n'.format(equiv_Hs))
        f.write('minbreak = {}\n'.format(minbreak))
        f.write('minform = {}\n'.format(minform))
        f.write('minchange = {}\n'.format(minchange))

    # Loop over Mothers
    for idx, mother in enumerate(glob.iglob(os.path.join(input_path, '**/*.xyz'), recursive=True)):
        xyz = read(mother)
        symbols, coords = xyz.get_chemical_symbols(), xyz.get_positions()
        mol = MolGraph(symbols=symbols, coords=coords)
        mol.infer_connections()
        name =  os.path.basename(os.path.dirname(mother))

        seeds = generate_driving_coords(
            mol,
            maxbreak=maxbreak,
            maxform=maxform,
            maxchange=maxchange,
            single_change=not ignore_single_change,
            equiv_Hs=equiv_Hs,
            minbreak=minbreak,
            minform=minform,
            minchange=minchange
        )
        print(f'{len(seeds)} Seeds were generated from {name}')

        output_path = os.path.join(pdir, '{}'.format(name))
        if not os.path.exists(output_path):
            os.mkdir(output_path)

        # Loop over seeds
        for idx, seed in enumerate(seeds):
            
            gsm_dir = os.path.join(output_path, f'gsm{idx:04}')
            if not os.path.exists(gsm_dir):
                os.mkdir(gsm_dir)
                
            isomers_file = os.path.join(gsm_dir, 'ISOMERS.txt')
            initial_file = os.path.join(gsm_dir, 'initial.xyz')
            bash_file = os.path.join(gsm_dir, 'gsm.sh')
            
            with open(bash_file, 'w') as f:
                f.write('''
gsm -xyzfile initial.xyz \\
    -mode SE_GSM \\
    -num_nodes 30 \\
    -package xTB_lot \\
    -isomers ISOMERS.txt \\
    -xyz_output_format multixyz \\
    -coordinate_type DLC > gsm_log 2>&1''')        
            
            with open(isomers_file, 'w') as f:
                f.write(str(seed))
            with open(initial_file, 'w') as f:
                f.write(str(len(symbols)) + '\n')
                f.write('\n')
                for symbol, xyz in zip(symbols, coords):
                    f.write('{0}  {1[0]: .10f}  {1[1]: .10f}  {1[2]: .10f}\n'.format(symbol, xyz))

    print('\nCreating GSM finished!')                    

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()
    
def get_parser():
    parser = argparse.ArgumentParser(description='Make GSM jobs from mother structures')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Input path of mother structures')    
    parser.add_argument('-o', '--output_path', required=True, 
                        help='Output path of gsm jobs')
    
    parser.add_argument('--maxbreak', type=int, default=2,
                        help='Maximum number of connections to break')
    parser.add_argument('--maxform', type=int, default=2,
                        help='Maximum number of connections to form')    
    parser.add_argument('--maxchange', type=int, default=3,
                        help='Maximum number of connections to change')
    
    parser.add_argument('--minbreak', type=int, default=0,
                        help='Minumum number of connections to break')
    parser.add_argument('--minform', type=int, default=0,
                        help='Minumum number of connections to form')    
    parser.add_argument('--minchange', type=int, default=1,
                        help='Minumum number of connections to change')
    
    parser.add_argument('--ignore_single_change', type=bool, default=True,
                        help='Do not consider single connection changes (e.g., nbreak=1, nform=0)')    
    parser.add_argument('--equiv_Hs', type=bool, default=False,
                        help='Create equivalent driving coordinates for the same reaction with different but\
                        equivalent hydrogens, i.e., hydrogens attached to non-cyclic tetrahedral carbons')
    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)



