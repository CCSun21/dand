import os
import sys
import glob
import argparse

from tqdm import tqdm 
from ase.db import connect

def main(args):
    
    print_args(args)
    
    input_path = args.input_path
    if not os.path.isdir(input_path):
        sys.exit(f"Error: '{input_path}' is not a directory.")
    output_path = args.output_path
    
    with connect(output_path) as db1:
        for f in glob.glob(os.path.join(input_path, '**/wb97x.db'), recursive=True):
            with connect(f) as db2:
                for row in tqdm(db2.select(), total=db2.count(), desc=f"{f}"):
                    db1.write(row.toatoms())

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()
    
def get_parser():
    parser = argparse.ArgumentParser(description='Merge db files in input directory')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Input path of directory containing db files to merge')    
    parser.add_argument('-o', '--output_path', required=True, 
                        help='Output path of the merged db file.')
    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)

