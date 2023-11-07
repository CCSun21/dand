#!/usr/bin/env python3
import os
import time
import argparse

from dandelion import __version__
from dandelion.refine.refine_forces import main as refine_forces, get_parser as refine_forces_parser
from dandelion.refine.compile_refined import main as compile_refined, get_parser as compile_refined_parser


def merge_args_with_defaults(module_parser, custom_args):
    """
    Merge custom arguments with module defaults.
    Args:
    - module_parser: the module parser function
    - custom_args: dictionary of custom arguments

    Returns:
    - argparse.Namespace: merged namespace of arguments
    """
    
    parser = module_parser()
    for action in parser._actions:
        if action.required:
            action.required = False

    defaults = vars(parser.parse_args([]))
    defaults.update(custom_args)

    for action in parser._actions:
        if not action.required and action.dest in custom_args:
            action.required = True

    return argparse.Namespace(**defaults)

def print_header(width=70):
    
    print(f'''                      

          ⢀⣀⣀⣀⣀⣀⡀       ⢀⢀⣀⢀⠞⠖⠁⠡⡂⡆ ⡠⢀⡀                     
         ⠺⢿⣿⣿⣿⣿⣿⣿⣷⣦⣠⣤⣤⣤⣄⣀⣀ ⡏⢸  ⢀ ⠣⠈ ⡠⡋⡨⡋⡂                  
           ⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣄⡀⡎⢀⡰⢀⢎⠌⢀⠔⣐⠠⣄⣀                
       ⢀ ⡔⢀⣴⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠿⠿⣿⣿⣷⣄⠂ ⢊⠎ ⠠⠂⡀⠕⠌⠌ ⡄⡠⢄            
    ⢀⡆⠄⠁⢈⢠⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣀   ⣀⣿⣿⣿⣆⠐    ⡨⠒⠁⡀⢠⣦⠍⠇⡀⢲⠂⡄⠄        
   ⠨⡀⠑⡈ ⢠⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡄   ⠈  ⣬⠠⣰⣿ ⢳⢹⡄⡆⠄⢀⢼       
 ⡄⠱⠈⠁⠑⢄⠐⣾⣿⣿⡿⠋⠁⣀⣠⣬⣽⣿⣿⣿⣿⣿⣿⠿⠿⠿⠿⠿⠿⠿⠿⠟⠁⡟⣅⡢⠁⠠⠜⡄⡑⢌⢧⡀ ⡀⣰⢁⡐⢁⢄⣡⣧⡤⠄   
⠠⡐⠓⠂⠌  ⢀⣿⣿⡏⢀⣴⣿⠿⠛⠉⠉⠶⢸⣿⣿⠿⠁⠢⠨⢀⣻⣿⣿⣿⣿⢟⣿⣝⠂  ⠠⡠⢆⠈⡂⠱⡇ ⣅⠫⠂⡠⢂⡪⠋  ⠁⡆  
⡶⠉ ⢀⡀⠁⡁⢸⣿⣿⢠⣾⡟⠁⣿⣿⡇ ⢀⠈⠉⠁    ⣀⠷⣹⣏⣷⢏⠹⠁    ⠈⢈ ⢇ ⢸⠱⢸⡏⡀⡶⡸⠎  ⠰⠁⡸   
⢈⡕⡈⠁⠐⠂⢀⢸⣿⣿⣾⠏⣿⣿⡿⣻⣿⢞⡢⠄ ⠈ ⡀⡤⠂⠁⠉⠌       ⢀⢀⠠⠐⢄ ⡀⢆⠎⢹⣶⣷⣧⡈⠈⠉⠤⠂⠉⢀⠱⡀ 
⢠⡊    ⠁⣸⣿⣿⣿⣀⠉⡻⡏⠋⠁ ⠁⠒⠒⡀⣍⠍⠁ ⡀ ⢠⠂     ⢀⠈⠄⢀⠄⡒⠅⠈⢄⢡ ⢿⣿⣷⣿⡄ ⠐⠄⠤ ⠜⢀ 
⠐⠁ ⠤⠒⢠⣾⣿⣿⣿⣿⣿⣷⣄⢄  ⢀ ⡏ ⢰⣃⠊⡐⠐⠁⢀⠈  ⣀ ⠰⠢⢀⠂⡰⠈⠂  ⡱⠂⢂⡇⡈⠻⢿⣿⠇   ⡤⠄⣀⡰⠁
    ⠁⣾⣿⣿⣿⣿⣿⣿⣿⣿⣦ ⠄ ⠉   ⠸⠫⢞⠈⣰⠈ ⡐⢲⣿⡏       ⢠⡾ ⣀⠊⢱ ⠠⡀    ⢈⢀⡐⠤⣕⡄
    ⢰⣿⡿⠛⠉   ⠈⠙⠛         ⠈⠈ ⠻⠔⠁⢸⡍⡇      ⢀⣏ ⢀⠠⠆ ⠣⡀⠈⡠⡀⠉⠢⡤⠢⣈⡡⣢⠦
⠈⠁           ⢻⣇               ⢸⡇⡇      ⣼⡿⠉  ⢀⡇ ⠑⡄⠑⣌⢄ ⠙⢄⠠⡪⣅ 
             ⠈⣾⡆              ⢸⣏⡇     ⢠⣿⠇   ⠸⢌⢢⢄⡠⠣⠈⠢⡁⡈⣎⢢⡬⠃ 

{"Energy refinement on samples using orca".center(width)}    
{("Ver. " + __version__  + " by mlee").center(width)}
''')

def print_separator(text, width=70):
    total_symbols_len = width - len(text) - 4  
    half_len = total_symbols_len // 2
    
    left_symbol = "║" + " " * (half_len - 1)
    right_symbol = " " * (total_symbols_len - half_len - 1) + "║"
    
    separator = left_symbol + '  ' + text + '  ' + right_symbol
    border = "╔" + "═" * (width-2) + "╗"
    end = "╚" + "═" * (width-2) + "╝"
    print("\n\n" + border)
    print(separator)
    print(end + "\n\n")


def main():
    args = parse_arguments()
    
    input_path = args.input_path
    max_workers = args.max_workers
    orcabinary = args.orca
    
    phases = [
        ("7. Refining forces", refine_forces, refine_forces_parser, {
            'input_path': os.path.join(input_path, 'xtb.h5'),
            'output_path': os.path.join(input_path, 'wb97x.db'),
            'orca' : orcabinary,
            'max_workers': max_workers            
        }),
        ("8. Compiling final samples", compile_refined, compile_refined_parser, {
            'input_path': os.path.join(input_path, 'wb97x.db'),
            'output_path': os.path.join(input_path, 'wb97x.h5')
        }),                  
    ]

    print_header()
    for title, function, parser, custom_args in phases:
        time.sleep(3)
        print_separator(title)
        merged_args = merge_args_with_defaults(parser, custom_args)
        function(merged_args)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Refine force on obtained samples,\
                                     Other parameters can be set in each modules')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Input path of working directory containing xtb.h5')    
    parser.add_argument('-n', '--max_workers', type=int, required=True, 
                        help='Number of worker processes')
    parser.add_argument('--orca', required=True, 
                        help="Path of the orca binary file")

    return parser.parse_args()


if __name__ == "__main__":
    main()
