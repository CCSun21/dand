import os
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm

# conda activate ts
# check whether gsm is killed when you interrupted
# use like "nohup python -u 2_run_gsm_jobs > gsm.out &"


def run_gsm_script(script_dir):
    #print(f"Executing in directory: {script_dir}")
    subprocess.run('bash gsm.sh', cwd=script_dir, capture_output=True, text=True, shell=True)

def main(args):
    
    print_args(args)
    
    input_path = args.input_path
    max_workers = args.max_workers

    # Find all directories containing gsm.sh scripts
    script_dirs = [dirpath for dirpath, _, filenames in os.walk(input_path) if "gsm.sh" in filenames]

    bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_gsm_script, script_dir) for script_dir in script_dirs] 
        
        for future in tqdm(as_completed(futures), desc='GSM on seeds', 
                           total=len(script_dirs), smoothing=0, bar_format=bar_format, ncols=70):
            pass # just update the tqdm

    print('GSM finished!')

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()
    
def get_parser():
    parser = argparse.ArgumentParser(description='Run GSM jobs concurrently')
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help='Base directory of mothers bearing seeds')
    parser.add_argument('-n', '--max_workers', type=int, default=1, 
                        help='Number of worker processes')

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)
