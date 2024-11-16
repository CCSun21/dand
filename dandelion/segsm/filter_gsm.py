import os
import sys
import glob
import shutil
import argparse

from rdkit import RDLogger
from ase.io import read, write
from openbabel import openbabel
from .ard_gsm.mol import MolGraph, SanitizationError
#from ard_gsm.mol import MolGraph, SanitizationError

# Suppress Noisy warning in the filter
RDLogger.logger().setLevel(RDLogger.CRITICAL)
openbabel.obErrorLog.SetOutputLevel(openbabel.obError)

'''
Faith of pyGSM run

1) png is not made
- xTB not converge
- pyGSM suicide on his criteria

2) png is made
- Exiting early -> should filter out
- Ran out of iterations -> also includes potential rxn
- Converged -> very rare
'''



def parse_gsm_log(keyword, content):
    """Find the value associated with a keyword in a text content."""
    # For TS_energy, we're expecting a float, so we use a different pattern
    if keyword == "TS energy:":
        pattern = f"{keyword} ([+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)"
    else:
        pattern = f"{keyword} (\d+)"
    
    import re
    matches = re.findall(pattern, content)
    
    # Return the matched value; assume there's only one match
    if matches:
        return matches[0][0]  # Due to group structures, we take the first element
    else:
        return None


def get_gsm_data(home, seed, string):
    try:
        with open(os.path.join(home, seed, string, 'gsm_log'), 'r') as f:
            content = f.read()
    except FileNotFoundError:
        return None

    nodes = []
    try:
        with open(os.path.join(home, seed, string, 'opt_converged_000.xyz'), 'r') as f:
            for i in range(30):
                try:
                    nodes.append(read(f, i))
                except:
                    break
    except FileNotFoundError:
        return None

    return {
        "TS_energy"         : float(parse_gsm_log("TS energy:", content)),
        "reactant_idx"      : int(parse_gsm_log("min reactant node:", content)),
        "product_idx"       : int(parse_gsm_log("min product node", content)),
        "TS_idx"            : int(parse_gsm_log("TS node is", content)),
        "nodes"             : nodes,
        'energies'          : [float(list(node.info.keys())[0]) for node in nodes]
    }



def profile_filter(strings, home, seed, barrier_max, barrier_min, delta_e_min):
    '''
    Given gsm success reactions,
    Filter strings by TS_index and Barrier height and delta_e.
    '''
    filtered = {}
    for string in strings:
        data = get_gsm_data(home, seed, string)
        if not data:
            continue

        if data["TS_idx"] >= data["product_idx"]: # wrong ts
            continue
        if (data["TS_energy"] > barrier_max) or (data["TS_energy"] < barrier_min): # too high or low barrier
            continue
        if abs(data['energies'][data['product_idx']]) * 627.503 < delta_e_min: # maybe reactant==product
            continue

        product_graph = MolGraph(symbols=data["nodes"][data["product_idx"]].get_chemical_symbols(),
                                  coords=data["nodes"][data["product_idx"]].get_positions(),
                                  energy=float(list(data["nodes"][data["product_idx"]].info.keys())[0]))
    
        filtered[string] = {
            'reactant': data["nodes"][data["reactant_idx"]],
            'product': data["nodes"][data["product_idx"]],
            'ts': data["nodes"][data["TS_idx"]],
            'product_graph': product_graph,
            'ts_energy': data["TS_energy"]
        }

    return filtered

def structure_filter(reactions):
    '''
    Chemically absurd products are filtered here. (graph->pybel->inchi->smiles)
    SMILES are constructed, and saved to the dictionary for the unique filter.
    '''

    filtered = {}
    
    for rxn, data in reactions.items():
        try:
            smiles = data['product_graph'].perceive_smiles()
            filtered[rxn] = data
            filtered[rxn]['product_smiles'] = smiles
        except SanitizationError:
            continue
    return filtered

def unique_filter(reactions):
    '''
    Duplicates are filtered based on SMILES. 
    If there are more than one of same SMILES, pick the lowest barrier reaction.
    '''
    unique = {}
    for rxn, data in reactions.items():
        smiles = data['product_smiles']
        ts_energy = data['ts_energy']
        if smiles not in unique or ts_energy < unique[smiles]['ts_energy']:
            unique[smiles] = {
                'reaction_key': rxn,
                'ts_energy': ts_energy,
                'reactant': data['reactant'],
                'product': data['product'],
                'ts': data['ts'],
            }
    return unique

def save_unique_reactions(home, output_path, seed, reactions):
    for smiles, data in reactions.items():
        reaction_dir = os.path.join(output_path, seed, data['reaction_key'])
        os.makedirs(reaction_dir, exist_ok=True)

        file_types = ["reactant", "ts", "product"]
        for f_type in file_types:
            write(os.path.join(reaction_dir, f"{f_type}.xyz"), data[f_type])
            write(os.path.join(reaction_dir, f"{f_type}.png"), data[f_type])

        shutil.copyfile(os.path.join(home, seed, data['reaction_key'], '0000_string.png'),
                        os.path.join(reaction_dir, 'string.png'))

        shutil.copyfile(os.path.join(home, seed, data['reaction_key'], 'opt_converged_000.xyz'),
                        os.path.join(reaction_dir, 'string.xyz'))
        
def main(args):
    
    print_args(args)
    
    input_path  = args.input_path
    if not os.path.isdir(input_path):
        sys.exit(f"Error: '{input_path}' is not a directory.")
    output_path = args.output_path
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    barrier_max = args.barrier_max
    barrier_min = args.barrier_min
    delta_e_min = args.delta_e_min
    
    mothers =  [d for d in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, d))]
    for mother in mothers:
        print('\n◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢◤◢')
        print(f'mother: {mother}'.center(35))
        driving_coordinates = list(glob.iglob(os.path.join(input_path, f'{mother}/*/gsm_log')))
        success_strings = [path.split('/')[-2] for path in glob.iglob(os.path.join(input_path, f'{mother}/*/0000_string.png'))]

        profile_filtered_strings = profile_filter(success_strings, input_path, mother, barrier_max, barrier_min, delta_e_min)
        structure_filtered_strings = structure_filter(profile_filtered_strings)
        unique_reactions = unique_filter(structure_filtered_strings)

        print(f'Initial seeds:                 {len(driving_coordinates):>5}')
        print(f'GSM success reactions:         {len(success_strings):>5}')
        print(f'Profile filtered reactions:    {len(profile_filtered_strings):>5}')
        print(f'Structure filtered reactions:  {len(structure_filtered_strings):>5}')
        print(f'Unique reactions:              {len(unique_reactions):>5}')

        save_unique_reactions(input_path, output_path, mother, unique_reactions)

    print('\nFiltering GSM finished!')

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
                        help='Input path of finished gsm jobs')    
    parser.add_argument('-o', '--output_path', required=True, 
                        help='Output path of filtered gsm jobs')
    
    parser.add_argument('--barrier_min', type=int, default=5)
    parser.add_argument('--barrier_max', type=int, default=200)    
    parser.add_argument('--delta_e_min', type=int, default=5)

    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)


