import os
import sys
import argparse
import subprocess

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def obabel_command(input_data, input_format, output_str, options=[], output_path=None):
    cmd = ['obabel', '-i', input_format] + input_data + ['-O', output_str] + options
    full_output_path = os.path.join(output_path, output_str) if output_path else output_str
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=output_path)
    return full_output_path

def obabel_from_smiles(smiles_str, output_str, options=[], output_path=None):
    cmd = ['obabel', '-ismi', '-', '-O', output_str] + options
    full_output_path = os.path.join(output_path, output_str) if output_path else output_str
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=output_path)
    process.communicate(input=smiles_str.encode())
    return full_output_path

def cleanup_files(output_path, files_to_remove):
    for file in files_to_remove:
        file_path = os.path.join(output_path, file)
        if os.path.exists(file_path):
            os.remove(file_path)

def main(args):
    print_args(args)
    
    input_path = os.path.abspath(args.input_path)
    if not os.path.isfile(input_path):
        sys.exit(f"Error: '{input_path}' is not a file.")
    output_path = os.path.abspath(args.output_path)
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    with open(input_path, 'r') as f:
        lines = f.readlines()
        lines = list(map(lambda s: s.strip(), lines))

    for m, mol_smi in enumerate(lines):
        print(f'==={m+1}th molecules : {mol_smi} ')
        mol = Chem.MolFromSmiles(mol_smi)
        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True)
        isomers = tuple(EnumerateStereoisomers(mol, options=opts))
        for i, isomer_smi in enumerate(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
            print(f'-{i+1}th isomer : {isomer_smi}')
            
            gen3d_file = obabel_from_smiles(isomer_smi, 'gen3d.xyz', ['--gen3d'], output_path=output_path)
            confab_file = obabel_command([os.path.basename(gen3d_file)], 'xyz', 'confab.sdf', ['--confab', '--rcutoff', '1.0'], output_path=output_path)
            obabel_command([os.path.basename(confab_file)], 'sdf', f'm{m+1}-i{i+1}-c.xyz', ['-m'], output_path=output_path)

    cleanup_files(output_path, ['confab.sdf', 'gen3d.xyz'])


def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()

def get_parser():
    parser = argparse.ArgumentParser(description="Generate Iso/Conformers from SMILES using RDkit and Obabel")
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help="Path of the input SMILES string file")
    parser.add_argument('-o', '--output_path', type=str, required=True,
                        help='Path of output directory to store Iso/Conformers.')
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)