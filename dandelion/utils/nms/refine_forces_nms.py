import os
import argparse
from itertools import repeat
from concurrent.futures import ProcessPoolExecutor, as_completed

import h5py
import logging
from tqdm import tqdm
from ase import Atoms
from ase.db import connect
from ase.calculators.orca import ORCA


# thank you for https://github.com/ZimmermanGroup/ORCA-Basis-Sets
custom_basis = '''
%basis
 newgto Br
 S   6
 1         0.1137182000E+06       0.1717696000E-02
 2         0.1707444000E+05       0.1316744000E-01
 3         0.3889576000E+04       0.6504553000E-01
 4         0.1097096000E+04       0.2269505000E+00
 5         0.3520624000E+03       0.4768357000E+00
 6         0.1207002000E+03       0.3583677000E+00
 S   6
 1         0.2471138000E+04       0.2243687000E-02
 2         0.5893838000E+03       0.2994853000E-01
 3         0.1918738000E+03       0.1256009000E+00
 4         0.7295339000E+02      -0.9832786000E-03
 5         0.3005839000E+02      -0.6013141000E+00
 6         0.1252927000E+02      -0.4913983000E+00
 P   6
 1         0.2471138000E+04       0.3790182000E-02
 2         0.5893838000E+03       0.2995979000E-01
 3         0.1918738000E+03       0.1318228000E+00
 4         0.7295339000E+02       0.3432708000E+00
 5         0.3005839000E+02       0.4642345000E+00
 6         0.1252927000E+02       0.2079387000E+00
 S   6
 1         0.1096411000E+03      -0.5975683000E-02
 2         0.3858948000E+02       0.5542122000E-01
 3         0.1637818000E+02       0.2681200000E+00
 4         0.7221836000E+01      -0.1543606000E+00
 5         0.3263697000E+01      -0.7206306000E+00
 6         0.1465499000E+01      -0.3316437000E+00
 P   6
 1         0.1096411000E+03      -0.6907483000E-02
 2         0.3858948000E+02      -0.3041432000E-01
 3         0.1637818000E+02       0.4602725000E-01
 4         0.7221836000E+01       0.3650689000E+00
 5         0.3263697000E+01       0.4949232000E+00
 6         0.1465499000E+01       0.2090394000E+00
 S   3
 1         0.2103651000E+01       0.3029029000E+00
 2         0.7547050000E+00      -0.2152659000E+00
 3         0.3005140000E+00      -0.9633941000E+00
 P   3
 1         0.2103651000E+01      -0.2826714000E-01
 2         0.7547050000E+00       0.3503065000E+00
 3         0.3005140000E+00       0.7182446000E+00
 S   1
 1         0.1090710000E+00       0.1000000000E+01 
 P   1
 1         0.1090710000E+00       0.1000000000E+01
 D   3
 1         0.6225514000E+02       0.7704229000E-01
 2         0.1731284000E+02       0.3707384000E+00
 3         0.5607915000E+01       0.7097628000E+00
 D   1
 1         0.1746486000E+01       1.0000000
 end
end
'''



class tqdm_hour(tqdm):
    """Provides an `hours per iteration` format parameter."""
    @property
    def format_dict(self):
        d = super(tqdm_hour, self).format_dict
        rate_hr = '{:.1f}'.format(1/d["rate"] / 3600) if d["rate"] else '?'
        d.update(rate_hr=(rate_hr + ' hour/' + d['unit']))
        return d

class tqdm_minute(tqdm):
    """Provides a `minutes per iteration` format parameter"""
    @property
    def format_dict(self):
        d = super(tqdm_minute, self).format_dict
        rate_min = '{:.0f}'.format(1/d["rate"] / 60) if d["rate"] else '?'
        d.update(rate_min=(rate_min + ' min/' + d['unit']))
        return d

bar_format_hr = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_hr}{postfix}]'
bar_format_min = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_min}{postfix}]'
bar_format_points = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'

def get_unique_ids_from_db(output_path):
    """Extract all unique IDs from the ASE database."""
    unique_ids = set()
    with connect(output_path) as db: 
        for row in db.select():
            data = row.data
            if "unique_id" in data:
                unique_ids.add(data['unique_id'])
    return unique_ids


def already_calculated(unique_id, unique_id_list):
    """Check if a unique ID has already been processed."""
    return unique_id in unique_id_list

def compute_force(coord, atomic_numbers, unique_id, output_path):
    """Compute forces using ORCA for a given set of coordinates."""
    atoms = Atoms(positions=coord, numbers=atomic_numbers)
    atoms.calc = ORCA(
        label=os.path.join(os.path.dirname(output_path), f"orca/{unique_id}/{unique_id}"),
        orcasimpleinput="wB97X 6-31G(d) NoTrah",
        orcablocks=custom_basis
    )
    try:
        atoms.get_forces()  # Forces and energy will be stored in the calculator of the Atoms object.
        return atoms
    except: # not converged
        return None

def delete_orca_results(unique_id, output_path, file_exts=['gbw', 'engrad', 'densities', 'ase']):
    """Delete specific file types based on the unique ID."""
    for ext in file_exts:
        file_path = os.path.join(os.path.dirname(output_path), f"orca/{unique_id}/{unique_id}.{ext}")
        if os.path.exists(file_path):
            os.remove(file_path)
            
def main(args):
    print_args(args)
    
    input_path = args.input_path
    if not os.path.isfile(input_path):
        sys.exit(f"Error: '{input_path}' is not a file.")
    output_path = args.output_path
    max_workers = args.max_workers
    orcabinary = args.orca

    os.environ["ASE_ORCA_COMMAND"] = f"{orcabinary} PREFIX.inp > PREFIX.out 2>&1"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    logging.basicConfig(filename=os.path.join(os.path.dirname(output_path), 'orca_errors.log'), 
                        level=logging.ERROR, 
                        format='%(asctime)s %(levelname)s: %(message)s', 
                        datefmt='%Y-%m-%d %H:%M:%S')

    unique_ids_from_db = get_unique_ids_from_db(output_path)

    tasks = []
    with connect(input_path) as input_db:
        for row in tqdm(input_db.select(), desc="Reading ASE DB"):
            unique_id = row.unique_id
            if already_calculated(unique_id, unique_ids_from_db):
                continue
            atoms = row.toatoms()
            tasks.append((atoms.positions, atoms.numbers, unique_id, output_path))

    results_batch = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_atoms = {executor.submit(compute_force, *task): task[2] for task in tasks}

        for future in tqdm(as_completed(future_to_atoms), total=len(future_to_atoms), desc="Calculating", leave=False):
            unique_id = future_to_atoms[future]
            atoms_result = future.result()
            if atoms_result is not None:
                results_batch.append((atoms_result, {'unique_id': unique_id}))

            if len(results_batch) >= max_workers:
                with connect(output_path) as output_db:
                    for atoms, data in results_batch:
                        output_db.write(atoms, data=data)
                results_batch.clear()

    # Write any remaining results
    with connect(output_path) as output_db:
        for atoms, data in results_batch:
            output_db.write(atoms, data=data)

    print('Calculations completed and saved to the output ASE database!')

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()

def get_parser():
    parser = argparse.ArgumentParser(description="Compute energies and forces and store in ASE database")
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help="Path of the input XTB HDF5 file")
    parser.add_argument('-o', '--output_path', required=True, 
                        help="Path of the output wB97X ASE database")
    parser.add_argument('-n', '--max_workers', type=int, default=1, 
                        help="Number of worker processes")
    parser.add_argument('--orca', required=True, 
                        help="Path of the orca binary file")

    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)

