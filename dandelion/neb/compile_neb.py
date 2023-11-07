import os
import sys
import json
import hashlib
import argparse
import itertools

import h5py
import ase.db
import numpy as np
from tqdm import tqdm
from ase.units import Hartree, Bohr


def get_hash(row):
    s = str(row.positions) + row.formula
    return int(hashlib.sha1(s.encode("utf-8")).hexdigest(), 16) % (10 ** 8)

def write_rxn(h5file, fmaxs_path, db_path, rxn):
    fmaxs = json.load(open(fmaxs_path))

    skip_next = False
    first = True
    cum_fmax = 0

    with ase.db.connect(db_path) as db:
        for i, (fmax, path) in enumerate(zip(fmaxs, sliced_it(10, db.select("")))):
            cum_fmax += fmax
            skip_this = skip_next
            skip_next = False
            last = i == len(fmaxs) - 1

            if last:
                skip_this = False

            if cum_fmax < 0.1:
                skip_next = True

            else:
                cum_fmax = 0

            if skip_this:
                continue

            if not first:
                path = path[1:-1]

            # reactant and product is sampled once
            # (all points -2) // 8 ==0
            
            forces_path = np.array([row.forces for row in path])
            positions_path = np.array([row.positions for row in path])
            energy_path = np.array([row.energy for row in path])

            if first:
                forces = forces_path
                positions = positions_path
                energy = energy_path
                reactant = path[0]  # pylint: disable=undefined-loop-variable
                product = path[-1]  # pylint: disable=undefined-loop-variable

            else:
                forces = np.concatenate((forces, forces_path), axis=0)
                positions = np.concatenate((positions, positions_path), axis=0)
                energy = np.concatenate((energy, energy_path), axis=0)

            first = False

    transition_state = path[  # pylint: disable=undefined-loop-variable
        np.argmax(energy_path)
    ]

    formula = reactant.formula
    atomic_numbers = reactant.numbers

    if formula in h5file:
        grp = h5file[formula]
    else:
        grp = h5file.create_group(formula)

    subgrp = grp.create_group(rxn)
    single_molecule(reactant, subgrp.create_group("reactant"))
    single_molecule(transition_state, subgrp.create_group("transition_state"))
    single_molecule(product, subgrp.create_group("product"))

    dict_ = {
        "forces": forces,
        "positions": positions,
        "energy": energy,
        "atomic_numbers": atomic_numbers,
    }
    write_group(dict_, subgrp)


def single_molecule(molecule, subgrp):
    dict_ = {
        "forces": np.expand_dims(molecule.forces, 0),
        "positions": np.expand_dims(molecule.positions, 0),
        "energy": np.expand_dims(molecule.energy, 0),
        "atomic_numbers": molecule.numbers,
        "hash": get_hash(molecule),
    }
    write_group(dict_, subgrp)


def write_group(dict_, grp):
    grp.create_dataset("atomic_numbers", data=dict_["atomic_numbers"])
    grp.create_dataset("GFN2-xTB.forces", data=dict_["forces"])
    grp.create_dataset("GFN2-xTB.energy", data=dict_["energy"])
    grp.create_dataset("positions", data=dict_["positions"])

    if "hash" in dict_:
        grp.create_dataset("hash", data=dict_["hash"])


def sliced_it(n, iterable):
    it = iter(iterable)
    while True:
        chunk = itertools.islice(it, n)
        yield list(chunk)


def main(args):  
    
    print_args(args)
    
    input_path = args.input_path
    output_path = args.output_path
    
    try:
        rxns = json.load(open(input_path))
    except IsADirectoryError:
        print('Input path should include reactions.json!')
        sys.exit(1)
    h5file = h5py.File(output_path, "w")

    data = h5file.create_group("data")
    indexfile = open(output_path + ".index.json", "w")
    index = {}

    bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
    for i, path in tqdm(enumerate(rxns), total=len(rxns), desc="Compiling reactions", bar_format=bar_format, ncols=70):

        fmaxs_path = os.path.join(path, "fmaxs.json")
        db_path = os.path.join(path, "neb.db")

        new_rxn_name = f"rxn{str(i).zfill(4)}"
        write_rxn(data, fmaxs_path, db_path, new_rxn_name)
        index[new_rxn_name] = os.path.basename(path)

    json.dump(index, indexfile, indent=4)

    print('Compiling finished!')

def print_args(args):
    print()
    print("Arguments provided:")
    arg_dict = vars(args)
    for key, value in arg_dict.items():
        print(f"  {key}: {value}")
    print()

def get_parser():
    parser = argparse.ArgumentParser(description="Compile filtered neb jobs to xtb h5 file.")
    
    parser.add_argument('-i', '--input_path', required=True, 
                        help="Path of reactions.json, contains all reactions that should be included in the dataset ")
    parser.add_argument('-o', '--output_path', required=True, 
                        help="Path to the h5 file to write to")
    
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args)

