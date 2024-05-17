import os
import sys
import json
import logging
import argparse
import shutil, errno
from ase.atoms import Atoms
from tabulate import tabulate
from ase.io import read, write
from ase.io.jsonio import decode
from maggma.stores import MongoStore
import numpy as np


def get_docs(launchpad_file, query, collections_name='fireworks'):
    tasks_store = MongoStore.from_launchpad_file(launchpad_file, collections_name)
    tasks_store.connect()
    docs = list(tasks_store.query(query))
    return docs


def atomic_number_to_symbol(atomic_number):
    periodic_table = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B',
        6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    }
    return periodic_table.get(atomic_number, 'X')  # 'X' for unknown element


def atoms_from_json(atoms_json):
    atoms_dict = json.loads(atoms_json)

    atomic_numbers = atoms_dict['numbers']['__ndarray__'][2]
    positions = atoms_dict['positions']['__ndarray__'][2]

    natoms = int(len(positions)/3)
    coords = np.reshape(positions, (natoms, 3))

    atomic_symbols = [atomic_number_to_symbol(i) for i in atomic_numbers]
    return Atoms(symbols=atomic_symbols, positions=coords)


def download_data(docs, base_dir):
    for doc in docs:
        tag = doc['metadata']['tag']
        index = tag.split('-')[-1].split('_')[0]
        # print('tag:', tag)
        if tag.startswith('TS0'):
            name = 'TS'
            structure = atoms_from_json(doc['output']['atoms']['atoms_json'])
            '''
            print(
                  name,
                  structure,
                  structure.get_positions()
                 )
            '''
        elif tag.startswith('irc'):
            name = tag.split('-')[1]
            structure = atoms_from_json(doc['output']['atoms']['atoms_json'])
            '''
            print(
                  name,
                  structure,
                  structure.get_positions()
                 )
            '''
        else:
            continue

        output_dir = os.path.join(base_dir, f"{index}")

        os.makedirs(output_dir, exist_ok=True)

        output_file = os.path.join(output_dir, f'{name}.xyz')
        write(output_file, structure, format='xyz')


def main():
    tag = 'nov15'
    noise_level = '00'

    query0 = {
              "metadata.tags.class": tag,
               "metadata.tag": {
                        '$regex': f'.*0.*_noise{noise_level}'
                       }
             }

    launchPad_file = os.path.join(os.environ["HOME"],"fw_config/my_launchpad.yaml")
    collections = 'quacc_results0'

    docs0 = get_docs(launchPad_file, query0, collections_name=collections)
    download_data(docs0, "/global/scratch/users/kumaranu/saved_files/neb_nn_inputs/")

    '''
    for doc in docs0:
        # print(doc.keys())
        print(doc['metadata']['tag'])
        #a = 'irc-reverse0-021_noise00'
        index = a.split('-')[-1].split('_')[0]
        if a.split('-')[1] == 'reverse0':
            print('reverse')

        #os.system(f'cp -rL {src} {dst}')
    print(len(docs0))
    '''


if __name__ == "__main__":
    main()

