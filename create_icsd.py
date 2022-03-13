#!/usr/bin/env python3

from multiprocessing import Pool
import numpy as np
import wyckoff_list
import pychemia
import os
import json
import argparse
import time
import string
from re import findall
from collections import OrderedDict
import pandas as pd
import pymatgen.io.cif as cif

class ICSDDatabase(pychemia.db.PyChemiaDB):
    def __init__(
            self,
            name='icsd',
            host='localhost',
            port=27020,
            user=None,
            passwd=None,
            ssl=False,
            replicaset=None):
        super().__init__(name, host, port, user, passwd, ssl,
                         replicaset)

        self._symprec = 1e-5
        self._verbose = True
        self.paths = []

    def symprec(self, symprec=1e-5):
        self._symprec = symprec

    def verbose(self, verbose=True):
        self._verbose=verbose
        
    def exists(self, path):
        if not self.entries.find_one({'properties.path': path}) is None:
            return True
        else:
            return False

    def check_dirs(self, path='.'):
        dirs = os.listdir(path)
        for idir in dirs:
            if os.path.isdir(path + os.sep + idir):
                if contains_cif(path + os.sep + idir):
                    self.paths.append(path + os.sep + idir)
                self.check_dirs(path + os.sep + idir)

    def run(self, nproc=1):
        ndirs = len(self.paths)
        if self._verbose:
            print(f"Number of directories to be analyzed: {ndirs}. Number of processors: {nproc}")
        with Pool(nproc) as p:
            data = p.map(extract_data, self.paths)
            for datum in data:
                self.update_path(datum)
        return


    def update_path(self, data):
        entry = self.entries.find_one({'$and': [{"properties.code_ICSD":data['properties']['code_ICSD']},
                                                {"properties.standardized":data['properties']["standardized"]},
                                                {"properties.theoretical":data['properties']["theoretical"]}
            ]})
        
        if entry is None:
            self.insert(**data)
        else:
            self.update(entry['_id'],
                        structure=data['structure'],
                        properties=data['properties'],)
        return

    def to_json(self, filename="xc_db_analysis.json"):
        ret = []
        for ientry in self.entries.find():
            ientry['_id'] = str(ientry['_id'])
            ret.append(ientry)
        with open(filename, 'w') as wf:
            json.dump(
                ret,
                wf,
                sort_keys=True,
                indent=4,
                separators=(
                    ',',
                    ': '))
        return


def contains_cif(path):
    dirs = os.listdir(path)
    if  any(['.cif' in x for x in dirs]):
        return True
    else:
        return False
    

def extract_data(path, symprec=1e-5):
    '''This extract assumes there is only one file in each directory
    needs to change
    '''
    path = os.getcwd() + os.sep + path
    path = path.replace("{}.{}".format(os.sep, os.sep), os.sep)
    dirs = os.listdir(path)
    for ifile in dirs:
        if ".cif" in ifile:
            break
    print(path+os.sep+ifile)
    cif_item = cif.CifParser(path+os.sep+ifile)
    structure = cif_item.get_structures()[0]
    symbols = [site.element.value for site in structure.species]
    cell = structure.lattice.matrix
    reduced = structure.frac_coords
    structure = pychemia.core.Structure(cell=cell, symbols=symbols, reduced=reduced)
    cif_dict = cif_item.as_dict()
    properties = {'path': path}
    properties['cif'] = cif_dict
    properties['bib'] = cif_item.get_bibtex_string()
    icsd = [cif_dict[keys]['_database_code_ICSD'] for keys in cif_dict][0]
    properties['code_ICSD'] = icsd
    properties['standardized'] = 'experimental' not in ifile
    properties['theoretical'] = 'theoritical' in ifile
    crystal = pychemia.crystal.CrystalSymmetry(structure)
    mp_id = get_mp_id(path)
    crystal_family = crystal.crystal_system(symprec)
    properties['mp_id'] = mp_id
    properties['crystal_family'] = crystal_family
    properties['crystal_symmetry'] = crystal.get_space_group_type(symprec)
    properties['wyckoff_analysis'] =  get_wyckoffs(
        structure,
        symprec=symprec)
    properties['lattice_degrees_of_freedom'] = get_lattice_degrees_of_freedom(
        crystal_family.lower(), structure.lattice)
    # print(f"mp_id : {mp_id:<30}| icsd code  : {icsd:>30}")
    return {'structure': structure,
        'properties': properties,
        }


def lattice_to_dict(lattice):
    return {'a': lattice.a,
            'b': lattice.b,
            'c': lattice.c,
            'alpha': lattice.alpha,
            'beta': lattice.beta,
            'gamma': lattice.gamma,
            }
    
    
def get_wyckoffs(structure, symprec=1e-5):
    crystal_symmetry = pychemia.crystal.CrystalSymmetry(structure)
    spg = crystal_symmetry.number(symprec=symprec)
    wyckoffs = np.array(
        crystal_symmetry.get_symmetry_dataset(
            symprec=symprec)['wyckoffs'],
        dtype='<U3')
    representation = wyckoff_list.get_wyckoff(spg)
    dof = np.zeros((structure.natom, 3), dtype=int)
    take_wyckoff = np.zeros((structure.natom))
    take_lattice = np.zeros((structure.natom))
    take_lattice[0] = 1
    wyck_idx = {}
    i = 1
    for letter in string.ascii_lowercase:
        wyck_idx[letter] = -1 * i
        i += 1

    for ispc in structure.species:
        idx = np.array(structure.symbols) == ispc
        for wyckoff_letter in np.unique(wyckoffs[idx]):
            idx2 = wyckoffs[idx] == wyckoff_letter
            take_wyckoff[np.where(idx)[0][idx2][0]] = 1
    current_representation = []
    mmm = []
    www = []
    for iatom in range(structure.natom):
        wyckoff_letter = wyckoffs[iatom]
        www.append(wyckoff_letter)
        rep = representation[wyck_idx[wyckoff_letter]]
        multplicity = len(rep)
        mmm.append(multplicity)
        wyckoffs[iatom] = str(multplicity) + wyckoff_letter
        if take_wyckoff[iatom] == 1:
            important = False
            if 'x' in str(rep):
                dof[iatom, 0] = 1
                important = True
            if 'y' in str(rep):
                dof[iatom, 1] = 1
                important = True
            if 'z' in str(rep):
                dof[iatom, 2] = 1
                important = True
            if important:
                current_representation.append(rep)

    return {'wyckoffs': wyckoffs[take_wyckoff.astype(np.bool_)].tolist(),
            'symmetry_species': np.array(structure.symbols)[take_wyckoff.astype(np.bool_)].tolist(),
            'inner_degrees_of_freedom': dof.tolist(),
            'representation': current_representation,
            'multiplicity': mmm,
            'wyckoff_letter': www}


def get_lattice_degrees_of_freedom(family, lattice):
    a, b, c = [round(x, 5) for x in [lattice.a, lattice.b, lattice.c]]
    if family in ['cubic', 'trigonal']:
        if a==b and b==c:
            return [1, 0, 0]
        else:
            return [1/3, 1/3, 1/3]
    elif family in ['hexagonal', 'tetragonal']:
        if a == b or b == c:
            return [1, 0, 1]
        elif a == c :
            return [1, 1, 0]
        else:
            return [1/3, 1/3, 1/3]
    elif family in ['triclinic', 'monoclinic', 'orthorhombic']:
        return [1, 1, 1]


def get_mp_id(path):
    mp_id = ''
    for x in path.split(os.sep):
        if "mp" in x:
            mp_id = findall("mp-[0-9]*", x)[0]
    return mp_id

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--path",
        dest="path",
        type=str,
        help='Path where you want the analysis to be done',
        default='.')
    parser.add_argument("--name",
                        dest="name",
                        type=str,
                        help="MongoDB name",
                        default=os.getcwd().split(os.sep)[-1])
    parser.add_argument(
        "--host",
        dest="host",
        type=str,
        help="MongoDB hostname",
        default='localhost')
    parser.add_argument(
        "--port",
        dest="port",
        type=int,
        help="MongoDB port",
        default=27020)
    parser.add_argument(
        "--user",
        dest='user',
        type=str,
        help="MongoDB username",
        default=None)
    parser.add_argument(
        "--passwd",
        dest='passwd',
        type=str,
        help="MongoDB password",
        default=None)
    parser.add_argument(
        "--symprec",
        dest='symprec',
        type=float,
        help="distance threshold for symmtery analysis in Angstrom",
        default=1e-5)
    parser.add_argument(
        "--mode",
        dest='mode',
        nargs='+',
        help="what is expected from the script select from:\n update, export",
        default=['export'])
    parser.add_argument(
        "-np",
        dest="nproc",
        help="Number of processors to be used",
        type=int,
        default=1
        )
    parser.add_argument(
        "--verbose",
        action=argparse.BooleanOptionalAction,
        help='To print the results of the on going analysis.',
        default="--verbose")
    args = parser.parse_args()
    db_settings = dict(name=args.name,
                       host=args.host,
                       port=args.port,
                       user=args.name,
                       passwd=args.passwd)
    db = ICSDDatabase(
        args.name,
        args.host,
        args.port,
        args.user,
        args.passwd)
    db.symprec(args.symprec)
    db.verbose(args.verbose)
    for mode in args.mode:
        if mode in ['update', 'initiate']:
            db.check_dirs(args.path)
            db.run(args.nproc)
        elif mode == 'export':
            db.to_json()



        
