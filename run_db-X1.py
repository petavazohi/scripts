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

class CalcDatabase(pychemia.db.PyChemiaDB):
    def __init__(
            self,
            name='calc_database',
            host='localhost',
            port=27020,
            user=None,
            passwd=None,
            ssl=False,
            replicaset=None,
            calc_type='relax'):
        super().__init__(name, host, port, user, passwd, ssl,
                         replicaset)

        self.calc_type = calc_type
        self._symprec = 1e-5
        self._verbose = True
        if self.calc_type == 'relax':
            self.strict = True
        else:
            self.strict = False
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
        if "DO_NOT_ANALYZE" in [x.upper() for x in dirs]:
            return
        for idir in dirs:
            if os.path.isdir(path + os.sep + idir):
                if contains_calculation(path + os.sep + idir, self.strict):
                    self.paths.append(path + os.sep + idir)
                    # self.update_path(path + os.sep + idir)
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
        entry = self.entries.find_one({'properties.path': data['properties']['path']})
        if entry is None:
            self.insert(**data)
        else:
            self.update(entry['_id'],
                        structure=data['structure'],
                        properties=data['properties'],
                        status=data['status'])
        return

    def to_json(self, filename="xc_db_analysis.json"):
        ret = []
        for ientry in self.entries.find({'status.status.relaxed': True}):

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


def contains_calculation(path, strict=True):
    dirs = os.listdir(path)
    if not strict:
        if "POSCAR" in dirs or any(['.vasp' in x for x in dirs]):
            return True
        else:
            return False
    else:
        if ("KPOINTS" in dirs and
            "INCAR" in dirs and
            "POSCAR" in dirs and
            "POTCAR" in dirs):
            return True
        else:
            return False

    

def extract_data(path, symprec=1e-5):
    path = os.getcwd() + os.sep + path
    path = path.replace("{}.{}".format(os.sep, os.sep), os.sep)
    dirs = os.listdir(path)
    if "DON_NOT_ANALYZE" in [x.upper() for x in dirs]:
        return

    init_poscar = ""
    for ifile in dirs:
        if 'init' in ifile:
            init_poscar = ifile
    if init_poscar == "" :
        init_poscar = "POSCAR"

        
    structure = pychemia.code.vasp.read_poscar(path + os.sep + init_poscar)
    crystal = pychemia.crystal.CrystalSymmetry(structure)
    properties = {'path': path}
    status = get_report(path)
    if os.path.exists(path + os.sep + "INCAR"):
        incar = pychemia.code.vasp.read_incar(
            path + os.sep + "INCAR").variables
        properties['incar'] = incar
    if os.path.exists(path + os.sep + "KPOINTS"):
        kpoints = pychemia.code.vasp.read_kpoints(
            path + os.sep + "KPOINTS").to_dict
        properties['kpoints'] = kpoints
    if os.path.exists(
            path +
            os.sep +
            "CONTCAR") and status['status']['relaxed']:
        final_structure = pychemia.code.vasp.read_poscar(
            path + os.sep + "CONTCAR")
        structure = final_structure
        final_structure.reduced[final_structure.reduced.round(3) >= 1] -= 1
        final_structure.reduced[final_structure.reduced.round(3) < 0] += 1
        final_crystal = pychemia.crystal.CrystalSymmetry(final_structure)
        final_crystal_family = final_crystal.crystal_system(symprec)
        properties['final'] = {
            'structure': final_structure.to_dict,
            'lattice': lattice_to_dict(final_structure.lattice),
            'crystal_family': final_crystal_family,
            'crystal_symmetry': final_crystal.get_space_group_type(symprec),
            'wyckoff_analysis': get_wyckoffs(
                structure,
                symprec=symprec),
            'lattice_degrees_of_freedom': get_lattice_degrees_of_freedom(
                final_crystal_family.lower(), final_structure.lattice),
            }
    if os.path.exists(
            path +
            os.sep +
            "vasprun.xml") and status['status']['relaxed']:
        vasprun = pychemia.code.vasp.VaspXML(path + os.sep + "vasprun.xml")
        if not vasprun.has_diverged:
            properties['band_gap'] = vasprun.band_gap
            gap = vasprun.band_gap['total']['gap']
            if  gap == 0 :
                properties['conductivity'] = "conductor"
            elif gap > 0 and gap <= 1 :
                properties['conductivity'] = "semiconductor"
            else:
                properties['conductivity'] = "insulator"
            properties['xc'] = get_xc(vasprun.potcar_info, incar)
        else :
            status['status']['relaxed'] = False
    mp_id = get_mp_id(path)
    properties['mp_id'] = mp_id
    properties['time'] = time.ctime(os.path.getmtime(path))
    if init_poscar:
        init_structure = pychemia.code.vasp.read_poscar(
            path + os.sep + init_poscar)
        init_structure.reduced[init_structure.reduced.round(3) >= 1] -= 1
        init_structure.reduced[init_structure.reduced.round(3) < 0] += 1
        init_crystal = pychemia.crystal.CrystalSymmetry(init_structure)
        init_crystal_family = init_crystal.crystal_system(symprec)
        properties['initial'] = {
            'structure': init_structure.to_dict,
            'lattice': lattice_to_dict(init_structure.lattice),
            'crystal_family': init_crystal_family,
            'crystal_symmetry': init_crystal.get_space_group_type(symprec)}
    return {'structure': structure,
            'properties': properties,
            'status': status['status']}


def check_finished(path):
    if os.path.exists(path + os.sep + 'OUTCAR'):
        with open(path + os.sep + 'OUTCAR') as rf:
            text = rf.read()
        if "I REFUSE TO CONTINUE WITH THIS SICK JOB" in text or "copy CONTCAR" in text:
            return False
        elif "reached required accuracy" in text and "Major page faults" in text:
            return True
    else:
        return False


def lattice_to_dict(lattice):
    return {'a': lattice.a,
            'b': lattice.b,
            'c': lattice.c,
            'alpha': lattice.alpha,
            'beta': lattice.beta,
            'gamma': lattice.gamma,
            }
    
def get_xc(potcar_info, incar):
    psp = [potcar_info[x]['pseudopotential'] for x in potcar_info]
    has_pbe = any(['PBE' in x for x in psp])
    if not has_pbe:
        return "LDA"
    else:
        if "GGA" not in incar and 'METAGGA' not in incar:
            return "PBE"
        elif "GGA" in incar:
            if incar['GGA'] == 'PS':
                return "PBEsol"
            elif incar['GGA'] == 'AM':
                return 'AM05'
        elif "METAGGA" in incar:
            return incar['METAGGA']

    
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


def get_report(path, verbose=True):
    # verbose=False
    report = {}
    report['status'] = {}
    if os.path.exists(path + os.sep + "kpoint_report.json"):
        # with open(path + os.sep + "kpoint_report.json", 'r') as rf:
        #     kp = json.load(rf)
        #     report['kpoint'] = kp["output"]
        report['status']['kpoint_converged'] = True
    else:
        report['status']['kpoint_converged'] = False
    if os.path.exists(path + os.sep + "encut_report.json"):
        # with open(path + os.sep + "encut_report.json", 'r') as rf:
        #     ecut = json.load(rf)
        #     report['ecut'] = ecut["output"]
        report['status']['encut_converged'] = True
    else:
        report['status']['encut_converged'] = False
    if os.path.exists(
            path +
            os.sep +
            "relax_report.json") and check_finished(path):
        # with open(path + os.sep + "relax_report.json", 'r') as rf:
        #     relax = json.load(rf)
        #     report['relax'] = relax["output"]
        report['status']['relaxed'] = True
    else:
        report['status']['relaxed'] = False
    if verbose:
        print(
            "path : {: <50}| kpoint : {: >5} | encut : {: >5}| relax : {: >5}|".format(
                (os.sep).join(path.split(os.sep)[-4:]),
                report['status']['kpoint_converged'],
                report['status']['encut_converged'],
                report['status']['relaxed']))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--path",
        dest="path",
        type=str,
        help='Path where you want the analysis to be done',
        default='.')
    parser.add_argument(
        "--job",
        dest="job_type",
        type=str,
        help="Type of the jobs needed to be calculated choose from:\n relax, ebs, dos, elastic, phonon",
        default='relax')
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
    db = CalcDatabase(
        args.name,
        args.host,
        args.port,
        args.user,
        args.passwd,
        calc_type=args.job_type)
    db.symprec(args.symprec)
    db.verbose(args.verbose)
    for mode in args.mode:
        if mode in ['update', 'initiate']:
            db.check_dirs(args.path)
            db.run(args.nproc)
        elif mode == 'export':
            db.to_json()

