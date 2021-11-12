#!/usr/bin/env python3

import numpy as np
import wyckoff_list
import pychemia
import os
import json
import argparse
import time
import string
from re import findall


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
        if self.calc_type == 'relax':
            self.strict = True
        else:
            self.strict = False

    def exists(self, path):
        if not self.entries.find_one({'properties.path': path}) is None:
            return True
        else:
            return False

    def contains_calculation(self, path):
        dirs = os.listdir(path)
        if not self.strict:
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

    def check_dirs(self, path='.'):
        dirs = os.listdir(path)
        if "DO_NOT_ANALYZE" in [x.upper() for x in dirs]:
            return
        for idir in dirs:
            if os.path.isdir(path + os.sep + idir):
                if self.contains_calculation(path + os.sep + idir):
                    self.update_path(path + os.sep + idir)
                self.check_dirs(path + os.sep + idir)
        return

    def extract_data(self, path, symprec=1e-5):
        path = os.getcwd() + os.sep + path
        path = path.replace("{}.{}".format(os.sep, os.sep))
        dirs = os.listdir(path)
        if "DON_NOT_ANALYZE" in [x.upper() for x in dirs]:
            return
        init_poscar = "POSCAR"
        for ifile in dirs:
            if '.vasp' in ifile:
                init_poscar = ifile
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
            contcar = pychemia.code.vasp.read_poscar(
                path + os.sep + "CONTCAR").to_dict
            properties['structure_relaxed'] = contcar
        if os.path.exists(
                path +
                os.sep +
                "vasprun.xml") and status['status']['relaxed']:
            vasprun = pychemia.code.vasp.VaspXML(path + os.sep + "vasprun.xml")
            properties['band_gap'] = vasprun.band_gap
        mp_id = get_mp_id(path)
        properties['mp_id'] = mp_id
        properties['time'] = time.ctime(os.path.getmtime(path))
        properties['crystal_family'] = crystal.crystal_system()
        properties['lattice_degrees_of_freedom'] = get_lattice_degrees_of_freedom(
            properties['crystal_family'].lower())
        properties['crystal_symmetry'] = crystal.get_space_group_type()
        properties['wyckoff_analysis'] = get_wyckoffs(
            structure, symprec=symprec)

        return {'structure': structure,
                'properties': properties,
                'status': status}

    def update_path(self, path, symprec=1e-5):
        data = self.extract_data(path, symprec)
        entry = self.entries.find_one({'properties.path': path})
        if entry is None:
            self.insert(**data)
        else:
            self.update(entry,
                        structure=data['structure'],
                        properties=data['properties'],
                        status=data['status'])


def check_finished(path):
    if os.path.exists(path + os.sep + 'OUTCAR'):
        with open(path + os.sep + 'OUTCAR') as rf:
            text = rf.read()
        if "I REFUSE TO CONTINUE WITH THIS SICK JOB" in text or "copy CONTCAR" in text:
            return False
        elif "reached required accuracy" in text and "Major page faults" in text:
            return True
    else :
        return False


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
    for iatom in range(structure.natom):
        wyckoff_letter = wyckoffs[iatom]
        rep = representation[wyck_idx[wyckoff_letter]]
        multplicity = len(rep)
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
            'representation': current_representation}


def get_lattice_degrees_of_freedom(family):
    if family in ['cubic', 'trigonal']:
        return 1
    elif family in ['hexagonal', 'tetragonal']:
        return 2
    elif family in ['triclinic', 'monoclinic', 'orthorhombic']:
        return 3


def get_mp_id(path):
    mp_id = ''
    for x in path.split(os.sep):
        if "mp" in x:
            mp_id = findall("mp-[0-9]*", x)[0]
    return mp_id


def get_report(path):
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
    print(
        "path : {: >30}| kpoint : {: >5} | encut : {: >5}| relax : {: >5}|".format(
            path.split('/./')[-1],
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
        help="Length threshold for symmtery analysis",
        default=1e-5)

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
    db.check_dirs(args.path)
            