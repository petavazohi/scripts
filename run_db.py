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
        # print(path)
        path = path.replace("{}.{}".format(os.sep, os.sep), os.sep)
        dirs = os.listdir(path)
        if "DON_NOT_ANALYZE" in [x.upper() for x in dirs]:
            return

        init_poscar = False
        for ifile in dirs:
            if 'init' in ifile:
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
            final_structure = pychemia.code.vasp.read_poscar(
                path + os.sep + "CONTCAR")
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
            properties['band_gap'] = vasprun.band_gap
            properties['xc'] = get_xc(vasprun.potcar_info, incar)
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
        return

    def to_excel(self,
                 filename="xc_db_analysis.xlsx",
                 query={'status.status.relaxed': True}):
        ret = OrderedDict({'Material ID': [],
                           'Exchange correlation':[],
                           'Space group name': [],
                           'Space group number': [],
                           'Family': [],
                           'No. lattice degrees of freedom': [],
                           'Bandgap': [], 'Direct gap': [],
                           'Atom symbol': [],
                           'No. Wyckoff degrees of freedom': [],
                           'No. species': [], 'No. atoms': [],
                           'a freedom': [], "b freedom": [], "c freedom": [],
                           'x freedom': [], "y freedom": [], "z freedom": [],
                           'Multiplicity': [],
                           'Wyckoff letter': [],
                           'a': [], 'b': [], 'c': [],
                           'x': [], 'y': [], 'z': [],
                           'alpha': [], 'beta': [], 'gamma': [],
                           'a ICSD': [], 'b ICSD': [], 'c ICSD': [],
                           'x ICSD': [], 'y ICSD': [], 'z ICSD': [],
                           'alpha ICSD': [], 'beta ICSD': [], 'gamma ICSD': [], 'Same spg': [],
                           })
        for entry in self.entries.find(query):
            _id = str(entry['_id'])
            natom = entry['properties']['final']['structure']['natom']
            nspecies = entry['properties']['final']['structure']['nspecies']
            reduced = entry['properties']['final']['structure']['reduced']
            reduced_initial = entry['properties']['initial']['structure']['reduced']
            symbols = entry['properties']['final']['structure']['symbols']
            inner = entry['properties']['final']['wyckoff_analysis']['inner_degrees_of_freedom']
            outer =  entry['properties']['final']['lattice_degrees_of_freedom']
            gap_spin_channel = min([x for x in entry['properties']['band_gap']])
            a_fin = entry['properties']['final']['lattice']['a']
            b_fin = entry['properties']['final']['lattice']['b']
            c_fin = entry['properties']['final']['lattice']['c']
            a_init = entry['properties']['initial']['lattice']['a']
            b_init = entry['properties']['initial']['lattice']['b']
            c_init = entry['properties']['initial']['lattice']['c']
            for iatom in range(natom):
                ret['Material ID'].append(entry['properties']['mp_id'])
                ret['Exchange correlation'].append(entry['properties']['xc'])
                ret['Space group number'].append(entry['properties']['final']['crystal_symmetry']['number'])
                ret['Space group name'].append(entry['properties']['final']['crystal_symmetry']['international_short'])
                ret['Family'].append(entry['properties']['final']['crystal_family'])
                if iatom == 0:
                    ret['No. lattice degrees of freedom'].append(sum(outer))
                else:
                    ret['No. lattice degrees of freedom'].append(0)
                ret['Bandgap'].append(entry['properties']['band_gap'][gap_spin_channel]['gap'])
                ret['Direct gap'].append(entry['properties']['band_gap'][gap_spin_channel]['direct'])
                ret['Atom symbol'].append(symbols[iatom])
                ret['No. Wyckoff degrees of freedom'].append(sum(inner[iatom]))
                ret['No. species'].append(nspecies)
                ret['No. atoms'].append(natom)
                
                ret['a freedom'].append(outer[0])
                ret["b freedom"].append(outer[1])
                ret["c freedom"].append(outer[2])
                
                ret['x freedom'].append(inner[iatom][0])
                ret["y freedom"].append(inner[iatom][1])
                ret["z freedom"].append(inner[iatom][2])
                ret['Multiplicity'].append(entry['properties']['final']['wyckoff_analysis']['multiplicity'][iatom])
                ret['Wyckoff letter'].append(entry['properties']['final']['wyckoff_analysis']['wyckoff_letter'][iatom])
                ret['a'].append(entry['properties']['final']['lattice']['a'])
                ret['b'].append(entry['properties']['final']['lattice']['b'])
                ret['c'].append(entry['properties']['final']['lattice']['c'])
                ret['x'].append(reduced[iatom][0])
                ret['y'].append(reduced[iatom][1])
                ret['z'].append(reduced[iatom][2])
                ret['alpha'].append(entry['properties']['final']['lattice']['alpha'])
                ret['beta'].append(entry['properties']['final']['lattice']['beta'])
                ret['gamma'].append(entry['properties']['final']['lattice']['gamma'])
                ret['a ICSD'].append(entry['properties']['initial']['lattice']['a'])
                ret['b ICSD'].append(entry['properties']['initial']['lattice']['b'])
                ret['c ICSD'].append(entry['properties']['initial']['lattice']['c'])
                ret['x ICSD'].append(reduced_initial[iatom][0])
                ret['y ICSD'].append(reduced_initial[iatom][1])
                ret['z ICSD'].append(reduced_initial[iatom][2])
                ret['alpha ICSD'].append(entry['properties']['initial']['lattice']['alpha'])
                ret['beta ICSD'].append(entry['properties']['initial']['lattice']['beta'])
                ret['gamma ICSD'].append(entry['properties']['initial']['lattice']['gamma'])
                ret['Same spg'].append(entry['properties']['initial']['crystal_symmetry']['number'] == entry['properties']['final']['crystal_symmetry']['number'])
        df = pd.DataFrame.from_dict(ret)
        
        df['x error'] = [diff_reduced_coordinate(x,x_prime) for x, x_prime in zip(df['x'], df['x ICSD'])]
        df['y error'] = [diff_reduced_coordinate(x,x_prime) for x, x_prime in zip(df['y'], df['y ICSD'])]
        df['z error'] = [diff_reduced_coordinate(x,x_prime) for x, x_prime in zip(df['z'], df['z ICSD'])]
        df['a error'] = (df['a'] - df['a ICSD'])/df['a']
        df['b error'] = (df['b'] - df['b ICSD'])/df['b']
        df['c error'] = (df['c'] - df['c ICSD'])/df['c']
        df['(RSME xyz).(DOF)'] = (((df['x error']*df['x freedom'])**2+
                                   (df['y error']*df['y freedom'])**2+
                                   (df['z error']*df['z freedom'])**2)/(df['No. Wyckoff degrees of freedom']))**(1/2)
        df['(MAE xyz).(DOF)'] = (df['x error']*df['x freedom']+
                                 df['y error']*df['y freedom']+
                                 df['z error']*df['z freedom'])/(df['No. Wyckoff degrees of freedom'])

        df['(RSME abc).(DOF)'] = (((df['a error']*df['a freedom'])**2+
                                   (df['b error']*df['b freedom'])**2+
                                   (df['c error']*df['c freedom'])**2)/(df['No. lattice degrees of freedom']))**(1/2)
        df['(MAE abc).(DOF)'] = (abs(df['a error'])*df['a freedom']+
                                 abs(df['b error'])*df['b freedom']+
                                 abs(df['c error'])*df['c freedom'])/(df['No. lattice degrees of freedom'])
        df['RSME abc'] = (((df['a error'])**2+
                           (df['b error'])**2+
                           (df['c error'])**2)/3)**(1/2)*df['No. lattice degrees of freedom']/df['No. lattice degrees of freedom']
        df['ME abc'] = (abs(df['a error'])+
                        abs(df['b error'])+
                        abs(df['c error']))/3*df['No. lattice degrees of freedom']/df['No. lattice degrees of freedom']
        df.to_excel(f'{self.name}-db.xlsx')
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



def diff_reduced_coordinate(x, x_prime):
    diff = abs(x - x_prime)
    if diff<=0.5:
        ret = diff
    elif diff>0.5:
        if x>=0.5:
            ret = abs(x-1-x_prime)
        elif x<0.5:
            ret = abs(x+1-x_prime)
    return ret
    
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
        help="Length threshold for symmtery analysis",
        default=1e-5)
    parser.add_argument(
        "--mode",
        dest='mode',
        nargs='+',
        help="what is expected from the script select from:\n update, export",
        default=['export'])

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
    for mode in args.mode:
        if mode in ['update', 'initiate']:
            db.check_dirs(args.path)
        elif mode == 'export_json':
            db.to_json()
        elif mode == 'export_excel':
            db.to_excel()
