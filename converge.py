#!/usr/bin/env python3                                                                                                                                                                                              
import sys
import os
from ase import *

def append_file(file_name,message):
    if not os.path.exists(file_name):
        wf = open(file_name,'w')
    else :
        wf = open(file_name,'a')
    wf.write(message)
    wf.close()
    return

def kpoints_max(structure, kppa=1000, return_kgrid=False):
    """Returns the max number of kpoints allowed during the
    kpoints convergence test. Based on the routine used in
    pymatgen.
    
    "struture" must be an ase.atoms object.

    """
    from math import sqrt, ceil

    lengths = [sqrt(sum(map(lambda y: y**2,structure.cell[i]))) \
                   for i in range(3)]
    ngrid = kppa / structure.get_number_of_atoms()
    mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1.0/3.0)
    num_div = [int(ceil(1.0 / lengths[i] * mult)) for i in range(3)]
    num_div = [i if i > 0 else 1 for i in num_div]
    if return_kgrid:
        return num_div[0],num_div[1],num_div[2]
    else:
        return num_div[0]*num_div[1]*num_div[2]
                            

def converge_kpoints(atoms, fh, tol=1e-3*27.211, vasp_params={}, kppa=1000):
    """converge the k-point sampling"""

    import os.path
    from math import sqrt, ceil
    from ase.calculators.vasp import Vasp

    if(os.path.exists('in.kpt')):    
        f = open('in.kpt', 'r')
        kptr = f.readline().split()
        f.close()
        return int(kptr[0]), int(kptr[1]), int(kptr[2])

    rcell = atoms.get_reciprocal_cell()
    vol   = 1.0/atoms.get_volume()
    nat   = atoms.get_number_of_atoms()
    
    # lets get the cell ratios
    a1 = sqrt(rcell[0][0]**2 + rcell[0][1]**2 + rcell[0][2]**2)
    a2 = sqrt(rcell[1][0]**2 + rcell[1][1]**2 + rcell[1][2]**2)
    a3 = sqrt(rcell[2][0]**2 + rcell[2][1]**2 + rcell[2][2]**2)

    if fh != None:
        append_file(fh,"==== Convergence with respect to k-points:\n")

    kpt = 1
    E1, E2 = 0.0, 0.0
    kpt_1  = [0, 0, 0]
    kpt_2  = [0, 0, 0]

    kpt_max = kpoints_max(atoms, kppa=kppa)
    
    while(True):
        if kpt == 1:
            kpt_new = [1, 1, 1]
        else:
            factor  = pow(kpt/vol, 1.0/3.0)
            kpt_new = [int(max(ceil(factor*a1), 1)), int(max(ceil(factor*a2), 1)), int(max(ceil(factor*a3), 1))]

        if kpt_new != kpt_1:
            calc = Vasp(kpts=(kpt_new[0], kpt_new[1], kpt_new[2]), ediff=1E-05, **vasp_params)
            print(kpt_new,vasp_params)
            atoms.set_calculator(calc)
            try:  # sometimes vasp crashes for some k-point sets
                print(atoms.get_potential_energy())
                E_new = atoms.get_potential_energy() / nat

                os.remove('WAVECAR')

                if fh != None:
                    append_file(fh,'target = %d \t test = %dx%dx%d = %d \tEnergy = %f \n' % (kpt,kpt_new[0],kpt_new[1],kpt_new[2],kpt_new[0]*kpt_new[1]*kpt_new[2],E_new))


                converged = ((E1 != 0.0) and (abs(E_new - E1) < tol) and (abs(E_new - E2) < tol))
                if converged or 2*kpt > kpt_max:
                    if not converged:
                        if abs(E_new - E1) < tol:
                            kpt_2 = kpt_1
                        else:
                            kpt_2 = kpt_new
                    if fh != None:
                        append_file(fh,"\nConverged : %dx%dx%d \n" % (kpt_2[0],kpt_2[1],kpt_2[2]))
                    break
            except:
                pass


            kpt_2 = kpt_1[:]
            kpt_1 = kpt_new[:]

            E2 = E1
            E1 = E_new
    
        kpt = kpt*2
        if(kpt > 10000):
            append_file("\nToo many k-points. Aborting")
            sys.exit("Too many k-points. Aborting")

    calc.clean()

    # make sure that we have at least the MP k-points
    kpt_min = kpoints_max(atoms, kppa=1000, return_kgrid=True)
    if(kpt_min[0] > kpt_2[0] or kpt_min[1] > kpt_2[1] or kpt_min[2] > kpt_2[2]):
        return kpt_min[0], kpt_min[1], kpt_min[2]
    else:
        return kpt_2[0], kpt_2[1], kpt_2[2]


def read_file(in_file):
    """reads a file and calculates the primitive unit cell"""

    import os
    import subprocess
    from copy import copy
    from ase import io
    import numpy as np

    basename, extension = os.path.splitext(in_file)
    basename = os.path.basename(basename)
    if (basename == 'POSCAR') or (basename == 'CONTCAR'):
        form = 'vasp'
    elif extension == '.ascii':
        form = 'v-sim'
    else:
        form = extension[1:]

    # generate the primitive unit cell
    if form == 'cif':
        cifReplace = {
            "data_ ": "data_",
            "_space_group_symop_operation_xyz": "_symmetry_equiv_pos_as_xyz",
            " 0.33333 ": " 0.333333333333333 ",
            " 0.66667 ": " 0.666666666666667 ",
            # now for the H-M symbols
            "P 2/m 2/m 2/m": "Pmmm", # 47
            "P 21/m 2/m 2/a": "Pmma", # 51
            "P 2/n 21/n 2/a": "Pnna", # 52
            "P 2/m 2/n 21/a": "Pmna", # 53
            "P 21/c 2/c 2/a": "Pcca", # 54
            "P 21/b 21/a 2/m": "Pbam", # 55
            "P 2/b 21/c 21/m": "Pbcm", # 57
            "P 21/n 21/n 2/m": "Pnnm", # 58
            "P 21/m 2/n 21/m (origin choice 1)": "Pmnm", # 59
            "P 21/b 2/c 21/n": "Pbcn", # 60
            "P 21/n 21/m 21/a": "Pnma", # 62
            "C 2/m 2/c 21/m": "Cmcm", # 63
            "C 2/m 2/c 21/a": "Cmca", # 64
            "C 2/m 2/m 2/m": "Cmmm", # 65
            "C 2/m 2/m 2/a":  "Cmma", # 67
            "F 2/m 2/m 2/m": "Fmmm", # 69
            "F 2/d 2/d 2/d (origin choice 1)": "Fddd", # 70
            "I 2/m 2/m 2/m": "Immm", # 71
            "I 21/m 21/m 21/a": "Imma", # 74
            "P 4/m 2/m 2/m": "P4/mmm", # 123
            "P 4/n 21/m 2/m (origin choice 1)": "P4/nmm", # 129
            "P 42/m 2/m 2/c": "P42/mmc", # 131
            "P 42/m 21/n 2/m": "P42/mnm", # 136
            "I 4/m 2/m 2/m": "I4/mmm", # 139
            "R 3 (hexagonal axes)": "R3", # 146
            "R 3 m (hexagonal axes)": "R3m", # 160
            "P -3 2/m 1": "P-3m1", # 164
            "R -3 2/m (hexagonal axes)": "R-3m", # 166
            "P 6/m 2/m 2/m": "P6/mmm", # 191
            "P 63/m 2/m 2/c":  "P63/mmc", # 194
            "P 42/m -3 2/n": "Pm-3n", # 223
        }    

        # we do not trust the reading of cif files of ase so we use cif2cell
        # to convert them to abinit format first
        pid  = str(os.getpid())

        # do a bit of search and replace so that cif2cell can read cif
        infile = open(in_file)
        outfile = open("/tmp/" + pid + ".cif", 'w')

        for line in infile:
            for src, target in cifReplace.iteritems():
                line = line.replace(src, target)
            outfile.write(line)

        infile.close()
        outfile.close()

        subprocess.call("cif2cell /tmp/" + pid + ".cif -p abinit -o /tmp/" + pid + ".abinit", shell=True)
        subprocess.call("python /scratch/petavazohi/MHM/refine/correct.acell.py /tmp/"+ pid + ".abinit > /tmp/"+ pid + ".corr", shell=True)
        subprocess.call("mv /tmp/"+ pid + ".corr /tmp/"+ pid + ".abinit", shell=True)
        cell = io.read("/tmp/" + pid + ".abinit", format="abinit")
        os.remove("/tmp/" + pid + ".cif")
        os.remove("/tmp/" + pid + ".abinit")
    else:
        cell = io.read(in_file, format=form)

    # check if volume is negative
    ucell = cell.get_cell()
    if np.linalg.det(ucell) < 0.0:
        tmp = copy(ucell[1, :])
        ucell[1, :] = copy(ucell[2, :])
        ucell[2, :] = copy(tmp)
        cell.set_cell(ucell)
        
        cpos = cell.get_positions()
        tmp = copy(cpos[:, 1])
        cpos[:, 1] = copy(cpos[:, 2])
        cpos[:, 2] = copy(tmp)
        cell.set_positions(cpos)

    return cell


def clean_vasp():
    """Delete unecessary vasp files"""

    import os
    
    files = ['CHG', 'CHGCAR', 'IBZKPT', 'OSZICAR',
         'PCDAT', 'WAVECAR', 'XDATCAR', 'PROCAR', 'ase-sort.dat',
         'LOCPOT', 'AECCAR0', 'AECCAR1', 'AECCAR2']
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass
