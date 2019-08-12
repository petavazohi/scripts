#!/usr/bin/env python3

import sys
import os
from ase import *
from ase import io
from ase import Atoms, Atom
from ase.calculators.vasp import Vasp
import numpy as np
import spglib
from converge import *
import re


vasp_params = {
    "xc": "PBE",
    "isym": 2, # run with (2) or without (0) symmetries
    "symprec": 1e-8,
#    "symprec": 1e-6,   ## SKS changed it later... may be better to switch it back to 1e-08
    "gamma": False,
    "algo": "Fast",
    'encut' : 800,
    #"algo": "Damped",
    #"time": 0.5,
    #"lhfcalc": True,
    #"hfscreen": 0.2,
    #"encutfock": 0,
    #"lsorbit" : True,

}

basename, extension = os.path.splitext(sys.argv[1])
vasp_params["pstress"] = 10.0*float(sys.argv[2])
vasp_params["ncore"] = 8
wf = 'convergence'

if os.path.exists(wf):
    rf = open(wf,'r')
    data = rf.read()
    rf.close()
    if 'Stress' in data :
        print('Calculation has finished')
        exit()
    stages  = [int(x) for x in re.findall('stage : ([0-9])',data)]
    if len(stages) == 0 :
        stage = 1
    else : 
        stage = stages[-1]


    if stage == 2 :
        kpt_1, kpt_2, kpt_3 =  [int(x) for x in re.findall('Converged\s*:\s*([0-9]*)x([0-9]*)x([0-9]*)',data)[0]]
        print('stage : %d' % stage)
        print('kpoint converged found : %dx%dx%d' %(kpt_1, kpt_2, kpt_3))
        if os.path.exists('CONTCAR'):
            print('found CONTCAR file, resuming the optimization')
            primitive_cell = read_file("CONTCAR")
        else :
            print('did not fine CONTCAR starting from %s'% sys.argv[1])
            primitive_cell = read_file(sys.argv[1])
    elif stage == 3 :
        print('stage : %d' % stage)
        print('reading CONTCAR')
        primitive_cell = read_file("CONTCAR")
    elif stage == 4 or stage ==5:
        kpt_1, kpt_2, kpt_3 =  [int(x) for x in re.findall('Converged\s*:\s*([0-9]*)x([0-9]*)x([0-9]*)',data)[1]]
        print('stage : %d' % stage)
        print('reading CONTCAR')
        primitive_cell = read_file("CONTCAR")
    else : 
        print('stage : 1')
        print('No %s file was found reading %s '% (wf,sys.argv[1]))
        stage = 1
        primitive_cell = read_file(sys.argv[1])
        
else : 
    print('stage : 1')
    print('No %s file was found reading %s '% (wf,sys.argv[1]))
    stage = 1
    primitive_cell = read_file(sys.argv[1])

print(primitive_cell)

vasp_params["ldau"]     = True
vasp_params["ldautype"] = 2
ldau_luj = {}
ldau_luj['Mn']= {'L':2,'U':3.9,'J':0}
ldau_luj['Bi']= {'L':-1,'U':0,'J':0}
vasp_params['ldau_luj']= ldau_luj
vasp_params["ldauprint"] = 2
vasp_params["lmaxmix"]   = 4 

nMn = sum(primitive_cell.get_atomic_numbers() == 25)
nBi = sum(primitive_cell.get_atomic_numbers() == 83)
magmom = []
for i in range(nMn) :
    magmom.append(0.5)
for i in range(nBi) :
    magmom.append(0.0)
   
primitive_cell.set_initial_magnetic_moments(magmom)


if stage == 1:
    if os.path.exists('convergence'):
        os.remove('convergence')
    append_file(wf,'stage : 1\n')
    # find converged k-point sampling
    kpt_1, kpt_2, kpt_3 = converge_kpoints(primitive_cell, wf, tol=0.002, vasp_params=vasp_params)
    # now we run the geometry optimization with vasp
    stage = 2

if stage == 2:

    append_file(wf,'stage : 2\n')
    append_file(wf,"\n==== Starting geometry:")
    append_file(wf, "Space group: "+ str(spglib.get_spacegroup(primitive_cell,symprec=1e-3))+'\n')
    append_file(wf, str(primitive_cell.get_cell())+ "\n"+ str(primitive_cell.get_scaled_positions())+'\n')
    calc = Vasp(kpts=(kpt_1, kpt_2, kpt_3), prec='high', ibrion=2, isif=3, ediff=1E-07, ediffg=-0.005, nsw=250, **vasp_params)
    calc.calculate(primitive_cell)
    clean_vasp()
    stage = 3 

if stage == 3:

    append_file(wf,'stage : 3\n')
    kpt_1, kpt_2, kpt_3 = converge_kpoints(primitive_cell, wf, tol=0.002, vasp_params=vasp_params)
    stage = 4

if stage == 4:

    append_file(wf,'stage : 4\n')
    calc = Vasp(kpts=(kpt_1, kpt_2, kpt_3), prec='high', ibrion=2, isif=3, ediff=1E-07, ediffg=-0.005, nsw=250, **vasp_params)
    calc.calculate(primitive_cell) # we run this twice
    primitive_cell.write(basename + '_opt.cif', format='cif')
    stage = 5

if stage == 5:
    # we now rerun vasp to get the ground-state
    # we use always the same cutoffs in order to ensure compatibility between the runs
    append_file(wf,'stage : 5\n')
    calc = Vasp(kpts=(kpt_1, kpt_2, kpt_3),  prec='high', ediff=1E-08, emin=-15, emax=15, **vasp_params)
    primitive_cell.set_calculator(calc)
    append_file(wf, "\n==== End geometry:")
    append_file(wf, "Space group: "+ str(spglib.get_spacegroup(primitive_cell,symprec=1e-3)))
    append_file(wf, str(primitive_cell.get_cell())+"\n"+ str(primitive_cell.get_scaled_positions()))
    append_file(wf, "\nForces:")
    append_file(wf, str(primitive_cell.get_forces()))
    append_file(wf,"\nStress:")
    append_file(wf, str(primitive_cell.get_stress()))
    clean_vasp()
