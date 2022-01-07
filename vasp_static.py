#!/usr/bin/env python3

import pychemia
import argparse
import json
import os
KPoints = pychemia.crystal.KPoints

parser = argparse.ArgumentParser()
parser.add_argument("--structure",dest="structure",type=str,help='poscar structure that you want to run',default = 'POSCAR')
parser.add_argument("-np" ,dest="np",type=int ,action="store", help="Number of MPI processes for the code",default = '4')
parser.add_argument("-xc",dest="xc",type=str,help='exchange correlation functional',default = 'PBE')
parser.add_argument('--incar_encut',dest='incar_encut',help="path to the extra variable you want to include in encut convergence",default ='none')
parser.add_argument('--incar_kpoint',dest='incar_kpoint',help="path to the extra variable you want to include in kpoint convergence",default='none')
parser.add_argument('--incar_static',dest='incar_static',help="path to the extra variable you want to include in statication",default='none')
parser.add_argument('--potcar_options',dest='psp_options',nargs='+',help="The options want to be used for pseudo potentials, eg: --potcar_options O_h Cu_pv")
parser.add_argument('--tags',dest='tags',nargs='+',help='If you want to have a tag for all calculations for example ENCUT 250', default=None)

args = parser.parse_args()

extra_vars = {}
extra_vars_encut = {}
extra_vars_kpoint = {}
extra_vars_static = {"EDIFF":1e-8}

if args.tags is not None:
    for itag in range(0,len(args.tags),2):
        key = args.tags[itag]
        value = args.tags[itag+1]
        extra_vars[key]=value
        
if args.xc == 'LDA' :
    pspdir= 'potpaw_LDA'
elif args.xc == 'PBE':
    pspdir = 'potpaw_PBE'
elif args.xc == 'PBEsol':
    pspdir = 'potpaw_PBE'
    extra_vars['GGA']="PS"
elif args.xc == 'AM05':
    pspdir = 'potpaw_PBE'
    extra_vars['GGA']="AM"
elif args.xc == 'SCAN':
    pspdir = 'potpaw_PBE'
    extra_vars['METAGGA']="SCAN"
    extra_vars['ALGO']='A'
    extra_vars['ISMEAR']=0
    extra_vars['LASPH'] =True
elif args.xc == 'MBJ':
    pspdir = 'potpaw_PBE'
    extra_vars['METAGGA']='MBJ'
    extra_vars['ALGO']='A'
    extra_vars['ISMEAR']=0
    extra_vars['LASPH'] =True
elif args.xc == "R2SCAN":
    pspdir = 'potpaw_PBE'
    extra_vars['METAGGA']='R2SCAN'
    extra_vars['ALGO']='A'
    extra_vars['ISMEAR']=0
    extra_vars['LASPH'] =True
elif args.xc == "RSCAN":
    pspdir = 'potpaw_PBE'
    extra_vars['METAGGA']='RSCAN'
    extra_vars['ALGO']='A'
    extra_vars['ISMEAR']=0
    extra_vars['LASPH'] =True

extra_vars['NBANDS'] = ' '

psp_options = {}
if args.psp_options:
    for option in args.psp_options:
        psp_options[option.split('_')[0]] = option.split('_')[1]


for ikey in extra_vars:
    extra_vars_encut[ikey] = extra_vars[ikey]
    extra_vars_kpoint[ikey] = extra_vars[ikey]
    extra_vars_static[ikey] = extra_vars[ikey]
    
if args.incar_encut != 'none' :
    incar_encut = pychemia.code.vasp.read_incar(str(args.incar_encut))
    extra_vars_encut = {}
    for tag in incar_encut:
        extra_vars_encut[tag] = incar_encut[tag]

if args.incar_kpoint != 'none' :
    incar_kpoint = pychemia.code.vasp.read_incar(str(args.incar_kpoint))
    extra_vars_kpoint = {}
    for tag in incar_kpoint:
        extra_vars_kpoint[tag] = incar_kpoint[tag]

if args.incar_static != 'none' :
    incar_static = pychemia.code.vasp.read_incar(str(args.incar_static))
    extra_vars_static = {}
    for tag in incar_static:
        extra_vars_static[tag] = incar_static[tag]



if 'ascii' in args.structure:
    st = pychemia.io.ascii.load(args.structure)
else : 
    st = pychemia.code.vasp.read_poscar(args.structure)



if 'ENCUT' in extra_vars:
    encut = int(extra_vars['ENCUT'])
    del extra_vars['ENCUT']
    del extra_vars_kpoint['ENCUT']
    del extra_vars_static['ENCUT']
                
elif not os.path.exists('encut_report.json'):
    encut_conv = pychemia.code.vasp.task.ConvergenceCutOffEnergy(structure=st,
                                                                 workdir='.',
                                                                 executable='vasp_std',
                                                                 pspdir=pspdir,
                                                                 extra_vars=extra_vars_encut,
                                                                 energy_tolerance=1E-3,
                                                                 psp_options=psp_options)
    encut_conv.run(nparal=args.np)
    encut_conv.save('encut_report.json')
    if not encut_conv.success:
        print("ENCUT convergence unsuccessful")
        exit()
    encut = encut_conv.best_encut
else :
    rf = open('encut_report.json','r')
    encut_data =json.load(rf)
    rf.close()
    encut = encut_data['output']['best_encut']

    
if not os.path.exists('kpoint_report.json'):
    kpt_conv = pychemia.code.vasp.task.ConvergenceKPointGrid(structure=st,
                                                             workdir='.',
                                                             executable='vasp_std',
                                                             pspdir=pspdir,
                                                             extra_vars=extra_vars_kpoint,
                                                             energy_tolerance=1E-3,
                                                             psp_options=psp_options)
    kpt_conv.run(args.np)
    kpt_conv.save('kpoint_report.json')
    if not kpt_conv.success:
        print('Kpoint convergence unsuccessful')
        exit()
    kp = kpt_conv.best_kpoints
    kgrid = kpt_conv.best_kpoints.grid

else :
    rf = open('kpoint_report.json','r')
    kpoint_data = json.load(rf)
    rf.close()
    kgrid = kpoint_data['output']['best_kp_grid']

static = pychemia.code.vasp.task.StaticCalculation(structure=st,
                                                   workdir='.',
                                                   executable='vasp_std',
                                                   encut=encut,
                                                   kpoints=KPoints(kmode='gamma', grid=kgrid),
                                                   extra_incar=extra_vars_static,
                                                   pspdir=pspdir,
                                                   psp_options=psp_options)
                                                   
static.run(args.np)
static.save('static_report.json')

