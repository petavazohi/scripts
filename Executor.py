#!/usr/bin/env python3

import os
import subprocess
import argparse
import time
import re
import numpy as np 
import math
import pychemia
import pymatgen.io.vasp as vasp

def create_kpoints(length):
    length = str(length)+'\n'
    comment = 'Automatic mesh'+'\n'
    nkp     = str(0)+'\n'
    method  = 'Auto'+'\n'
    wf = open("KPOINTS",'w')
    wf.write(comment)
    wf.write(nkp    )
    wf.write(method )
    wf.write(length )
    wf.close()
    return

def kpoint_manual(e_threshold,start,end,step,executable,nparal):
    # manual kmesh generation
    # in this case start = total number of kpoints to start with
    #              end   = maximum number of kpoints
    #              step  = jumps from
    address = os.getcwd()
    # I am adding this section as selective dynamics is still not implemented in pychemia
    try :
        structure = pychemia.code.vasp.read_poscar("POSCAR")
        cell = structure.cell
    except :
        structure = vasp.Poscar.from_file("POSCAR").structure
        cell = structure.lattice.matrix
    rec_cell = np.linalg.inv(cell).T * 2*np.pi
    b1 = np.linalg.norm(rec_cell[0,:])
    b2 = np.linalg.norm(rec_cell[1,:])
    b3 = np.linalg.norm(rec_cell[2,:])
    # rlv = Reciprocal Lattice Vector
    rlv = np.array([b1,b2,b3])
    kpnts = [] 
    toten= [] 
    for i in range(start,end,step): 
        j = np.arange(0,200,1e-3) 
        counter = 0 
        while math.ceil(rlv[np.argmax(rlv)]*j[counter]) != i : 
            counter += 1
        kpnts.append([math.ceil(b1*j[counter]),math.ceil(b2*j[counter]),math.ceil(b3*j[counter])]) 
    for ikpoint in kpnts:
        print("===================================================================================")
        print("===================================================================================")
        print("===================================================================================")
        print('Running VASP for kmesh %i %i %i'% (tuple(ikpoint)))
        kp = pychemia.crystal.KPoints(kmode='Monkhorst-pack')
        kp.set_grid(ikpoint)
        pychemia.code.vasp.kpoints.write_kpoints(kp=kp,filepath='KPOINTS')
        if not os.path.exists('INCAR'):
            incar = pychemia.code.vasp.VaspInput()
            incar.set_encut(1.4,POTCAR='POTCAR')
            incar['EDIFF' ] = 1e-4
            incar['NWRITE'] = 2
            incar['PREC'  ] = 'Accurate'
            incar['NCORE' ] = 4
            incar['SYSTEM'] = '-'.join(address.split('/')[-2:])
            incar.write("INCAR")
        runtime = execute(nparal,address,executable)
        rf = open("OUTCAR",'r')
        data = rf.read()
        rf.close()
        toten.append(float(re.findall("TOTEN\s*=\s*([-+0-9.]*)\s*eV",data)[-1]))
        wf = open("kpoint_convergence",'a')
        wf.write("kmesh = %i %i %i , TOTEN =%f \n" %(ikpoint[0],ikpoint[1],ikpoint[2],toten[-1]))
        wf.close()
    kpnts = np.array(kpnts)
    toten = np.array(toten)
    kpnt_idx = abs(toten - toten[-1]) < e_threshold
    best_kpnt = list(kpnts[kpnt_idx][0])
    wf = open('best_kpnt','w')
    wf.write('kpnt_cut = %i %i %i ' % tuple(best_kpnt))
    wf.close()

    return
    


def kpoint_convergence(e_threshold,start,end,step,executable,nparal):
    address = os.getcwd() 
    toten = []
    kmesh = []
    if not os.path.exists('INCAR'):
        incar = pychemia.code.vasp.VaspInput()
        incar.set_encut(1.4,POTCAR='POTCAR')
        incar['EDIFF' ] = 1e-4
        incar['NWRITE'] = 2
        incar['PREC'  ] = 'Accurate'
        incar['NCORE' ] = 4
        incar['SYSTEM'] = '-'.join(address.split('/')[-2:])
        incar.write("INCAR")
    # create potcar here 
    klengths = np.arange(start,end,step)
    for klength in klengths:
        create_kpoints(klength)
        print("===================================================================================")
        print("===================================================================================")
        print("===================================================================================")
        runtime = execute(nparal,address,executable)
        rf = open("OUTCAR",'r')
        data = rf.read()
        rf.close()
        toten.append(float(re.findall("TOTEN\s*=\s*([-+0-9.]*)\s*eV",data)[-1]))
        kmesh.append(re.findall("generate k-points for.*",data)[0])
        wf = open("kpoint_convergence",'a')        
        wf.write("kpoint length = %i ,kmesh = %s, TOTEN =%f \n" %(klength,kmesh[-1],toten[-1]))
        wf.close()
        if klength != start : # this is to check if it's not doing the calculation for the first time
            change = abs(toten[-2]-toten[-1])
            if change < e_threshold : # need to add more comparision here
                print("VASP calculations converged with k points length %i and kmesh %s " % (klength,kmesh[-1]))
                break
    print("VASP calculations converged with k points length %i and kmesh %s " % (klengths[conv_idx],kmesh[conv_idx]))
    return 

def encut_convergence(e_threshold,start,end,step,executable,nparal):
    address = os.getcwd()
    toten = []
    if not os.path.exists('KPOINTS'):
        create_kpoints(10)
    rf = open('POTCAR')
    potcar = rf.read()
    rf.close()
    encut_init = round(max([float(x) for x in re.findall('ENMAX\s*=\s*([0-9.]*);',potcar)])*1.3)
    if start < encut_init :
        print('Initial value provided for ENUCT is less than 1.3*ENMAX pseudo potential, replacing Estart with %f'% encut_init)
        start = round(encut_init,-2)
    encuts = np.arange(start,end,step)
    for iencut in encuts:
        incar = pychemia.code.vasp.VaspInput()
        incar['ENCUT']  = iencut
        incar['SYSTEM'] = '-'.join(address.split('/')[-2:])
        incar['EDIFF' ] = 1e-5
        incar['NWRITE'] = 2
        incar['PREC'  ] = 'Accurate'
        incar['NCORE' ] = 4

        incar.write("INCAR")
        print("===================================================================================")
        print("===================================================================================")
        print("===================================================================================")
        print('Running vasp with ENCUT = {}'.format(iencut))
        runtime = execute(nparal,address,executable)
        rf = open("OUTCAR",'r')
        data = rf.read()
        rf.close()
        toten.append(float(re.findall("TOTEN\s*=\s*([-+0-9.]*)\s*eV",data)[-1]))
        wf = open("encut_convergence",'a')
        wf.write("encut = %i , TOTEN =%f \n" %(iencut,toten[-1]))
        wf.close()
    toten = np.array(toten)
    encut_idx = abs(toten - toten[-1]) < e_threshold
    best_encut = encuts[encut_idx][0]
    wf = open('best_encut','w')
    wf.write('best_cut = {}'.format(best_encut))
    wf.close()
    return

def relax_structure(encut,kgrid,kmode,ismear,executable,nparal):
    address = os.getcwd()
    if not os.path.exists('INCAR'):
        incar = pychemia.code.vasp.VaspInput()
        if encut == None:
            incar.set_encut(1.4,POTCAR='POTCAR')
        else :
            incar.set_encut(encut)
            incar['SYSTEM']  = '-'.join(address.split('/')[-2:])
            incar['NWRITE']  = 3        # how much info to write in outcar, long MD 0,1 short MD 2,3 debug 4
            incar['PREC']    = 'Accurate'
            incar['ADDGRID'] = True    # additional support grid for augmentation charges
            incar['ISMEAR']  = ismear   # tetrahedron method with Blöchl corrections 
            incar['ISTART']  = 0        # does not read WAVECAR
            incar['LWAVE']   = False    # does not write WAVECAR
            incar['EDIFF']   = 1e-06
            incar['LREAL']   = 'a'      # projection in real space or reciprocal, a is automatic
            incar['NELMIN'] = 6         # minimum number of electronic steps
            incar['EDIFFG'] = -1E-03
            incar['NSW'] = 100
            incar['IBRION'] = 2         # how atoms are updated to move 
            incar['ISIF'] = 3
            incar['ISYM'] = 2
            incar['SIGMA'] = 0.05
            incar['NPAR' ] = int(nparal/16)
            incar.write("INCAR")
    if kmode == None and kgrid == None:
        create_kpoints(30)
    else : 
        kp = pychemia.crystal.KPoints(kmode=kmode)
        kp.set_grid(kgrid)
        pychemia.code.vasp.kpoints.write_kpoints(kp=kp,filepath='KPOINTS')
    runtime = execute(nparal,address,executable)
    return

def SCF(encut,kgrid,kmode,ismear,executable,nparal):
    address = os.getcwd()
    magmom = '' 
    try :
        structure = pychemia.code.vasp.read_poscar("POSCAR")
        comp = structure.get_composition() 
    except :
        structure = vasp.Poscar.from_file("POSCAR").structure
        comp = structure.composition
    for x in comp: 
        magmom += '%i*0.5 ' % comp[x] 
    incar = pychemia.code.vasp.VaspInput()
    if encut == None:
        incar.set_encut(1.4,POTCAR='POTCAR')
    else :
        incar.set_encut(encut)
    incar['SYSTEM']  = '-'.join(address.split('/')[-2:])
    incar['ISMEAR']  = ismear   # tetrahedron method with Blöchl corrections 
    incar['ISTART']  = 0        # does not read WAVECAR
    incar['EDIFF']   = 1e-06
    incar['LREAL']   = 'a'      # projection in real space or reciprocal, a is automatic
    incar['NELMIN'] = 6         # minimum number of electronic steps
    incar['IBRION'] = -1         # how atoms are updated to move 
    incar['ISPIN'] = 2
    incar['NCORE' ] = 4
    incar['MAGMOM'] = magmom
    incar['SIGMA'] = 0.05
    incar.write("INCAR")
    if kmode == None and kgrid == None:
        create_kpoints(30)
    else : 
        kp = pychemia.crystal.KPoints(kmode=kmode)
        kp.set_grid(kgrid)
        pychemia.code.vasp.kpoints.write_kpoints(kp=kp,filepath='KPOINTS')
    runtime = execute(nparal,address,executable)
    return


def execute(nparal,work_dir,vasp_exe):

     wf=open(work_dir+os.sep+'RUNNING','w')
     wf.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))
     wf.close()
     start_time=time.time()
     status = subprocess.call("mpirun -np {} {}".format(nparal,vasp_exe),cwd=work_dir ,shell=True)
     status = 0
     end_time=time.time()
     runtime=end_time-start_time
     if status== 0:
         print("VASP execution completed with returcode: %d runtime: %d secs" % (status, runtime))
     else:
         print("VASP execution failed with returcode: %d runtime: %d secs" % (status, runtime))
     os.remove(work_dir+os.sep+'RUNNING')
     return runtime

if __name__ == "__main__" : 
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure",dest="structure",type=str,help='poscar structure that you want to run',default = 'POSCAR')
    parser.add_argument("-np" ,dest="np",type=int ,action="store", help="Number of MPI processes for the code",default = '1')
    subparsers = parser.add_subparsers(dest='calc')
    parser_kpnt = subparsers.add_parser('kpoint_convergence')
    parser_kpnt.add_argument('--mode',type=str,default='auto')
    parser_kpnt.add_argument('--Kstart',type=int,default=10)
    parser_kpnt.add_argument('--Kend',type=int,default=100)
    parser_kpnt.add_argument('--Kstep',type=int,default=10)
    parser_kpnt.add_argument('--Ethreshold',type=float,default=1e-3)
    parser_encut = subparsers.add_parser('encut_convergence')
    parser_encut.add_argument('--Estart',type=float,default=100)
    parser_encut.add_argument('--Eend',type=float,default=900)
    parser_encut.add_argument('--Estep',type=float,default=50)
    parser_encut.add_argument('--Ethreshold',type=float,default=1e-3)
    parser.add_argument("--executable",dest="executable",type=str,action="store",help="vasp executable",default="vasp_std")
    parser_rlx = subparsers.add_parser('structure_relax')
    parser_rlx.add_argument('--Kmode',type=str,default='Monkhorst-pack')
    parser_rlx.add_argument('--Kgrid',type=int,nargs=3)
    parser_rlx.add_argument('--encut',type=float)
    parser_rlx.add_argument('--ismear',type=int,default=0)
    parser_scf = subparsers.add_parser('scf')
    parser_scf.add_argument('--Kmode',type=str,default='Monkhorst-pack')
    parser_scf.add_argument('--Kgrid',type=int,nargs=3)
    parser_scf.add_argument('--encut',type=float)
    parser_scf.add_argument('--ismear',type=int,default=0)
    args = parser.parse_args()
    if  args.calc == 'kpoint_convergence':
        if args.mode == 'auto':
            kpoint_convergence(e_threshold=args.Ethreshold,start=args.Kstart,end=args.Kend,step=args.Kstep,executable=args.executable,nparal=args.np)
        elif args.mode == 'manual':
            kpoint_manual(e_threshold=args.Ethreshold,start=args.Kstart,end=args.Kend,step=args.Kstep,executable=args.executable,nparal=args.np)
    elif args.calc == 'encut_convergence':
        encut_convergence(e_threshold=args.Ethreshold,start=args.Estart,end=args.Eend,step=args.Estep,executable=args.executable,nparal=args.np)
    elif args.calc == 'structure_relax' :
        relax_structure(encut=args.encut,kgrid=args.Kgrid,kmode=args.Kmode,ismear=args.ismear,executable=args.executable,nparal=args.np)
    elif args.calc == 'scf':
        SCF(encut=args.encut,kgrid=args.Kgrid,kmode=args.Kmode,ismear=args.ismear,executable=args.executable,nparal=args.np)
