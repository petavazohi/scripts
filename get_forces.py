#!/usr/bin/env python3

import re
import numpy as np

rf = open("OUTCAR",'r')
data = rf.read()
rf.close()
pos_force = []
for iforces in re.findall('POSITION.*TOTAL-FORCE.*\n[-\s]*\n([\s0-9.\n-]*)\n\s----[-\s]*total',data):
    lines = iforces.split('\n')
    istep = []
    for iline in lines:
        elements = [float(x) for x in iline.split()]
        istep.append(elements)
    pos_force.append(istep)

lattice_params = []
for ilattice in re.findall('length of vectors.*\n\s*([0-9.\s]*)',data):
    lattice = [float(x) for x in ilattice.split()]
    lattice = lattice[0:3]
    lattice_params.append(lattice)
lattice_params = np.array(lattice_params)

press = []
for ipress in re.findall('external pressure\s=*\s*([-0-9.]*)',data):
    press.append(float(ipress))
pos_force = np.array(pos_force)
pos = pos_force[:,:,0:3]
force = pos_force[:,:,3:6]
press = np.array(press)
wf = open('FORCE_POS_by_atm','w')

for iatom in range(force.shape[1]):    
    wf.write("atom  %i\n"% (iatom+1))
    positions = pos[:,iatom,:]
    forces = force[:,iatom,:]
    for istep in range(len(positions[:,0])):
        wf.write("{: > 3} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12}\n".format(
            istep,positions[istep,0],positions[istep,1],positions[istep,2],forces[istep,0],forces[istep,1],forces[istep,2]))
wf.close()

wf = open('FORCE_POS_by_step','w')
for istep in range(force.shape[0]):
    wf.write("step %i:\n"%(istep+1))
    positions = pos[istep,:,:]
    forces = force[istep,:,:]
    wf.write('external pressure = {} kB\n'.format(press[istep]))
    wf.write('a = {: >12} b = {: >12} c = {: >12}\n'.format(lattice_params[istep,0],lattice_params[istep,1],lattice_params[istep,2]))
    for iatom in range(len(positions)):

        wf.write("{: >3} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12}\n".format(
            iatom,positions[iatom,0],positions[iatom,1],positions[iatom,2],forces[iatom,0],forces[iatom,1],forces[iatom,2]))
        
wf.close()
