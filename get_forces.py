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
pos_force = np.array(pos_force)
pos = pos_force[:,:,0:3]
force = pos_force[:,:,3:6]
wf = open('FORCE_POS_by_atm','w')

for iatom in range(force.shape[1]):
    wf.write("atom  %i\n"% iatom)
    positions = pos[:,iatom,:]
    forces = force[:,iatom,:]
    for istep in range(len(positions)):
                       wf.write("{: > 3} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12}\n".format(
                           istep,positions[istep,0],positions[istep,1],positions[istep,2],forces[istep,0],forces[istep,1],forces[istep,2]))
wf.close()

wf = open('FORCE_POS_by_step','w')
for istep in range(force.shape[0]):
    wf.write("step %i:\n"%istep)
    positions = pos[istep,:,:]
    forces = force[istep,:,:]
    for iatom in range(len(positions)):

        wf.write("{: >3} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12}\n".format(
            iatom,positions[iatom,0],positions[iatom,1],positions[iatom,2],forces[iatom,0],forces[iatom,1],forces[iatom,2]))
        
wf.close()
