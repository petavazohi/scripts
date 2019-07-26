#!/usr/bin/env python3

import pymatgen.io.vasp as vasp
import pymatgen.electronic_structure.plotter as plotter
import argparse
import sys

def normalize(dos,factor):
    for spin in dos.densities:
        dos.densities[spin] /=  factor
    return dos
    

parser = argparse.ArgumentParser()
parser.add_argument("--xmin" ,dest="xmin",type=float ,action="store", help="low boundary of energy")
parser.add_argument("--xmax" ,dest="xmax",type=float ,action="store", help="up boundary of energy")
parser.add_argument("--ymin" ,dest="ymin",type=float ,action="store", help="low boundary of intensity")
parser.add_argument("--ymax" ,dest="ymax",type=float ,action="store", help="up boundary of intensity")
parser.add_argument("--show" ,dest="show",type=bool ,action="store", help="show plot",default= False)
parser.add_argument("--patoms",dest="patoms",nargs='+',action="store",help="atoms you want to projected on the dos")
parser.add_argument("--norm",dest='norm',type=float,action="store",help="normalization factor for the density of states to be devided by",default=1)
args = parser.parse_args()
outcar = vasp.Outcar('OUTCAR')

vasprun = vasp.Vasprun('vasprun.xml')
pl = plotter.DosPlotter()
dos = vasprun.complete_dos

# for spin in dos.densities:
#     integral = 0
#     ndos = dos.densities[spin].shape[0]
#     for i in range(1,ndos):
#         if dos.energies[i-1] < vasprun.efermi:
#             integral += (dos.densities[spin][i]+ dos.densities[spin][i-1])/2*(dos.energies[i]-dos.energies[i-1])
#     print(int(round(outcar.nelect)),integral)
#     dos.densities[spin] *=  int(outcar.nelect)/integral

dos = normalize(dos,args.norm)


pl.add_dos(dos=dos,label='Total')
if args.patoms != None :
    for element in dos.get_element_dos():
        if element.name in args.patoms :
            par_dos = normalize(dos.get_element_dos()[element],args.norm)
            pl.add_dos(dos=par_dos,label=element)

plt = pl.get_plot()
plt.axvline(x=0,linestyle='--',c='k')
if args.xmin != None and args.xmax != None :
    plt.xlim(args.xmin,args.xmax)
if args.ymin != None and args.ymax != None :
    plt.ylim(args.ymin,args.ymax)
if args.show :
    plt.show()
else :
    plt.savefig('dos.pdf')

