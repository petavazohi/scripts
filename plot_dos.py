#!/usr/bin/env python

import pymatgen.io.vasp as vasp
import pymatgen.electronic_structure.plotter as plotter
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("--xmin" ,dest="xmin",type=float ,action="store", help="low boundary of energy")
parser.add_argument("--xmax" ,dest="xmax",type=float ,action="store", help="up boundary of energy")
parser.add_argument("--ymin" ,dest="ymin",type=float ,action="store", help="low boundary of intensity")
parser.add_argument("--ymax" ,dest="ymax",type=float ,action="store", help="up boundary of intensity")
parser.add_argument("--show" ,dest="show",type=bool ,action="store", help="show plot",default= False)


args = parser.parse_args()

vasprun = vasp.Vasprun('vasprun.xml')
pl = plotter.DosPlotter()
dos = vasprun.complete_dos
pl.add_dos(dos=dos,label='Total')
for element in dos.get_element_dos():
    pl.add_dos(dos=dos.get_element_dos()[element],label=element)

plt = pl.get_plot()
if args.xmin != None and args.xmax != None :
    plt.xlim(args.xmin,args.xmax)
if args.ymin != None and args.ymax != None :
    plt.ylim(args.ymin,args.ymax)
if args.show :
    plt.show()
else :
    plt.savefig('dos.pdf')

