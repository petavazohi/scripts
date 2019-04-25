#!/usr/bin/env python

import pymatgen.io.vasp as vasp
import pymatgen.electronic_structure.plotter as plotter
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--efermi",dest="efermi",type=float,help='fermi energy from SCF')
parser.add_argument("--Emin" ,dest="Emin",type=float ,action="store", help="low boundary of energy")
parser.add_argument("--Emax" ,dest="Emax",type=float ,action="store", help="up boundary of energy")
parser.add_argument("--show" ,dest="show",type=bool ,action="store", help="show plot",default= False)


args = parser.parse_args()

vasprun = vasp.Vasprun('vasprun.xml')
bs = vasprun.get_band_structure()
pl = plotter.BSPlotter(bs)

if args.efermi != None:
    vasprun.efermi = args.efermi
if args.Emin != None and args.Emax !=None:
    plt = pl.get_plot(ylim=[args.Emin,args.Emax]).axhline(y=0,linestyle='--',linewidth=1)
else:
    plt = pl.get_plot().axhline(y=0,linestyle='--',linewidth=1)
if args.show :
    plt.figure.show()
else :
    plt.figure.savefig('bands.pdf')
