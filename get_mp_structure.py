#!/usr/bin/env python3

from pymatgen import MPRester, Composition
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--MPkey" ,dest="MPkey",type=str ,action="store", help="materials project key")
parser.add_argument("--mp_id",dest='mp_id',type=str,action='store',help='materials project id')

args = parser.parse_args()
mpr = MPRester(args.MPkey)
structure = mpr.get_structure_by_material_id(args.mp_id)
structure.to (fmt='POSCAR',filename="POSCAR")
