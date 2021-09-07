#!/usr/bin/env python3
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Composition
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--key" ,dest="MPkey",type=str ,action="store", help="materials project key")
parser.add_argument("--id",dest='mp_id',type=str,action='store',help='materials project id')

args = parser.parse_args()
mpr = MPRester(args.MPkey)
structure = mpr.get_structure_by_material_id(args.mp_id)
structure.to (fmt='POSCAR',filename="{}-{}.vasp".format(args.mp_id,structure.composition.reduced_formula))

