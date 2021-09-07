#!/usr/bin/env python3

from pymatgen.io import vasp
import numpy as np
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
import pymatgen

run = vasp.Vasprun('vasprun.xml')
outcar = vasp.Outcar('OUTCAR')
st = run.final_structure
symbols = np.array([x.value for x in st.species])
sa = SpacegroupAnalyzer(st)
print('**************')
print("conventional cell")
print(sa.get_conventional_standard_structure())
print(sa.get_conventional_standard_structure().volume)
print(sa.get_space_group_symbol(),sa.get_space_group_number())
print("******************")
print("Final structure")
print(run.final_structure)
print("vol={}".format(run.final_structure.volume))
print("******************")
print("Bandgap information")
print(run.get_band_structure().get_band_gap())
print("******************")
print("Magnetzation per site")
for ielement in range(len(symbols)):
    print("{}-{} {}".format(symbols[ielement], ielement,
                            outcar.magnetization[ielement]['tot']))
print("******************")
print("Magnetzation per species")
for element in st.symbol_set:
    temp = np.array(outcar.magnetization)[np.where(symbols == element)[0]]
    print("{} {}".format(element, np.average([abs(x['tot']) for x in temp])))
