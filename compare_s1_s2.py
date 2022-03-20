#!/usr/bin/env python3
import pychemia
query, 

from icsd_compare import sort_sites = {'$and': [{"properties.mp_id":"mp-1"},
                                                {"properties.standardized":True},
                                                {"properties.theoretical":False},
                                                {'properties.xc':'PBE'},
                                                {'properties.code_ICSD'},
                                                ]
                                       }




def match_structures(st1, st2):
    if st1.natom != st2.natom:
        print("not same size, compare later")
        return

mp_id = "mp-27507"
code_icsd = "1004"
icsd_db = pychemia.db.get_database(dict(name='ICSD',
                                        port=27020,))
com_db = pychemia.db.get_database(dict(name="XC_",
                                        port=27020,))
sigfigs = 2


entry_calc = com_db.entries.find_one({'$and': [{"properties.mp_id":mp_id},
                                               {'properties.xc':'PBE'},
                                               ]
                                      }
                                     )
print(entry_calc['properties']['path'])
_id = entry_calc['_id']
st1 = com_db.get_structure(_id)
cs1 = pychemia.crystal.CrystalSymmetry(st1)
st1 = cs1.refine_cell(1e-2)
st1 = sort_sites(st1)
# print(st1)
reduced = st1.reduced.round(sigfigs)


entry_icsd = icsd_db.entries.find_one({'$and': [{"properties.standardized":True},
                                                {"properties.theoretical":False},
                                                {'properties.code_ICSD':code_icsd},
                                                ]
                                       }
                                      )
_id = entry_icsd['_id']
st2 = icsd_db.get_structure(_id)
cs2 = pychemia.crystal.CrystalSymmetry(st2)
st2 = cs2.refine_cell(1e-2)
st2 = sort_sites(st2)
reduced_icsd = st2.reduced.round(sigfigs)
print(entry_icsd['properties']['path'])
print(mp_id, "|", code_icsd)

for iatom in range(st1.natom):
    print(f"{reduced[iatom][0]}, {reduced[iatom][1]}, {reduced[iatom][2]} | {reduced_icsd[iatom][0]}, {reduced_icsd[iatom][1]}, {reduced_icsd[iatom][2]}")
#print(st2)


