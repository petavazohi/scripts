#!/usr/bin/env python3

import numpy as np
import xlsxwriter

import bibtexparser
import argparse
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# from pymatgen.core.structure import Structure
import pychemia
import spglib
from pymatgen.core import Element
import wyckoff_list
import string
import networkx as nx
periodic = pychemia.utils.periodic

def to_excel(db_calc,
             db_icsd,
             filename="xc_db_analysis.xlsx",
             symprec=1e-5,
             query={'status.relaxed': True}):

    sigfigs = int(-np.log10(symprec))
    # only lattice parameters
    head_l = ['XC',
              'Chemical formula',
              'SPG symbol',
              'SPG number',
              'Family',
              'Lattice DOF',
              'Bandgap', 'Direct gap',
              'nspecies', 'natoms',
              'a DOF', "b DOF", "c DOF",
              'a', 'a ICSD', 'a error',
              'b', 'b ICSD', 'b error',
              'c', 'c ICSD',  'c error',
               '(RSME abc).(DOF)', '(MAE abc).(DOF)','(ME abc).(DOF)',
              'alpha', 'beta', 'gamma',
              'alpha ICSD', 'beta ICSD', 'gamma ICSD',
              'Same SPG',
              'MP ID',"ICSD code"]
    
    # lattice + species
    head_ls = ['XC',
               'Chemical formula',
               'SPG symbol',
               'SPG number',
               'Family',
               'Lattice DOF',
               'Bandgap', 'Direct gap',
               'Species symbols',
               "Element class",
               'Element block',
               "Element group",
               'Element period',
               'Element electronegativity',
               'nspecies', 'natoms',
               'a DOF', "b DOF", "c DOF",
               'a', 'a ICSD', 'a error',
               'b', 'b ICSD', 'b error',
               'c', 'c ICSD',  'c error',
               '(RSME abc).(DOF)', '(MAE abc).(DOF)','(ME abc).(DOF)',
               'alpha', 'beta', 'gamma',
               'alpha ICSD', 'beta ICSD', 'gamma ICSD',
               'Same SPG',
               'MP ID',"ICSD code",]
    
    # only Wyckoff
    head_w = ['XC',
              'Atom symbol',
              'SPG symbol',
              'SPG number',
              'Family',
              'Bandgap', 'Direct gap',
              "Element class",
              'Element block',
              "Element group",
              'Element electronegativity',
              'Element period',
              'Wyckoff DOF',
              'nspecies', 'natoms',
              'x DOF', "y DOF", "z DOF",
              'Multiplicity',
              'Wyckoff letter',
              "Representation",
              'x', 'x ICSD', 'x error',
              'y', 'y ICSD', 'y error',
              'z', 'z ICSD', 'z error',
              '(RSME xyz).(DOF)', '(MAE xyz).(DOF)',               
              'Same SPG',
              "Chemical formula",
              'MP ID',"ICSD code", 
              ]
    # references
    head_r = ['MP ID',"ICSD code", "Chemical formula", "Chemical name"]

    # metadata about the studied structures
    head_m = ['SPG symbol',
              'SPG number',
              'Family',
              'Species symbols',
              "Element class",
              'Element block',
              "Element group",
              'Element electronegativity',
              'Element period',
              'nspecies', 'natoms',
              'Chemical formula',
              'MP ID',"ICSD code"]
    
    
    workbook = xlsxwriter.Workbook(filename)
    worksheet_l = workbook.add_worksheet("lattice")
    worksheet_ls = workbook.add_worksheet("lattice+species")
    worksheet_w = workbook.add_worksheet("wyckoff")
    worksheet_m = workbook.add_worksheet("metadata")
    worksheet_r = workbook.add_worksheet("references")
    
    worksheet_l.freeze_panes(1, 2)  # Freeze the first row.worksheet.
    worksheet_ls.freeze_panes(1, 2)  # Freeze the first row.worksheet.
    worksheet_w.freeze_panes(1, 2)  # Freeze the first row.worksheet.
    worksheet_r.freeze_panes(1, 0)  # Freeze the first row.worksheet.
    worksheet_m.freeze_panes(1, 0)  # Freeze the first row.worksheet.
    # create formatings
    header_format = workbook.add_format({'bold': True,
                                         'align': 'center',
                                         'valign': 'vcenter',
                                         })
    header_format.set_text_wrap()
    header_format.set_border()
    center_format_nf = workbook.add_format({'align': 'center'})
    center_format_nf.set_border()
    float_format_nf = workbook.add_format({'num_format': '#,##0.00', 'align': 'center'})
    float_format_nf.set_border()
    int_format_nf = workbook.add_format({'num_format': '#,##0', 'align': 'center'})
    int_format_nf.set_border()
    percent_format_nf = workbook.add_format({'num_format': '0.00%', 'align': 'center'})
    percent_format_nf.set_border()


    center_format_c = workbook.add_format({'align': 'center'})    
    center_format_c.set_bg_color('#C7E4EE')
    center_format_c.set_border()
    float_format_c = workbook.add_format({'num_format': '#,##0.00', 'align': 'center'})
    float_format_c.set_bg_color('#C7E4EE')
    float_format_c.set_border()
    int_format_c = workbook.add_format({'num_format': '#,##0', 'align': 'center'})
    int_format_c.set_bg_color('#C7E4EE')
    int_format_c.set_border()
    percent_format_c = workbook.add_format({'num_format': '0.00%', 'align': 'center'})
    percent_format_c.set_bg_color('#C7E4EE')
    percent_format_c.set_border()
    float_format_red = workbook.add_format({'num_format': '#,##0.0', 'align': 'center'})
    float_format_red.set_bg_color('#D9ECDB')
    float_format_red.set_border()

    
    col_number_l = {}
    for i, title in enumerate(head_l):
        worksheet_l.write_string(0, i, title, header_format)
        col_number_l[title] = i
    
    col_number_ls = {}
    for i, title in enumerate(head_ls):
        worksheet_ls.write_string(0, i, title, header_format)
        col_number_ls[title] = i

    col_number_w = {}
    for i, title in enumerate(head_w):
        worksheet_w.write_string(0, i, title, header_format)
        col_number_w[title] = i
        
    col_number_r = {}
    for i, title in enumerate(head_r):
        worksheet_r.write_string(0, i, title, header_format)
        col_number_r[title] = i

    col_number_m = {}
    for i, title in enumerate(head_m):
        worksheet_m.write_string(0, i, title, header_format)
        col_number_m[title] = i

    added_materials = []
    row_l = 1
    row_w = 1
    row_ls = 1
    row_r = 1
    row_m = 1
    
    for entry_calc in db_calc.entries.find(query):
        _id = str(entry_calc['_id'])
        mp_id = entry_calc['properties']['mp_id']
        xc = entry_calc['properties']['xc']
        # if mp_id != 'mp-23018' or xc != 'PBE':
        #     continue
        # spg_no = entry_calc['properties']['final']['crystal_symmetry']['number']
        # spg_symbol = entry_calc['properties']['final']['crystal_symmetry']['international_short']
        # family = entry_calc['properties']['final']['crystal_family']
        # gap = entry['properties']['band_gap']['total']['gap']
        gap_spin_channel = min([x for x in entry_calc['properties']['band_gap']])
        bandgap = entry_calc['properties']['band_gap'][gap_spin_channel]['gap']
        direct_gap = entry_calc['properties']['band_gap'][gap_spin_channel]['direct']
        entry_icsd = icsd_db.entries.find_one({'$and': [{"properties.mp_id":mp_id},
                                                        {"properties.standardized":True},
                                                        {"properties.theoretical":False}
                                                        ]})
        if entry_icsd:
            code_icsd = entry_icsd['properties']['code_ICSD']
            formula = entry_icsd['properties']['cif'][f'{code_icsd}-ICSD']["_chemical_formula_structural"]
            name = entry_icsd['properties']['cif'][f'{code_icsd}-ICSD']["_chemical_name_common"]
            st_icsd = icsd_db.get_structure(entry_icsd['_id'])
            st_icsd = shift(st_icsd, symprec)
            st_icsd = get_conventional(st_icsd, symprec)
            for i in range(st_icsd.natom):
                st_icsd = shift(st_icsd, symprec, np.random.randint(st_icsd.natom), True)
                st_icsd = sort_sites(st_icsd, symprec)
                st_icsd = shift(st_icsd, symprec, np.random.randint(st_icsd.natom), True)
                st_icsd = get_primitive(st_icsd, symprec)
                st_icsd = shift(st_icsd, symprec, np.random.randint(st_icsd.natom), True)
                st_icsd = get_conventional(st_icsd, symprec)
            st_icsd = shift(st_icsd, symprec)
            st_icsd = get_standardized(st_icsd, symprec)
            st_icsd = sort_sites(st_icsd, symprec)
            # st_icsd.sort_axes()
            reduced_icsd = st_icsd.reduced
            cs_icsd = pychemia.crystal.CrystalSymmetry(st_icsd)
            spg_no = cs_icsd.number(symprec)
            spg_symbol = cs_icsd.symbol(symprec)
            family = cs_icsd.crystal_system().lower()
            a_icsd = st_icsd.lattice.a
            b_icsd = st_icsd.lattice.b
            c_icsd = st_icsd.lattice.c
            alpha_icsd = st_icsd.lattice.alpha
            beta_icsd = st_icsd.lattice.beta
            gamma_icsd = st_icsd.lattice.gamma
            spg_no_icsd = spg_no # entry_icsd['properties']['crystal_symmetry']['number']
            
            st = db_calc.get_structure(entry_calc['_id'])
            st = shift(st, symprec)
            st = get_conventional(st, symprec)
            for i in range(st.natom):
                st = shift(st, symprec, np.random.randint(st.natom), True)
                st = sort_sites(st, symprec)
                st = shift(st, symprec, np.random.randint(st.natom), True)
                st = get_primitive(st, symprec)
                st = shift(st, symprec, np.random.randint(st.natom), True)
                st = get_conventional(st, symprec)
            st = shift(st, symprec)
            st = get_standardized(st, symprec)
            st = sort_sites(st, symprec)
            # st.sort_axes()
            cs = pychemia.crystal.CrystalSymmetry(st)
            if st.natom != st_icsd.natom:
                print("---------------------------------------------------------")
                print(st.composition,"|", st_icsd.composition)
                print(mp_id, spg_no, "|", code_icsd, spg_no_icsd)
                print("---------------------------------------------------------")
                continue
            a = st.lattice.a
            b = st.lattice.b
            c = st.lattice.c
            reduced = st.reduced
            alpha = st.lattice.alpha
            beta = st.lattice.beta
            gamma = st.lattice.gamma
            # print(st_icsd)
            # print(st)
            

            # print(st.composition,"|", st_icsd.composition)
            print(mp_id, spg_no, "|", code_icsd, spg_no_icsd)
            natom = st.natom
            nspecies = st.nspecies
            symbols = st.symbols
            wyckoff_analysis = get_wyckoffs(st, symprec, True)
            crystal = pychemia.crystal.CrystalSymmetry(st)
            crystal_family = crystal.crystal_system(symprec)
            inner = wyckoff_analysis['inner_degrees_of_freedom']
            lattice_degrees_of_freedom = get_lattice_degrees_of_freedom(crystal_family.lower(), st.lattice, sigfigs)
            outer = lattice_degrees_of_freedom
            multiplicity = wyckoff_analysis['multiplicity']
            wyckoff_letter = wyckoff_analysis['wyckoff_letter']
            wyckoff_rep = wyckoff_analysis['representation']
            dof_format = [int_format_nf, float_format_red]

            da = (a-a_icsd)/a
            db = (b-b_icsd)/b
            dc = (c-c_icsd)/c
            if sum(outer) != 0:
                rmse_abc = ((da*outer[0])**2+
                            (db*outer[1])**2+
                            (dc*outer[2])**2/sum(outer))**1/2
                mae_abc = (abs(da*outer[0])+
                           abs(db*outer[1])+
                           abs(dc*outer[2]))/sum(outer)
                me_abc = (da*outer[0]+
                          db*outer[1]+
                          dc*outer[2])/sum(outer)
            else :
                rmse_abc = None
                mae_abc = None
                me_abc= None
            
            worksheet_l.write_number(row_l, col_number_l['Lattice DOF'], sum(outer), int_format_nf)
            worksheet_l.write_number(row_l, col_number_l['a DOF'], outer[0], dof_format[int(outer[0])])
            worksheet_l.write_number(row_l, col_number_l['b DOF'], outer[1], dof_format[int(outer[1])])
            worksheet_l.write_number(row_l, col_number_l['c DOF'], outer[2], dof_format[int(outer[2])])
            worksheet_l.write_number(row_l, col_number_l['a'], a, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['b'], b, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['c'], c, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['alpha'], alpha, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['beta'], beta, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['gamma'], gamma, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['a ICSD'], a_icsd, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['b ICSD'], b_icsd, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['c ICSD'], c_icsd, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['alpha ICSD'], alpha_icsd, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['beta ICSD'], beta_icsd, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['gamma ICSD'], gamma_icsd, float_format_nf)
            worksheet_l.write_string(row_l, col_number_l["MP ID"], mp_id, center_format_nf)
            worksheet_l.write(row_l, col_number_l["ICSD code"], code_icsd, center_format_nf)
            worksheet_l.write_string(row_l, col_number_l["Chemical formula"], formula, center_format_nf)
            worksheet_l.write_string(row_l, col_number_l['XC'], xc, center_format_nf)
            worksheet_l.write_number(row_l, col_number_l['SPG number'], spg_no, int_format_nf)
            worksheet_l.write_string(row_l, col_number_l["SPG symbol"], spg_symbol, center_format_nf)
            worksheet_l.write_string(row_l, col_number_l["Family"], family, center_format_nf)
            worksheet_l.write_number(row_l, col_number_l['Bandgap'], bandgap, float_format_nf)
            worksheet_l.write_boolean(row_l, col_number_l['Direct gap'], direct_gap, center_format_nf)
            worksheet_l.write_number(row_l, col_number_l['nspecies'], nspecies, int_format_nf)
            worksheet_l.write_number(row_l, col_number_l['natoms'], natom, int_format_nf)
            worksheet_l.write_boolean(row_l, col_number_l['Same SPG'], spg_no == spg_no_icsd, center_format_nf)
            worksheet_l.write_number(row_l, col_number_l['a error'], da, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['b error'], db, float_format_nf)
            worksheet_l.write_number(row_l, col_number_l['c error'], dc, float_format_nf)
            worksheet_l.write(row_l, col_number_l['(RSME abc).(DOF)'], rmse_abc, percent_format_nf)
            worksheet_l.write(row_l, col_number_l['(MAE abc).(DOF)'], mae_abc, percent_format_nf)
            worksheet_l.write(row_l, col_number_l['(ME abc).(DOF)'], me_abc, percent_format_nf)
            row_l += 1

            for ispc in range(nspecies):
                if ispc == 0:
                    int_format = int_format_c
                    dof_format = [int_format_c, float_format_red]
                    float_format = float_format_c
                    center_format = center_format_c
                    percent_format = percent_format_c
                else:
                    int_format = int_format_nf
                    dof_format = [int_format_nf, float_format_red]
                    float_format = float_format_nf
                    center_format = center_format_nf
                    percent_format = percent_format_nf
                sym = st.species[ispc]
                worksheet_ls.write_string(row_ls, col_number_ls["Species symbols"], sym, center_format)                
                worksheet_ls.write_string(row_ls, col_number_ls['Element class'], get_element_class(sym), center_format)
                worksheet_ls.write_string(row_ls, col_number_ls["Element block"], periodic.block(sym), center_format)
                worksheet_ls.write_number(row_ls, col_number_ls['Element group'], periodic.group(sym), int_format)
                worksheet_ls.write_number(row_ls, col_number_ls['Element period'], periodic.period(sym), int_format)
                worksheet_ls.write_number(row_ls, col_number_ls['Element electronegativity'], periodic.electronegativity(sym), float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['Lattice DOF'], sum(outer), int_format)
                worksheet_ls.write_number(row_ls, col_number_ls['a DOF'], outer[0], dof_format[int(outer[0])])
                worksheet_ls.write_number(row_ls, col_number_ls['b DOF'], outer[1], dof_format[int(outer[1])])
                worksheet_ls.write_number(row_ls, col_number_ls['c DOF'], outer[2], dof_format[int(outer[2])])
                worksheet_ls.write_number(row_ls, col_number_ls['a'], a, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['b'], b, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['c'], c, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['alpha'], alpha, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['beta'], beta, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['gamma'], gamma, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['a ICSD'], a_icsd, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['b ICSD'], b_icsd, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['c ICSD'], c_icsd, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['alpha ICSD'], alpha_icsd, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['beta ICSD'], beta_icsd, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['gamma ICSD'], gamma_icsd, float_format)
                worksheet_ls.write_string(row_ls, col_number_ls["MP ID"], mp_id, center_format)
                worksheet_ls.write(row_ls, col_number_ls["ICSD code"], code_icsd, center_format)
                worksheet_ls.write_string(row_ls, col_number_ls["Chemical formula"], formula, center_format_nf)
                worksheet_ls.write_string(row_ls, col_number_ls['XC'], xc, center_format)
                worksheet_ls.write_number(row_ls, col_number_ls['SPG number'], spg_no, int_format)
                worksheet_ls.write_string(row_ls, col_number_ls["SPG symbol"], spg_symbol, center_format)
                worksheet_ls.write_string(row_ls, col_number_ls["Family"], family, center_format)
                worksheet_ls.write_number(row_ls, col_number_ls['Bandgap'], bandgap, float_format)
                worksheet_ls.write_boolean(row_ls, col_number_ls['Direct gap'], direct_gap, center_format)
                worksheet_ls.write_number(row_ls, col_number_ls['nspecies'], nspecies, int_format)
                worksheet_ls.write_number(row_ls, col_number_ls['natoms'], natom, int_format)
                worksheet_ls.write_boolean(row_ls, col_number_ls['Same SPG'], spg_no == spg_no_icsd, center_format)
                worksheet_ls.write_number(row_ls, col_number_ls['a error'], da, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['b error'], db, float_format)
                worksheet_ls.write_number(row_ls, col_number_ls['c error'], dc, float_format)
                worksheet_ls.write(row_ls, col_number_ls['(RSME abc).(DOF)'], rmse_abc, percent_format)
                worksheet_ls.write(row_ls, col_number_ls['(MAE abc).(DOF)'], mae_abc, percent_format)
                worksheet_ls.write(row_ls, col_number_ls['(ME abc).(DOF)'], me_abc, percent_format)

                row_ls += 1
                
            
            for iatom in range(natom):
                if iatom == 0:
                    dof_format = [int_format_c, float_format_red]
                    int_format = int_format_c
                    float_format = float_format_c
                    center_format = center_format_c
                    percent_format = percent_format_c
                else:
                    int_format = int_format_nf
                    dof_format = [int_format_nf, float_format_red]
                    float_format = float_format_nf
                    center_format = center_format_nf
                    percent_format = percent_format_nf
                sym = symbols[iatom]
                dx = diff_reduced_coordinate(reduced[iatom][0], reduced_icsd[iatom][0])
                dy = diff_reduced_coordinate(reduced[iatom][1], reduced_icsd[iatom][1])
                dz = diff_reduced_coordinate(reduced[iatom][2], reduced_icsd[iatom][2])
                if sum(inner[iatom]) != 0:
                    rmse_xyz = ((dx*inner[iatom][0])**2+
                                (dy*inner[iatom][1])**2+
                                (dz*inner[iatom][2])**2/sum(inner[iatom]))**1/2
                    mae_xyz = (abs(dx*inner[iatom][0])+
                               abs(dy*inner[iatom][1])+
                               abs(dz*inner[iatom][2]))/sum(inner[iatom])
                else:
                    rmse_xyz = None
                    mae_xyz = None
                
                worksheet_w.write_string(row_w, col_number_w['Element class'], get_element_class(sym), center_format)
                worksheet_w.write_string(row_w, col_number_w["Element block"], periodic.block(sym), center_format)
                worksheet_w.write_number(row_w, col_number_w['Element group'], periodic.group(sym), int_format)
                worksheet_w.write_number(row_w, col_number_w['Element period'], periodic.period(sym), int_format)
                worksheet_w.write_number(row_w, col_number_w['Element electronegativity'], periodic.electronegativity(sym), float_format)
                
                worksheet_w.write_string(row_w, col_number_w['XC'], xc, center_format)
                worksheet_w.write_number(row_w, col_number_w['SPG number'], spg_no, int_format)
                worksheet_w.write_string(row_w, col_number_w["SPG symbol"], spg_symbol, center_format)
                worksheet_w.write_string(row_w, col_number_w["Family"], family, center_format)
                worksheet_w.write_number(row_w, col_number_w['Bandgap'], bandgap, float_format)
                worksheet_w.write_boolean(row_w, col_number_w['Direct gap'], direct_gap, center_format)
                worksheet_w.write_string(row_w, col_number_w["Atom symbol"], symbols[iatom], center_format)
                worksheet_w.write_number(row_w, col_number_w['Wyckoff DOF'], sum(inner[iatom]), int_format)
                worksheet_w.write_number(row_w, col_number_w['nspecies'], nspecies, int_format)
                worksheet_w.write_number(row_w, col_number_w['natoms'], natom, int_format)
                worksheet_w.write_string(row_w, col_number_w["MP ID"], mp_id, center_format)
                worksheet_w.write(row_w, col_number_w["ICSD code"], code_icsd, center_format)
                worksheet_w.write_string(row_w, col_number_w["Chemical formula"], formula, center_format_nf)
                worksheet_w.write_number(row_w, col_number_w['x DOF'], inner[iatom][0], dof_format[int(inner[iatom][0])])
                worksheet_w.write_number(row_w, col_number_w['y DOF'], inner[iatom][1], dof_format[int(inner[iatom][1])])
                worksheet_w.write_number(row_w, col_number_w['z DOF'], inner[iatom][2], dof_format[int(inner[iatom][2])])

                worksheet_w.write_number(row_w, col_number_w['Multiplicity'], multiplicity[iatom], int_format)
                worksheet_w.write_string(row_w, col_number_w['Wyckoff letter'], wyckoff_letter[iatom], center_format)
                worksheet_w.write_string(row_w, col_number_w['Representation'], str(wyckoff_rep[iatom]), center_format)
                
                worksheet_w.write_number(row_w, col_number_w['x'], reduced[iatom][0], float_format)
                worksheet_w.write_number(row_w, col_number_w['y'], reduced[iatom][1], float_format)
                worksheet_w.write_number(row_w, col_number_w['z'], reduced[iatom][2], float_format)                
                
                worksheet_w.write_number(row_w, col_number_w['x ICSD'], reduced_icsd[iatom][0], float_format)
                worksheet_w.write_number(row_w, col_number_w['y ICSD'], reduced_icsd[iatom][1], float_format)
                worksheet_w.write_number(row_w, col_number_w['z ICSD'], reduced_icsd[iatom][2], float_format)
                worksheet_w.write_boolean(row_w, col_number_w['Same SPG'], spg_no == spg_no_icsd, center_format)
                worksheet_w.write_number(row_w, col_number_w['x error'], dx, float_format)
                worksheet_w.write_number(row_w, col_number_w['y error'], dy, float_format)
                worksheet_w.write_number(row_w, col_number_w['z error'], dz, float_format)
                worksheet_w.write(row_w, col_number_w['(RSME xyz).(DOF)'], rmse_xyz, percent_format)
                worksheet_w.write(row_w, col_number_w['(MAE xyz).(DOF)'], mae_xyz, percent_format)

                row_w += 1
                
            if mp_id not in added_materials:
                bib = bibtexparser.loads(entry_icsd['properties']['bib'])
                worksheet_r.write_string(row_r, col_number_r["MP ID"], mp_id, center_format_nf)
                worksheet_r.write_string(row_r, col_number_r["ICSD code"], code_icsd, center_format_nf)
                worksheet_r.write_string(row_r, col_number_r["Chemical formula"], formula, center_format_nf)
                worksheet_r.write_string(row_r, col_number_r["Chemical name"], name, center_format_nf)
                for x in bib.entries[0]:
                    if x in col_number_r:
                        worksheet_r.write(row_r, col_number_r[x], bib.entries[0][x], center_format_nf)
                    else:
                        i = len(col_number_r)
                        worksheet_r.write_string(0, i, x, header_format)
                        col_number_r[x] = i
                        worksheet_r.write(row_r, col_number_r[x], bib.entries[0][x], center_format_nf)
                row_r+=1
                for ispc in range(nspecies):
                    if ispc == 0:
                        int_format = int_format_c
                        float_format = float_format_c
                        center_format = center_format_c
                        percent_format = percent_format_c
                    else:
                        int_format = int_format_nf
                        float_format = float_format_nf
                        center_format = center_format_nf
                        percent_format = percent_format_nf
                    sym = st.species[ispc]
                    worksheet_m.write_string(row_m, col_number_m['Element class'], get_element_class(sym), center_format)
                    worksheet_m.write_string(row_m, col_number_m["Element block"], periodic.block(sym), center_format)
                    worksheet_m.write_number(row_m, col_number_m['Element group'], periodic.group(sym), int_format)
                    worksheet_m.write_number(row_m, col_number_m['Element period'], periodic.period(sym), int_format)
                    worksheet_m.write_number(row_m, col_number_m['Element electronegativity'], periodic.electronegativity(sym), float_format)
                    worksheet_m.write_string(row_m, col_number_m["Species symbols"], sym, center_format)
                    worksheet_m.write_string(row_m, col_number_m["MP ID"], mp_id, center_format)
                    worksheet_m.write(row_m, col_number_m["ICSD code"], code_icsd, center_format)
                    worksheet_m.write_number(row_m, col_number_m['SPG number'], spg_no, int_format)
                    worksheet_m.write_string(row_m, col_number_m["SPG symbol"], spg_symbol, center_format)
                    worksheet_m.write_string(row_m, col_number_m["Family"], family, center_format)
                    worksheet_m.write_number(row_m, col_number_m['nspecies'], nspecies, int_format)
                    worksheet_m.write_number(row_m, col_number_m['natoms'], natom, int_format)
                    worksheet_m.write_string(row_m, col_number_m['Chemical formula'], formula, center_format)
                    row_m += 1
                added_materials.append(mp_id)
    workbook.close()
    return

def get_wyckoffs(structure, symprec=1e-2, order=False):
    sigfigs = int(-np.log10(symprec))
    structure.reduced[structure.reduced.round(sigfigs) >= 1] -= 1
    structure.reduced[structure.reduced.round(sigfigs) < 0] += 1
    crystal_symmetry = pychemia.crystal.CrystalSymmetry(structure)
    spg = crystal_symmetry.number(symprec=symprec)
    wyckoffs = np.array(
        crystal_symmetry.get_symmetry_dataset(
            symprec=symprec)['wyckoffs'],
        dtype='<U3')
    equi_atoms= crystal_symmetry.get_symmetry_dataset(
            symprec=symprec)['equivalent_atoms']
    representation = wyckoff_list.get_wyckoff(spg)
    dof = np.zeros((structure.natom, 3), dtype=int)
    take_wyckoff = np.zeros((structure.natom))
    take_lattice = np.zeros((structure.natom))
    take_lattice[0] = 1
    wyck_idx = {}
    i = 1
    for letter in string.ascii_lowercase:
        wyck_idx[letter] = -1 * i
        i += 1
    

    for ieq in set(equi_atoms):
        idx = equi_atoms == ieq
        take_wyckoff[np.where(idx)[0][0]] = 1
    wyckoff_lookup = {}
    for i, rep in enumerate(representation[::-1]):
        m = len(rep)
        l = string.ascii_lowercase[i]
        wyckoff_lookup[f"{m}{l}"] = rep
    current_representation = []
    if order:
        identifiers = []
        for iatom, sym in enumerate(structure.symbols):
            l = wyckoffs[iatom]
            m = len(representation[wyck_idx[l]])
            # wyckoffs[iatom] = str(m) + l
            e = equi_atoms[iatom]
            identifiers.append(f"{m}{l}-{sym}-{e}")
        identifiers = np.array(identifiers)
        current_representation = np.empty_like(structure.reduced, 
                                               dtype=(np.str_,10))
        for _id in np.unique(identifiers):
            idx = identifiers == _id
            frac_coords = structure.reduced[idx].round(sigfigs)
            
            ml = _id.split('-')[0]

            # maps coordinates to the wyckoff_list
            mapper = match(frac_coords, wyckoff_lookup[ml], symprec)
            # print(mapper)
            for ix, iy in enumerate(np.where(idx)[0]):
                current_representation[iy] = [x.strip() for x in wyckoff_lookup[ml][mapper[ix]]]
    mmm = []
    www = []
    for iatom in range(structure.natom):
        wyckoff_letter = wyckoffs[iatom]
        www.append(wyckoff_letter)
        rep = representation[wyck_idx[wyckoff_letter]]
        multplicity = len(rep)
        mmm.append(multplicity)
        wyckoffs[iatom] = str(multplicity) + wyckoff_letter
        if len(current_representation)!=0:
            if take_wyckoff[iatom] == 1:
                for ix in range(3):
                    rep = current_representation[iatom][ix]
                    if 'x' in rep or 'y' in rep or 'z' in rep:
                        dof[iatom, ix] = 1

    return {'wyckoffs': wyckoffs[take_wyckoff.astype(np.bool_)].tolist(),
            'symmetry_species': np.array(structure.symbols)[take_wyckoff.astype(np.bool_)].tolist(),
            'inner_degrees_of_freedom': dof.tolist(),
            'representation': current_representation,
            'multiplicity': mmm,
            'wyckoff_letter': www,
            "equivalent_atoms": equi_atoms}

def get_primitive(st, symprec=1e-2):
    cs = pychemia.crystal.CrystalSymmetry(st)
    return cs.find_primitive(symprec)

def get_conventional(st, symprec=1e-2):
    cs = pychemia.crystal.CrystalSymmetry(st)
    return cs.refine_cell(symprec)

def get_standardized(st, symprec=1e-2, primitive=False):
    cs = pychemia.crystal.CrystalSymmetry(st)
    return cs.get_new_structure(spglib.standardize_cell(cs.spglib_cell,symprec=symprec))


def get_lattice_degrees_of_freedom(family, lattice, sigfigs):
    a, b, c = [round(x, sigfigs) for x in [lattice.a, lattice.b, lattice.c]]
    if family in ['cubic', 'trigonal']:
        if a==b and b==c:
            return [1, 0, 0]
        else:
            return [1/3, 1/3, 1/3]
    elif family in ['hexagonal', 'tetragonal']:
        if a == b or b == c:
            return [1, 0, 1]
        elif a == c :
            return [1, 1, 0]
        else:
            return [1/3, 1/3, 1/3]
    elif family in ['triclinic', 'monoclinic', 'orthorhombic']:
        return [1, 1, 1]

def sort_sites(st, symprec=1e-2, ascending=True):
    # First: Sort sites using the distance to the origin
    sigfigs = int(-np.log10(symprec))
    # st.reduced = st.reduced.round(sigfigs)
    st.reduced[st.reduced >= 1] -= 1
    st.reduced[st.reduced < 0] += 1
    reduced = st.reduced.copy()
    reduced = reduced.round(1)
    norms = [np.linalg.norm(x).round(sigfigs) for x in st.reduced]
    x = reduced[:, 0]
    y = reduced[:, 1]
    z = reduced[:, 2]
    wyckoff_analysis = get_wyckoffs(st, symprec, False)
    equi_atoms= wyckoff_analysis['equivalent_atoms']
    symbols = []
    multiplicity = wyckoff_analysis['multiplicity']
    wyckoff_letter = wyckoff_analysis['wyckoff_letter']
    for iatom in range(st.natom):
        symbols.append(f"{multiplicity[iatom]}{wyckoff_letter[iatom]}{st.symbols[iatom]}{equi_atoms[iatom]}")
    
    # print sorted_indices
    to_sort = np.array(list(zip(symbols, x, y, z, norms)),
                       dtype=[('symbols', np.unicode_, 16),
                              ('x', np.float64),
                              ('y', np.float64),
                              ('z', np.float64),
                              ('norm', np.float64),])
    idx = np.argsort(to_sort, order=('symbols', 'x', 'y', 'z', 'norm'), )
    if ascending:
        idx = np.flip(idx)
    # Second: Sort again using the atomic number
    st.sort_sites_using_list(idx)
    return st.copy()


def shift(st, symprec=1e-2, which=0, random=False):
    st = sort_sites(st, symprec, False)
    shift_selections = [-3/4, -2/3, -1/2, -1/3, -1/4, 1/4, 1/3, 1/2, 2/3, 3/4]
    if random:
        a1 = np.random.choice(shift_selections)
        a2 = np.random.choice(shift_selections)
        a3 = np.random.choice(shift_selections)
    else:
        a1 = 0.0
        a2 = 0.0
        a3 = 0.0
    reduced = np.zeros_like(st.reduced)
    reduced[:, 0] = st.reduced[:, 0] - st.reduced[which,0] + a1
    reduced[:, 1] = st.reduced[:, 1] - st.reduced[which,1] + a2
    reduced[:, 2] = st.reduced[:, 2] - st.reduced[which,2] + a3
    st = pychemia.core.Structure(symbols=st.symbols, cell=st.cell, reduced=reduced)
    
    return st

def match(frac_coords, wyckoffs, symprec=1e-2):
    sigfigs = int(-np.log10(symprec))
    G = nx.Graph()
    def conv(xyz):
        if ('x' not in xyz) and ('y' not in xyz) and ('z' not in xyz):
            
            fraction = xyz.strip().split('/')
            
            if len(fraction) ==2:
                ret=float(fraction[0])/float(fraction[1])
            else:
                ret= float(fraction[0])
            if ret>=1:
                ret-=1
            elif ret<0:
                ret+=1
            return round(ret , sigfigs)
        else:
            return xyz
    def is_number(x):
        return type(x) is float
    mapper = {x:[] for x,_ in enumerate(wyckoffs)}
    for ipos, rep in enumerate(wyckoffs):
        x = conv(rep[0])
        y = conv(rep[1])
        z = conv(rep[2])
        G.add_node(f"{ipos}-w", coords=[x,y,z])
    for jpos, fc in enumerate(frac_coords):
        G.add_node(f"{jpos}-c", coords=fc)
    

    for ipos, rep in enumerate(wyckoffs):
        x = conv(rep[0])
        y = conv(rep[1])
        z = conv(rep[2])
        for jpos, fc in enumerate(frac_coords):
            xp = fc[0]
            yp = fc[1]
            zp = fc[2]
            
            if not is_number(x):
                if not is_number(y):
                    if not is_number(z): 
                        continue
                    else : 
                        if z==zp:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=1)
                else:
                    if not is_number(z):
                        if y==yp:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=1)
                    else :
                        if y==yp and z==zp:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=2)
            else:
                if not is_number(y):
                    if not is_number(z):
                        if x==xp:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=1)
                    else :
                        if x==xp and z==zp:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=2)
                else:
                    if not is_number(z):
                        if x==xp and y==yp:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=2)
                    else:
                        if x==xp and y==yp and z==zp:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=3)
    for ipos, rep in enumerate(wyckoffs):
        x = conv(rep[0])
        y = conv(rep[1])
        z = conv(rep[2])
        if not is_number(x):
            if not is_number(y):
                if not is_number(z): 
                    for jpos, fc in enumerate(frac_coords):
                        xp = fc[0]
                        yp = fc[1]
                        zp = fc[2]
                        if len(G.edges(f"{jpos}-c")) == 0:
                            G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=1)
                            break            

    for inode in G.nodes():
        edges = np.array([e for e in G.edges(inode, data='weight')])
        if len(edges) in [0,1]:
            continue
        weights = edges[:, 2].astype(int)
        idx = np.flip(np.argsort(weights))
        for i in range(1, len(idx)):
            u = edges[idx][i][0]
            v = edges[idx][i][1]
            G.remove_edge(u, v)
        jnode = edges[idx][0][1]
        edges = np.array([e for e in G.edges(jnode, data='weight')])
        for e in edges:
            u = e[0]
            v = e[1]
            if u==jnode and v==inode:
                continue
            else:
                G.remove_edge(u, v)
    for ipos, rep in enumerate(wyckoffs):
        if len(G.edges(f"{ipos}-w")) == 0:
            for jpos, fc in enumerate(frac_coords):
                if len(G.edges(f"{jpos}-c")) == 0:
                    G.add_edge(f"{ipos}-w", f"{jpos}-c", weight=1)
                    break    
    mapper = {}
    for e in G.edges():
        u = int(e[0].split("-")[0])
        v = int(e[1].split("-")[0])
        mapper[v]=u
    return mapper

                    
def get_element_class(sym):
    ele = Element(sym)
    if ele.is_actinoid:
        return "Actinoid"
    elif ele.is_alkali:
        return "Alkali metal"
    elif ele.is_alkaline:
        return "Alkaline metal"
    elif sym in ["C", "N", 'P', 'O', 'S', 'Se', 'F', 'Cl', 'Br', 'I', 'H']:
        return "Reactive nonmetals"
    elif ele.is_lanthanoid:
        return "Lanthanoid"
    elif ele.is_metalloid:
        return "Metalloid"
    elif ele.is_noble_gas:
        return "Noble gas"
    elif ele.is_post_transition_metal:
        return "Post-transition metal"
    elif ele.is_transition_metal:
        return "Transition metal"
    

    
def diff_reduced_coordinate(x, x_prime):
    diff = abs(x - x_prime)
    if diff<=0.5:
        ret = diff
    elif diff>0.5:
        if x>=0.5:
            ret = abs(x-1-x_prime)
        elif x<0.5:
            ret = abs(x+1-x_prime)
    return ret
    
if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--icsd_db",
        dest="icsd_db_name",
        type=str,
        help='Name of the ICSD MongoDB database')
    parser.add_argument(
        "--compare_db",
        dest="comp_db_name",
        type=str,
        help="Name of the comparison MongoDB database")
    parser.add_argument(
        "--host",
        dest="host",
        type=str,
        help="MongoDB hostname",
        default='localhost')
    parser.add_argument(
        "--port",
        dest="port",
        type=int,
        help="MongoDB port",
        default=27020)
    parser.add_argument(
        "--user",
        dest='user',
        type=str,
        help="MongoDB username",
        default=None)
    parser.add_argument(
        "--passwd",
        dest='passwd',
        type=str,
        help="MongoDB password",
        default=None)
    parser.add_argument(
        "--output",
        dest='output',
        type=str,
        help="Output file name",
        default="output.xlsx")
    parser.add_argument(
        "--verbose",
        action=argparse.BooleanOptionalAction,
        help='To print the results of the on going analysis.',
        default="--verbose")
    args = parser.parse_args()
    icsd_db = pychemia.db.get_database(dict(name=args.icsd_db_name,
                                            host=args.host,
                                            port=args.port,
                                            user=args.user,
                                            passwd=args.passwd))
    com_db = pychemia.db.get_database(dict(name=args.comp_db_name,
                                            host=args.host,
                                            port=args.port,
                                            user=args.user,
                                            passwd=args.passwd))
    to_excel(com_db, icsd_db, args.output, symprec=1e-3)
    
