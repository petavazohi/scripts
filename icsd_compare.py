#!/usr/bin/env python3
import pychemia
from collections import OrderedDict
import numpy as np
import xlsxwriter
import bibtexparser
import argparse
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# from pymatgen.core.structure import Structure
from pymatgen.core import Element
import wyckoff_list
import string
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
              'a', 'b', 'c',
              'alpha', 'beta', 'gamma',
              'a ICSD', 'b ICSD', 'c ICSD',
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
               'a', 'b', 'c',
               'alpha', 'beta', 'gamma',
               'a ICSD', 'b ICSD', 'c ICSD',
               'alpha ICSD', 'beta ICSD', 'gamma ICSD',
               'Same SPG',
               'MP ID',"ICSD code"]
    
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
              'x', 'y', 'z',
              'x ICSD', 'y ICSD', 'z ICSD',
              'Same SPG',
              "Chemical formula",
              'MP ID',"ICSD code"]
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
    precent_format_nf = workbook.add_format({'num_format': '0.00%', 'align': 'center'})
    precent_format_nf.set_border()


    center_format_c = workbook.add_format({'align': 'center'})    
    center_format_c.set_bg_color('#C7E4EE')
    center_format_c.set_border()
    float_format_c = workbook.add_format({'num_format': '#,##0.00', 'align': 'center'})
    float_format_c.set_bg_color('#C7E4EE')
    float_format_c.set_border()
    int_format_c = workbook.add_format({'num_format': '#,##0', 'align': 'center'})
    int_format_c.set_bg_color('#C7E4EE')
    int_format_c.set_border()
    precent_format_c = workbook.add_format({'num_format': '0.00%', 'align': 'center'})
    precent_format_c.set_bg_color('#C7E4EE')
    precent_format_c.set_border()
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

        if mp_id != 'mp-27507' or xc != 'PBE':
            continue
        print(entry_calc['properties']['path'])
        spg_no = entry_calc['properties']['final']['crystal_symmetry']['number']
        spg_symbol = entry_calc['properties']['final']['crystal_symmetry']['international_short']
        family = entry_calc['properties']['final']['crystal_family']
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
            print(entry_icsd['properties']['path'])
            formula = entry_icsd['properties']['cif'][f'{code_icsd}-ICSD']["_chemical_formula_structural"]
            name = entry_icsd['properties']['cif'][f'{code_icsd}-ICSD']["_chemical_name_common"]
            st_icsd = icsd_db.get_structure(entry_icsd['_id'])
            # st_icsd.sort_axes()
            cs_icsd = pychemia.crystal.CrystalSymmetry(st_icsd)
            st_icsd = cs_icsd.refine_cell(symprec)
            st_icsd = sort_sites(st_icsd.copy())
            a_icsd = st_icsd.lattice.a
            b_icsd = st_icsd.lattice.b
            c_icsd = st_icsd.lattice.c
            reduced_icsd = st_icsd.reduced.round(sigfigs)
            reduced_icsd[reduced_icsd >= 1] -= 1
            reduced_icsd[reduced_icsd < 0] += 1
            alpha_icsd = st_icsd.lattice.alpha
            beta_icsd = st_icsd.lattice.beta
            gamma_icsd = st_icsd.lattice.gamma
            spg_no_icsd = entry_icsd['properties']['crystal_symmetry']['number']
            
            st = db_calc.get_structure(entry_calc['_id'])
            # st.sort_axes()
            cs = pychemia.crystal.CrystalSymmetry(st)
            st = cs.refine_cell(symprec)
            st = sort_sites(st.copy())
            # st = standardize(st, symprec)
            if st.natom != st_icsd.natom:
                print("---------------------------------------------------------")
                print(st.composition,"|", st_icsd.composition)
                print(mp_id, spg_no, "|", code_icsd, spg_no_icsd)
                print("---------------------------------------------------------")
                continue
            a = st.lattice.a
            b = st.lattice.b
            c = st.lattice.c
            reduced = st.reduced.round(sigfigs)
            reduced[reduced >= 1] -= 1
            reduced[reduced < 0] += 1
            alpha = st.lattice.alpha
            beta = st.lattice.beta
            gamma = st.lattice.gamma
            # print(reduced)
            # print(reduced_icsd)
            

            # print(st.composition,"|", st_icsd.composition)
            print(mp_id, spg_no, "|", code_icsd, spg_no_icsd)
            natom = st.natom
            nspecies = st.nspecies
            symbols = st.symbols
            wyckoff_analysis = get_wyckoffs(st, symprec)
            crystal = pychemia.crystal.CrystalSymmetry(st)
            crystal_family = crystal.crystal_system(symprec)
            inner = wyckoff_analysis['inner_degrees_of_freedom']
            lattice_degrees_of_freedom = get_lattice_degrees_of_freedom(crystal_family.lower(), st.lattice, sigfigs)
            outer = lattice_degrees_of_freedom
            multiplicity = wyckoff_analysis['multiplicity']
            wyckoff_letter = wyckoff_analysis['wyckoff_letter']
            dof_format = [int_format_nf, float_format_red]
            
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
            worksheet_l.write_string(row_l, col_number_l["ICSD code"], code_icsd, center_format_nf)
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
            row_l += 1

            for ispc in range(nspecies):
                if ispc == 0:
                    int_format = int_format_c
                    dof_format = [int_format_c, float_format_red]
                    float_format = float_format_c
                    center_format = center_format_c
                else:
                    int_format = int_format_nf
                    dof_format = [int_format_nf, float_format_red]
                    float_format = float_format_nf
                    center_format = center_format_nf
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
                worksheet_ls.write_string(row_ls, col_number_ls["ICSD code"], code_icsd, center_format)
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
                row_ls += 1
                
            
            for iatom in range(natom):
                if iatom == 0:
                    dof_format = [int_format_c, float_format_red]
                    int_format = int_format_c
                    float_format = float_format_c
                    center_format = center_format_c
                else:
                    int_format = int_format_nf
                    dof_format = [int_format_nf, float_format_red]
                    float_format = float_format_nf
                    center_format = center_format_nf
                sym = symbols[iatom]
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
                worksheet_w.write_string(row_w, col_number_w["ICSD code"], code_icsd, center_format)
                worksheet_w.write_string(row_w, col_number_w["Chemical formula"], formula, center_format_nf)
                worksheet_w.write_number(row_w, col_number_w['x DOF'], inner[iatom][0], dof_format[int(inner[iatom][0])])
                worksheet_w.write_number(row_w, col_number_w['y DOF'], inner[iatom][1], dof_format[int(inner[iatom][1])])
                worksheet_w.write_number(row_w, col_number_w['z DOF'], inner[iatom][2], dof_format[int(inner[iatom][2])])

                worksheet_w.write_number(row_w, col_number_w['Multiplicity'], multiplicity[iatom], int_format)
                worksheet_w.write_string(row_w, col_number_w['Wyckoff letter'], wyckoff_letter[iatom], center_format)
                                
                worksheet_w.write_number(row_w, col_number_w['x'], reduced[iatom][0], float_format)
                worksheet_w.write_number(row_w, col_number_w['y'], reduced[iatom][1], float_format)
                worksheet_w.write_number(row_w, col_number_w['z'], reduced[iatom][2], float_format)                
                
                worksheet_w.write_number(row_w, col_number_w['x ICSD'], reduced_icsd[iatom][0], float_format)
                worksheet_w.write_number(row_w, col_number_w['y ICSD'], reduced_icsd[iatom][1], float_format)
                worksheet_w.write_number(row_w, col_number_w['z ICSD'], reduced_icsd[iatom][2], float_format)
                print(f"{reduced[iatom][0]}, {reduced[iatom][1]}, {reduced[iatom][2]} | {reduced_icsd[iatom][0]}, {reduced_icsd[iatom][1]}, {reduced_icsd[iatom][2]}")
                worksheet_w.write_boolean(row_w, col_number_w['Same SPG'], spg_no == spg_no_icsd, center_format)
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
                    else:
                        int_format = int_format_nf
                        float_format = float_format_nf
                        center_format = center_format_nf
                    sym = st.species[ispc]
                    worksheet_m.write_string(row_m, col_number_m['Element class'], get_element_class(sym), center_format)
                    worksheet_m.write_string(row_m, col_number_m["Element block"], periodic.block(sym), center_format)
                    worksheet_m.write_number(row_m, col_number_m['Element group'], periodic.group(sym), int_format)
                    worksheet_m.write_number(row_m, col_number_m['Element period'], periodic.period(sym), int_format)
                    worksheet_m.write_number(row_m, col_number_m['Element electronegativity'], periodic.electronegativity(sym), float_format)
                    worksheet_m.write_string(row_m, col_number_m["Species symbols"], sym, center_format)
                    worksheet_m.write_string(row_m, col_number_m["MP ID"], mp_id, center_format)
                    worksheet_m.write_string(row_m, col_number_m["ICSD code"], code_icsd, center_format)
                    worksheet_m.write_number(row_m, col_number_m['SPG number'], spg_no, int_format)
                    worksheet_m.write_string(row_m, col_number_m["SPG symbol"], spg_symbol, center_format)
                    worksheet_m.write_string(row_m, col_number_m["Family"], family, center_format)
                    worksheet_m.write_number(row_m, col_number_m['nspecies'], nspecies, int_format)
                    worksheet_m.write_number(row_m, col_number_m['natoms'], natom, int_format)
                    worksheet_m.write_string(row_m, col_number_m['Chemical formula'], formula, center_format)
                    row_m += 1
                added_materials.append(mp_id)
        # df = pd.DataFrame.from_dict(ret)

        # df['x error'] = [diff_reduced_coordinate(x, x_prime) for x, x_prime in zip(df['x'], df['x ICSD'])]
        # df['y error'] = [diff_reduced_coordinate(x, x_prime) for x, x_prime in zip(df['y'], df['y ICSD'])]
        # df['z error'] = [diff_reduced_coordinate(x, x_prime) for x, x_prime in zip(df['z'], df['z ICSD'])]
        # df['a error'] = (df['a'] - df['a ICSD'])/df['a']
        # df['b error'] = (df['b'] - df['b ICSD'])/df['b']
        # df['c error'] = (df['c'] - df['c ICSD'])/df['c']
        # df['(RSME xyz).(DOF)'] = (((df['x error']*df['x freedom'])**2+
        #                            (df['y error']*df['y freedom'])**2+
        #                            (df['z error']*df['z freedom'])**2)/(df['No. Wyckoff degrees of freedom']))**(1/2)
        # df['(MAE xyz).(DOF)'] = (df['x error']*df['x freedom']+
        #                          df['y error']*df['y freedom']+
        #                          df['z error']*df['z freedom'])/(df['No. Wyckoff degrees of freedom'])

        # df['(RSME abc).(DOF)'] = (((df['a error']*df['a freedom'])**2+
        #                            (df['b error']*df['b freedom'])**2+
        #                            (df['c error']*df['c freedom'])**2)/(df['No. lattice degrees of freedom']))**(1/2)
        # df['(MAE abc).(DOF)'] = (abs(df['a error'])*df['a freedom']+
        #                          abs(df['b error'])*df['b freedom']+
        #                          abs(df['c error'])*df['c freedom'])/(df['No. lattice degrees of freedom'])
        # df['RSME abc'] = (((df['a error'])**2+
        #                    (df['b error'])**2+
        #                    (df['c error'])**2)/3)**(1/2)*df['No. lattice degrees of freedom']/df['No. lattice degrees of freedom']
        # df['ME abc'] = (abs(df['a error'])+
        #                 abs(df['b error'])+
        #                 abs(df['c error']))/3*df['No. lattice degrees of freedom']/df['No. lattice degrees of freedom']
        # df.to_excel(f'{self.name}-db.xlsx')
    workbook.close()
    return

# def standardize(structure, symprec=1e-5):
#     st = Structure(lattice = structure.cell, species=structure.symbols, coords=structure.reduced)
#     st = get_primitive(st, symprec)
#     # st = get_conventioanl(st, symprec)
#     symbols = [site.symbol for site in st.species]
#     cell = st.lattice.matrix
#     reduced = st.frac_coords
#     return pychemia.core.Structure(cell=cell, symbols=symbols, reduced=reduced)


# def get_primitive(structure, symprec=1e-5):
#     """The input and output are pymatgen structures
#     """
#     spg_an = SpacegroupAnalyzer(structure, symprec=symprec)
#     return spg_an.get_primitive_standard_structure()

# def get_conventioanl(structure, symprec=1e-5):
#     """The input and output are pymatgen structures
#     """
#     spg_an = SpacegroupAnalyzer(structure, symprec=symprec)
#     return spg_an.get_conventional_standard_structure()


def get_wyckoffs(structure, symprec=1e-5):
    crystal_symmetry = pychemia.crystal.CrystalSymmetry(structure)
    spg = crystal_symmetry.number(symprec=symprec)
    wyckoffs = np.array(
        crystal_symmetry.get_symmetry_dataset(
            symprec=symprec)['wyckoffs'],
        dtype='<U3')
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

    for ispc in structure.species:
        idx = np.array(structure.symbols) == ispc
        for wyckoff_letter in np.unique(wyckoffs[idx]):
            idx2 = wyckoffs[idx] == wyckoff_letter
            take_wyckoff[np.where(idx)[0][idx2][0]] = 1
    current_representation = []
    mmm = []
    www = []
    for iatom in range(structure.natom):
        wyckoff_letter = wyckoffs[iatom]
        www.append(wyckoff_letter)
        rep = representation[wyck_idx[wyckoff_letter]]
        multplicity = len(rep)
        mmm.append(multplicity)
        wyckoffs[iatom] = str(multplicity) + wyckoff_letter
        if take_wyckoff[iatom] == 1:
            important = False
            if 'x' in str(rep):
                dof[iatom, 0] = 1
                important = True
            if 'y' in str(rep):
                dof[iatom, 1] = 1
                important = True
            if 'z' in str(rep):
                dof[iatom, 2] = 1
                important = True
            if important:
                current_representation.append(rep)

    return {'wyckoffs': wyckoffs[take_wyckoff.astype(np.bool_)].tolist(),
            'symmetry_species': np.array(structure.symbols)[take_wyckoff.astype(np.bool_)].tolist(),
            'inner_degrees_of_freedom': dof.tolist(),
            'representation': current_representation,
            'multiplicity': mmm,
            'wyckoff_letter': www}


def get_primitive(st, symprec=1e-5):
    cs = pychemia.crystal.CrystalSymmetry(st)
    return cs.find_primitive(symprec)

def get_conventional(st, symprec=1e-5):
    cs = pychemia.crystal.CrystalSymmetry(st)
    return cs.refine_cell(symprec)


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

def sort_sites(st):
    # First: Sort sites using the distance to the origin
    norms = [np.linalg.norm(x).round(2) for x in st.reduced]
    x = st.reduced[:, 0].round(2)
    y = st.reduced[:, 1].round(2)
    z = st.reduced[:, 2].round(2)
    symbols = st.symbols
    # print sorted_indices
    to_sort = np.array(list(zip(symbols, norms, x, y, z)),
                       dtype=[('symbols', np.unicode_, 16),
                              ('norm', np.float64),
                              ('x', np.float64),
                              ('y', np.float64),
                              ('z', np.float64)])
    print(to_sort)
    idx = np.argsort(to_sort, order=('symbols', 'norm', 'x', 'y', 'z'))
    # Second: Sort again using the atomic number
    st.sort_sites_using_list(idx)

    return st.copy()
                                                                        



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
    to_excel(com_db, icsd_db, args.output, symprec=1e-2)
    
