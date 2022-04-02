#!/usr/bin/env python3

import numpy as np
import xlsxwriter
from multiprocessing import Pool
import bibtexparser
import argparse
import string
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# from pymatgen.core.structure import Structure
import pychemia
import spglib
from pymatgen.core import Element
import wyckoff_list
import string
import networkx as nx
import time

periodic = pychemia.utils.periodic
x_list = []
for x in string.ascii_uppercase:
    x_list.append(x)
for x in string.ascii_uppercase:
    for y in string.ascii_uppercase:
        x_list.append(x+y)


# def get_structure(args):
#     st_calc = args[0]['structure']
#     st_icsd = args[1]['structure']
#     symprec = args[2]
#     icsd_id =  str(args[0]['_id'])
#     print(f'{st_calc.formula:10}  |  {icsd_id:30}')
#     # print(st, st_icsd, symprec)
#     sigfigs = int(-np.log10(symprec))
#     try:
#         for i in range(st_icsd.natom):
#             st_icsd = shift(st_icsd, symprec, np.random.randint(st_icsd.natom), True)
#             st_icsd = sort_sites(st_icsd, symprec)
#             st_icsd = shift(st_icsd, symprec, np.random.randint(st_icsd.natom), True)
#             st_icsd = get_primitive(st_icsd, symprec)
#             st_icsd = shift(st_icsd, symprec, np.random.randint(st_icsd.natom), True)
#             st_icsd = get_conventional(st_icsd, symprec)
#             # time.sleep(0.2 + 1 * np.random.random())
#         st_icsd = shift(st_icsd, symprec)
#         st_icsd = get_standardized(st_icsd, symprec)
#         st_icsd = sort_sites(st_icsd, symprec)

#         st_calc = shift(st_calc, symprec)
#         st_calc = get_conventional(st_calc, symprec)
#         for i in range(st_calc.natom):
#             st_calc = shift(st_calc, symprec, np.random.randint(st_calc.natom), True)
#             st_calc = sort_sites(st_calc, symprec)
#             st_calc = shift(st_calc, symprec, np.random.randint(st_calc.natom), True)
#             st_calc = get_primitive(st_calc, symprec)
#             st_calc = shift(st_calc, symprec, np.random.randint(st_calc.natom), True)
#             st_calc = get_conventional(st_calc, symprec)
#         st_calc = shift(st_calc, symprec)
#         st_calc = get_standardized(st_calc, symprec)
#         st_calc = sort_sites(st_calc, symprec)
#     except:
#         st_calc = args[0]['structure']
#         st_icsd = args[1]['structure']
#         st_calc.symbols = ["Error" for x in st_calc.symbols]

#     return [{'structure':st_calc, '_id':args[0]['_id']}, {'structure':st_icsd, '_id':args[1]['_id']}]


def get_structure(args):
    st_calc = args[0]['structure']
    st_icsd = args[1]['structure']
    symprec = args[2]
    icsd_id =  str(args[0]['_id'])
    print(f'{st_calc.formula:10}  |  {icsd_id:30}')
    # print(st, st_icsd, symprec)
    sigfigs = int(-np.log10(symprec))
    
    try:
        if (st_calc.natom != st_icsd.natom):
            st_calc = rattle(st_calc)
            st_icsd = rattle(st_icsd)
        st_calc, st_icsd = match_structures(st_calc, st_icsd, count=0)
    
    except:
        st_calc = args[0]['structure']
        st_icsd = args[1]['structure']
        st_calc.symbols = ["Error" for x in st_calc.symbols]

    return [{'structure':st_calc, '_id':args[0]['_id']}, {'structure':st_icsd, '_id':args[1]['_id']}]

def parallel_query(db_calc, db_icsd,
                   symprec=1e-5, nproc=2,
                   query={'status.relaxed': True}):
    to_standardize = []
    print("preparing data for parallel")
    sigfigs = int(-np.log10(symprec))
    for i, entry_calc in enumerate(db_calc.entries.find(query)):
        mp_id = entry_calc['properties']['mp_id']
        entry_icsd = db_icsd.entries.find_one({'$and': [{"properties.mp_id":mp_id},
                                                        {"properties.standardized":True},
                                                        {"properties.theoretical":False}
                                                        ]})
        if entry_icsd:
            _id_icsd = entry_icsd['_id']
            _id = entry_calc['_id']
            st_icsd = db_icsd.get_structure(entry_icsd['_id'])
            st = db_calc.get_structure(entry_calc['_id'])
            to_standardize.append([{'structure':st, '_id':_id}, {'structure':st_icsd, '_id':_id_icsd}, symprec])
    with Pool(nproc) as p:
        results = p.map(get_structure, to_standardize)
    return results


def to_excel(db_calc, db_icsd,
             filename="xc_db_analysis.xlsx",
             symprec=1e-5,
             nproc=1,
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
              'a ICSD', 'b ICSD', 'c ICSD',
              'a error', 'b error','c error',
              'Ea', 'Eb', 'Ec',
              'Isodirectional error',
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
               'a', 'b', 'c',
               'a ICSD', 'b ICSD', 'c ICSD',
               'a error', 'b error','c error',
               'Isodirectional error',
               'Ea', 'Eb', 'Ec',
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
              'x', 'y', 'z',
              'x ICSD', 'y ICSD', 'z ICSD',
              'a ICSD', 'b ICSD', 'c ICSD',
              'Ex', 'Ey', 'Ez',
              'x error','y error', 'z error',
              '(RSME xyz).(DOF)', '(MAE xyz).(DOF)',
              'skip',
              'Same SPG',
              "Chemical formula",
              'MP ID',"ICSD code",
              'Threshold',
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

    # worksheet.write_rich_string(0,
    #                             col_number_l[],
    #                             'j = k',
    #                             superscript, '(n-1)',
    #                             center)
    
    col_number_ls = {}
    for i, title in enumerate(head_ls):
        worksheet_ls.write_string(0, i, title, header_format)
        col_number_ls[title] = i

    col_number_w = {}
    for i, title in enumerate(head_w):
        worksheet_w.write_string(0, i, title, header_format)
        col_number_w[title] = i

    threshold = 0.3
    worksheet_w.write_number(1, i, threshold, float_format_nf)
    cn_thr = "$"+x_list[i]+"$"+str(2)
        
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

    entries = parallel_query(db_calc, db_icsd,
                             symprec=symprec, nproc=nproc,
                             query=query)
    wf = open("errors.txt", 'w')
    for ientry in entries:
        _id_calc = ientry[0]["_id"]
        _id_icsd = ientry[1]["_id"]
        st = ientry[0]['structure']
        st_icsd = ientry[1]['structure']
        entry_calc = db_calc.entries.find_one({'_id':_id_calc})
        entry_icsd = db_icsd.entries.find_one({'_id':_id_icsd})
        if "Error" in st.symbols:
            print("++++++++++++++++++++++++++")
            print("++++++++++++++++++++++++++")
            print("++++++++++++++++++++++++++")
            print(st.formula, _id_calc)
            print("++++++++++++++++++++++++++")
            print("++++++++++++++++++++++++++")
            print("++++++++++++++++++++++++++")
            wf.write(st.formula+str(_id_calc)+"\n")
            continue
        mp_id = entry_calc['properties']['mp_id']
        xc = entry_calc['properties']['xc']
        gap_spin_channel = min([x for x in entry_calc['properties']['band_gap']])
        bandgap = entry_calc['properties']['band_gap'][gap_spin_channel]['gap']
        direct_gap = entry_calc['properties']['band_gap'][gap_spin_channel]['direct']
        code_icsd = entry_icsd['properties']['code_ICSD']
        formula = entry_icsd['properties']['cif'][f'{code_icsd}-ICSD']["_chemical_formula_structural"]
        name = entry_icsd['properties']['cif'][f'{code_icsd}-ICSD']["_chemical_name_common"]
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
        print(f'{mp_id:10}  |  {code_icsd:10}')
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

        cn_a_err = x_list[col_number_l['a error']]+str(row_l+1)

        cn_a = x_list[col_number_l['a']]+str(row_l+1)
        cn_a_icsd = x_list[col_number_l['a ICSD']]+str(row_l+1)
        worksheet_l.write_formula(cn_a_err, f'={cn_a}-{cn_a_icsd}', float_format_nf)

        cn_b_err = x_list[col_number_l['b error']]+str(row_l+1)
        cn_b = x_list[col_number_l['b']]+str(row_l+1)
        cn_b_icsd = x_list[col_number_l['b ICSD']]+str(row_l+1)
        worksheet_l.write_formula(cn_b_err, f'={cn_b}-{cn_b_icsd}', float_format_nf)

        cn_c_err = x_list[col_number_l['c error']]+str(row_l+1)
        cn_c = x_list[col_number_l['c']]+str(row_l+1)
        cn_c_icsd = x_list[col_number_l['c ICSD']]+str(row_l+1)
        worksheet_l.write_formula(cn_c_err, f'={cn_c}-{cn_c_icsd}', float_format_nf)

        cn_isodir = x_list[col_number_l['Isodirectional error']]+str(row_l+1)
        worksheet_l.write_formula(cn_isodir, f'=IF(ABS(SUM(SIGN({cn_a_err})+SIGN({cn_b_err})+SIGN({cn_c_err})))=3, TRUE, FALSE)', float_format_nf)
        
        # relative errors
        cn_Ea = x_list[col_number_l['Ea']]+str(row_l+1)
        worksheet_l.write_formula(cn_Ea, f'=({cn_a}-{cn_a_icsd})/{cn_a_icsd}', float_format_nf)

        cn_Eb = x_list[col_number_l['Eb']]+str(row_l+1)
        worksheet_l.write_formula(cn_Eb, f'=({cn_b}-{cn_b_icsd})/{cn_b_icsd}', float_format_nf)

        cn_Ec = x_list[col_number_l['Ec']]+str(row_l+1)
        worksheet_l.write_formula(cn_Ec, f'=({cn_c}-{cn_c_icsd})/{cn_c_icsd}', float_format_nf)

        if sum(outer) != 0:
            cn_out_a = x_list[col_number_l['a DOF']]+str(row_l+1)
            cn_out_b = x_list[col_number_l['b DOF']]+str(row_l+1)
            cn_out_c = x_list[col_number_l['c DOF']]+str(row_l+1)

            cn_rsme_abc = x_list[col_number_l['(RSME abc).(DOF)']]+str(row_l+1)
            worksheet_l.write(cn_rsme_abc, f'=(({cn_a_err}*{cn_out_a})^2+({cn_b_err}*{cn_out_b})^2+({cn_c_err}*{cn_out_c})^2)/SUM({cn_out_a}, {cn_out_b}, {cn_out_c})',  float_format_nf)

            cn_mae_abc = x_list[col_number_l['(MAE abc).(DOF)']]+str(row_l+1)
            worksheet_l.write(cn_mae_abc, f'=(ABS({cn_a_err}*{cn_out_a})+ABS({cn_b_err}*{cn_out_b})+ABS({cn_c_err}*{cn_out_c}))/SUM({cn_out_a}, {cn_out_b}, {cn_out_c})',  float_format_nf)
            
            cn_me_abc = x_list[col_number_l['(ME abc).(DOF)']]+str(row_l+1)
            worksheet_l.write(cn_me_abc, f'=({cn_a_err}*{cn_out_a}+{cn_b_err}*{cn_out_b}+{cn_c_err}*{cn_out_c})/SUM({cn_out_a}, {cn_out_b}, {cn_out_c})',  float_format_nf)
            
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

            cn_a_err = x_list[col_number_ls['a error']]+str(row_ls+1)
            cn_a = x_list[col_number_ls['a']]+str(row_ls+1)
            cn_a_icsd = x_list[col_number_ls['a ICSD']]+str(row_ls+1)
            worksheet_ls.write_formula(cn_a_err, f'={cn_a}-{cn_a_icsd}', float_format_nf)
            
            cn_b_err = x_list[col_number_ls['b error']]+str(row_ls+1)
            cn_b = x_list[col_number_ls['b']]+str(row_ls+1)
            cn_b_icsd = x_list[col_number_ls['b ICSD']]+str(row_ls+1)
            worksheet_ls.write_formula(cn_b_err, f'={cn_b}-{cn_b_icsd}', float_format_nf)

            cn_c_err = x_list[col_number_ls['c error']]+str(row_ls+1)
            cn_c = x_list[col_number_ls['c']]+str(row_ls+1)
            cn_c_icsd = x_list[col_number_ls['c ICSD']]+str(row_ls+1)
            worksheet_ls.write_formula(cn_c_err, f'={cn_c}-{cn_c_icsd}', float_format_nf)

            cn_isodir = x_list[col_number_ls['Isodirectional error']]+str(row_ls+1)
            worksheet_ls.write_formula(cn_isodir, f'=IF(ABS(SUM(SIGN({cn_a_err})+SIGN({cn_b_err})+SIGN({cn_c_err})))=3, TRUE, FALSE)', float_format_nf)

            
            # relative errors
            cn_Ea = x_list[col_number_ls['Ea']]+str(row_ls+1)
            worksheet_ls.write_formula(cn_Ea, f'=({cn_a}-{cn_a_icsd})/{cn_a_icsd}', float_format_nf)

            cn_Eb = x_list[col_number_ls['Eb']]+str(row_ls+1)
            worksheet_ls.write_formula(cn_Eb, f'=({cn_b}-{cn_b_icsd})/{cn_b_icsd}', float_format_nf)

            cn_Ec = x_list[col_number_ls['Ec']]+str(row_ls+1)
            worksheet_ls.write_formula(cn_Ec, f'=({cn_c}-{cn_c_icsd})/{cn_c_icsd}', float_format_nf)

            if sum(outer) != 0:
                cn_out_a = x_list[col_number_ls['a DOF']]+str(row_ls+1)
                cn_out_b = x_list[col_number_ls['b DOF']]+str(row_ls+1)
                cn_out_c = x_list[col_number_ls['c DOF']]+str(row_ls+1)

                cn_rsme_abc = x_list[col_number_ls['(RSME abc).(DOF)']]+str(row_ls+1)
                worksheet_ls.write(cn_rsme_abc, f'=SQRT((({cn_a_err}*{cn_out_a})^2+({cn_b_err}*{cn_out_b})^2+({cn_c_err}*{cn_out_c})^2)/SUM({cn_out_a}, {cn_out_b}, {cn_out_c}))',  float_format_nf)

                cn_mae_abc = x_list[col_number_ls['(MAE abc).(DOF)']]+str(row_ls+1)
                worksheet_ls.write(cn_mae_abc, f'=(ABS({cn_a_err}*{cn_out_a})+ABS({cn_b_err}*{cn_out_b})+ABS({cn_c_err}*{cn_out_c}))/SUM({cn_out_a}, {cn_out_b}, {cn_out_c})',  float_format_nf)

                cn_me_abc = x_list[col_number_ls['(ME abc).(DOF)']]+str(row_ls+1)
                worksheet_ls.write(cn_me_abc, f'=({cn_a_err}*{cn_out_a}+{cn_b_err}*{cn_out_b}+{cn_c_err}*{cn_out_c})/SUM({cn_out_a}, {cn_out_b}, {cn_out_c})',  float_format_nf)

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
            worksheet_w.write_string(row_w, col_number_w["Chemical formula"], formula, center_format)
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

            worksheet_w.write_number(row_w, col_number_w['a ICSD'], a_icsd, float_format)
            worksheet_w.write_number(row_w, col_number_w['b ICSD'], b_icsd, float_format)
            worksheet_w.write_number(row_w, col_number_w['c ICSD'], c_icsd, float_format)

            worksheet_w.write_boolean(row_w, col_number_w['Same SPG'], spg_no == spg_no_icsd, center_format)


            cn_Ex = x_list[col_number_w['Ex']]+str(row_w+1)
            cn_x = x_list[col_number_w['x']]+str(row_w+1)
            cn_x_icsd = x_list[col_number_w['x ICSD']]+str(row_w+1)
            worksheet_w.write_formula(cn_Ex, f'=MIN(ABS({cn_x}-{cn_x_icsd}), ABS({cn_x}-1-{cn_x_icsd}), ABS({cn_x}+1-{cn_x_icsd}))', float_format)


            cn_Ey = x_list[col_number_w['Ey']]+str(row_w+1)
            cn_y = x_list[col_number_w['y']]+str(row_w+1)
            cn_y_icsd = x_list[col_number_w['y ICSD']]+str(row_w+1)
            worksheet_w.write_formula(cn_Ey, f'=MIN(ABS({cn_y}-{cn_y_icsd}), ABS({cn_y}-1-{cn_y_icsd}), ABS({cn_y}+1-{cn_y_icsd}))', float_format)


            cn_Ez = x_list[col_number_w['Ez']]+str(row_w+1)
            cn_z = x_list[col_number_w['z']]+str(row_w+1)
            cn_z_icsd = x_list[col_number_w['z ICSD']]+str(row_w+1)
            worksheet_w.write_formula(cn_Ez, f'=MIN(ABS({cn_z}-{cn_z_icsd}), ABS({cn_z}-1-{cn_z_icsd}), ABS({cn_z}+1-{cn_z_icsd}))', float_format)

            cn_x_err = x_list[col_number_w['x error']]+str(row_w+1)
            cn_a_icsd = x_list[col_number_w['a ICSD']]+str(row_w+1)
            worksheet_w.write_formula(cn_x_err, f'=MIN(ABS({cn_x}-{cn_x_icsd}), ABS({cn_x}-1-{cn_x_icsd}), ABS({cn_x}+1-{cn_x_icsd}))*{cn_a_icsd}', float_format)
            
            cn_y_err = x_list[col_number_w['y error']]+str(row_w+1)
            cn_b_icsd = x_list[col_number_w['b ICSD']]+str(row_w+1)
            worksheet_w.write_formula(cn_y_err, f'=MIN(ABS({cn_y}-{cn_y_icsd}), ABS({cn_y}-1-{cn_y_icsd}), ABS({cn_y}+1-{cn_y_icsd}))*{cn_b_icsd}', float_format)

            cn_z_err = x_list[col_number_w['z error']]+str(row_w+1)
            cn_c_icsd = x_list[col_number_w['c ICSD']]+str(row_w+1)
            worksheet_w.write_formula(cn_z_err, f'=MIN(ABS({cn_z}-{cn_z_icsd}), ABS({cn_z}-1-{cn_z_icsd}), ABS({cn_z}+1-{cn_z_icsd}))*{cn_c_icsd}', float_format)

            if sum(inner[iatom]) != 0:
                cn_in_x = x_list[col_number_w['x DOF']]+str(row_w+1)
                cn_in_y = x_list[col_number_w['y DOF']]+str(row_w+1)
                cn_in_z = x_list[col_number_w['z DOF']]+str(row_w+1)

                cn_rsme_xyz = x_list[col_number_w['(RSME xyz).(DOF)']]+str(row_w+1)
                worksheet_w.write(cn_rsme_xyz, f'=SQRT((({cn_x_err}*{cn_in_x})^2+({cn_y_err}*{cn_in_y})^2+({cn_z_err}*{cn_in_z})^2)/SUM({cn_in_x}, {cn_in_y}, {cn_in_z}))',  float_format)

                cn_mae_xyz = x_list[col_number_w['(MAE xyz).(DOF)']]+str(row_w+1)
                worksheet_w.write(cn_mae_xyz, f'=(ABS({cn_x_err}*{cn_in_x})+ABS({cn_y_err}*{cn_in_y})+ABS({cn_z_err}*{cn_in_z}))/SUM({cn_in_x}, {cn_in_y}, {cn_in_z})',  float_format)

                cn_skip = x_list[col_number_w['skip']]+str(row_w+1)
                worksheet_w.write_formula(cn_skip, f'=OR(IF({cn_x_err}*{cn_in_x}>{cn_thr}, TRUE, FALSE), IF({cn_y_err}*{cn_in_y}>{cn_thr}, TRUE, FALSE), IF({cn_z_err}*{cn_in_z}>{cn_thr}, TRUE, FALSE))', center_format)
                

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
    cn_x_err_start = x_list[col_number_w['x error']]+str(2)
    cn_z_err_end = x_list[col_number_w['z error']]+str(row_w)


    float_format_red = workbook.add_format({'num_format': '#,##0.0', 'align': 'center', 'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
    float_format_red.set_border()

    worksheet_w.conditional_format(f'{cn_x_err_start}:{cn_z_err_end}', {'type':     'cell',
                                                                        'criteria': '>=',
                                                                        'value':    cn_thr,
                                                                        'format':   float_format_red})
    workbook.close()
    wf.close()
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



def match_structures(st1, st2, symprec=1e-2, count=0):
    sigfigs = int(-np.log10(symprec))
    st1.reduced[st1.reduced.round(sigfigs) >= 1] -= 1
    st2.reduced[st2.reduced.round(sigfigs) < 0] += 1
    if count>100:
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print(f"reached max recurssion {count}")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        
        return st1, st2
    
    G = nx.Graph()
    mapper = {x:[] for x,_ in enumerate(st1.reduced)}
    for ipos, fc1 in enumerate(st1.reduced):
        G.add_node(f"{ipos}-s1", coords=fc1)
    for jpos, fc2 in enumerate(st2.reduced):
        G.add_node(f"{jpos}-s2", coords=fc2)
    

    for ipos, fc1 in enumerate(st1.reduced):
        for jpos, fc2 in enumerate(st2.reduced):
            diff = np.array([min(abs(x1-x2), abs(x1-1-x2), abs(x1+1-x2)) for x1, x2 in zip(fc1, fc2)])
            weight = (diff.round() <= [0.1, 0.1, 0.1]).sum()
            if weight==3 :
                G.add_edge(f"{ipos}-s1", f"{jpos}-s2", weight=1/max(1e-5,np.linalg.norm(diff)))
    for inode in G.nodes():
        edges = np.array([e for e in G.edges(inode, data='weight')])
        if len(edges) in [0,1]:
            continue
        
        weights = edges[:, 2].astype(float)
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
    mapper = {}
    perfect_match = True
    
    for inode in G.nodes():
        edges = np.array([e for e in G.edges(inode)])
        if len(edges) != 1:
            perfect_match = False
        
    if perfect_match:
        for e in G.edges():
            u = int(e[0].split("-")[0])
            v = int(e[1].split("-")[0])
            mapper[v]=u
        idx = np.zeros(st1.natom).astype(int)
        for iatom in range(st1.natom):
            idx[iatom] = mapper[iatom]
        st1.sort_sites_using_list(idx)

    else:
        st1 = rattle(st1, symprec)
        st2 = rattle(st2, symprec)
        
        st1, st2 = match_structures(st1, st2, symprec, count+1)
        
    return st1, st2


def rattle(st, symprec=1e-2):
    st = get_primitive(st, symprec)
    st = shift(st, symprec, -1, True)
    st = get_conventional(st, symprec)
    st = get_standardized(st, symprec)
    st = shift(st, symprec, -1, True)
    st = sort_sites(st, symprec)
    st = shift(st, symprec, -1, True)
    st = get_primitive(st, symprec)
    st = shift(st, symprec, -1, True)
    st = get_conventional(st, symprec)
    st = shift(st, symprec, -1, True)
    st = sort_sites(st, symprec)
    st = shift(st, symprec, -1, True)
    st = get_primitive(st, symprec)
    st = sort_sites(st, symprec)
    st = shift(st, symprec, -1, True)
    st = get_conventional(st, symprec)
    # st = get_primitive(st, symprec)
    # st = get_standardized(st, symprec)
    st = sort_sites(st, symprec, False)
    # st = st_sort_axes(st)
    return st

                    
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
    return (min(abs(x-x_prime), abs(x-1-x_prime), abs(x+1-x_prime)))
 
if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db_icsd",
        dest="db_icsd_name",
        type=str,
        help='Name of the ICSD MongoDB database')
    parser.add_argument(
        "--compare_db",
        dest="db_calc_name",
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
        "--output","-o",
        dest='output',
        type=str,
        help="Output file name",
        default="output.xlsx")
    parser.add_argument(
        "--nprocessors","-np",
        dest="nproc",
        type=int,
        help="Number of processors to be used",
        default=1)
    parser.add_argument(
        "--verbose","-v",
        action=argparse.BooleanOptionalAction,
        help='To print the results of the on going analysis.',
        default="--verbose")
    args = parser.parse_args()
    db_icsd = pychemia.db.get_database(dict(name=args.db_icsd_name,
                                            host=args.host,
                                            port=args.port,
                                            user=args.user,
                                            passwd=args.passwd))
    db_calc = pychemia.db.get_database(dict(name=args.db_calc_name,
                                            host=args.host,
                                            port=args.port,
                                            user=args.user,
                                            passwd=args.passwd))
    to_excel(db_calc, db_icsd, args.output, symprec=1e-3, nproc=args.nproc)
    
