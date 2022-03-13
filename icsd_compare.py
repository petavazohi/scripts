#!/usr/bin/env python3
import pychemia
from collections import OrderedDict
import pandas as pd

    def to_excel(self,
                 filename="xc_db_analysis.xlsx",
                 query={'status.relaxed': True}):
        ret = OrderedDict({'Material ID': [],
                           'Exchange correlation':[],
                           'Space group name': [],
                           'Space group number': [],
                           'Family': [],
                           'No. lattice degrees of freedom': [],
                           'Bandgap': [], 'Direct gap': [],
                           'Atom symbol': [],
                           'No. Wyckoff degrees of freedom': [],
                           'No. species': [], 'No. atoms': [],
                           'a freedom': [], "b freedom": [], "c freedom": [],
                           'x freedom': [], "y freedom": [], "z freedom": [],
                           'Multiplicity': [],
                           'Wyckoff letter': [],
                           'a': [], 'b': [], 'c': [],
                           'x': [], 'y': [], 'z': [],
                           'alpha': [], 'beta': [], 'gamma': [],
                           'a ICSD': [], 'b ICSD': [], 'c ICSD': [],
                           'x ICSD': [], 'y ICSD': [], 'z ICSD': [],
                           'alpha ICSD': [], 'beta ICSD': [], 'gamma ICSD': [], 'Same spg': [],
                           })
        for entry in self.entries.find(query):
            _id = str(entry['_id'])
            print(_id)
            natom = entry['properties']['final']['structure']['natom']
            nspecies = entry['properties']['final']['structure']['nspecies']
            reduced = entry['properties']['final']['structure']['reduced']
            reduced_initial = entry['properties']['initial']['structure']['reduced']
            symbols = entry['properties']['final']['structure']['symbols']
            inner = entry['properties']['final']['wyckoff_analysis']['inner_degrees_of_freedom']
            outer =  entry['properties']['final']['lattice_degrees_of_freedom']
            # gap = entry['properties']['band_gap']['total']['gap']
            gap_spin_channel = min([x for x in entry['properties']['band_gap']])
            a_fin = entry['properties']['final']['lattice']['a']
            b_fin = entry['properties']['final']['lattice']['b']
            c_fin = entry['properties']['final']['lattice']['c']
            a_init = entry['properties']['initial']['lattice']['a']
            b_init = entry['properties']['initial']['lattice']['b']
            c_init = entry['properties']['initial']['lattice']['c']
            for iatom in range(natom):
                ret['Material ID'].append(entry['properties']['mp_id'])
                ret['Exchange correlation'].append(entry['properties']['xc'])
                ret['Space group number'].append(entry['properties']['final']['crystal_symmetry']['number'])
                ret['Space group name'].append(entry['properties']['final']['crystal_symmetry']['international_short'])
                ret['Family'].append(entry['properties']['final']['crystal_family'])
                if iatom == 0:
                    ret['No. lattice degrees of freedom'].append(sum(outer))
                else:
                    ret['No. lattice degrees of freedom'].append(0)
                ret['Bandgap'].append(entry['properties']['band_gap'][gap_spin_channel]['gap'])
                ret['Direct gap'].append(entry['properties']['band_gap'][gap_spin_channel]['direct'])
                ret['Atom symbol'].append(symbols[iatom])
                ret['No. Wyckoff degrees of freedom'].append(sum(inner[iatom]))
                ret['No. species'].append(nspecies)
                ret['No. atoms'].append(natom)
                
                ret['a freedom'].append(outer[0])
                ret["b freedom"].append(outer[1])
                ret["c freedom"].append(outer[2])
                
                ret['x freedom'].append(inner[iatom][0])
                ret["y freedom"].append(inner[iatom][1])
                ret["z freedom"].append(inner[iatom][2])
                ret['Multiplicity'].append(entry['properties']['final']['wyckoff_analysis']['multiplicity'][iatom])
                ret['Wyckoff letter'].append(entry['properties']['final']['wyckoff_analysis']['wyckoff_letter'][iatom])
                ret['a'].append(entry['properties']['final']['lattice']['a'])
                ret['b'].append(entry['properties']['final']['lattice']['b'])
                ret['c'].append(entry['properties']['final']['lattice']['c'])
                ret['x'].append(reduced[iatom][0])
                ret['y'].append(reduced[iatom][1])
                ret['z'].append(reduced[iatom][2])
                ret['alpha'].append(entry['properties']['final']['lattice']['alpha'])
                ret['beta'].append(entry['properties']['final']['lattice']['beta'])
                ret['gamma'].append(entry['properties']['final']['lattice']['gamma'])
                ret['a ICSD'].append(entry['properties']['initial']['lattice']['a'])
                ret['b ICSD'].append(entry['properties']['initial']['lattice']['b'])
                ret['c ICSD'].append(entry['properties']['initial']['lattice']['c'])
                ret['x ICSD'].append(reduced_initial[iatom][0])
                ret['y ICSD'].append(reduced_initial[iatom][1])
                ret['z ICSD'].append(reduced_initial[iatom][2])
                ret['alpha ICSD'].append(entry['properties']['initial']['lattice']['alpha'])
                ret['beta ICSD'].append(entry['properties']['initial']['lattice']['beta'])
                ret['gamma ICSD'].append(entry['properties']['initial']['lattice']['gamma'])
                ret['Same spg'].append(entry['properties']['initial']['crystal_symmetry']['number'] == entry['properties']['final']['crystal_symmetry']['number'])
        df = pd.DataFrame.from_dict(ret)
        
        df['x error'] = [diff_reduced_coordinate(x, x_prime) for x, x_prime in zip(df['x'], df['x ICSD'])]
        df['y error'] = [diff_reduced_coordinate(x, x_prime) for x, x_prime in zip(df['y'], df['y ICSD'])]
        df['z error'] = [diff_reduced_coordinate(x, x_prime) for x, x_prime in zip(df['z'], df['z ICSD'])]
        df['a error'] = (df['a'] - df['a ICSD'])/df['a']
        df['b error'] = (df['b'] - df['b ICSD'])/df['b']
        df['c error'] = (df['c'] - df['c ICSD'])/df['c']
        df['(RSME xyz).(DOF)'] = (((df['x error']*df['x freedom'])**2+
                                   (df['y error']*df['y freedom'])**2+
                                   (df['z error']*df['z freedom'])**2)/(df['No. Wyckoff degrees of freedom']))**(1/2)
        df['(MAE xyz).(DOF)'] = (df['x error']*df['x freedom']+
                                 df['y error']*df['y freedom']+
                                 df['z error']*df['z freedom'])/(df['No. Wyckoff degrees of freedom'])

        df['(RSME abc).(DOF)'] = (((df['a error']*df['a freedom'])**2+
                                   (df['b error']*df['b freedom'])**2+
                                   (df['c error']*df['c freedom'])**2)/(df['No. lattice degrees of freedom']))**(1/2)
        df['(MAE abc).(DOF)'] = (abs(df['a error'])*df['a freedom']+
                                 abs(df['b error'])*df['b freedom']+
                                 abs(df['c error'])*df['c freedom'])/(df['No. lattice degrees of freedom'])
        df['RSME abc'] = (((df['a error'])**2+
                           (df['b error'])**2+
                           (df['c error'])**2)/3)**(1/2)*df['No. lattice degrees of freedom']/df['No. lattice degrees of freedom']
        df['ME abc'] = (abs(df['a error'])+
                        abs(df['b error'])+
                        abs(df['c error']))/3*df['No. lattice degrees of freedom']/df['No. lattice degrees of freedom']
        df.to_excel(f'{self.name}-db.xlsx')
        return



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
    
