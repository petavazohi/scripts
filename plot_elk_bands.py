import argparse
import re
import numpy as np
import matplotlib.pyplot as plt



class data:
    def __init__(self,kvector,band,spd):
        self.kvector = kvector
        self.band    = band
        self.spd     = spd


parser = argparse.ArgumentParser()

parser.add_argument("--mode",dest="mode",type=str,help='Mode of bandsplot, plain, parametric',choices=['plain','parametric'],default='plain')
parser.add_argument('--orbitals',nargs='+',type=int,default=[0,1,2,3],help='orbitals you want to project')
parser.add_argument('--atoms',nargs='+',type=int,help='atoms you want to project')
parser.add_argument('--infile',type=str,help='Input file',default='BANDS.OUT')
parser.add_argument('--elimit',type=float,nargs=2)
parser.add_argument('--cmap',type=str,default='jet')

args = parser.parse_args()

rf = open('elk.in')
elkIn = rf.read()
rf.close()


nvortex,nkpoint = [int(x) for x in re.findall('plot1d\n\s*([0-9]*)\s*([0-9]*)',elkIn)[0]]
ticks = ["$%s$" % (x.replace(',','').replace('vlvp1d','').replace(' ','')) for x in re.findall('plot1d\n\s*[0-9]*\s*[0-9]*.*\n'+nvortex*'.*:(.*)\n',elkIn)[0]]
nspecies = int(re.findall('atoms\n\s*([0-9]*)',elkIn)[0])
species = {}
for ispc in re.findall('\'([A-Za-z]*.in)\'.*\n\s*([0-9]*)',elkIn):
    species[ispc[0]] =  int(ispc[1])


tick_pos = np.zeros(shape=(len(ticks,)))
rf = open('BANDLINES.OUT','r')
bandLines = rf.readlines()
rf.close()
itick = 0
for iline in range(0,len(bandLines),3):
    tick_pos[itick] = float(bandLines[iline].split()[0])
    itick += 1

if args.mode == 'plain' :
    rf = open(args.infile,'r')
    lines = rf.readlines()
    rf.close()
    nband = int(len(lines)/(nkpoint+1))
    bands = np.zeros(shape=(nband,nkpoint))
    kposition = np.zeros(shape=(nkpoint,))
    iline = 0        
    for iband in range(nband):
        for ikpoint in range(nkpoint):
            bands[iband,ikpoint] =  float(lines[iline].split()[1])
            iline += 1
            if ikpoint == nkpoint -1 :
                iline +=1
    bands *= 27.21138386
    iline = 0
    for iband in range(1):
        for ikpoint in range(nkpoint):
            kposition[ikpoint] =  float(lines[iline].split()[0])
            iline += 1
            if ikpoint == nkpoint -1 :
                iline +=1
                
    for iband in range(nband):
        plt.plot(kposition,bands[iband,:],'b')



if args.mode == 'parametric':
    from matplotlib.collections import LineCollection
    import matplotlib
    file_names = [] 
    ispc = 1 
    for spc in species: 
        for iatom in range(species[spc]): 
            file_names.append('BAND_S{:02d}_A{:04d}.OUT'.format(ispc,iatom+1))
        ispc += 1 
    file_names = np.array(file_names)
    files_to_read = file_names[args.atoms]
    
    
    rf = open(files_to_read[0],'r')
    lines = rf.readlines()
    rf.close()
    nband = int(len(lines)/(nkpoint+1))
    bands = np.zeros(shape=(nband,nkpoint))
    kposition = np.zeros(shape=(nkpoint,))
    iline = 0
    for iband in range(nband):
        for ikpoint in range(nkpoint):
            bands[iband,ikpoint] =  float(lines[iline].split()[1])
            iline += 1
            if ikpoint == nkpoint -1 :
                iline +=1
    bands *= 27.21138386
    iline = 0
    for iband in range(1):
        for ikpoint in range(nkpoint):
            kposition[ikpoint] =  float(lines[iline].split()[0])
            iline += 1
            if ikpoint == nkpoint -1 :
                iline +=1
    color = np.zeros(shape=(nband,nkpoint))
    orbitals = np.array(args.orbitals)
    for ifile in files_to_read :
        rf = open(ifile,'r')
        lines = rf.readlines()
        rf.close()
        iline = 0
        print(ifile)
        for iband in range(nband):
            for ikpoint in range(nkpoint):
                temp = np.array([float(x) for x in lines[iline].split()])
                color[iband,ikpoint] = color[iband,ikpoint]+temp[orbitals+3].sum()
                iline += 1
                if ikpoint == nkpoint -1 :
                    iline +=1

    
    norm = matplotlib.colors.Normalize(color.min(), color.max())
    fig = plt.figure()
    gca = fig.gca()

    for iband in range(nband):
        points = np.array([kposition,bands[iband,:]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=plt.get_cmap(args.cmap), norm=norm)
        lc.set_array(color[iband,:])
        lc.set_linewidth(0.5)
        gca.add_collection(lc)
    cb = plt.colorbar(lc)
        
                
if args.elimit != None :
    plt.ylim(args.elimit[0],args.elimit[1])
else :
    plt.ylim(bands.min()-0.1*abs(bands.min()),bands.max()+0.1*abs(bands.max()))
for itick in tick_pos :
    plt.axvline(itick,linewidth=2,color='k')
plt.xlim(kposition.min(),kposition.max())
plt.xticks(tick_pos,ticks)
plt.show()
