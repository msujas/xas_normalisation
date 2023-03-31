import os
from glob import glob
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

direc = r'Z:\visitor\ch6617\bm31\20230329\CMS029/'
os.chdir(direc)
thetaOffset = 0


dspacing = 3.133789
planck = 6.62607015e-34
charge = 1.60217663e-19
speedOfLight = 299792458

fluoCounter = 'xmap_roi00'
counterNames = ['ZapEnergy','TwoTheta','mon_3','mon_4','mon_1','ion_1_2','ion_1_3','ion_1_1',fluoCounter]
counterNames_NF = [c for c in counterNames if c != fluoCounter] #NF - no fluorescence

def angle_to_kev(angle): #NB the TwoTheta data in the .dat files is really theta
    wavelength = 2*dspacing*np.sin(angle*np.pi/(180))
    wavelength_m = wavelength*10**(-10)
    energy_kev = planck*speedOfLight/(wavelength_m*charge*1000)
    return np.round(energy_kev,6)

for root, dirs, files in os.walk(os.getcwd()):
    if 'columns' in root:
        continue
    print(root)
    currentdir = root + '/'

    os.chdir(currentdir)
    datfiles = glob('*.dat')
    datfiles = [currentdir + file for file in datfiles]

    if len(datfiles) == 0:
        continue

    print(os.getcwd())
    for file in datfiles:
        basename = os.path.splitext(os.path.basename(file))[0]
        if not os.path.exists(currentdir+'columns/'):
            os.makedirs(currentdir+'columns/')
        newdir = f'columns/{basename}/'
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        if not os.path.exists(f'{newdir}regrid/'):
            os.makedirs(f'{newdir}regrid/')
        spectrum_count = -1


        f = open(file,'r')
        lines = f.readlines()
        f.close()
        onscan = False
        dfFilteredDct = {}
        for c,line in enumerate(lines):
            if '#S' in line and 'zapline' in line:
                newstring = ''
                newstring += line
                spectrum_count += 1
                onscan = True

            elif '#T' in line and onscan:
                timeStep = int(line.split()[1])/1000
            elif '#D' in line and onscan:
                newstring += line
            elif '#L' in line:
                dfstart = c +1
                columns = line.replace('#L ','').split()
            elif '#C' in line and onscan:
                dfend = c
                onscan = False
                if dfend-dfstart <= 1:
                    continue
                
                df = pd.read_csv(file, skiprows = dfstart, nrows = dfend - dfstart, delim_whitespace = True, header = None)
                df.columns = columns
                if fluoCounter in columns:
                    filtCols = counterNames
                else:
                    filtCols = counterNames_NF
                dfFilteredDct[spectrum_count] = df[filtCols].copy(deep=True)
                #dfFilteredDct[spectrum_count].set_index('ZapEnergy',inplace = True)
                dfFilteredDct[spectrum_count]['Theta_offset'] =  dfFilteredDct[spectrum_count]['TwoTheta'].apply(lambda x: np.round(x + thetaOffset,7))
                dfFilteredDct[spectrum_count]['ZapEnergy_offset'] = dfFilteredDct[spectrum_count]['Theta_offset'].apply(angle_to_kev)
                dfFilteredDct[spectrum_count].set_index('ZapEnergy_offset',inplace = True)

                for counter in filtCols: #removing unused counters
                    if 'mon_' in counter or 'ion_1' in counter:
                        if np.max(dfFilteredDct[spectrum_count][counter].values) < 1000*timeStep:
                            dfFilteredDct[spectrum_count].drop(counter,axis = 1,inplace = True)
    

                newfile = f'{newdir}/{basename}_{spectrum_count:02d}.dat'
                newfilerg = f'{newdir}/regrid/{basename}_{spectrum_count:02d}.dat'
                if len([col for col in dfFilteredDct[spectrum_count].columns if 'mon_' in col]) == 0:
                    if os.path.exists(newfilerg):
                        os.remove(newfilerg)
                    if os.path.exists(newfile):
                        os.remove(newfile)
                    continue
                f = open(newfile,'w')
                f.write(newstring)
                f.close()
                dfFilteredDct[spectrum_count].to_csv(newfile,sep = ' ',mode = 'a')
                f2 = open(newfilerg,'w')
                f2.write(newstring)
                f2.close()
        if len(dfFilteredDct) == 0:
            continue
        ZElens = [len(dfFilteredDct[c].index.values) for c in dfFilteredDct]

        ZEindex = ZElens.index(max(ZElens))
        ZEkey = list(dfFilteredDct.keys())[ZEindex]
        ZE = dfFilteredDct[ZEkey].index.values


        diffs = np.array([ZE[i]-ZE[i-1] for i in range(1,len(ZE))])
        spacing = np.mean(diffs).round(4)


        ZEmin = 0
        ZEmax = -1
        #os.chdir(f'{newdir}regrid/')

        for c in dfFilteredDct:

            newfilerg = f'{newdir}regrid/{basename}_{c:02d}.dat'
            regridDF = pd.DataFrame()
            if len([col for col in dfFilteredDct[c].columns if 'mon_' in col]) == 0:
                continue
            for n in range(5):
                
                try:
                    grid = np.arange(ZE[ZEmin].round(4),ZE[ZEmax].round(5),spacing)
                    grid = grid.round(5)
                    if len(dfFilteredDct[c].index.values) < len(grid) - 5:
                        print(f'{file} spectrum {c} too short, couldn\'t be regridded')
                        if os.path.exists(newfilerg):
                            os.remove(newfilerg)
                        break
                    for counter in dfFilteredDct[c].columns:
                        if 'mon' in counter or 'ion_1' in counter or fluoCounter in counter:
                            gridfunc = interp1d(dfFilteredDct[c].index.values,dfFilteredDct[c][counter].values)
                            regridDF[counter] = gridfunc(grid).round(1)

                    
                    regridDF.index = grid
                    regridDF.index.name = 'energy_offset(keV)'
                    if n != 0:
                        print(f'ZEmin = {ZEmin}, ZEmax {ZEmax}')
                    break
                except ValueError as e:
                    if 'below' in str(e):
                        ZEmin +=1
                    elif 'above' in str(e):
                        ZEmax -= 1
            if len(regridDF.columns) == 0:
                continue
            regridDF.to_csv(newfilerg,sep = ' ',mode = 'a')
