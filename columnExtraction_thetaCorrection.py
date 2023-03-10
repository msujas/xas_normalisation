import os
from glob import glob
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

direc = r'C:\Users\kenneth1a\Documents\beamlineData\ch6616/'
os.chdir(direc)
thetaOffset = 0.018991970


dspacing = 3.133789
planck = 6.62607015e-34
charge = 1.60217663e-19
speedOfLight = 299792458
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
        onscan = False
        dfFilteredDct = {}
        for c,line in enumerate(lines):
            if '#S' in line and 'zapline' in line:
                newstring = ''
                newstring += line
                spectrum_count += 1
                onscan = True


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
                f.close()
                df = pd.read_csv(file, skiprows = dfstart, nrows = dfend - dfstart, delim_whitespace = True, header = None)
                df.columns = columns
                if 'xmap_roi00' in columns:
                    filtCols = ['ZapEnergy','TwoTheta','mon_3','mon_4','ion_1_2','ion_1_3','xmap_roi00']
                else:
                    filtCols = ['ZapEnergy','TwoTheta','mon_3','mon_4','ion_1_2','ion_1_3']
                dfFilteredDct[spectrum_count] = df[filtCols].copy(deep=True)
                #dfFilteredDct[spectrum_count].set_index('ZapEnergy',inplace = True)
                dfFilteredDct[spectrum_count]['Theta_offset'] =  dfFilteredDct[spectrum_count]['TwoTheta'].apply(lambda x: np.round(x + thetaOffset,7))
                dfFilteredDct[spectrum_count]['ZapEnergy_offset'] = dfFilteredDct[spectrum_count]['Theta_offset'].apply(angle_to_kev)
                dfFilteredDct[spectrum_count].set_index('ZapEnergy_offset',inplace = True)

                newfile = f'{newdir}{basename}_{spectrum_count:02d}.dat'
                newfilerg = f'{newdir}regrid/{basename}_{spectrum_count:02d}.dat'

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
            for n in range(5):

                try:
                    grid = np.arange(ZE[ZEmin].round(4),ZE[ZEmax].round(5),spacing)
                    grid = grid.round(5)
                    gridfuncMon = interp1d(dfFilteredDct[c].index.values,dfFilteredDct[c]['mon_4'].values)
                    regridDF['mon_4'] = gridfuncMon(grid).round(1)
                    gridfuncMon3 = interp1d(dfFilteredDct[c].index.values,dfFilteredDct[c]['mon_3'].values)
                    regridDF['mon_3'] = gridfuncMon3(grid).round(1)
                    grindfuncI12 = interp1d(dfFilteredDct[c].index.values,dfFilteredDct[c]['ion_1_2'].values)
                    regridDF['ion_1_2'] = grindfuncI12(grid).round(1)
                    grindfuncI13 = interp1d(dfFilteredDct[c].index.values,dfFilteredDct[c]['ion_1_3'].values)
                    regridDF['ion_1_3'] = grindfuncI13(grid).round(1)
                    if 'xmap_roi00' in dfFilteredDct[0].columns:
                        grindfuncxmap = interp1d(dfFilteredDct[c].index.values,dfFilteredDct[c]['xmap_roi00'].values)
                        regridDF['xmap_roi00'] = grindfuncxmap(grid).round(1)
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
                    if n == 4:
                        print(len(ZE))
                        print(len(dfFilteredDct[c].index.values))
                        print(f'{file} spectrum {c} couldn\'t be regridded')
                        print(str(e))
                        os.remove(newfilerg)
            if len(regridDF.columns) == 0:
                continue
            regridDF.to_csv(newfilerg,sep = ' ',mode = 'a')
