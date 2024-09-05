import os, re
from glob import glob
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


direc = r'C:\Users\kenneth1a\Documents\beamlineData\May2024carlos'
thetaOffset = 0


dspacing = 3.13379
planck = 6.62607015e-34
charge = 1.60217663e-19
speedOfLight = 299792458

digits = 4
fluoCounter = 'xmap_roi00'
monPattern = 'mon_'
ion1Pattern = 'ion_1'
counterNames = ['ZapEnergy','TwoTheta', 'mon_2','mon_3','mon_4','mon_1','ion_1_2','ion_1_3','ion_1_1', 'Det_1', 'Det_2', 'Det_3',fluoCounter]
counterNames_NF = [c for c in counterNames if c != fluoCounter] #NF - no fluorescence
xColumns = ['ZapEnergy','TwoTheta']
monCounters = ['mon_1', 'mon_2', 'mon_3', 'mon_4']
i1counters = ['ion_1_1', 'ion_1_2', 'ion_1_3', 'Det_1', 'Det_2', 'Det_3']
i2name = 'I2'

def angle_to_kev(angle): #NB the TwoTheta data in the .dat files is really theta
    wavelength = 2*dspacing*np.sin(angle*np.pi/(180))
    wavelength_m = wavelength*10**(-10)
    energy_kev = planck*speedOfLight/(wavelength_m*charge*1000)
    return np.round(energy_kev,6)

def processFile(file, fileDct, currentdir, thetaOffset, startSpectrum = 0):
    f = open(file,'r')
    data = f.read()
    f.close()
    filemtime = os.path.getmtime(file)
    if not 'zapline' in data:
        fileDct[file] = [filemtime,-1]
        return
    
    basename = os.path.splitext(os.path.basename(file))[0]
    if not os.path.exists(currentdir+'columns/'):
        os.makedirs(currentdir+'columns/')
    newdir = f'{currentdir}/columns/{basename}/'
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    if not os.path.exists(f'{newdir}/regrid/'):
        os.makedirs(f'{newdir}/regrid/')
    spectrum_count = -1


    f = open(file,'r')
    lines = f.readlines()
    f.close()
    print(file)
    scanStart = False
    onscan = False
    for c,line in enumerate(lines):
        if '#S' in line and 'zapline' in line:
            newstring = ''
            newstring += line
            spectrum_count += 1
            if spectrum_count < startSpectrum:
                continue
            onscan = True

        elif '#T' in line and onscan:
            timeStep = int(line.split()[1])/1000
        elif '#D' in line and onscan:
            newstring += line
        elif '#L' in line and onscan:
            dfstart = c +1
            columns = line.replace('#L ','').split()
            scanStart = True
            lineno = 0
        elif scanStart and '#' not in line:
            lineSplit = np.array([np.fromstring(line,sep = ' ')])
            if lineno == 0:
                array = lineSplit

                lineno += 1
            else:
                array = np.append(array,lineSplit,axis = 0)
            
        elif '#C' in line and onscan:
            dfend = c
            scanStart = False
            onscan = False
            if dfend-dfstart <= 1:
                continue
            df = pd.DataFrame(data=array,columns=columns)

            dfFiltered = df[xColumns].copy(deep=True)
            dfFiltered['Theta_offset'] =  dfFiltered['TwoTheta'].apply(lambda x: np.round(x + thetaOffset,7))
            dfFiltered['ZapEnergy_offset'] = dfFiltered['Theta_offset'].apply(angle_to_kev)
            dfFiltered.set_index('ZapEnergy_offset',inplace = True)

            usedMon = df[monCounters].max().idxmax()
            
            dfFiltered[usedMon] = df[usedMon].values
            newfile = f'{newdir}/{basename}_{spectrum_count:0{digits}d}.dat'
            if np.max(dfFiltered[usedMon].values) < timeStep*10000 or np.min(dfFiltered[usedMon].values) < 10*timeStep:
                if os.path.exists(newfile):
                    os.remove(newfile)
                continue
            usedI1s = [col for col in i1counters if np.max(df[col].values) > 10000*timeStep]
            if usedI1s:
                usedI1 = usedI1s[0] #df[i1counters].max().idxmax()
                dfFiltered[usedI1] = df[usedI1].values
            usedI2 = ""
            if len(usedI1s) > 1:
                usedI2 = usedI1s[1]
                dfFiltered[i2name] = df[usedI2].values
            if fluoCounter in columns:
                if df[fluoCounter].values.max() > 300*timeStep:
                    dfFiltered[fluoCounter] = df[fluoCounter].values
            if not usedI1s and fluoCounter not in dfFiltered.columns:
                if os.path.exists(newfile):
                    os.remove(newfile)
                continue
            
            f2 = open(newfile,'w')
            f2.write(newstring)
            f2.close()
            dfFiltered.to_csv(newfile,sep = ' ',mode = 'a')
            print(newfile)

    fileDct[file] = [filemtime,spectrum_count]
    

def regrid(coldir, unit = 'keV'):
    if not os.path.exists(coldir):
        return
    print(f'regridding {coldir}')
    print(coldir)
    files = glob(f'{coldir}/*.dat')
    files.sort()
    if len(files) == 0:
        return
    dfFilteredDct = {}
    headers = []
    if unit == 'keV':
        escale = 1
    elif unit == 'eV':
        escale = 1000
    else:
        escale = 1
        unit = 'keV'
    for c,file in enumerate(files):
        basefile = os.path.basename(file)
        df = pd.read_csv(file,index_col=0, sep = ' ', comment='#')
        minindex = np.argmin(df.index.values)
        dfFilteredDct[basefile]= df.iloc[minindex:]
        f = open(file,'r')
        lines = f.readlines()
        f.close()
        header = ''.join([line for line in lines if '#' in line])
        headers.append(header)
    ZElens = [len(dfFilteredDct[basefile].index.values) for basefile in dfFilteredDct]
    ZEmins = np.array([np.min(dfFilteredDct[file].index.values) for file in dfFilteredDct])
    greatestMin = np.max(ZEmins)
    ZEindex = ZElens.index(max(ZElens))
    ZEkey = list(dfFilteredDct.keys())[ZEindex]
    ZE = dfFilteredDct[ZEkey].index.values
    spacing = np.round((ZE[-1] - ZE[0])/(len(ZE)-1),6)
    
    no_tries = 30
    grid = np.round(np.arange((greatestMin+spacing),ZE[-1],spacing),5)
    
    for c, file in enumerate(dfFilteredDct):
        ZEmin = 0
        ZEmax = -1
        Emin = grid[0]
        Emax = grid[-1]
        newfilerg = f'{coldir}/regrid/{file}'
        if len(dfFilteredDct[file].index.values) < len(grid) - no_tries:
            print(f'{file} too short, couldn\'t be regridded')
            if os.path.exists(newfilerg):
                os.remove(newfilerg)
            continue
        f = open(newfilerg,'w')
        f.write(headers[c])
        f.close()
        regridDF = pd.DataFrame()
        if len([col for col in dfFilteredDct[file].columns if monPattern in col]) == 0:
            continue
        while Emin < dfFilteredDct[file].index.values[0]:
            ZEmin += 1
            Emin = grid[ZEmin]
        while Emax > dfFilteredDct[file].index.values[-1]:
            ZEmax -= 1
            Emax = grid[ZEmax]
        newgrid = grid[ZEmin:ZEmax]
        newgrid = newgrid.round(5)
        monCounter = [c for c in dfFilteredDct[file].columns if monPattern in c][0]
        usedi1counters = [c for c in dfFilteredDct[file].columns if c in i1counters]
        if usedi1counters:
            i1counter = usedi1counters[0]
            muT = np.log(dfFilteredDct[file][monCounter].values/dfFilteredDct[file][i1counter].values)
            gridfunc = interp1d(dfFilteredDct[file].index.values,muT)
            muTregrid = gridfunc(newgrid)
            regridDF['muT'] = muTregrid
        if i2name in dfFilteredDct[file].columns:
            mu2 = np.log(dfFilteredDct[file][i1counter].values/dfFilteredDct[file][i2name].values)
            gridfunc = interp1d(dfFilteredDct[file].index.values,mu2)
            mu2regrid  = gridfunc(newgrid)
            regridDF['mu2'] = mu2regrid
        if fluoCounter in dfFilteredDct[file].columns:
            muF = dfFilteredDct[file][fluoCounter]/dfFilteredDct[file][monCounter]
            gridfunc = interp1d(dfFilteredDct[file].index.values,muF)
            muFregrid = gridfunc(newgrid)
            regridDF['muF'] = muFregrid
        for counter in dfFilteredDct[file].columns:
            if monPattern in counter or counter in i1counters or fluoCounter in counter or counter == i2name:
                gridfunc = interp1d(dfFilteredDct[file].index.values,dfFilteredDct[file][counter].values)
                regridDF[counter] = gridfunc(newgrid).round(1)
        if unit == 'eV':
            newgrid = (newgrid*escale).round(2)
        regridDF.index = newgrid
        regridDF.index.name = f'energy_offset({unit})'
        if len(regridDF.columns) == 0:
            continue
        regridDF.to_csv(newfilerg,sep = ' ',mode = 'a')



def run(direc,thetaOffset=0, unit = 'keV'):

    os.chdir(direc)
    fileDct = {} #dictionary with files as keys, values: [modified time, last spectrum]

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
            processFile(file, fileDct, currentdir,  thetaOffset)
            basename = os.path.splitext(os.path.basename(file))[0]

            regrid(f'{currentdir}columns/{basename}', unit = unit)
            
    return fileDct

if __name__ == '__main__':
    fileDct = run(direc=direc, thetaOffset=thetaOffset)
