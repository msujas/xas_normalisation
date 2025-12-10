import os, re
from glob import glob
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from functools import partial

thetaOffset = 0

dspacing = 3.13429 #3.13379 old value, before 8/2025
planck = 6.62607015e-34
charge = 1.60217663e-19
speedOfLight = 299792458

digits = 4
fluoCounter = 'xmap_roi00'
fluoCounters = ['xmap_roi00', 'Det_5']
monPattern = 'mon_'
ion1Pattern = 'ion_1'
counterNames = ['ZapEnergy','TwoTheta', 'mon_2','mon_3','mon_4','mon_1','ion_1_2','ion_1_3','ion_1_1', 'Det_1', 'Det_2', 'Det_3'] + fluoCounters
counterNames_NF = [c for c in counterNames if c != fluoCounter] #NF - no fluorescence
xColumns = ['ZapEnergy','TwoTheta']
monCounters = ['mon_1', 'mon_2', 'mon_3', 'mon_4']
i1counters = ['ion_1_1', 'ion_1_2', 'ion_1_3', 'Det_1', 'Det_2', 'Det_3']
i2name = 'I2'

eList = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 
         'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Cs_K', 'Ba', 'Ba_K', 
         'La', 'La_K', 'Ce', 'Ce_K', 'Pr', 'Pr_K', 'Nd', 'Nd_K', 'Pm', 'Pm_K', 'Sm', 'Sm_K', 'Eu', 'Eu_K', 'Gd', 'Gd_K', 
         'Tb', 'Tb_K', 'Dy', 'Dy_K', 'Ho', 'Ho_K', 'Er', 'Er_K', 'Tm', 'Tm_K', 'Yb', 'Yb_K', 'Lu', 'Lu_K', 'Hf', 'Ta', 
         'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']

def angle_to_kev(angle, dspacing = dspacing): #NB the TwoTheta data in the .dat files is really theta
    #n lam = 2d sin(theta)
    #E = hc/lam
    #V = E/qe
    wavelength = 2*dspacing*np.sin(angle*np.pi/(180))
    wavelength_m = wavelength*10**(-10)
    energy_kev = planck*speedOfLight/(wavelength_m*charge*1000)
    return np.round(energy_kev,6)


def getoutdir(file, coldir, subdir):
    basename = os.path.splitext(os.path.basename(file))[0]
    filesplit = basename.split('_')
    #edge = '_'.join(basename.split('_')[-2:])
    method = filesplit[-1]
    fileStart = '_'.join(filesplit[:-1])
    #print(fileStart)
    element = [e for e in eList if fileStart.endswith(e)][0]


    edge = f'{element}_{method}'
    match subdir:
        case 'edge': newdir = f'{coldir}/{edge}/'
        case 'file': newdir = f'{coldir}/{basename}/'
        case _: raise ValueError('subdir must be "edge" or "file"')
    return newdir

def processFile(file, fileDct, currentdir, thetaOffset, startSpectrum = 0, subdir = 'edge', dspacing = dspacing):
    f = open(file,'r')
    data = f.read()
    f.close()
    filemtime = os.path.getmtime(file)
    if not 'zapline' in data:
        fileDct[file] = [filemtime,-1]
        return
    angle_to_kev_func = partial(angle_to_kev, dspacing=dspacing)
    basename = os.path.splitext(os.path.basename(file))[0]
    coldir = currentdir+'columns/'
    if thetaOffset != 0:
        coldir = currentdir+f'columns{thetaOffset:.3f}/'
    if not os.path.exists(coldir):
        os.makedirs(coldir)
    try:
        newdir = getoutdir(file, coldir, subdir)
    except IndexError:
        print(f'{file} does not have element information')
        return
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
            newstring += line
        elif '#D' in line and onscan:
            newstring += line
        elif '#L' in line and onscan:
            dfstart = c +1
            columns = line.replace('#L ','').split()
            scanStart = True
            lineno = 0
        elif scanStart and '#' not in line and line:
            lineSplit = np.array([np.fromstring(line,sep = ' ')])
            if lineno == 0:
                array = lineSplit
                lineShape = lineSplit.shape
                lineno += 1
            elif lineSplit.shape == lineShape:
                array = np.append(array,lineSplit,axis = 0)
            
        elif ('#C' in line or not line) and onscan:
            dfend = c
            scanStart = False
            onscan = False
            if dfend-dfstart <= 1:
                continue
            df = pd.DataFrame(data=array,columns=columns)

            dfFiltered = df[xColumns].copy(deep=True)
            dfFiltered['Theta_offset'] =  dfFiltered['TwoTheta'].apply(lambda x: np.round(x + thetaOffset,7))
            dfFiltered['ZapEnergy_offset'] = dfFiltered['Theta_offset'].apply(angle_to_kev_func)
            dfFiltered.set_index('ZapEnergy_offset',inplace = True)
            dfFiltered.index.name = '#ZapEnergy_offset'
            energy = dfFiltered.index.values

            usedMon = df[monCounters].sum().idxmax()
            
            dfFiltered[usedMon] = df[usedMon].values
            newfile = f'{newdir}/{basename}_{spectrum_count:0{digits}d}.dat'
            if np.min(dfFiltered[usedMon].values) < 1000*timeStep: #check if beam was off during scan
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

            usedFluos = [col for col in fluoCounters if col in df.columns and np.max(df[col].values) > 50]
            for fluoCounter in usedFluos:
                dfFiltered[fluoCounter] = df[fluoCounter].values

            if (not usedI1s and not usedFluos) or np.max(energy) - np.min(energy) < 0.1:
                if os.path.exists(newfile):
                    os.remove(newfile)
                continue
            
            f2 = open(newfile,'w')
            f2.write(newstring)
            f2.close()
            dfFiltered.to_csv(newfile,sep = ' ',mode = 'a')
            print(newfile)

    fileDct[file] = [filemtime,spectrum_count]
    

def merge(regriddir, unit = 'keV'):
    os.makedirs(f'{regriddir}/merge/',exist_ok=True)
    mergedct = {}
    files = set(['_'.join(file.split('_')[:-1]) for file in glob(f'{regriddir}/*.dat')])
    print(files)
    energycol = f'energy_offset({unit})'
    for file in files:
        basefile = os.path.basename(file)
        files2 = glob(f'{file}*.dat')
        print(files2)
        mergedct[file] = {}
        e0 = 0
        eend = 100000
        for f in files2:
            print(f)
            fr = open(f,'r')
            headcol = [line for line in fr.readlines() if line.startswith('#')][-1].replace('\n','').replace('#','')
            fr.close()
            cols = headcol.split()
            df = pd.read_csv(f, comment='#', sep = ' ', header=None)
            df.columns = cols
            mergedct[file][f] = df
            energy = df[energycol].values

            e0t = round(energy[0],5)
            eendt = round(energy[-1],5)
            if e0t > e0:
                e0 = e0t
            if eendt < eend:
                eend = eendt
        muFsum = 0
        muTsum = 0
        muTcount = 0
        muFcount = 0
        minindex = np.abs(energy-e0).argmin()
        maxindex = np.abs(energy-eend).argmin()
        energyAxis = energy[minindex:maxindex+1]
        for f in mergedct[file]:
            energy = mergedct[file][f][energycol].values
            minindex = np.abs(energy-e0).argmin()
            maxindex = np.abs(energy-eend).argmin()
            df = mergedct[file][f].iloc[minindex:maxindex+1]
            if 'muT' in df.columns:
                muTsum += df['muT'].values[minindex:maxindex+1]
                muTcount += 1
            if 'muF1' in df.columns:
                muFsum += df['muF1'].values[minindex:maxindex+1]
                muFcount += 1
        if muTcount:
            muTsum = muTsum/muTcount
            np.savetxt(f'{regriddir}/merge/{basefile}_T_merge.dat',np.array([energyAxis,muTsum]).transpose(),fmt = '%.5f', 
                       header=f'{energycol} muT')
        if muFcount:
            muFsum = muFsum/muFcount
            np.savetxt(f'{regriddir}/merge/{basefile}_F_merge.dat',np.array([energyAxis,muFsum]).transpose(),fmt = '%.5f',
                       header=f'{energycol} muF')
        

def regrid(coldir, unit = 'keV', averaging = 1, i1countersRG = None, monCountersRG = None):
    if i1countersRG == None:
        i1countersRG = i1counters
    if monCountersRG == None:
        monCountersRG = monCounters
    if not os.path.exists(coldir):
        return
    print(f'regridding {coldir}')
    print(coldir)
    files = glob(f'{coldir}/*.dat')
    files.sort()


    if len(files) == 0:
        return
    if not os.path.exists(f'{coldir}/regrid'):
        os.makedirs(f'{coldir}/regrid')
    avdir = f'{coldir}/regridAv{averaging}/'
    if averaging > 1 and not os.path.exists(avdir):
        os.makedirs(avdir)

    dfFilteredDct = {}
    headers = []
    if unit == 'keV':
        escale = 1
    elif unit == 'eV':
        escale = 1000
    else:
        escale = 1
        unit = 'keV'
    for i,file in enumerate(files):
        f = open(file,'r')
        lines = f.readlines()
        f.close()
        header = [line for line in lines if '#' in line]
        colnames = header[-1].split()
        header = ''.join(header[:-1])
        headers.append(header)
        

        basefile = os.path.basename(file)
        df = pd.read_csv(file, sep = ' ', comment='#', names = colnames)
        df = df.set_index(colnames[0])
        minindex = np.argmin(df.index.values)
        dfFilteredDct[basefile]= df.iloc[minindex:]

    
    ZElens = [len(dfFilteredDct[basefile].index.values) for basefile in dfFilteredDct]
    ZEmins = np.array([np.min(dfFilteredDct[file].index.values) for file in dfFilteredDct])
    greatestMin = np.max(ZEmins)
    ZEindex = ZElens.index(max(ZElens))
    ZEkey = list(dfFilteredDct.keys())[ZEindex]
    ZE = dfFilteredDct[ZEkey].index.values
    spacing = np.round((ZE[-1] - ZE[0])/(len(ZE)-1),6)
    
    no_tries = 30
    grid = np.round(np.arange((greatestMin+spacing),ZE[-1],spacing),5)
    averagingCount = 0
    averagingDct = {}
    avheader = ''
    fluoAv = []
    transAv = []
    oldbasefile = ''
    for n, file in enumerate(dfFilteredDct):
        basefile = '_'.join(file.split('_')[:-1])
        if basefile != oldbasefile:
            averagingCount = 0
            avheader = ''
            fluoAv = []
            transAv = []
        oldbasefile = basefile
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
        f.write(headers[n])
        f.close()
        regridDF = pd.DataFrame()
        if len([col for col in dfFilteredDct[file].columns if col in monCountersRG]) == 0:
            continue
        while Emin < dfFilteredDct[file].index.values[0]:
            ZEmin += 1
            Emin = grid[ZEmin]
        while Emax > dfFilteredDct[file].index.values[-1]:
            ZEmax -= 1
            Emax = grid[ZEmax]
        newgrid = grid[ZEmin:ZEmax]
        newgrid = newgrid.round(5)
        monCounter = [c for c in dfFilteredDct[file].columns if c in monCountersRG][0]
        usedi1counters = [c for c in dfFilteredDct[file].columns if c in i1countersRG]
        if usedi1counters:
            trans = True
        else:
            trans = False
        usedFluos = [col for col in dfFilteredDct[file].columns if col in fluoCounters]
        if usedFluos:
            fluo = True
        else:
            fluo=False

        i2 = i2name in dfFilteredDct[file].columns
        fluoAv.append(fluo)
        transAv.append(trans)
        if trans:
            i1counter = usedi1counters[0]
            muT = np.log(dfFilteredDct[file][monCounter].values/dfFilteredDct[file][i1counter].values)
            gridfunc = interp1d(dfFilteredDct[file].index.values,muT)
            muTregrid = gridfunc(newgrid)
            regridDF['muT'] = muTregrid
        if i2:
            mu2 = np.log(dfFilteredDct[file][i1counter].values/dfFilteredDct[file][i2name].values)
            gridfunc = interp1d(dfFilteredDct[file].index.values,mu2)
            mu2regrid  = gridfunc(newgrid)
            regridDF['mu2'] = mu2regrid

        for c2,fluoCounter in enumerate(usedFluos):
            muF = dfFilteredDct[file][fluoCounter]/dfFilteredDct[file][monCounter]
            gridfunc = interp1d(dfFilteredDct[file].index.values,muF)
            muFregrid = gridfunc(newgrid)
            regridDF[f'muF{c2+1}'] = muFregrid

        for counter in dfFilteredDct[file].columns:
            if counter in monCountersRG or counter in i1countersRG or counter in fluoCounters or counter == i2name:
                gridfunc = interp1d(dfFilteredDct[file].index.values,dfFilteredDct[file][counter].values)
                regridDF[counter] = gridfunc(newgrid).round(1)
        if unit == 'eV':
            newgrid = (newgrid*escale).round(2)
        regridDF.index = newgrid
        regridDF.index.name = f'#energy_offset({unit})'

        if len(regridDF.columns) == 0:
            continue
        regridDF.to_csv(newfilerg,sep = ' ',mode = 'a')
        
        if averaging <= 1:
            continue

        averagingDct[averagingCount] = {}
        averagingDct[averagingCount]['ZEmin'] = ZEmin
        averagingDct[averagingCount]['ZEmax'] = ZEmax
        averagingDct[averagingCount]['ZapEnergy'] = newgrid
        avheader += headers[n]
        if usedi1counters:
            averagingDct[averagingCount]['muT'] = regridDF['muT'].values

        for c,fluoCounter in enumerate(usedFluos):
            averagingDct[averagingCount][f'muF{c+1}'] = regridDF[f'muF{c+1}'].values
            averagingDct[averagingCount][fluoCounter] = regridDF[fluoCounter].values
        if averagingCount == averaging-1:
            averagingCount = -1
            zemins = np.array([averagingDct[i]['ZapEnergy'][0] for i in np.arange(averaging, dtype = np.uint8)])
            zemaxs = np.array([averagingDct[i]['ZapEnergy'][-1] for i in np.arange(averaging, dtype = np.uint8)])
            zeminIs = np.array([averagingDct[i]['ZEmin'] for i in np.arange(averaging, dtype = np.uint8)])
            zemaxIs = np.array([averagingDct[i]['ZEmax'] for i in np.arange(averaging, dtype = np.uint8)])
            zeminmax = np.max(zemins)
            zemaxmin = np.min(zemaxs)
            zeminI = np.max(zeminIs)
            zemaxI = np.min(zemaxIs)
            avGrid = grid[zeminI:zemaxI]
            avMuT = np.empty(shape = (len(avGrid),averaging))
            avMuF = {}
            avFluoRaw = {}
            for i in range(len(usedFluos)):
                avMuF[i] =  np.empty(shape = (0,len(avGrid)))
                avFluoRaw[i] = np.empty(shape = (0,len(avGrid)))
            for i in range(averaging):
                zapenergy = averagingDct[i]['ZapEnergy']
                minindex = np.argmin(np.abs(zapenergy - zeminmax))
                maxindex = np.argmin(np.abs(zapenergy - zemaxmin))
                if transAv[i]:
                    avMuT[:,i] = averagingDct[i]['muT'][minindex:maxindex+1]
                for c,fluoCounter in enumerate(usedFluos):
                    if not f'muF{c+1}':
                        continue
                    avMuF[c] = np.append(avMuF[c] ,[averagingDct[i][f'muF{c+1}'][minindex:maxindex+1]], axis = 0)
                    avFluoRaw[c] = np.append(avFluoRaw[c],[averagingDct[i][fluoCounter][minindex:maxindex+1]], axis = 0)
            avDf = pd.DataFrame()
            avDf[f'#ZapEnergy({unit})'] = avGrid
            if transAv[i]:
                avMuT = np.average(avMuT,axis = 1)
                avDf['muT'] = avMuT
            for c,fluoCounter in enumerate(usedFluos):
                avMuF[c] = np.average(avMuF[c], axis = 0)
                avFluoRaw[c] = np.average(avFluoRaw[c],axis = 0)
                avDf[f'muF{c+1}'] = avMuF[c]
                avDf[fluoCounter] = avFluoRaw[c]
            
            fname = f'{avdir}/{file}'
            f = open(fname,'w')
            f.write(avheader)
            f.close()
            avDf.to_csv(fname, index = False, mode = 'a', sep = ' ')
            avheader = ''
            fluoAv = []
            transAv = []
        averagingCount += 1
    merge(f'{coldir}/regrid',unit=unit)

def run(direc,thetaOffset=0, unit = 'keV', averaging = 1, elements = None, excludeElements = None, subdir = 'edge', dspacing=dspacing):

    os.chdir(direc)
    fileDct = {} #dictionary with files as keys, values: [modified time, last spectrum]

    for root, dirs, files in os.walk(os.getcwd()):
        if 'columns' in root:
            continue
        print(root)
        currentdir = root + '/'
        coldir = f'{currentdir}columns'
        if thetaOffset != 0:
            coldir = f'{currentdir}columns{thetaOffset:.3f}'
        os.chdir(currentdir)

        if elements:
            datfiles = []
            for e in elements:
                datfiles += glob(f'*_{e}_*.dat')
        elif excludeElements:
            datfiles = glob('*.dat')
            for file in files:
                for e in excludeElements:
                    if f'_{e}_' in file:
                        datfiles.remove(file)
                

        else:
            datfiles = glob('*.dat')
        
        datfiles = [currentdir + file for file in datfiles]

        if len(datfiles) == 0:
            continue
        print(os.getcwd())
        for file in datfiles:
            processFile(file, fileDct, currentdir,  thetaOffset, subdir = subdir, dspacing=dspacing)
            try:
                outdir = getoutdir(file, coldir, subdir)           
            except IndexError:
                print(f'{file} doesn\'t have correct name format')
                continue

            regrid(f'{outdir}', unit = unit, averaging=averaging)
            
    return fileDct

