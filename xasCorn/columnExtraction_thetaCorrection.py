import os
from glob import glob
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from functools import partial
import xasCorn.xasNormalisation as xasn
import logging
import pathlib

logger = logging.getLogger()
home = pathlib.Path.home()
logdir = f'{home}/.log/xascorn'
os.makedirs(logdir,exist_ok=True)
logfile = f'{logdir}.log'
logging.basicConfig(filename=logfile, level = logging.INFO, format = '%(asctime)s %(levelname)-8s %(message)s',
                        datefmt = '%Y/%m/%d_%H:%M:%S')
thetaOffset = 0

dspacing = 3.13439 #3.13429 before 6/2026, 3.13379 before 8/2025
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


class FileInfo():
    def __init__(self, mtime, scanno=-1):
        self.mtime = mtime
        self.scanno = scanno

class XasProcessor():
    def __init__(self,unit = 'keV', thetaOffset = 0 , dspacing=dspacing, averaging = 1, elements:list = None, 
                 excludeElements:list = None, subdir = 'edge', cpsThreshold = 10000):
        self.fileDct: dict[str, FileInfo] = {}
        self.unit = unit
        self.thetaOffset = thetaOffset
        self.dspacing = dspacing
        self.averaging = averaging
        self.elements = elements
        self.excludeElements = excludeElements
        self.subdir = subdir
        self.columnsubdir = 'columns'
        self.cpsThreshold = cpsThreshold
        if thetaOffset != 0:
            self.columnsubdir += f'{thetaOffset:.3f}'
        self.angle_to_kev_func = partial(angle_to_kev, dspacing=self.dspacing)
        self.runNormalisation = partial(xasn.run,unit = self.unit, coldirname=self.columnsubdir, elements=self.elements, 
                                        excludeElements=self.excludeElements, averaging= self.averaging)
    def getoutdir(self,file):
        coldir = f'{os.path.dirname(file)}/{self.columnsubdir}'
        basename = os.path.splitext(os.path.basename(file))[0]
        filesplit = basename.split('_')
        method = filesplit[-1]
        fileStart = '_'.join(filesplit[:-1])
        element = [e for e in eList if fileStart.endswith(e)][0]
        edge = f'{element}_{method}'
        match self.subdir:
            case 'edge': newdir = f'{coldir}/{edge}/'
            case 'file': newdir = f'{coldir}/{basename}/'
            case _: raise ValueError('subdir must be "edge" or "file"')
        return newdir
    def processFile(self,file, startSpectrum = 0, savefiles = True) -> pd.DataFrame:
        dfFiltered = pd.DataFrame()
        currentdir = os.path.dirname(file)
        f = open(file,'r')
        data = f.read()
        f.close()
        filemtime = os.path.getmtime(file)
        if not 'zapline' in data:
            self.fileDct[file] = FileInfo(filemtime,-1)
            return
        if self.elements:
            for e in self.elements:
                if not f'_{e}_exafs.dat' in file and not f'_{e}_xanes.dat' in file:
                    self.fileDct[file] = FileInfo(filemtime,-1)
                    return
        elif self.excludeElements:
            for e in self.excludeElements:
                if f'_{e}_exafs.dat' in file or f'_{e}_xanes.dat' in file:
                    self.fileDct[file] = FileInfo(filemtime,-1)
                    return
        
        basename = os.path.splitext(os.path.basename(file))[0]
        coldir = f'{currentdir}/{self.columnsubdir}/'

        if not os.path.exists(coldir):
            os.makedirs(coldir)
        try:
            newdir = self.getoutdir(file)
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
                dfFiltered['Theta_offset'] =  dfFiltered['TwoTheta'].apply(lambda x: np.round(x + self.thetaOffset,7))
                dfFiltered['ZapEnergy_offset'] = dfFiltered['Theta_offset'].apply(self.angle_to_kev_func)
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
                usedI1s = [col for col in i1counters if np.max(df[col].values) > self.cpsThreshold*timeStep and np.min(df[col].values) > 1]
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
                if not savefiles:
                    continue
                f2 = open(newfile,'w')
                f2.write(newstring)
                f2.close()
                dfFiltered.to_csv(newfile,sep = ' ',mode = 'a')
                print(newfile)
        self.fileDct[file] = FileInfo(filemtime,spectrum_count)
        return dfFiltered #returns DF of last scan
        
    def merge(self,regriddir):
        if not os.path.exists(regriddir):
            return
        os.makedirs(f'{regriddir}/merge/',exist_ok=True)
        print(f'merging {regriddir}')
        mergedct = {}
        files = set(['_'.join(file.split('_')[:-1]) for file in glob(f'{regriddir}/*.dat')])
        energycol = f'energy_offset({self.unit})'
        match os.path.split(regriddir)[-1]:
            case 'trans':
                mucol = 'muT'
                filepart = 'T'
                head = 'muT'
            case 'fluo':
                mucol = 'muF1'
                filepart = 'F'
                head = 'muF'
            case _:
                print("not valid regrid directory, skipping")
                return
        for file in files:
            basefile = os.path.basename(file)
            files2 = glob(f'{file}*.dat')
            mergedct[file] = {}
            e0 = 0
            eend = 100000
            musum = 0
            count = 0
            for f in files2:
                fr = open(f,'r')
                headcol = [line for line in fr.readlines() if line.startswith('#')][-1].replace('\n','').replace('#','')
                fr.close()
                cols = headcol.split()
                df = pd.read_csv(f, comment='#', sep = ' ', header=None)
                df.columns = cols
                mergedct[file][f] = df
                energy = df[energycol].values
                musum += df[mucol].values
                count += 1

            musum = musum/count
            np.savetxt(f'{regriddir}/merge/{basefile}_{filepart}_merge.dat',np.array([energy,musum]).transpose(),fmt = '%.5f', 
                    header=f'{energycol} {head}')

    def getElement(self,coldir):
        for e in eList:
            if f'{e}_exafs' or f'{e}_xanes' in coldir:
                return e
        
    def regrid(self,coldir,   i1countersRG = None, monCountersRG = None):
        if i1countersRG == None:
            i1countersRG = i1counters
        if monCountersRG == None:
            monCountersRG = monCounters
        if not os.path.exists(coldir):
            return
        element = self.getElement(coldir)
        if self.elements:
            if element not in self.elements:
                return
        elif self.excludeElements:
            if element in self.excludeElements:
                return
        print(f'regridding {coldir}')
        print(coldir)
        files = glob(f'{coldir}/*.dat')
        files.sort()


        if len(files) == 0:
            return
        transdir = f'{coldir}/regrid/trans'
        fluodir = f'{coldir}/regrid/fluo'


        dfFilteredDct:dict[str, pd.DataFrame] = {}
        headers = []
        if self.unit == 'keV':
            escale = 1
        elif self.unit == 'eV':
            escale = 1000
        else:
            escale = 1
            self.unit = 'keV'
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
        maxlen = max(ZElens)
        length_tolerance = 30
        dellist = []
        for l, item in zip(ZElens, dfFilteredDct):
            if l < maxlen - length_tolerance:
                dellist.append(item)
        for d in dellist:
            print(f'{d} too short, not regridding')
            outfileT = f'{transdir}/{d}'
            outfileF = f'{fluodir}/{d}'
            if os.path.exists(outfileT): os.remove(outfileT)
            if os.path.exists(outfileF): os.remove(outfileF)
            dfFilteredDct.pop(d,None)
        ZElens = [len(dfFilteredDct[basefile].index.values) for basefile in dfFilteredDct] # regenerating due to deleted values
        ZEmins = np.array([np.min(dfFilteredDct[file].index.values) for file in dfFilteredDct])
        ZEmaxs = np.array([np.max(dfFilteredDct[file].index.values) for file in dfFilteredDct])
        greatestMin = np.max(ZEmins)
        smallestMax = np.min(ZEmaxs)
        ZEindex = ZElens.index(max(ZElens))
        ZEkey = list(dfFilteredDct.keys())[ZEindex]
        ZE = dfFilteredDct[ZEkey].index.values
        spacing = (ZE[-1] - ZE[0])/(len(ZE)-1)
        
        
        grid = np.round(np.arange((greatestMin+spacing),smallestMax,spacing),6)
        fluoAv = []
        transAv = []
        oldbasefile = ''
        for n, file in enumerate(dfFilteredDct):
            basefile = '_'.join(file.split('_')[:-1])
            if basefile != oldbasefile:
                fluoAv = []
                transAv = []
            oldbasefile = basefile
            newfilergT = f'{transdir}/{file}'
            newfilergF = f'{fluodir}/{file}'

            regridDF = pd.DataFrame()
            if len([col for col in dfFilteredDct[file].columns if col in monCountersRG]) == 0:
                continue
            monCounter = [c for c in dfFilteredDct[file].columns if c in monCountersRG][0]
            usedi1counters = [c for c in dfFilteredDct[file].columns if c in i1countersRG]
            if usedi1counters:
                trans = True
                os.makedirs(transdir,exist_ok=True)
                f = open(newfilergT,'w')
                f.write(headers[n])
                f.close()
            else:
                trans = False
            usedFluos = [col for col in dfFilteredDct[file].columns if col in fluoCounters]
            if usedFluos:
                fluo = True
                os.makedirs(fluodir,exist_ok=True)
                f = open(newfilergF,'w')
                f.write(headers[n])
                f.close()
            else:
                fluo=False

            i2 = i2name in dfFilteredDct[file].columns
            fluoAv.append(fluo)
            transAv.append(trans)
            def tryregrid(mu, colname):
                try:
                    gridfunc = interp1d(dfFilteredDct[file].index.values,mu)
                    muregrid = gridfunc(grid)
                    regridDF[colname] = muregrid
                    return muregrid
                except ValueError as e:
                    logging.error(f'problem regridding {file}\n{e}')
                    raise e

            if trans:
                i1counter = usedi1counters[0]
                muT = np.log(dfFilteredDct[file][monCounter].values/dfFilteredDct[file][i1counter].values)
                tryregrid(muT, 'muT')

            if i2:
                mu2 = np.log(dfFilteredDct[file][i1counter].values/dfFilteredDct[file][i2name].values)
                tryregrid(mu2, 'mu2')

            for c2,fluoCounter in enumerate(usedFluos):
                muF = dfFilteredDct[file][fluoCounter]/dfFilteredDct[file][monCounter]
                tryregrid(muF, f'muF{c2+1}')

            for counter in dfFilteredDct[file].columns: #saving regrid of original counters
                if counter in monCountersRG or counter in i1countersRG or counter in fluoCounters or counter == i2name:
                    countervalues = dfFilteredDct[file][counter].values
                    try:
                        gridfunc = interp1d(dfFilteredDct[file].index.values,countervalues)
                        regridDF[counter] = gridfunc(grid).round(1)
                    except ValueError as e:
                        logging.error(f'problem regridding counters in {file}\n{e}')
                        raise e

            if self.unit == 'eV':
                grid = (grid*escale).round(2)
            regridDF.index = grid
            regridDF.index.name = f'#energy_offset({self.unit})'

            if len(regridDF.columns) == 0:
                continue
            if trans:
                regridDF[['muT', monCounter, i1counter]].to_csv(newfilergT,sep=' ', mode='a')
            if fluo:
                regridDF[['muF1', monCounter, *usedFluos]].to_csv(newfilergF,sep = ' ',mode = 'a')
            

        self.average(transdir,averaging=self.averaging)
        self.average(fluodir, averaging=self.averaging)
        self.merge(transdir)
        self.merge(fluodir)

    def average(self, regriddir:str, averaging:int ):
        if averaging <=1:
            return
        files = glob(f'{regriddir}/*.dat')
        if regriddir[-1] == '\\' or regriddir[-1] == '/':
            regriddir = regriddir[:-1]
        technique = os.path.basename(regriddir)
        avdir = f'{regriddir}/../../regridAv{averaging}/{technique}'
        os.makedirs(avdir, exist_ok=True)
        print(f'averaging {averaging} for {regriddir}')
        def getbasename(file:str):
            basefile = os.path.basename(file)
            return '_'.join(basefile.split('_')[:-1])
        
        basenames = []
        for file in files:
            basename = getbasename(file)
            basenames.append(basename)
        basenames = set(basenames)
        for basename in basenames:
            i = 0
            for file in sorted(glob(f'{regriddir}/{basename}*.dat')):
                energy, mu = np.loadtxt(file, usecols=(0,1), unpack=True, comments='#')
                f = open(file,'r')
                header = ''.join([line for line in f.readlines() if line[0]=='#'][:-1])
                f.close()
                if i == 0:
                    avlist = []
                    fullheader = ''
                    minenergy = energy[0]
                    maxenergy = energy[-1]
                if np.min(energy) > minenergy:
                    minenergy = energy[0]
                if np.max(energy) < maxenergy:
                    maxenergy = energy[-1]
                    
                fullheader += header
                avlist.append(np.array([energy,mu]))
                if i == averaging-1:
                    techstring = technique[0].upper()
                    if techstring == 'F':
                        techstring += '1'
                    fullheader += f'#energy_offset(keV) mu{techstring}'
                    minindex =  np.argmin(np.abs(energy-minenergy))
                    maxindex = np.argmin(np.abs(energy-maxenergy))
                    energyaxis = energy[minindex:maxindex+1]
                    avmu = np.empty(shape = (len(energyaxis), averaging))
                    for c,av in enumerate(avlist):
                        minindex = np.argmin(np.abs(av[0]-minenergy))
                        maxindex = np.argmin(np.abs(av[0]-maxenergy))
                        try:
                            avtrunc = av[1][minindex:maxindex+1]
                            avmu[:,c] = avtrunc
                        except Exception as e:
                            print(minenergy,maxenergy, c)
                            print(energyaxis[0], energyaxis[-1])
                            print(avtrunc.shape)
                            print(avtrunc)
                            print(av.shape)
                            print(len(energyaxis), len(avtrunc[0]) )
                            raise e
                    avmu = np.mean(avmu, axis = 1)
                    avdata = np.array([energyaxis,avmu])
                    basefile = os.path.basename(file)
                    newfile = f'{avdir}/{basefile}'
                    print(f'saving {newfile}')
                    np.savetxt(newfile, avdata.transpose(), fmt = "%.6f", header=fullheader, comments='')
                    i = -1
                i +=1
            
    def run(self,direc):
        for root, _dirs, _files in os.walk(direc):
            if 'columns' in root:
                continue
            print(root)
            datfiles = glob(f'{root}/*.dat')
            if len(datfiles) == 0:
                continue
            for file in datfiles:
                self.processFile(file)
                try:
                    outdir = self.getoutdir(file)           
                except IndexError:
                    print(f'{file} doesn\'t have correct name format')
                    continue
                if self.fileDct[file].scanno == -1:
                    continue
                self.regrid(outdir)


def getDFedgeStep(dfFiltered:pd.DataFrame):
    '''
    get edge step for a scan
    '''
    cols = dfFiltered.columns
    moncounter = [col for col in cols if col in monCounters][0]
    i1s = [col for col in cols if col in i1counters]
    fcs = [col for col in cols if col in fluoCounters]
    tstep = 0
    fstep = 0
    mon = dfFiltered[moncounter]
    e = dfFiltered.index.values
    if i1s:
        i1c = i1s[0]
        muT = np.log10(mon/dfFiltered[i1c])
        ds = pd.Series(data = muT, index = e)
        gT = xasn.normalise(ds)
        tstep = gT.edge_step
    if fcs:
        fc = fcs[0]
        muF = dfFiltered[fc]/mon
        ds = pd.Series(data= muF, index = e)
        gF = xasn.normalise(ds)
        fstep = gF.edge_step
    return tstep, fstep

def getlastscan(file):
    f = open(file,'r')
    s = 0
    for line in f:
        if line.startswith('#S'):
            s = int(line.split()[1])
    f.close()
    return s

def filegetedgestep(file):
    lastscan = getlastscan(file)
    df = XasProcessor().processFile(file, startSpectrum=lastscan-1,savefiles=False)
    return getDFedgeStep(df)


