import os, re
from glob import glob
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


direc = r'C:\Users\kenneth1a\Documents\beamlineData\a311222'
thetaOffset = 0


dspacing = 3.133789
planck = 6.62607015e-34
charge = 1.60217663e-19
speedOfLight = 299792458

digits = 4
fluoCounter = 'xmap_roi00'
monPattern = 'mon_'
ion1Pattern = 'ion_1'
counterNames = ['ZapEnergy','TwoTheta','mon_3','mon_4','mon_1','ion_1_2','ion_1_3','ion_1_1',fluoCounter]
counterNames_NF = [c for c in counterNames if c != fluoCounter] #NF - no fluorescence

def angle_to_kev(angle): #NB the TwoTheta data in the .dat files is really theta
    wavelength = 2*dspacing*np.sin(angle*np.pi/(180))
    wavelength_m = wavelength*10**(-10)
    energy_kev = planck*speedOfLight/(wavelength_m*charge*1000)
    return np.round(energy_kev,6)

def processFile(file, fileDct, currentdir, thetaOffset, startSpectrum = 0):
    f = open(file,'r')
    data = f.read()
    f.close()
    if not 'zapline' in data:
        return
    filemtime = os.path.getmtime(file)
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
    onscan = False
    dfFilteredDct = {}
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
                if monPattern in counter or ion1Pattern in counter:
                    if np.max(dfFilteredDct[spectrum_count][counter].values) < 1000*timeStep:
                        dfFilteredDct[spectrum_count].drop(counter,axis = 1,inplace = True)
            if np.max(dfFilteredDct[spectrum_count][fluoCounter].values) < 10:
                dfFilteredDct[spectrum_count].drop(fluoCounter,axis = 1,inplace = True)


            newfile = f'{newdir}/{basename}_{spectrum_count:0{digits}d}.dat'

            if len([col for col in dfFilteredDct[spectrum_count].columns if monPattern in col]) == 0:

                if os.path.exists(newfile):
                    os.remove(newfile)
                continue
            f = open(newfile,'w')
            f.write(newstring)
            f.close()
            dfFilteredDct[spectrum_count].to_csv(newfile,sep = ' ',mode = 'a')
            print(newfile)

    fileDct[file] = [filemtime,spectrum_count]
    

def regrid(coldir):
    if not os.path.exists(coldir):
        return
    print(coldir)
    files = glob(f'{coldir}/*.dat')
    files.sort()
    if len(files) == 0:
        return
    pattern = '_'
    for n in range(digits):
        pattern += '[0-9]'
    pattern += '.dat'
    basename = re.sub(pattern,'',os.path.basename(files[0]))
    dfFilteredDct = {}
    headers = []
    for c,file in enumerate(files):
        basefile = os.path.basename(file)
        dfFilteredDct[basefile] = pd.read_csv(file,index_col=0, sep = ' ', comment='#')
        f = open(file,'r')
        lines = f.readlines()
        f.close()
        header = ''.join([line for line in lines if '#' in line])
        headers.append(header)
    ZElens = [len(dfFilteredDct[basefile].index.values) for basefile in dfFilteredDct]

    ZEindex = ZElens.index(max(ZElens))
    ZEkey = list(dfFilteredDct.keys())[ZEindex]
    ZE = dfFilteredDct[ZEkey].index.values



    spacing = np.round((ZE[-1] - ZE[0])/(len(ZE)-1),6)


    ZEmin = 0
    ZEmax = -1
    #os.chdir(f'{newdir}regrid/')
    no_tries = 50
    for c, file in enumerate(dfFilteredDct):

        newfilerg = f'{coldir}/regrid/{file}'

        f = open(newfilerg,'w')
        f.write(headers[c])
        f.close()
        regridDF = pd.DataFrame()
        if len([col for col in dfFilteredDct[file].columns if monPattern in col]) == 0:
            continue
        for n in range(no_tries):
            
            try:
                grid = np.arange(ZE[ZEmin].round(4),ZE[ZEmax].round(5),spacing)
                grid = grid.round(5)
                if len(dfFilteredDct[file].index.values) < len(grid) - no_tries + 1:
                    print(f'{basename} spectrum {c} too short, couldn\'t be regridded')
                    if os.path.exists(newfilerg):
                        os.remove(newfilerg)
                    break
                for counter in dfFilteredDct[file].columns:
                    if monPattern in counter or ion1Pattern in counter or fluoCounter in counter:
                        gridfunc = interp1d(dfFilteredDct[file].index.values,dfFilteredDct[file][counter].values)
                        regridDF[counter] = gridfunc(grid).round(1)

                
                regridDF.index = grid
                regridDF.index.name = 'energy_offset(keV)'
                if n != 0:
                    print(f'{basename} spectrum {c} ZEmin = {ZEmin}, ZEmax {ZEmax}')
                break
            except ValueError as e:
                if 'below' in str(e):
                    ZEmin +=1
                elif 'above' in str(e):
                    ZEmax -= 1
                if n == no_tries - 1:
                    print(f'{basename} spectrum {c} too short, couldn\'t be regridded')
                    if os.path.exists(newfilerg):
                        os.remove(newfilerg)
        if len(regridDF.columns) == 0:
            continue
        regridDF.to_csv(newfilerg,sep = ' ',mode = 'a')



def run(direc,thetaOffset=0):

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
        #digits = math.ceil(math.log10(len(datfiles)))
        print(os.getcwd())
        for file in datfiles:
            processFile(file, fileDct, currentdir,  thetaOffset)
            basename = os.path.splitext(os.path.basename(file))[0]

            regrid(f'{currentdir}columns/{basename}')
            
    return fileDct

if __name__ == '__main__':
    fileDct = run(direc=direc, thetaOffset=thetaOffset)
