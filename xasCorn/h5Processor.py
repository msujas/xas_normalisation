import h5py
import numpy as np
import os, time
import pandas as pd
from glob import glob
import re
from .columnExtraction_thetaCorrection import regrid

def removeByteLabel(string):
    string = string.replace('b\'','')
    string = string.replace('\'','')
    return string

def h5ToDat(hfile, filingIndex = 0):
    braggAxis = 'dcmbragg_trig'
    energyAxis = 'energy_enc'
    currentdir = os.path.dirname(os.path.realpath(hfile))
    basefile = os.path.basename(hfile)
    sampleName = re.sub('_[0-9][0-9][0-9][0-9]','',basefile).replace('.h5','')
    os.makedirs(f'{currentdir}/columns/{sampleName}',exist_ok=True)
    file = h5py.File(hfile,'r')
    keys = list(file.keys())
    print(keys)
    for i in range(len(keys)):
        key = keys[i]
        scanIndex = filingIndex + i
        newSampleName = f'{sampleName}_{scanIndex:04d}'
        x = file[key]
        if not 'measurement' in list(x.keys()):
            continue
        meas = x['measurement']
        dt = removeByteLabel(str(x['start_time'].__array__()))
        title = removeByteLabel(str(x['title'].__array__()))
        endReason = removeByteLabel(str(x['end_reason'].__array__()))
        sampleMeta = removeByteLabel(str(x['sample']['name'].__array__()))
        columns = list(meas.keys())
        
        if not 'SUCCESS' in endReason:
            print(f'{key}: {endReason}')

        df = pd.DataFrame()
        for col in columns:
            try:
                df[col] = meas[col]
            except ValueError:
                print(f'{col} has wrong data length, skipping')
        if 'trigscan' in title:
            try:
                dfEnergy = df.pop(energyAxis)
                df.insert(0,energyAxis,dfEnergy)
                
            except:
                print(f'no energy axis in {hfile}, {key}')
        columns = df.columns
        columnString = '#'+' '.join(columns)
        header = f'#S {scanIndex}\n'
        header += f'#title {title}\n'
        header += f'#sample {sampleMeta}\n'
        header += f'#dt {dt}\n'
        header += f'#er {endReason}\n'
        header += f'{columnString}\n'
        fname = f'{currentdir}/columns/{sampleName}/{newSampleName}.dat'
        f = open(fname,'w')
        f.write(header)
        f.close()
        df.to_csv(fname,sep = ' ', index = False,mode='a', header=False)
        print(fname)
        
        '''
        appendFile = f'{currentdir}/{sampleName}.dat'
        if scanIndex == 0 and os.path.exists(appendFile):
            os.remove(appendFile)
        f = open(appendFile,'a')
        f.write(header)
        f.close()
        df.to_csv(appendFile,sep=' ',index= False,mode = 'a', header = False)
        f = open(appendFile,'a')
        f.write('\n')
        f.close()
        '''
    return filingIndex+i

def runLoop():
    mtimedct = {}
    fileSearch = True
    while True:
        fileIndexDct = {}
        for root, dirs, files in os.walk('.'):
            h5files = glob(f'{root}/*.h5')
            for file in h5files:
                basefile = os.path.basename(file)
                sampleName = re.sub('_[0-9][0-9][0-9][0-9]','',basefile).replace('.h5','')
                sampleNamedir = f'{root}/{sampleName}'
                if not sampleNamedir in fileIndexDct:
                    fileIndexDct[sampleNamedir] = 0
                else:
                    fileIndexDct[sampleNamedir] += 1
                mtime = os.path.getmtime(file)
                if file in mtimedct and mtimedct[file] == mtime:
                    continue
                
                mtimedct[file] = mtime
                print(file)
                if fileSearch:
                    newi = h5ToDat(file, fileIndexDct[sampleNamedir])
                else:
                    newi = 0
                    subfiles = glob(f'{sampleNamedir}*.h5')
                    for sfile in subfiles:
                        newi = h5ToDat(sfile, newi)

                fileSearch = True
                #regrid(f'{root}/columns/{sampleName}',monCountersRG=['ct01','ct02','ct03','ct04'], i1countersRG=['ct05','ct06','ct07','ct08'])
                fileIndexDct[sampleNamedir] = newi
        if fileSearch:
            print('looking for new files')
            fileSearch = False
        time.sleep(1)

if __name__ == '__main__':
    for root,dirs,files in os.walk(r'C:\Users\kenneth1a\Documents\beamlineData\Dec2024_h5'):
        h5files = glob(f'{root}/*.h5')
        for file in h5files:
            print(file)
            h5ToDat(file)