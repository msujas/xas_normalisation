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
        header = f'#S {scanIndex} {title}\n'
        header += f'#NSAMPLE {sampleMeta}\n'
        header += f'#DT {dt}\n'
        header += f'#ER {endReason}\n'
        header += f'{columnString}\n'
        fname = f'{currentdir}/columns/{sampleName}/{newSampleName}.dat'
        f = open(fname,'w')
        f.write(header)
        f.close()
        df.to_csv(fname,sep = ' ', index = False,mode='a', header=False)
        print(fname)
        

    return filingIndex+i

def appendFile(coldir):
        files = glob(f'{coldir}/*.dat')
        if not files:
            return
        sampleName = re.sub('_[0-9][0-9][0-9][0-9].dat','',files[0])
        sampleName = os.path.basename(sampleName)
        string = ''
        for file in files:
            f = open(file,'r')
            string += f.read()
            string += '\n'
            f.close()
        
        appendFile = f'{coldir}/../../{sampleName}.dat'
        f = open(appendFile,'w')
        f.write(string)
        
def runLoop():
    mtimedct = {}
    fileSearch = False
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
                if not fileSearch:
                    newi = h5ToDat(file, fileIndexDct[sampleNamedir])
                else:
                    newi = 0
                    subfiles = glob(f'{sampleNamedir}.h5')
                    subfiles += glob(f'{sampleNamedir}_*.h5')
                    subfiles.sort()
                    for sfile in subfiles:
                        newi = h5ToDat(sfile, newi)
                        newi += 1

                fileSearch = False
                #regrid(f'{root}/columns/{sampleName}',monCountersRG=['ct01','ct02','ct03','ct04'], i1countersRG=['ct05','ct06','ct07','ct08'])
                fileIndexDct[sampleNamedir] = newi
        if not fileSearch:
            print('looking for new files')
            fileSearch = True
        time.sleep(1)
