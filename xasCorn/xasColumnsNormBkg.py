'''
script for running the column extraction and xas normalisation in the background,
constantly looking for new files or if files have been modified and rerunning in
a subdirectory if it finds something there
'''
if not __name__ == '__main__':
    from . import columnExtraction_thetaCorrection, xasNormalisation
else:
    import columnExtraction_thetaCorrection, xasNormalisation
import os
from glob import glob
import time
import argparse

direc = r'C:\Users\kenneth1a\Documents\beamlineData\xasTest'
thetaOffset = 0
waitTime = 1

def main(direc = os.path.curdir,thetaOffset = 0, waitTime = 1):
    '''
    arguments: 
    -to - theta offset to apply to mono to correct energy range (default 0)
    directory - the directory to run in (default current directory in cmd)
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('directory',nargs='?', default=direc, help = 'the directory to run in. Default is current directory')
    parser.add_argument('-to','--thetaOffset', type=float, default=thetaOffset, 
                        help = 'theta offset to apply to monochromator for energy correction')
    parser.add_argument('-u', '--unit', default='keV', type = str, help='keV or eV (default keV) for the regrid and normalised files')
    parser.add_argument('-av', '--averaging', default=1, type = int, help='number of spectra to average over, default 1 (no averaging)')
    parser.add_argument('-e', '--elements',default=None, type = str, help = 'comma separated list of elements, e.g. "-e Fe,Cu", if nothing given will search all elements')
    parser.add_argument('-ee', '--excludeElements',default=None, type = str, help = 'comma separated list of elements to exclude, e.g. "-ee Fe,Cu", if nothing given will not exculde any (doesn\'t do anything if used with -e)')
    args = parser.parse_args()
    thetaOffset = args.thetaOffset
    elements = args.elements
    excludeElements = args.excludeElements
    if elements:
        elements = elements.split(',')
    if excludeElements:
        excludeElements = excludeElements.split(',')
    
    direc = os.path.realpath(args.directory)
    if direc[-1] == '/' or direc[-1] == '\\':
        direc = direc[:-1]
    unit = args.unit
    averaging = args.averaging
    if unit != 'keV' and unit != 'eV':
        print('unit must be "keV" or "eV"')
        return
    print(direc)
    columndirname = 'columns'
    if thetaOffset != 0:
        columndirname = f'columns{thetaOffset:.3f}'

    print('running column extraction')
    fileDct = columnExtraction_thetaCorrection.run(direc,thetaOffset, unit = unit, averaging = averaging, elements=elements, excludeElements=excludeElements)
    print('running normalisation')
    xasNormalisation.run(direc, unit = unit, coldirname=columndirname, elements=elements, excludeElements=excludeElements, averaging=averaging)
    repeat = True
    while True:
        if repeat == True:
            print('looking for new files')
            repeat = False
        for root, dirs, files in os.walk(direc):
            currentdir = root + '/'
            if 'columns' in root:
                continue
            os.chdir(currentdir)
            
            
            if elements:
                datfiles = []
                for e in elements:
                    datfiles += glob(f'*_{e}_*.dat')
            elif excludeElements:
                datfiles = glob(f'*.dat')
                for file in datfiles:
                    for e in excludeElements:
                        if f'_{e}_' in file:
                            datfiles.remove(file)
            else:
                datfiles = glob(f'*.dat')
            datfiles = [currentdir + file for file in datfiles]
            for file in datfiles:
                if 'tempLog' in file:
                    continue
                if file not in list(fileDct.keys()) or os.path.getmtime(file) != fileDct[file][0]:
                    if file not in list(fileDct.keys()):
                        startSpectrum = 0
                    else:
                        startSpectrum = fileDct[file][1]-1
                    repeat = True
                    print(f'running column extraction on {file}')
                    columnExtraction_thetaCorrection.processFile(file, fileDct, currentdir, thetaOffset, startSpectrum=startSpectrum)
                    basename = os.path.splitext(os.path.basename(file))[0]
                    columnExtraction_thetaCorrection.regrid(f'{currentdir}{columndirname}/{basename}', unit = unit, averaging=averaging)
                    columnDir = os.path.basename(file).replace('.dat','')
                    normdir = f'{root}/{columndirname}/{columnDir}'
                    if os.path.exists(normdir):
                        print(f'running normalisation in {normdir}')
                        xasNormalisation.run(normdir, unit = unit, coldirname=columndirname, elements=elements, excludeElements=excludeElements, averaging=averaging)

        time.sleep(waitTime)

if __name__ == '__main__':
    main(direc = direc,thetaOffset=thetaOffset,waitTime=waitTime)