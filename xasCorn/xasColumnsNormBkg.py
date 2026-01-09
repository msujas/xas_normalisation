'''
script for running the column extraction and xas normalisation in the background,
constantly looking for new files or if files have been modified and rerunning in
a subdirectory if it finds something there
'''
if not __name__ == '__main__':
    from . import columnExtraction_thetaCorrection
else:
    import columnExtraction_thetaCorrection
import os
from glob import glob
import time
import argparse

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
                        help = ('theta offset to apply to monochromator for energy correction. For side crystal '
                                'approximately -0.103 with default d-spacing 3.13429'))
    parser.add_argument('-u', '--unit', default='keV', type = str, help='keV or eV (default keV) for the regrid and normalised files')
    parser.add_argument('-av', '--averaging', default=1, type = int, help='number of spectra to average over, default 1 (no averaging)')
    parser.add_argument('-e', '--elements',default=None, type = str, help = 'comma separated list of elements, e.g. "-e Fe,Cu", '
    'if nothing given will search all elements')
    parser.add_argument('-ee', '--excludeElements',default=None, type = str, help = 'comma separated list of elements to exclude, e.g. ' \
    '"-ee Fe,Cu", if nothing given will not exculde any (doesn\'t do anything if used with -e)')
    parser.add_argument('-s','--subdir',type = str, default = 'edge', help= 'how to arrange output folders. "edge" - everything with ' \
    'the same element and type (EXAFS or XANES) will go in the same folder, or "file" - every file gets its own folder. '\
        'Default "edge"')
    parser.add_argument('-d','--dspacing',type=float, default=3.13429, help='monochromator d-spacing (default 3.13429). If data from '\
                        'before 9/2025, should be 3.13379 (or use the new value and apply a theta offset, try -0.005198). ' \
                        'Recalibrated to 3.13429 9/2025, so use this after 8/2025')

    args = parser.parse_args()
    thetaOffset = args.thetaOffset
    elements = args.elements
    excludeElements = args.excludeElements
    subdir = args.subdir
    dspacing = args.dspacing
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

    print('running column extraction')
    fileprocessor = columnExtraction_thetaCorrection.XasProcessor(unit = unit, thetaOffset=thetaOffset, dspacing=dspacing, averaging=averaging,
                                                                  elements=elements, excludeElements=excludeElements, subdir=subdir)
    fileprocessor.run(direc)
    print('running normalisation')
    fileprocessor.runNormalisation(direc)
    repeat = True
    while True:
        if repeat == True:
            print('looking for new files')
            repeat = False
        for root, dirs, files in os.walk(direc):
            if 'columns' in root:
                continue
            datfiles = glob(f'{root}/*.dat')
            for file in datfiles:
                if 'tempLog' in file:
                    continue
                if file not in list(fileprocessor.fileDct.keys()) or os.path.getmtime(file) != fileprocessor.fileDct[file].mtime:
                    if file not in list(fileprocessor.fileDct.keys()):
                        startSpectrum = 0
                    else:
                        startSpectrum = fileprocessor.fileDct[file].scanno
                    repeat = True
                    print(f'running column extraction on {file}')
                    fileprocessor.processFile(file,  startSpectrum=startSpectrum)
                    try:
                        outdir = fileprocessor.getoutdir(file)
                    except IndexError:
                        print(f'{file} doesn\'t have correct name format')
                        continue
                    fileprocessor.regrid(outdir)

                    if os.path.exists(outdir):
                        print(f'running normalisation in {outdir}')
                        fileprocessor.runNormalisation(outdir)

        time.sleep(waitTime)