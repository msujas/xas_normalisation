'''
script for running the column extraction and xas normalisation in the background,
constantly looking for new files or if files have been modified and rerunning in
a subdirectory if it finds something there
'''
try:
    from . import columnExtraction_thetaCorrection, xasNormalisation
except ImportError:
    import columnExtraction_thetaCorrection, xasNormalisation
import os
from glob import glob
import time
import argparse

direc = r'C:\Users\kenneth1a\Documents\beamlineData\May2024carlos'
thetaOffset = 0
waitTime = 1

def main(direc = '.',thetaOffset = 0, waitTime = 1):
    '''
    arguments: 
    -to - theta offset to apply to mono to correct energy range (default 0)
    directory - the directory to run in (default current directory in cmd)
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('directory',nargs='?', default=direc, help = 'the directory to run in')
    parser.add_argument('-to','--thetaOffset', type=float, default=thetaOffset, 
                        help = 'theta offset to apply to monochromator for energy correction')
    args = parser.parse_args()

    direc = args.directory
    if direc[-1] == '/' or direc[-1] == '\\':
        direc = direc[:-1]

    
    
    print('running column extraction')
    fileDct = columnExtraction_thetaCorrection.run(direc,thetaOffset)
    print('running normalisation')
    xasNormalisation.run(direc)
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
            datfiles = glob('*.dat')
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
                    columnExtraction_thetaCorrection.regrid(f'{currentdir}columns/{basename}')
                    columnDir = os.path.basename(file).replace('.dat','')
                    normdir = f'{root}/columns/{columnDir}'
                    if os.path.exists(normdir):
                        print(f'running normalisation in {root}/columns/{columnDir}')
                        xasNormalisation.run(f'{root}/columns/{columnDir}')

        time.sleep(waitTime)

if __name__ == '__main__':
    main(direc = direc,thetaOffset=thetaOffset,waitTime=waitTime)