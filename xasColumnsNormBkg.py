'''
script for running the column extraction and xas normalisation in the background,
constantly looking for new files or if files have been modified and rerunning in
a subdirectory if it finds something there
'''
import columnExtraction_thetaCorrection, xasNormalisation
import os
from glob import glob
import time

direc = r'C:\Users\kenneth1a\Documents\beamlineData\a311222/'
thetaOffset = 0
waitTime = 1

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
            if file not in list(fileDct.keys()) or os.path.getmtime(file) != fileDct[file]:
                repeat = True
                print(f'running column extraction on {file}')
                #newfileDct = columnExtraction_thetaCorrection.run(root,thetaOffset)
                columnExtraction_thetaCorrection.processFile(file, fileDct, currentdir, columnExtraction_thetaCorrection.fluoCounter,
                                                             columnExtraction_thetaCorrection.counterNames, columnExtraction_thetaCorrection.counterNames_NF,
                                                             columnExtraction_thetaCorrection.monPattern, columnExtraction_thetaCorrection.ion1Pattern,
                                                             thetaOffset,digits=4)
                columnDir = os.path.basename(file).replace('.dat','')
                print(f'running normalisation in {root}/columns/{columnDir}')
                xasNormalisation.run(f'{root}/columns/{columnDir}')

    time.sleep(waitTime)