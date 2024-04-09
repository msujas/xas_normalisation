from larch import Group
from larch.math import lincombo_fit
import numpy as np
import os
from glob import glob
import pandas as pd

fluoCounter = 'xmap_roi00'
monPattern = 'mon_'
i1Pattern = 'ion_1_'
direc=  r'C:\Users\kenneth1a\Documents\beamlineData\a311222\Ex_Situ_231105\columns\CHE34_fresh_Re_xanes\regrid\norm'

def run(direc, components = None):
    os.chdir(direc)
    filesT = glob('*T.nor')
    filesT.sort()
    filesF = glob('*F.nor')
    filesF.sort()
    groups = []
    if components == None:
        components = [0,-1]
    for c,file in enumerate(filesT):
        e, mu = np.loadtxt(file,unpack = True,comments='#')
        if np.max(e) < 1000:
            scale = 1000
        else:
            scale = 1

        group = Group(name=str(c), energy=e*scale,norm=mu)
        groups.append(group)
    componentGroups = [groups[i] for i in components]
    fitGroups = []
    for g in groups:
        fitGroup = lincombo_fit(g, componentGroups,arrayname ='norm')
        fitGroups.append(fitGroup)
    return fitGroups
    
if __name__ == '__main__':
    fits = run(direc)
    os.chdir(direc)
    for c,fit in enumerate(fits):
        compVals = list(fit.weights.values())
        if c == 0:
            compNames = list(fit.weights.keys())
            df = pd.DataFrame(columns=compNames)
            df.index.name = 'scan'
        df.loc[c] = compVals
    print(df)
    df.to_csv('linCombo.txt',sep='\t',float_format='%.5f')
    #print(fits)

        


            
        
