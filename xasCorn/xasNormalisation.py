from larch.xafs import find_e0, pre_edge
from larch import Group
import numpy as np
import os, re
import pandas as pd

direc = r''



def normalise(ds, exafsnorm = 3, xanesnorm = 1):
    '''
    Takes a pandas Series with values as mu and index as energy as an argument.
    The normalisation orders seem to be different between Athena and Larch. 
    1 and 2 in Larch seem to correspond to 2 and 3 in Athena (linear and quadratic), respectively. 0 Seems not to correspond to 1, however.
    0 appears linear, but with a shallower gradient than 1. The documentation recommends 0 if < 50 eV used to fit post-edge.
    The values used in this are roughly the same as the defaults
    Some distributions of Larch don't normalise data properly in keV, so data is converted to eV for normalising.
    Returns a Larch Group object with e0 and normalisation calculated.
    '''
    energy = ds.index.values
    if np.max(energy) < 1000:
        kev = True
    else:
        kev = False
    if kev: #converting axis to eV
        scale = 1000
    else:
        scale = 1
    group = Group(energy = energy*scale, mu = ds.values)
    
    find_e0(group = group, energy = group.energy, mu = group.mu)
    pre1 = group.energy[0] - group.e0
    pre2 = -30
    if group.energy[-1] - group.energy[0] > 500: #EXAFS
        post1 = 150
        nnorm = exafsnorm
    else: #XANES
        post1 = 65
        nnorm = xanesnorm
    post2 = group.energy[-1] - group.e0
    pre_edge(group = group,energy = group.energy, mu = group.mu, e0 = group.e0, pre1=pre1,pre2=pre2,
             norm1 = post1, norm2=post2, nnorm = nnorm)
    return group

def run(direc):
    if not os.path.exists(direc):
        return
    os.chdir(direc)

    fluorescenceCounter = 'xmap_roi00'
    monPattern = 'mon_'
    ion1Pattern = 'ion_1'
    muFheader = 'muF'
    muTheader = 'muT'
    
    for root,dirs,files in os.walk(os.getcwd()):
        if not 'regrid' in root or 'merge' in root or 'norm' in root:
            continue
        os.chdir(root)
        if not os.path.exists('merge/'):
            os.makedirs('merge/')
        if not os.path.exists('norm/trans'):
            os.makedirs('norm/trans')
        if not os.path.exists('norm/fluo'):
            os.makedirs('norm/fluo')
        datfiles = [file for file in files if file.endswith('.dat')]
        datfiles.sort()
        if len(datfiles) == 0:
            continue
        print(root)
        Emins = np.array([])
        Emaxs = np.array([])
        dfmergedct = {}
        fluoList = []
        for c,file in enumerate(datfiles):
            
            f = open(file,'r')
            header = ''.join(f.readlines()[:2]).replace('#','')
            f.close()
            if c == 0:
                df0 = pd.read_csv(file,sep = ' ',comment = '#',index_col = 0)

                if muTheader in df0.columns:
                    transmission = True
                else:
                    transmission = False

            df = pd.read_csv(file,sep = ' ',comment = '#',index_col = 0, header = 0)
            if muFheader in df.columns:
                fluorescence = True
            else:
                fluorescence = False
            fluoList.append(fluorescence)
            dfmergedct[file] = pd.DataFrame()
            Emins = np.append(Emins,df.index.values[0])
            Emaxs = np.append(Emaxs,df.index.values[-1])
            dfmergedct[file]['energy_offset(keV)'] = df.index.values
            E = df.index.values
            if transmission:
                muT = df[muTheader].values
                dfmergedct[file]['muT'] = muT
            if fluorescence:
                muFluo = df[muFheader].values
                dfmergedct[file]['muF'] = muFluo
      
        E0merge = np.max(Emins)
        EendMerge = np.min(Emaxs)
        dfMergeT = pd.DataFrame()
        dfMergeF = pd.DataFrame()
        for c,file in enumerate(dfmergedct):
            basefileT = file.replace('.dat','T')
            basefileF = file.replace('.dat','F')
            E = dfmergedct[file]['energy_offset(keV)'].values
            minindex = np.abs(E - E0merge).argmin()
            maxindex = np.abs(E - EendMerge).argmin()

            if fluoList[c] == True:
                muFluo = dfmergedct[file]['muF'].loc[minindex:].values #making individual files start with same E value to make plotting easier
                ds = pd.Series(index = E[minindex:],data = muFluo)
                groupF = normalise(ds)
                fileF = f'norm/fluo/{basefileF}.nor'
                print(fileF)
                np.savetxt(fileF,np.array([E[minindex:],groupF.flat]).transpose(),header = f'{header}Energy(keV) mu_norm',fmt = '%.5f')
                if c == 0:
                    dfMergeF['energy_offset(keV)'] = dfmergedct[file]['energy_offset(keV)'].loc[minindex:maxindex].values
                    dfMergeF = dfMergeF.set_index('energy_offset(keV)')
                dfMergeF[c] = dfmergedct[file]['muF'].loc[minindex:maxindex].values
            if transmission:
                muT =  dfmergedct[file]['muT'].loc[minindex:].values
                ds = pd.Series(index = E[minindex:],data = muT)
                groupT = normalise(ds)
                fileT = f'norm/trans/{basefileT}.nor'
                print(fileT)
                np.savetxt(fileT,np.array([E[minindex:],groupT.flat]).transpose(),header = f'{header}Energy(keV) mu_norm',fmt = '%.5f')
                if c == 0:
                    dfMergeT['energy_offset(keV)'] = dfmergedct[file]['energy_offset(keV)'].loc[minindex:maxindex].values
                    dfMergeT = dfMergeT.set_index('energy_offset(keV)')
                dfMergeT[c] = dfmergedct[file]['muT'].loc[minindex:maxindex].values

        if transmission:
            dfmergeTrans = dfMergeT.mean(axis = 1)
            dfmergeTrans.name = 'mu'
            basefileTmerge = re.sub('[0-9][0-9][0-9][0-9].dat','T',file)
            dfmergeTrans.to_csv(f'merge/{basefileTmerge}_merge.dat',sep = ' ')
            if np.max(dfmergeTrans.values) - np.min(dfmergeTrans.values) > 0.1: #checking to see if data is real
                groupTmerge = normalise(dfmergeTrans)
                fileTmerge = f'merge/{basefileTmerge}.nor'
                np.savetxt(fileTmerge,np.array([groupTmerge.energy,groupTmerge.flat]).transpose(),header = '#Energy(keV) mu_norm',fmt = '%.5f')
        if True in fluoList:
            dfmergeFluo = dfMergeF.mean(axis = 1)
            dfmergeFluo.name = 'mu'
            basefileFmerge = re.sub('[0-9][0-9][0-9][0-9].dat','F',file)
            dfmergeFluo.to_csv(f'merge/{basefileFmerge}_merge.dat', sep = ' ')
            if np.max(dfmergeFluo.values) - np.min(dfmergeFluo.values) > 0.1:
                groupFmerge = normalise(dfmergeFluo)
                fileFmerge = f'merge/{basefileFmerge}.nor'
                np.savetxt(fileFmerge,np.array([groupFmerge.energy,groupFmerge.flat]).transpose(),header = '#Energy(keV) mu_norm',fmt = '%.5f')
if __name__ == '__main__':
    run(direc = direc)