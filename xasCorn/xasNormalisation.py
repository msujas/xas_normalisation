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
    pre2Function = -(group.e0*0.0033 + 4.04) #-30
    if pre2Function - pre1 > 30:
        pre2 = pre2Function #-30
    else:
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


def run(direc, unit = 'keV', coldirname = 'columns'):
    if not os.path.exists(direc):
        return
    os.chdir(direc)

    fluorescenceCounter = 'xmap_roi00'
    monPattern = 'mon_'
    ion1Pattern = 'ion_1'
    muFheader = 'muF'
    muTheader = 'muT'
    
    for root,dirs,files in os.walk(os.getcwd()):
        if not 'regrid' in root or coldirname not in root or 'merge' in root or 'norm' in root:
            continue
        if coldirname == 'columns' and ('columns-' in root or re.search('colummns[0-9]',root)):
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
        transmissionList = []
        headers = []
        for c,file in enumerate(datfiles):
            
            f = open(file,'r')
            header = ''.join(f.readlines()[:2]).replace('#','')
            f.close()
            headers.append(header)
            fluorescence = False
            transmission = False
            df = pd.read_csv(file,sep = ' ',comment = '#',index_col = 0, header = 0)
            
            if muFheader in df.columns:
                values = df[fluorescenceCounter].values
                fluorescence = not(np.inf in values or np.max(values) < 100)

            if muTheader in df.columns:
                values = df[muTheader].values
                if not np.inf in values:
                    transmission = True
                
            fluoList.append(fluorescence)
            transmissionList.append(transmission)
            dfmergedct[file] = pd.DataFrame()
            Emins = np.append(Emins,df.index.values[0])
            Emaxs = np.append(Emaxs,df.index.values[-1])
            dfmergedct[file][f'energy_offset({unit})'] = df.index.values
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
            E = dfmergedct[file][f'energy_offset({unit})'].values
            minindex = np.abs(E - E0merge).argmin()
            maxindex = np.abs(E - EendMerge).argmin()

            if fluoList[c]:
                muFluo = dfmergedct[file]['muF'].loc[minindex:].values #making individual files start with same E value to make plotting easier
                ds = pd.Series(index = E[minindex:],data = muFluo)
                groupF = normalise(ds)
                e0 = groupF.e0
                edgeStep = groupF.edge_step
                fileF = f'norm/fluo/{basefileF}.nor'
                print(fileF)
                spectrumHeader = f'{headers[c]}edge: {e0} eV\nedge step: {edgeStep}\nEnergy({unit}) mu_norm'
                np.savetxt(fileF,np.array([E[minindex:],groupF.flat]).transpose(),header = spectrumHeader,fmt = '%.5f')
                if c == 0:
                    dfMergeF[f'energy_offset({unit})'] = dfmergedct[file][f'energy_offset({unit})'].loc[minindex:maxindex].values
                    dfMergeF = dfMergeF.set_index(f'energy_offset({unit})')
                dfMergeF[c] = dfmergedct[file]['muF'].loc[minindex:maxindex].values
            if transmissionList[c]:
                muT =  dfmergedct[file]['muT'].loc[minindex:].values
                ds = pd.Series(index = E[minindex:],data = muT)
                groupT = normalise(ds)
                fileT = f'norm/trans/{basefileT}.nor'
                if len(groupT.flat) != len(E[minindex:]):
                    print(f'normalisation for {fileT} cut some data, skipping')
                    continue
                e0 = groupT.e0
                edgeStep = groupT.edge_step
                
                spectrumHeader = f'{headers[c]}edge: {e0} eV\nedge step: {edgeStep}\nEnergy({unit}) mu_norm'
                print(fileT)
                np.savetxt(fileT,np.array([E[minindex:],groupT.flat]).transpose(),header = spectrumHeader,fmt = '%.5f')
                if c == 0:
                    dfMergeT[f'energy_offset({unit})'] = dfmergedct[file][f'energy_offset({unit})'].loc[minindex:maxindex].values
                    dfMergeT = dfMergeT.set_index(f'energy_offset({unit})')
                dfMergeT[c] = dfmergedct[file]['muT'].loc[minindex:maxindex].values

        if True in transmissionList:
            dfmergeTrans = dfMergeT.mean(axis = 1)
            dfmergeTrans.name = 'mu'
            basefileTmerge = re.sub('[0-9][0-9][0-9][0-9].dat','T',file)
            dfmergeTrans.to_csv(f'merge/{basefileTmerge}_merge.dat',sep = ' ')
            groupTmerge = normalise(dfmergeTrans)
            e0 = groupTmerge.e0
            edgeStep = groupTmerge.edge_step
            fileTmerge = f'merge/{basefileTmerge}.nor'
            spectrumHeader = f'edge: {e0} eV\nedge step: {edgeStep}\nEnergy({unit}) mu_norm'
            np.savetxt(fileTmerge,np.array([dfmergeTrans.index.values,groupTmerge.flat]).transpose(),header = spectrumHeader,fmt = '%.5f')

        if True in fluoList:
            dfmergeFluo = dfMergeF.mean(axis = 1)
            dfmergeFluo.name = 'mu'
            basefileFmerge = re.sub('[0-9][0-9][0-9][0-9].dat','F',file)
            dfmergeFluo.to_csv(f'merge/{basefileFmerge}_merge.dat', sep = ' ')
            groupFmerge = normalise(dfmergeFluo)
            fileFmerge = f'merge/{basefileFmerge}.nor'
            e0 = groupFmerge.e0
            edgeStep = groupFmerge.edge_step
            spectrumHeader = f'edge: {e0} eV\nedge step: {edgeStep}\nEnergy({unit}) mu_norm'
            np.savetxt(fileFmerge,np.array([dfmergeFluo.index.values,groupFmerge.flat]).transpose(),header = spectrumHeader,fmt = '%.5f')
            
if __name__ == '__main__':
    run(direc = direc)