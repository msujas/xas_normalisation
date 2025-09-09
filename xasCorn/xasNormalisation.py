from larch.xafs import find_e0, pre_edge
from larch import Group
import numpy as np
import os, re
import pandas as pd
from glob import glob

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
    try:
        find_e0(group = group, energy = group.energy, mu = group.mu)
    except ValueError as e:
        print('couldn\'t normalise file')
        raise ValueError
        

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


def normaliseRG(regriddir, unit = 'keV'):
    if not 'regrid' in regriddir or 'norm' in regriddir:
        return
    os.makedirs(f'{regriddir}/norm',exist_ok=True)
    transdir = f'{regriddir}/norm/trans/'
    fluodir = f'{regriddir}/norm/fluo'
    files = glob(f'{regriddir}/*.dat')
    energycol = f'#energy_offset({unit})'
    for file in files:
        print(file)
        f = open(file,'r')
        header = [line for line in f.readlines() if line.startswith('#')]
        f.close()
        columns = header[-1].replace('\n','').split()
        header = ''.join(header[:-1])
        df = pd.read_csv(file,sep = ' ', header= None, comment='#')
        df.columns = columns
        energy = df[energycol].values
        filen = os.path.basename(file.replace('.dat','.nor'))
        if 'muT' in columns:
            if not os.path.exists(transdir):
                os.makedirs(transdir)
            muT = df['muT']
            muT.index = energy
            groupT = normalise(muT)
            headerT = header + f'\n#edge: {groupT.e0}\n'
            headerT += f'#edge step: {groupT.edge_step}\n'
            headerT += ' '.join(columns)     
            np.savetxt(f'{regriddir}/norm/trans/{filen}', np.array([energy,groupT.flat]).transpose(),header=headerT, fmt = '%.5f')
        if 'muF1' in columns:
            if not os.path.exists(fluodir):
                os.makedirs(fluodir)
            muF = df['muF1']
            muF.index = energy
            groupF = normalise(muF)
            headerF = f'{header}\n#edge: {groupF.e0}\n'
            headerF += f'#edge step: {groupF.edge_step}\n'
            headerF += ' '.join(columns)
            np.savetxt(f'{regriddir}/norm/fluo/{filen}',np.array([energy,groupF.flat]).transpose(),header=headerF, fmt = '%.5f')
    mergeFiles = glob(f'{regriddir}/merge/*.dat')
    for file in mergeFiles:
        data = np.loadtxt(file,unpack=True,comments='#')
        energy = data[0]
        mu = data[1]
        ds = pd.Series(data = mu, index = energy)
        groupMerge = normalise(ds)
        header = f'edge: {groupMerge.e0}\n'
        header += f'edge step: {groupMerge.edge_step}\n'
        header += f'energy({unit}) mu_norm'
        np.savetxt(file.replace('.dat','.nor'), np.array([energy,groupMerge.flat]).transpose(), header=header, fmt = '%.5f')

def run(direc, unit = 'keV', coldirname = 'columns', elements = None, excludeElements = None, averaging = 1):
    if not os.path.exists(direc):
        return
    os.chdir(direc)
    for root,dirs,files in os.walk(os.getcwd()):
        if not 'regrid' in root or coldirname not in root or 'merge' in root or 'norm' in root:
            continue
        if averaging < 2 and 'regridAv' in root:
            continue
        elif averaging > 1 and 'regridAv' in root and f'regridAv{averaging}' not in root:
            continue 
        if coldirname == 'columns' and ('columns-' in root or re.search('colummns[0-9]',root)):
            continue
        if elements:
            skip = True
            for e in elements:
                if f'_{e}_' in root:
                    skip = False
                    break 
            if skip:
                continue
        elif excludeElements:
            skip = False
            for e in excludeElements:
                if f'_{e}_' in root:
                    skip = True
                    break
            if skip:
                continue
        print(root)
        normaliseRG(root, unit)

def run2(direc, unit = 'keV', coldirname = 'columns', elements = None, excludeElements = None, averaging = 1):
    if not os.path.exists(direc):
        return
    os.chdir(direc)

    fluorescenceCounter = 'xmap_roi00'
    monPattern = 'mon_'
    ion1Pattern = 'ion_1'
    muFheaders = ['muF1','muF2']
    muTheader = 'muT'
    mu2header = 'mu2'
    
    for root,dirs,files in os.walk(os.getcwd()):
        if not 'regrid' in root or coldirname not in root or 'merge' in root or 'norm' in root:
            continue
        if averaging < 2 and 'regridAv' in root:
            continue
        elif averaging > 1 and 'regridAv' in root and f'regridAv{averaging}' not in root:
            continue 

        if coldirname == 'columns' and ('columns-' in root or re.search('colummns[0-9]',root)):
            continue
        if elements:
            skip = True
            for e in elements:
                if f'_{e}_' in root:
                    skip = False
                    break 
            if skip:
                continue
        elif excludeElements:
            skip = False
            for e in excludeElements:
                if f'_{e}_' in root:
                    skip = True
                    break
            if skip:
                continue
        
        os.chdir(root)
        if not os.path.exists('merge/'):
            os.makedirs('merge/')

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
        mu2list = []
        headers = []
        usedMuFList = []
        for c,file in enumerate(datfiles):
            
            f = open(file,'r')
            header = [line for line in f.readlines() if '#' in line]
            colnames = header[-1].split()
            header = header[:-1]
            header = ''.join(header).replace('#','')
            f.close()
            headers.append(header)
            fluorescence = False
            transmission = False
            i2 = False
            df = pd.read_csv(file,sep = ' ',comment = '#',index_col = 0, names = colnames)
            usedMuF = [col for col in df.columns if col in muFheaders]
            usedMuFList.append(usedMuF)
            if usedMuF:
                fluorescence = True
        
            '''
            for c,muFheader in enumerate(usedMuF):
                if muFheader in df.columns:
                    #values = df[fluorescenceCounter].values
                    #fluorescence[c] = not(np.inf in values or np.max(values) < 100)
                    fluorescence[c] = True
            '''
        
            if muTheader in df.columns:
                values = df[muTheader].values
                if not np.inf in values:
                    transmission = True
            if mu2header in df.columns:
                values = df[mu2header].values
                if not np.inf in values:
                    i2=True
                
            fluoList.append(fluorescence)
            transmissionList.append(transmission)
            mu2list.append(i2)
            dfmergedct[file] = pd.DataFrame()
            Emins = np.append(Emins,df.index.values[0])
            Emaxs = np.append(Emaxs,df.index.values[-1])
            dfmergedct[file][f'energy_offset({unit})'] = df.index.values
            E = df.index.values
            if transmission:
                muT = df[muTheader].values
                dfmergedct[file]['muT'] = muT
            if fluorescence:
                muFluo = {}
                for c,muFheader in enumerate(usedMuF):
                    muFluo[c] = df[muFheader].values
                    dfmergedct[file][muFheader] = muFluo[c]
            if i2:
                mu2 =  df[mu2header].values
                dfmergedct[file]['mu2'] = mu2
      
        E0merge = np.max(Emins)
        EendMerge = np.min(Emaxs)
        dfMergeT = pd.DataFrame()
        dfMergeF = pd.DataFrame()
        dfMergeMu2 = pd.DataFrame()
        doFluo = True in fluoList
        if True in transmissionList:
            if not os.path.exists('norm/trans'):
                os.makedirs('norm/trans')
        if doFluo:
            if not os.path.exists('norm/fluo'):
                os.makedirs('norm/fluo')
        
        for c,file in enumerate(dfmergedct):
            basefileT = file.replace('.dat','T')
            basefileF = file.replace('.dat','F')
            E = dfmergedct[file][f'energy_offset({unit})'].values
            minindex = np.abs(E - E0merge).argmin()
            maxindex = np.abs(E - EendMerge).argmin()

            for n,fluoName in enumerate(usedMuFList[c]):
                muFluo = dfmergedct[file][fluoName].loc[minindex:].values #making individual files start with same E value to make plotting easier
                ds = pd.Series(index = E[minindex:],data = muFluo)
                try:
                    groupF = normalise(ds)
                except ValueError:
                    print(f'couldn\'t normalise {basefileF}.nor')
                    continue
                e0 = groupF.e0
                edgeStep = groupF.edge_step
                fileF = f'norm/fluo/{basefileF}.nor'
                print(fileF)
                if n == 0:
                    spectrumHeader = f'{headers[c]}edge: {e0} eV\nedge step: {edgeStep}\nEnergy({unit})'
                    array = np.array([E[minindex:]])
                array = np.append(array,[groupF.flat], axis = 0)
                spectrumHeader += f' {fluoName}_norm'
                if n == len(usedMuFList[c]) -1:
                    np.savetxt(fileF,np.array([E[minindex:],groupF.flat]).transpose(),header = spectrumHeader,fmt = '%.5f')
                if c == 0:
                    dfMergeF[f'energy_offset({unit})'] = dfmergedct[file][f'energy_offset({unit})'].loc[minindex:maxindex].values
                    dfMergeF = dfMergeF.set_index(f'energy_offset({unit})')
                dfMergeF[c] = dfmergedct[file][fluoName].loc[minindex:maxindex].values
            if transmissionList[c]:
                muT =  dfmergedct[file]['muT'].loc[minindex:].values
                ds = pd.Series(index = E[minindex:],data = muT)
                try:
                    groupT = normalise(ds)
                except ValueError:
                    print(f'couldn\'t normalise {basefileT}.nor')
                    continue
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
            if mu2list[c]:
                if c == 0:
                    dfMergeMu2[f'energy_offset({unit})'] = dfmergedct[file][f'energy_offset({unit})'].loc[minindex:maxindex].values 
                    dfMergeMu2 = dfMergeMu2.set_index(f'energy_offset({unit})')
                dfMergeMu2[c] = dfmergedct[file]['mu2'].loc[minindex:maxindex].values
        
        if True in transmissionList:
            dfmergeTrans = dfMergeT.mean(axis = 1)
            dfmergeTrans.name = 'mu'
            basefileTmerge = re.sub('[0-9][0-9][0-9][0-9].dat','T',file)
            #dfmergeTrans.to_csv(f'merge/{basefileTmerge}_merge.dat',sep = ' ')
            groupTmerge = normalise(dfmergeTrans)
            e0 = groupTmerge.e0
            edgeStep = groupTmerge.edge_step
            fileTmerge = f'merge/{basefileTmerge}.nor'
            spectrumHeader = f'edge: {e0} eV\nedge step: {edgeStep}\nEnergy({unit}) mu_norm'
            np.savetxt(fileTmerge,np.array([dfmergeTrans.index.values,groupTmerge.flat]).transpose(),header = spectrumHeader,fmt = '%.5f')
        
        if doFluo:
            dfmergeFluo = dfMergeF.mean(axis = 1)
            dfmergeFluo.name = 'mu'
            basefileFmerge = re.sub('[0-9][0-9][0-9][0-9].dat','F',file)
            #dfmergeFluo.to_csv(f'merge/{basefileFmerge}_merge.dat', sep = ' ')
            groupFmerge = normalise(dfmergeFluo)
            fileFmerge = f'merge/{basefileFmerge}.nor'
            e0 = groupFmerge.e0
            edgeStep = groupFmerge.edge_step
            spectrumHeader = f'edge: {e0} eV\nedge step: {edgeStep}\nEnergy({unit}) mu_norm'
            np.savetxt(fileFmerge,np.array([dfmergeFluo.index.values,groupFmerge.flat]).transpose(),header = spectrumHeader,fmt = '%.5f')
        if True in mu2list:
            dfMergeMu2 = dfMergeMu2.mean(axis=1)
            dfMergeMu2.name = 'mu2'
            basefileMu2 = re.sub('[0-9][0-9][0-9][0-9].dat','mu2',file)
            dfMergeMu2.to_csv(f'merge/{basefileMu2}_merge.dat', sep = ' ')
        
if __name__ == '__main__':
    run(direc = direc)