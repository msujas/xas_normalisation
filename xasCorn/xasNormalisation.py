from larch.xafs import find_e0, pre_edge
from larch import Group
import numpy as np
import os, re
import pandas as pd
from glob import glob
from functools import partial


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

savenorm = partial(np.savetxt, fmt = '%.5f', comments = '#')

def normaliseRG(regriddir, unit = 'keV'):
    if not 'regrid' in regriddir or 'norm' in regriddir:
        return
    normdir = f'{regriddir}/norm'
    files = glob(f'{regriddir}/*.dat')
    if not files:
        return
    os.makedirs(normdir,exist_ok=True)
    for file in files:
        print(file)
        f = open(file,'r')
        header = [line.replace('#','') for line in f.readlines() if line.startswith('#')]
        f.close()
        columns = header[-1].replace('\n','').split()
        energycol = columns[0]
        header = ''.join(header[:-1])
        df = pd.read_csv(file,sep = ' ', header= None, comment='#')
        df.columns = columns
        energy = df[energycol].values
        filen = os.path.basename(file.replace('.dat','.nor'))
        mucol = [col for col in columns if col == 'muT' or col == 'muF1'][0]

        mu = df[mucol]
        mu.index = energy
        try:
            group = normalise(mu)
            header = header + f'edge: {group.e0}\n'
            header += f'edge step: {group.edge_step}\n'
            header += f'{columns[0]} {mucol}norm'    
            savenorm(f'{normdir}/{filen}', np.array([energy,group.flat]).transpose(),header=header)
        except AttributeError:
            print(f'couldn\'t normalise {mucol} for {file}')

    mergeFiles = glob(f'{regriddir}/merge/*.dat')
    for file in mergeFiles:
        data = np.loadtxt(file,unpack=True,comments='#')
        energy = data[0]
        mu = data[1]
        ds = pd.Series(data = mu, index = energy)
        try:
            groupMerge = normalise(ds)
            header = f'edge: {groupMerge.e0}\n'
            header += f'edge step: {groupMerge.edge_step}\n'
            header += f'energy({unit}) mu_norm'
            savenorm(file.replace('.dat','.nor'), np.array([energy,groupMerge.flat]).transpose(), 
                    header=header)
        except AttributeError:
            print(f'couldn\'t normalise {file}')


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
                if f'{e}_xanes' in root or f'{e}_exafs' in root:
                    skip = False
                    break 
            if skip:
                continue
        elif excludeElements:
            skip = False
            for e in excludeElements:
                if f'{e}_xanes' in root or f'{e}_exafs' in root:
                    skip = True
                    break
            if skip:
                continue
        print(root)
        normaliseRG(root, unit)

      
