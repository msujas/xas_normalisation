from larch.xafs import find_e0, pre_edge
from larch import Group
import numpy as np
import os, re
import pandas as pd

direc = r'C:\Users\kenneth1a\Documents\beamlineData\a311231\Fl-CrFoil_Cr_exafs\regrid'



def normalise(ds,group, kev = True):
    '''
    The normalisation orders seem to be different between Athena and Larch. 
    1 and 2 in Larch seem to correspond to 2 and 3 in Athena (linear and quadratic), respectively. 0 Seems not to correspond to 1, however.
    0 appears linear, but with a shallower gradient than 1. The documentation recommends 0 if < 50 eV used to fit post-edge.
    The values used in this are roughly the same as the defaults
    Some distributions of Larch don't normalise data properly in keV, so data is converted to eV for normalising
    '''
    group.energy = ds.index.values
    group.mu = ds.values

    find_e0(group = group, energy = group.energy, mu = group.mu)
    pre1 = group.energy[0] - group.e0
    pre2 = -0.02
    if group.energy[-1] - group.energy[0] > 0.5: #EXAFS
        post1 = 0.15
        nnorm = 2 
    else: #XANES
        post1 = 0.065
        nnorm = 1
    post2 = group.energy[-1] - group.e0
    if kev: #converting axis to eV
        scale = 1000
    else:
        scale = 1
    pre_edge(group = group,energy = group.energy*scale, mu = group.mu, e0 = group.e0*scale, pre1=pre1*scale,pre2=pre2*scale,
             norm1 = post1*scale, norm2=post2*scale, nnorm = nnorm)

def run(direc):
    if not os.path.exists(direc):
        return
    os.chdir(direc)

    fluorescenceCounter = 'xmap_roi00'
    monPattern = 'mon_'
    ion1Pattern = 'ion_1'
    
    for root,dirs,files in os.walk(os.getcwd()):
        if not 'regrid' in root or 'merge' in root or 'norm' in root:
            continue
        os.chdir(root)
        if not os.path.exists('merge/'):
            os.makedirs('merge/')
        if not os.path.exists('norm/'):
            os.makedirs('norm/')

        datfiles = [file for file in files if file.endswith('.dat')]
        datfiles.sort()
        if len(datfiles) == 0:
            continue
        print(root)
        for c,file in enumerate(datfiles):
            print(file)
            f = open(file,'r')
            header = ''.join(f.readlines()[:2]).replace('#','')
            f.close()
            if c == 0:
                df0 = pd.read_csv(file,sep = '\s',comment = '#',index_col = 0)
                if fluorescenceCounter in df0.columns:
                    dfFluo = pd.DataFrame()
                    dfFluo.index = df0.index
                if len([col for col in df0.columns if ion1Pattern in col]) == 0:
                    transmission = False
                else:
                    transmission = True
                    i1_counter = [col for col in df0.columns if ion1Pattern in col][0]
                dfTrans = pd.DataFrame()
                dfTrans.index = df0.index
            df = pd.read_csv(file,sep = '\s',comment = '#',index_col = 0, header = 0)
            mon_counter = [col for col in df.columns if monPattern in col][0]

                
            E = df.index.values

            if c != 0 and len(df.index.values) < len(dfTrans.index.values): #removing data points if data sizes don't match for averaging

                E0 = dfTrans.index.values
                mismatchsize = len(E0) - len(E)
                print(f'mismatched data size')
                print(f'cutting data by {mismatchsize} points to fit')
                for n in range(mismatchsize):
                    if E[-1] != E0[-1]:
                        dfTrans.drop(E0[-1],axis = 0,inplace = True)
                        if fluorescenceCounter in df0.columns:
                            dfFluo.drop(E0[-1],axis = 0,inplace = True)
                    elif E[0] != E0[0]:
                        dfTrans.drop(E0[0],axis = 0,inplace=True)
                        if fluorescenceCounter in df0.columns:
                            dfFluo.drop(E0[0],axis = 0,inplace = True)
                    E0 = dfTrans.index.values
            elif c != 0 and len(df.index.values) > len(dfTrans.index.values): #removing data points if data sizes don't match
                
                E0 = dfTrans.index.values
                mismatchsize = len(E) - len(E0)
                print(f'mismatched data size')
                print(f'cutting data by {mismatchsize} points to fit')
                for n in range(mismatchsize):
                    if E[-1] != E0[-1] and len(E0):
                        df.drop(E[-1],axis = 0,inplace = True)

                    elif E[0] != E0[0] and len(E0):
                        df.drop(E[0],axis = 0,inplace=True)

                    E = df.index.values
    
            if fluorescenceCounter in df.columns:
        
                muFluo = df[fluorescenceCounter].values/df[mon_counter].values
                dfFluo[c] = muFluo
                ds = pd.Series(index = dfFluo.index,data = muFluo)

                if np.max(df[mon_counter].values) > 500 and np.max(df[fluorescenceCounter].values) > 10: #filtering out bad data
                    groupF = Group()
                    normalise(ds,groupF)
                    basefileF = file.replace('.dat','F.nor')
                    fileF = f'norm/{basefileF}'

                    np.savetxt(fileF,np.array([E,groupF.flat]).transpose(),header = f'{header}Energy(keV) mu_norm',fmt = '%.5f')

            if transmission:
                try:
                    muTrans = np.log(df[mon_counter].values/df[i1_counter].values)
                    dfTrans[c] = muTrans

                    if np.max(df[mon_counter].values) < 500:
                        continue
                    ds = pd.Series(index = dfTrans.index,data = muTrans)
                    groupT = Group()
                    normalise(ds,groupT)
                    basefileT = file.replace('.dat','T.nor')
                    fileT = f'norm/{basefileT}'
                    np.savetxt(fileT,np.array([E,groupT.flat]).transpose(),header = f'{header}Energy(keV) mu_norm',fmt = '%.5f')
                except KeyError:
                    print(f'transmission too low in {file}')
        if transmission:
            dfmergeTrans = dfTrans.mean(axis = 1)
            dfmergeTrans.name = 'mu'
            dfmergeTrans.to_csv('merge/transMerge.dat',sep = ' ')
            if np.max(dfmergeTrans.values) - np.min(dfmergeTrans.values) > 0.1: #checking to see if data is real
                groupTmerge = Group()
                normalise(dfmergeTrans,groupTmerge)
                basefileTmerge = re.sub('[0-9][0-9][0-9][0-9].dat','T.nor',file)
                fileTmerge = f'merge/{basefileTmerge}'
                np.savetxt(fileTmerge,np.array([groupTmerge.energy,groupTmerge.flat]).transpose(),header = '#Energy(keV) mu_norm',fmt = '%.5f')
        if fluorescenceCounter in df0.columns:
            dfmergeFluo = dfFluo.mean(axis = 1)
            dfmergeFluo.name = 'mu'
            dfmergeFluo.to_csv('merge/fluoMerge.dat', sep = ' ')
            if np.max(dfmergeFluo.values) - np.min(dfmergeFluo.values) > 0.1:
                groupFmerge = Group()
                normalise(dfmergeFluo,groupFmerge)
                basefileFmerge = re.sub('[0-9][0-9][0-9][0-9].dat','F.nor',file)
                fileFmerge = f'merge/{basefileFmerge}'
                np.savetxt(fileFmerge,np.array([groupFmerge.energy,groupFmerge.flat]).transpose(),header = '#Energy(keV) mu_norm',fmt = '%.5f')
if __name__ == '__main__':
    run(direc = direc)