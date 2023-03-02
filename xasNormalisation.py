from larch.xafs import autobk, find_e0, pre_edge
from larch import Group
import numpy as np
import os
from glob import glob
import pandas as pd

direc = r'Z:\visitor\ch6616\bm31\20230228/'
os.chdir(direc)

highEmono = 'mon_3'
lowEmono = 'mon_4'
highEi1 = 'ion_1_3'
lowEi1 = 'ion_1_2'
highE = 17.5
def normalise(ds,group):
    group.energy = ds.index.values
    group.mu = ds.values

    find_e0(group = group, energy = group.energy, mu = group.mu)
    pre1 = group.energy[0] - group.e0
    pre2 = -0.03
    if group.energy[-1] - group.energy[0] > 0.5: #EXAFS
        post1 = 0.15
        nnorm = 2
    else: #XANES
        post1 = 0.065
        nnorm = 1
    post2 = group.energy[-1] - group.e0
    pre_edge(group = group,energy = group.energy, mu = group.mu, e0 = group.e0, pre1=pre1,pre2=pre2,norm1 = post1, norm2=post2, nnorm = nnorm)


for root,dirs,files in os.walk(os.getcwd()):
    if 'regrid' in root and not 'merge' in root and not 'norm' in root:

        os.chdir(root)
        if not os.path.exists('merge/'):
            os.makedirs('merge/')
        if not os.path.exists('norm/'):
            os.makedirs('norm/')

        datfiles = [file for file in files if file.endswith('.dat')]
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

                if 'xmap_roi00' in df0.columns:
                    dfFluo = pd.DataFrame()
                    dfFluo.index = df0.index

                dfTrans = pd.DataFrame()
                dfTrans.index = df0.index
            df = pd.read_csv(file,sep = '\s',comment = '#',index_col = 0)

            E = df.index.values
            if E[0] >= highE:
                mon_counter = highEmono
                i1_counter = highEi1
            else:
                mon_counter = lowEmono
                i1_counter = lowEi1
            if c != 0 and len(df.index.values) < len(dfTrans.index.values): #removing data points if data sizes don't match for averaging

                E0 = dfTrans.index.values
                mismatchsize = len(E0) - len(E)
                print(f'{root} - mismatched data size')
                print(f'cutting data by {mismatchsize} points to fit')
                for n in range(mismatchsize):
                    if E[-1] != E0[-1]:
                        dfTrans.drop(E0[-1],axis = 0,inplace = True)
                        if 'xmap_roi00' in df0.columns:
                            dfFluo.drop(E0[-1],axis = 0,inplace = True)
                    elif E[0] != E0[0]:
                        dfTrans.drop(E0[0],axis = 0,inplace=True)
                        if 'xmap_roi00' in df0.columns:
                            dfFluo.drop(E0[0],axis = 0,inplace = True)
                    E0 = dfTrans.index.values
            elif c != 0 and len(df.index.values) > len(dfTrans.index.values): #removing data points if data sizes don't match
                
                E0 = dfTrans.index.values
                mismatchsize = len(E) - len(E0)
                print(f'{root} - mismatched data size')
                print(f'cutting data by {mismatchsize} points to fit')
                for n in range(mismatchsize):
                    if E[-1] != E0[-1]:
                        df.drop(E0[-1],axis = 0,inplace = True)

                    elif E[0] != E0[0]:
                        df.drop(E0[0],axis = 0,inplace=True)
                    E = df.index.values
      
            if 'xmap_roi00' in df0.columns:
        
                muFluo = df['xmap_roi00'].values/df[mon_counter].values
                dfFluo[c] = muFluo
                ds = pd.Series(index = dfFluo.index,data = muFluo)

                if np.max(df[mon_counter].values) > 500 and np.max(df['xmap_roi00'].values) > 10:
                    groupF = Group()
                    normalise(ds,groupF)
                    basefileF = file.replace('.dat','F.nor')
                    fileF = f'norm/{basefileF}'

                    np.savetxt(fileF,np.array([E,groupF.flat]).transpose(),header = f'{header}Energy(keV) mu_norm',fmt = '%.5f')


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
        dfmergeTrans = dfTrans.mean(axis = 1)
        dfmergeTrans.name = 'mu'
        dfmergeTrans.to_csv('merge/transMerge.dat',sep = ' ')
        if np.max(dfmergeTrans.values) - np.min(dfmergeTrans.values) > 0.1: #checking to see if data is real
            groupTmerge = Group()
            normalise(dfmergeTrans,groupTmerge)
            basefileTmerge = file.replace('.dat','T.nor')
            fileTmerge = f'merge/{basefileTmerge}'
            np.savetxt(fileTmerge,np.array([groupTmerge.energy,groupTmerge.flat]).transpose(),header = '#Energy(keV) mu_norm',fmt = '%.5f')
        if 'xmap_roi00' in df0.columns:
            dfmergeFluo = dfFluo.mean(axis = 1)
            dfmergeFluo.name = 'mu'
            dfmergeFluo.to_csv('merge/fluoMerge.dat', sep = ' ')
            if np.max(dfmergeFluo.values) - np.min(dfmergeFluo.values) > 0.1:
                groupFmerge = Group()
                normalise(dfmergeFluo,groupFmerge)
                basefileFmerge = file.replace('.dat','F.nor')
                fileFmerge = f'merge/{basefileFmerge}'
                np.savetxt(fileFmerge,np.array([groupFmerge.energy,groupFmerge.flat]).transpose(),header = '#Energy(keV) mu_norm',fmt = '%.5f')
