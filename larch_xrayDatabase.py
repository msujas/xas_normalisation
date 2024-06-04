from larch.xray import xray_edges, atomic_symbol
import pandas as pd
import numpy as np
import os
#os.chdir(r'C:\Users\kenneth1a\Documents\beamlineData/')
edges = ['K','L1','L2','L3','main']
edgesdf = pd.DataFrame(columns = edges)
for n in range(22,84):
    element = atomic_symbol(n)
    elementEdges = []
    mainEdge = 0
    for e in edges[:-1]:
        edgeValue = np.round(xray_edges(n)[e][0]/1000,4)
        if e == 'K' and edgeValue < 35:
            mainEdge = edgeValue
        if edgeValue > 5 and mainEdge == 0 and e == 'L3':
            mainEdge = edgeValue
        
        elementEdges.append(edgeValue)
    elementEdges.append(mainEdge)
    edgesdf.loc[element] = elementEdges
print(edgesdf)
edgesdf.to_csv('xrayEdges.dat', sep = '\t', header=True, index=True)

xanes_programDF = pd.DataFrame(columns = ['start','stop','step'])
exafs_programDF = pd.DataFrame(columns = ['start','stop','step'])
xanes_programDFside = pd.DataFrame(columns = ['start','stop','step'])
exafs_programDFside = pd.DataFrame(columns = ['start','stop','step'])
def xanes_step(startE):
    return (0.0186*startE + 0.2145).round(1)
def exafs_step(startE):
    return (0.0234*startE + 0.3754).round(1)
for e in edgesdf.index.values:

    edgeValue = edgesdf.loc[e]['main']
    if edgeValue > 25:
        startDiff = 0.2
    else:
        startDiff = 0.1
    start_xanes = (edgeValue - startDiff).round(2)
    stop_xanes = (edgeValue + 0.15).round(2)
    start_exafs = (edgeValue - startDiff).round(2)
    stop_exafs = (edgeValue + 1).round(2)
    step_xanes = xanes_step(start_xanes)
    step_exafs = exafs_step(start_exafs)
    '''
    if edgeValue < 10:
        step_xanes = 0.3
        step_exafs = 0.5
    elif edgeValue >= 10 and edgeValue < 17:
        step_xanes = 0.5
        step_exafs = 0.7
    elif edgeValue >= 17:
        step_xanes = 0.7
        step_exafs = 1
    '''
    xanes_programDF.loc[e] = [start_xanes,stop_xanes,step_xanes]
    exafs_programDF.loc[e] = [start_exafs,stop_exafs,step_exafs]
    if edgesdf.loc[e]['K'] > 35 and edgesdf.loc[e]['K'] < 41:
        element_name = f'{e}_K'
        startDiff = 0.2
        edgeValue = edgesdf.loc[e]['K']
        start_xanes = (edgeValue - startDiff).round(2)
        stop_xanes = (edgeValue + 0.15).round(2)
        start_exafs = (edgeValue - startDiff).round(2)
        stop_exafs = (edgeValue + 1).round(2)
        step_xanes = xanes_step(start_xanes)
        step_exafs = exafs_step(start_exafs)
        xanes_programDF.loc[element_name] = [start_xanes,stop_xanes,step_xanes]
        exafs_programDF.loc[element_name] = [start_exafs,stop_exafs,step_exafs]
print(exafs_programDF)
exafs_programDF.to_csv('exafs_programs.dat',sep = ',')
xanes_programDF.to_csv('xanes_programs.dat',sep = ',')
