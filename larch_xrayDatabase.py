from larch.xray import xray_edges, atomic_symbol
import pandas as pd
import numpy as np
import os
from xasCorn.columnExtraction_thetaCorrection import angle_to_kev, dspacing, speedOfLight, planck, charge

#os.chdir(r'C:\Users\kenneth1a\Documents\beamlineData/')




def kev_to_angle(energykeV):
    Ej = charge*energykeV*1000
    freq = Ej/planck
    wavelength = (speedOfLight/freq)*10**10
    theta = np.arcsin(wavelength/(2*dspacing))*180/np.pi
    return theta

    
edges = ['K','L1','L2','L3','main']
edgesdf = pd.DataFrame(columns = edges)
thetaOffset =  0.103 # 0.115 calculated with old d-spacing #degrees
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
#0.0047*edge + 0.0623 (edge in keV)
def edgeStart(edgeValue):
    return 0.0045*edgeValue + 0.0685
for e in edgesdf.index.values:

    edgeValue = edgesdf.loc[e]['main']
    startDiff = edgeStart(edgeValue)
    start_xanes = (edgeValue - startDiff).round(2)
    stop_xanes = (edgeValue + 0.15).round(2)
    start_exafs = (edgeValue - startDiff).round(2)
    stop_exafs = (edgeValue + 1).round(2)
    step_xanes = xanes_step(start_xanes)
    step_exafs = exafs_step(start_exafs)

    thetaValue = kev_to_angle(edgeValue)
    thetaSide = thetaValue + thetaOffset
    edgeSide = angle_to_kev(thetaSide)
    startXanesSide = (edgeSide - startDiff).round(2)
    stopXanesSide = (edgeSide+0.15).round(2)
    startExafsSide = (edgeSide-startDiff).round(2)
    stopExafsSide = (edgeSide+1).round(2)
    
    xanes_programDF.loc[e] = [start_xanes,stop_xanes,step_xanes]
    exafs_programDF.loc[e] = [start_exafs,stop_exafs,step_exafs]

    xanes_programDFside.loc[e] = [startXanesSide,stopXanesSide,step_xanes]
    exafs_programDFside.loc[e] = [startExafsSide,stopExafsSide,step_exafs]

    if edgesdf.loc[e]['K'] > 35 and edgesdf.loc[e]['K'] < 65:
        element_name = f'{e}_K'
        edgeValue = edgesdf.loc[e]['K']
        startDiff = edgeStart(edgeValue)
        start_xanes = (edgeValue - startDiff).round(2)
        stop_xanes = (edgeValue + 0.15).round(2)
        start_exafs = (edgeValue - startDiff).round(2)
        stop_exafs = (edgeValue + 1).round(2)
        step_xanes = xanes_step(start_xanes)
        step_exafs = exafs_step(start_exafs)
        xanes_programDF.loc[element_name] = [start_xanes,stop_xanes,step_xanes]
        exafs_programDF.loc[element_name] = [start_exafs,stop_exafs,step_exafs]

        thetaValue = kev_to_angle(edgeValue)
        thetaSide = thetaValue + thetaOffset
        edgeSide = angle_to_kev(thetaSide)
        startXanesSide = (edgeSide - startDiff).round(2)
        stopXanesSide = (edgeSide+0.15).round(2)
        startExafsSide = (edgeSide-startDiff).round(2)
        stopExafsSide = (edgeSide+1).round(2)
        xanes_programDFside.loc[element_name] = [startXanesSide,stopXanesSide,step_xanes]
        exafs_programDFside.loc[element_name] = [startExafsSide,stopExafsSide,step_exafs]

print(exafs_programDF)
exafs_programDF.to_csv('exafs_programs.dat',sep = ',')
xanes_programDF.to_csv('xanes_programs.dat',sep = ',')

exafs_programDFside.to_csv('exafs_programs_side.dat',sep = ',')
xanes_programDFside.to_csv('xanes_programs_side.dat',sep = ',')

edgeList = xanes_programDF.index.tolist()
edgestr = '\n'.join(edgeList)
f = open('xasCorn/edgeList.txt','w')
f.write(edgestr)
f.close()