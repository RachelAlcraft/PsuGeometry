'''
This script creates a file for hand chosen geos, so it is slower as they are no serialised.
This particular report is designed to look at hydrogen bonding on the carbonyl oxygen
'''

import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

geos = ['TAU','N:CA','CA:C','C:O','C:N+1']

print('### LOADING csv files ###')
dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bbden_adjusted.csv")
dataPdbLap = pd.read_csv(help.loadPath + "bblap_adjusted.csv")
# ensure data is correctly restricted
dataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True,True)
dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,False,True)
dataPdbLap = help.applyRestrictions(dataPdbLap,True,True,True,False,True)

locsMax = []
locsLap = []
locsMax.append('ID')
locsLap.append('ID')
for geo in geos:
    dataPdbCut[geo +'_Orig'] = dataPdbCut[geo]
    dataPdbAdj[geo + '_MaxAdj'] = dataPdbAdj[geo]
    dataPdbLap[geo + '_LapAdj'] = dataPdbLap[geo]
    locsMax.append(geo + '_MaxAdj')
    locsLap.append(geo + '_LapAdj')

dataPdbAdj = dataPdbAdj[locsMax]
dataPdbLap = dataPdbLap[locsLap]

mergedDataSet = dataPdbCut
mergedDataSet['ID2'] =mergedDataSet['ID']
mergedDataSet = mergedDataSet.set_index('ID').join(dataPdbAdj.set_index('ID'))
for geo in geos:
    mergedDataSet[geo +'_MaxDiff'] = mergedDataSet[geo +'_Orig'] - mergedDataSet[geo +'_MaxAdj']

mergedDataSet = mergedDataSet.set_index('ID2').join(dataPdbLap.set_index('ID'))
for geo in geos:
    mergedDataSet[geo +'_LapDiff'] = mergedDataSet[geo +'_Orig'] - mergedDataSet[geo +'_LapAdj']
    mergedDataSet[geo + '_MaxLapDiff'] = mergedDataSet[geo + '_MaxDiff'] - mergedDataSet[geo + '_LapAdj']


mergedDataSet.to_csv(help.loadPath + "MergedEvidenced.csv", index=False)