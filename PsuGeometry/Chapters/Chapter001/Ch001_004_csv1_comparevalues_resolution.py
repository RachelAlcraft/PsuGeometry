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
dataPdbAdj = pd.read_csv(help.loadPath + "bb_adjusted.csv")

# ensure data is correctly restricted
dataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True)
dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,False)

# embellish with contacts


locs = []
locs.append('ID')
for geo in geos:
    dataPdbCut[geo +'_Orig'] = dataPdbCut[geo]
    dataPdbAdj[geo + '_Adj'] = dataPdbAdj[geo]
    locs.append(geo + '_Adj')

dataPdbAdj = dataPdbAdj[locs]

mergedDataSet = dataPdbCut
mergedDataSet = mergedDataSet.set_index('ID').join(dataPdbAdj.set_index('ID'))
for geo in geos:
    mergedDataSet[geo +'_Diff'] = mergedDataSet[geo +'_Orig'] - mergedDataSet[geo +'_Adj']



mergedDataSet.to_csv(help.loadPath + "MergedEvidenced.csv", index=False)