
'''
In this file we load all the csv files and look at some distributions of geos
per pdb to check if any pdbs stand out
'''


import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

print('### LOADING csv files ###') # bit rubbish but we didn;t change the object references with dssp
dataPdbUn = pd.read_csv(help.loadPath + "bb_unrestricted.csv")
dataPdbRes = pd.read_csv(help.loadPath + "bb_restricted.csv")
dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bbden_adjusted.csv")
dataPdbLap = pd.read_csv(help.loadPath + "bblap_adjusted.csv")
# ensure data is correctly restricted
dataPdbUn = help.applyRestrictions(dataPdbUn,True,False,False,False,False)
dataPdbRes = help.applyRestrictions(dataPdbRes,True,True,True,True,False)
dataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True,True)
dataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,False,True)
dataPdbLap = help.applyRestrictions(dataPdbLap,True,True,True,False,True)

pdbListIn = dataPdbUn["pdbCode"].unique()
georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

print('### Shorten csv files ###')#Turn the data into a describe for the volumns we are interested in
cutdataPdbUn = dataPdbUn[['pdbCode','C:O','N:CA','CA:C','C:N+1','TAU']]
cutdataPdbRes = dataPdbRes[['pdbCode','C:O','N:CA','CA:C','C:N+1','TAU']]
cutdataPdbCut = dataPdbCut[['pdbCode','C:O','N:CA','CA:C','C:N+1','TAU']]
cutdataPdbAdj = dataPdbAdj[['pdbCode','C:O','N:CA','CA:C','C:N+1','TAU']]
cutdataPdbLap = dataPdbLap[['pdbCode','C:O','N:CA','CA:C','C:N+1','TAU']]
#cutdataPdbCut01 = dataPdbCut01[['pdbCode','C:O','N:CA','CA:C','C:N+1','TAU']]
#cutdataPdbAdj01 = dataPdbAdj01[['pdbCode','C:O','N:CA','CA:C','C:N+1','TAU']]

print('### Describe csv files ###')#Turn the data into a describe for the volumns we are interested in
descdataPdbUn = cutdataPdbUn.groupby('pdbCode').describe()
descdataPdbRes = cutdataPdbRes.groupby('pdbCode').describe()
descdataPdbCut = cutdataPdbCut.groupby('pdbCode').describe()
descdataPdbAdj = cutdataPdbAdj.groupby('pdbCode').describe()
descdataPdbLap = cutdataPdbLap.groupby('pdbCode').describe()
#descdataPdbCut01 = cutdataPdbCut01.groupby('pdbCode').describe()
#descdataPdbAdj01 = cutdataPdbAdj01.groupby('pdbCode').describe()

print('### Change column names ###')#Turn the data into a describe for the volumns we are interested in
descdataPdbUn['pdbCode'] = descdataPdbUn.index
descdataPdbRes['pdbCode'] = descdataPdbRes.index
descdataPdbCut['pdbCode'] = descdataPdbCut.index
descdataPdbAdj['pdbCode'] = descdataPdbAdj.index
descdataPdbLap['pdbCode'] = descdataPdbLap.index
#descdataPdbCut01['pdbCode'] = descdataPdbCut01.index
#descdataPdbAdj01['pdbCode'] = descdataPdbAdj01.index


descdataPdbUn.columns = [['C:O count','C:O mean','C:O std','C:O min','C:O 25%','C:O 50%','C:O 75%','C:O max',
'N:CA count','N:CA mean','N:CA std','N:CA min','N:CA 25%','N:CA 50%','N:CA 75%','N:CA max',
'CA:C count','CA:C mean','CA:C std','CA:C min','CA:C 25%','CA:C 50%','CA:C 75%','CA:C max',
'C:N+1 count','C:N+1 mean','C:N+1 std','C:N+1 min','C:N+1 25%','C:N+1 50%','C:N+1 75%','C:N+1 max',
'TAU count','TAU mean','TAU std','TAU min','TAU 25%','TAU 50%','TAU 75%','TAU max','pdbCode']]

descdataPdbRes.columns = [['C:O count','C:O mean','C:O std','C:O min','C:O 25%','C:O 50%','C:O 75%','C:O max',
'N:CA count','N:CA mean','N:CA std','N:CA min','N:CA 25%','N:CA 50%','N:CA 75%','N:CA max',
'CA:C count','CA:C mean','CA:C std','CA:C min','CA:C 25%','CA:C 50%','CA:C 75%','CA:C max',
'C:N+1 count','C:N+1 mean','C:N+1 std','C:N+1 min','C:N+1 25%','C:N+1 50%','C:N+1 75%','C:N+1 max',
'TAU count','TAU mean','TAU std','TAU min','TAU 25%','TAU 50%','TAU 75%','TAU max','pdbCode']]

descdataPdbCut.columns = [['C:O count','C:O mean','C:O std','C:O min','C:O 25%','C:O 50%','C:O 75%','C:O max',
'N:CA count','N:CA mean','N:CA std','N:CA min','N:CA 25%','N:CA 50%','N:CA 75%','N:CA max',
'CA:C count','CA:C mean','CA:C std','CA:C min','CA:C 25%','CA:C 50%','CA:C 75%','CA:C max',
'C:N+1 count','C:N+1 mean','C:N+1 std','C:N+1 min','C:N+1 25%','C:N+1 50%','C:N+1 75%','C:N+1 max',
'TAU count','TAU mean','TAU std','TAU min','TAU 25%','TAU 50%','TAU 75%','TAU max','pdbCode']]

descdataPdbAdj.columns = [['C:O count','C:O mean','C:O std','C:O min','C:O 25%','C:O 50%','C:O 75%','C:O max',
'N:CA count','N:CA mean','N:CA std','N:CA min','N:CA 25%','N:CA 50%','N:CA 75%','N:CA max',
'CA:C count','CA:C mean','CA:C std','CA:C min','CA:C 25%','CA:C 50%','CA:C 75%','CA:C max',
'C:N+1 count','C:N+1 mean','C:N+1 std','C:N+1 min','C:N+1 25%','C:N+1 50%','C:N+1 75%','C:N+1 max',
'TAU count','TAU mean','TAU std','TAU min','TAU 25%','TAU 50%','TAU 75%','TAU max','pdbCode']]

descdataPdbLap.columns = [['C:O count','C:O mean','C:O std','C:O min','C:O 25%','C:O 50%','C:O 75%','C:O max',
'N:CA count','N:CA mean','N:CA std','N:CA min','N:CA 25%','N:CA 50%','N:CA 75%','N:CA max',
'CA:C count','CA:C mean','CA:C std','CA:C min','CA:C 25%','CA:C 50%','CA:C 75%','CA:C max',
'C:N+1 count','C:N+1 mean','C:N+1 std','C:N+1 min','C:N+1 25%','C:N+1 50%','C:N+1 75%','C:N+1 max',
'TAU count','TAU mean','TAU std','TAU min','TAU 25%','TAU 50%','TAU 75%','TAU max','pdbCode']]

descdataPdbUn['ID'] = descdataPdbUn['pdbCode']
descdataPdbRes['ID'] = descdataPdbRes['pdbCode']
descdataPdbCut['ID'] = descdataPdbCut['pdbCode']
descdataPdbAdj['ID'] = descdataPdbAdj['pdbCode']
descdataPdbLap['ID'] = descdataPdbLap['pdbCode']


print('### Save csv files ###')#Turn the data into a describe for the volumns we are interested in
descdataPdbUn.to_csv(help.loadPath + "DescribeGeos_Unrestricted.csv", index=False)
descdataPdbRes.to_csv(help.loadPath + "DescribeGeos_Restricted.csv", index=False)
descdataPdbCut.to_csv(help.loadPath + "DescribeGeos_Cut.csv", index=False)
descdataPdbAdj.to_csv(help.loadPath + "DescribeGeos_AdjustedMax.csv", index=False)
descdataPdbLap.to_csv(help.loadPath + "DescribeGeos_AdjustedLap.csv", index=False)



