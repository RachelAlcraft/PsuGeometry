
'''
In this file we compare individual geos
to see if any pdbs are problematic
'''

import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

print('### LOADING csv files ###')
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

pdbListIn = dataPdbCut["pdbCode"].unique()
pdbListIn.sort()

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

print('### Describing all ###')
allpdbUndesc = dataPdbUn[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()
allpdbResdesc = dataPdbRes[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()
allpdbCutdesc = dataPdbCut[['C:O','N:CA','CA:C','C:N+1','TAU']].describe()
allpdbAdjdesc = dataPdbAdj[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()
allpdbLapdesc = dataPdbLap[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()

georep.addCsv(data=allpdbUndesc, title='All Unrestricted')
georep.addCsv(data=allpdbResdesc, title='All Restricted')
georep.addCsv(data=allpdbCutdesc, title='All Cut')
georep.addCsv(data=allpdbAdjdesc, title='All Max Adjusted')
georep.addCsv(data=allpdbLapdesc, title='All Lap Adjusted')

for pdb in pdbListIn:
    print('### ',pdb,' ###')
    onepdbUn = dataPdbUn.query("pdbCode == '" + pdb + "'")
    onepdbRes = dataPdbRes.query("pdbCode == '" + pdb + "'")
    onepdbCut = dataPdbCut.query("pdbCode == '" + pdb + "'")
    onepdbAdj = dataPdbAdj.query("pdbCode == '" + pdb + "'")
    onepdbLap = dataPdbLap.query("pdbCode == '" + pdb + "'")

    onepdbUndesc = onepdbUn[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()
    onepdbResdesc = onepdbRes[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()
    onepdbCutdesc = onepdbCut[['C:O','N:CA','CA:C','C:N+1','TAU']].describe()
    onepdbAdjdesc = onepdbAdj[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()
    onepdbLapdesc = onepdbLap[['C:O', 'N:CA', 'CA:C', 'C:N+1', 'TAU']].describe()

    georep.addCsv(data=onepdbUndesc, title=pdb + ' Unrestricted')
    georep.addCsv(data=onepdbResdesc, title=pdb + ' Restricted')
    georep.addCsv(data=onepdbCutdesc,title=pdb + ' Cut')
    georep.addCsv(data=onepdbAdjdesc,title=pdb + ' Max Adjusted')
    georep.addCsv(data=onepdbLapdesc, title=pdb + ' Lap Adjusted')


georep.printToHtml('PDB Geos', 5, 'geos_per_pdb')


