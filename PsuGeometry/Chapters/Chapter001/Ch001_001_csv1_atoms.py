'''
In this file we load all the maxima differernce files to report on how much
each atom type varies
and to decide on the evidential atoms for each pdb
'''


import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

#1. get the pdb list
print('### CREATING csv files ###')
pdbListIn = help.getPDBList()

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
pdbdata = []
pdbbaddata = []
pdboccdata = []


for pdb in pdbListIn:
    pdbPath = help.filesADJRoot + 'MaximaDifferences_' + pdb + '.csv'
    onepdbmax = pd.read_csv(pdbPath)
    onepdbmax['ID'] = onepdbmax['pdbCode']
    pdbdata.append(onepdbmax)

    pdbBadPath = help.filesADJRoot + 'MaximaDifferences_' + pdb + '_Bad_.csv'
    onepdbmaxbad = pd.read_csv(pdbBadPath)
    onepdbmaxbad['ResNo'] = onepdbmaxbad['ResNo'].astype(str)
    onepdbmaxbad['BAD'] = onepdbmaxbad['pdbCode'] +onepdbmaxbad['Chain'] + onepdbmaxbad['ResNo'] + onepdbmaxbad['AtomType']
    pdbbaddata.append(onepdbmaxbad)

    pdbOccPath = help.filesADJRoot + 'OccupancyMaxima_' + pdb + '.csv'
    onepdbmaxocc = pd.read_csv(pdbOccPath)
    onepdbmaxocc['ResNo'] = onepdbmaxocc['ResNo'].astype(str)
    onepdbmaxocc = onepdbmaxocc.query('Fraction < 1')
    onepdbmaxocc['BAD1'] = onepdbmaxocc['pdbCode']
    onepdbmaxocc['BAD'] = onepdbmaxocc['pdbCode'] + onepdbmaxocc['Chain'] + onepdbmaxocc['ResNo'] + onepdbmaxocc['AtomType']
    pdboccdata.append(onepdbmaxocc)




pdbcsv = pd.concat(pdbdata)
pdbbadcsv = pd.concat(pdbbaddata)
pdbocccsv = pd.concat(pdboccdata)

pdbcsv.to_csv(help.loadPath + "AllAtoms_Maxima.csv", index=False)

badonly =pdbbadcsv['BAD']
badonly.to_csv(help.loadPath + "AllBadAtoms_Maxima.csv", index=False)

occonly =pdbocccsv['BAD']
occonly.to_csv(help.loadPath + "AllBadAtoms_Occupant.csv", index=False)



