
'''
In this file we load all the maxima differernce files to report on how much
each atom type varies
and to decide on the evidential atoms for each pdb
'''


import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

pdbPath = help.loadPath + 'AllAtoms_Maxima.csv'
atomsdata = pd.read_csv(pdbPath)
atomsdataBB = atomsdata.query("AtomType == 'N' or AtomType == 'CA' or AtomType == 'C' or AtomType == 'O'")
atomsdataSC = atomsdata.query("AtomType == 'NE2' or AtomType == 'CB' or AtomType == 'CD' or AtomType == 'OE1'")
atomsdataSC = atomsdataSC.query("AA == 'GLN'")
atomsdataBF = atomsdata.query("BFactor <= 15")

pdbListIn = atomsdata["pdbCode"].unique()
pdbListIn.sort()
pdbListIn = pdbListIn[:40]
#pdbListIn = ['1ix9']

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

atomsdataN = atomsdata.query("AtomType == 'N'")
atomsdataCA = atomsdata.query("AtomType == 'CA'")
atomsdataC = atomsdata.query("AtomType == 'C'")
atomsdataO = atomsdata.query("AtomType == 'O'")

georep.addHistogram(data=atomsdataN,geoX='LapDiff', title='N all',hue='ID')
georep.addHistogram(data=atomsdataCA,geoX='LapDiff',title='CA all',hue='ID')
georep.addHistogram(data=atomsdataC,geoX='LapDiff',title='C all',hue='ID')
georep.addHistogram(data=atomsdataO,geoX='LapDiff',title='O all')

georep.addScatter(data=atomsdataN, geoX='LapDiff',geoY='BFactor', title='N All', hue='Difference',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataCA, geoX='LapDiff',geoY='BFactor', title='CA All', hue='Difference',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataC, geoX='LapDiff',geoY='BFactor', title='C All', hue='Difference',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataO, geoX='LapDiff',geoY='BFactor', title='O All', hue='Difference',categorical=False,sort='RAND',palette='jet')

georep.addScatter(data=atomsdataBB, geoX='Difference', geoY='LapDiff', title='Backbone', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
georep.addScatter(data=atomsdataSC, geoX='Difference', geoY='LapDiff', title='GLN', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
georep.addScatter(data=atomsdataBB, geoX='Vadj', geoY='Laplacian', title='Backbone', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
georep.addScatter(data=atomsdataSC, geoX='Vadj', geoY='Laplacian', title='GLN', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')

georep.addScatter(data=atomsdata, geoX='Difference', geoY='LapDiff', title='All', hue='BFactor', categorical=False, sort='RAND', palette='jet')
georep.addScatter(data=atomsdata, geoX='Vadj', geoY='Laplacian', title='All', hue='Difference', categorical=False, sort='ASC', palette='jet')
georep.addScatter(data=atomsdataBF, geoX='Difference', geoY='LapDiff', title='All BF<=15', hue='BFactor', categorical=False, sort='RAND', palette='jet')
georep.addScatter(data=atomsdataBF, geoX='Vadj', geoY='Laplacian', title='All BF<=15', hue='Difference', categorical=False, sort='ASC', palette='jet')


for pdb in pdbListIn:
    onepdbmax = atomsdata.query("pdbCode == '" + pdb + "'")
    onepdbmaxBB = onepdbmax.query("AtomType == 'N' or AtomType == 'CA' or AtomType == 'C' or AtomType == 'O'")
    onepdbmaxSC = onepdbmax.query("AA == 'GLN'")
    onepdbmaxSC = onepdbmaxSC.query("AtomType == 'NE2' or AtomType == 'CB' or AtomType == 'CD' or AtomType == 'OE1'")


    georep.addScatter(data=onepdbmaxBB, geoX='Vadj', geoY='LapDiff', title='BB ' + pdb, hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
    georep.addScatter(data=onepdbmaxSC, geoX='Vadj', geoY='LapDiff', title='SC ' + pdb, hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
    georep.addScatter(data=onepdbmax, geoX='Vadj', geoY='LapDiff', title='' + pdb, hue='BFactor', categorical=False, sort='RAND', palette='jet')
    georep.addScatter(data=onepdbmax, geoX='Vadj', geoY='Laplacian', title='' + pdb, hue='Difference', categorical=False, sort='ASC', palette='jet')

georep.printToHtml('Laplacian Analysis', 4, 'laplacian')