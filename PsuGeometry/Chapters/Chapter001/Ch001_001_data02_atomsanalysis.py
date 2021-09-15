
'''
In this file we load all the maxima differernce files to report on how much
each atom type varies
and to decide on the evidential atoms for each pdb
'''


import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

allPath = help.loadPath + 'AllAtoms_Maxima.csv'
allAdjusteddata = pd.read_csv(allPath)

allAdjusteddata['DiffOrigLap'] =allAdjusteddata['DenOrig'] - allAdjusteddata['DenLap']
allAdjusteddata['DiffOrigDen'] =allAdjusteddata['DenOrig'] - allAdjusteddata['DenDen']

atomsdataBB = allAdjusteddata.query("AtomType == 'N' or AtomType == 'CA' or AtomType == 'C' or AtomType == 'O'")
atomsdataBF = allAdjusteddata.query("BFactor <= 15")


pdbListIn = allAdjusteddata["pdbCode"].unique()
pdbListIn.sort()
pdbListIn = pdbListIn[:40]

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
georepLap = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
georepDen = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

atomsdataN = allAdjusteddata.query("AtomType == 'N'")
atomsdataCA = allAdjusteddata.query("AtomType == 'CA'")
atomsdataC = allAdjusteddata.query("AtomType == 'C'")
atomsdataO = allAdjusteddata.query("AtomType == 'O'")

georep.addScatter(data=allAdjusteddata, geoX='DenDen',geoY='DenLap', title='Compare Differences', hue='DenOrig',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=allAdjusteddata, geoX='Dist_Den',geoY='Dist_Lap', title='Compare Differences', hue='Dist_Den_Lap',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=allAdjusteddata, geoX='Dist_Den_Lap',geoY='BFactor', title='Compare Differences', hue='LapLap',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=allAdjusteddata, geoX='DenDen',geoY='DenLap', title='Compare Differences', hue='AtomType',categorical=True,sort='RAND',palette='jet_r')

georep.addHistogram(data=atomsdataN,geoX='Dist_Den_Lap', title='N Lap-Den Diff',hue='ID')
georep.addHistogram(data=atomsdataCA,geoX='Dist_Den_Lap',title='CA Lap-Den Diff',hue='ID')
georep.addHistogram(data=atomsdataC,geoX='Dist_Den_Lap',title='C Lap-Den Diff',hue='ID')
georep.addHistogram(data=atomsdataO,geoX='Dist_Den_Lap',title='O Lap-Den Diff')

georep.addHistogram(data=atomsdataN,geoX='DiffOrigDen', title='N Density Diff',hue='ID')
georep.addHistogram(data=atomsdataCA,geoX='DiffOrigDen',title='CA Density Diff',hue='ID')
georep.addHistogram(data=atomsdataC,geoX='DiffOrigDen',title='C Density Diff',hue='ID')
georep.addHistogram(data=atomsdataO,geoX='DiffOrigDen',title='O Density Diff')

georep.addHistogram(data=atomsdataN,geoX='DiffOrigLap', title='N Laplacian Diff',hue='ID')
georep.addHistogram(data=atomsdataCA,geoX='DiffOrigLap',title='CA Laplacian Diff',hue='ID')
georep.addHistogram(data=atomsdataC,geoX='DiffOrigLap',title='C Laplacian Diff',hue='ID')
georep.addHistogram(data=atomsdataO,geoX='DiffOrigLap',title='O Laplacian Diff')

georep.addScatter(data=atomsdataN, geoX='Dist_Den_Lap',geoY='BFactor', title='N Density Adjusted', hue='DiffOrigDen',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataCA, geoX='Dist_Den_Lap',geoY='BFactor', title='CA Density Adjusted', hue='DiffOrigDen',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataC, geoX='Dist_Den_Lap',geoY='BFactor', title='C Density Adjusted', hue='DiffOrigDen',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataO, geoX='Dist_Den_Lap',geoY='BFactor', title='O Density Adjusted', hue='DiffOrigDen',categorical=False,sort='RAND',palette='jet')

georep.addScatter(data=atomsdataN, geoX='Dist_Den_Lap',geoY='BFactor', title='N Laplacian Adjusted', hue='DiffOrigLap',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataCA, geoX='Dist_Den_Lap',geoY='BFactor', title='CA Laplacian Adjusted', hue='DiffOrigLap',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataC, geoX='Dist_Den_Lap',geoY='BFactor', title='C Laplacian Adjusted', hue='DiffOrigLap',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=atomsdataO, geoX='Dist_Den_Lap',geoY='BFactor', title='O Laplacian Adjusted', hue='DiffOrigLap',categorical=False,sort='RAND',palette='jet')

georep.addScatter(data=atomsdataBB, geoX='DenDen', geoY='LapDen', title='Backbone Density', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
georep.addScatter(data=atomsdataBB, geoX='DenLap', geoY='LapLap', title='Backbone Laplacian', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
georep.addScatter(data=atomsdataBF, geoX='DenDen', geoY='LapDen', title='Backbone Density<=15', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
georep.addScatter(data=atomsdataBF, geoX='DenLap', geoY='LapLap', title='Backbone Laplacian<=15', hue='AtomType', categorical=True, sort='RAND', palette='jet_r')

georep.addScatter(data=atomsdataBB, geoX='DiffOrigDen', geoY='Dist_Den_Lap', title='Density Adj', hue='BFactor', categorical=False, sort='RAND', palette='jet')
georep.addScatter(data=atomsdataBB, geoX='DiffOrigLap', geoY='Dist_Den_Lap', title='Laplacian Adj', hue='BFactor', categorical=False, sort='RAND', palette='jet')
georep.addScatter(data=atomsdataBF, geoX='DiffOrigDen', geoY='Dist_Den_Lap', title='Density Adj<=15', hue='BFactor', categorical=False, sort='RAND', palette='jet')
georep.addScatter(data=atomsdataBF, geoX='DiffOrigLap', geoY='Dist_Den_Lap', title='Laplacian Adj<=15', hue='BFactor', categorical=False, sort='RAND', palette='jet')

for pdb in pdbListIn:
    onepdbmax = allAdjusteddata.query("pdbCode == '" + pdb + "'")
    onepdbmaxBB = onepdbmax.query("AtomType == 'N' or AtomType == 'CA' or AtomType == 'C' or AtomType == 'O'")
    onepdbmaxSC = onepdbmax.query("AA == 'GLN'")
    onepdbmaxSC = onepdbmaxSC.query("AtomType == 'NE2' or AtomType == 'CB' or AtomType == 'CD' or AtomType == 'OE1'")

    georepDen.addScatter(data=onepdbmaxBB, geoX='DenDen', geoY='Dist_Den_Lap', title='BB ' + pdb, hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
    georepDen.addScatter(data=onepdbmaxSC, geoX='DenDen', geoY='Dist_Den_Lap', title='SC ' + pdb, hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
    georepDen.addScatter(data=onepdbmax, geoX='DenDen', geoY='Dist_Den_Lap', title='' + pdb, hue='BFactor', categorical=False, sort='RAND', palette='jet')
    georepDen.addScatter(data=onepdbmax, geoX='DenDen', geoY='LapDen', title='' + pdb, hue='DiffOrigDen', categorical=False, sort='ASC', palette='jet')

    georepLap.addScatter(data=onepdbmaxBB, geoX='DenLap', geoY='Dist_Den_Lap', title='BB ' + pdb, hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
    georepLap.addScatter(data=onepdbmaxSC, geoX='DenLap', geoY='Dist_Den_Lap', title='SC ' + pdb, hue='AtomType', categorical=True, sort='RAND', palette='jet_r')
    georepLap.addScatter(data=onepdbmax, geoX='DenLap', geoY='Dist_Den_Lap', title='' + pdb, hue='BFactor', categorical=False,sort='RAND', palette='jet')
    georepLap.addScatter(data=onepdbmax, geoX='DenLap', geoY='LapLap', title='' + pdb, hue='DiffOrigLap', categorical=False, sort='ASC', palette='jet')

georep.printToHtml('Topology Analysis', 4, 'topology')
georepLap.printToHtml('Topology Analysis, Laplacian Adjusted', 4, 'topology_laplacian_adjusted')
georepDen.printToHtml('Topology Analysis, Density Adjusted', 4, 'topology_density_adjusted')