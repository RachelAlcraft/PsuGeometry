'''
This script creates a file for hand chosen geos, so it is slower as they are no serialised.
This particular report is designed to look at hydrogen bonding on the carbonyl oxygen
'''

import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

geos = ['TAU','TAU+1','TAU-1','CA:C:O','O:C:N+1','CA:C:N+1','CA-1:CA:CA+1',
        'O-1:C-1','C-1:N','N:CA','CA:C','C:O','C:N+1','N+1:CA+1','CA+1:C+1','C+1:O+1',
        'PHI','PSI','OMEGA','CA-1:C-1:N:CA',
        'CA-1:CA','CA:CA+1','C-1:C','C:C+1','N-1:N','N:N+1',
        'CA-1:N','CA-1:O-1','O-1:N','C-1:CA','N:C','CA:O','CA:N+1','O:N+1','C:CA+1','N+1:C+1',
        'O-1:CA','N:O','O:CA+1','N+1:O+1','N-1:O-1']


title='Finding evidential residues'
fileName = 'evidential'


print('### LOADING csv files ###') # bit rubbish but we didn;t change the object references with dssp
dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bbden_adjusted.csv")
dataPdbLap = pd.read_csv(help.loadPath + "bblap_adjusted.csv")
# Find restrictions
#Reduced On lap-diff <0.02
ev1DataPdbCut = help.applyRestrictions(dataPdbCut,True,True,True,True,True)
ev1DataPdbAdj = help.applyRestrictions(dataPdbAdj,True,True,True,False,True)
ev1DataPdbLap = help.applyRestrictions(dataPdbLap,True,True,True,False,True)
#Reduced on resolution
ev2DataPdbCut =ev1DataPdbCut.query('RES < 0.85')
ev2DataPdbAdj =ev1DataPdbAdj.query('RES < 0.85')
ev2DataPdbLap =ev1DataPdbLap.query('RES < 0.85')

print('### Creating scatter files ###')

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

georep.addComment('1) First get an idea of how the data compares with all resolution and topology')
georep.addComment('')
georep.addComment('')
georep.addScatter(data=dataPdbCut, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Original coords',categorical=False,palette='viridis_r', sort='DESC')
georep.addScatter(data=dataPdbAdj, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Den Adj coords',categorical=False,palette='viridis_r', sort='DESC')
georep.addScatter(data=dataPdbLap, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Lap Adj coords',categorical=False,palette='viridis_r', sort='DESC')

georep.addComment('The probability and resolution do not correlate')
georep.addComment('So we are not dealing with the rarity effect')
georep.addComment('')
georep.addHexBins(data=dataPdbCut, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')
georep.addHexBins(data=dataPdbAdj, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')
georep.addHexBins(data=dataPdbLap, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')

georep.addComment('2) It looks like lower Dist_Den_Lap values are more similar')
georep.addComment('this is because we would expect the maximum value and the minimum 2nd derivative to be in the same place')
georep.addComment('so reduce to <0.02 Dist_Lap_Den')
georep.addScatter(data=ev1DataPdbCut, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Original coords',categorical=False,palette='viridis_r', sort='DESC')
georep.addScatter(data=ev1DataPdbAdj, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Den Adj coords',categorical=False,palette='viridis_r', sort='DESC')
georep.addScatter(data=ev1DataPdbLap, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Lap Adj coords',categorical=False,palette='viridis_r', sort='DESC')

georep.addHexBins(data=ev1DataPdbCut, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')
georep.addHexBins(data=ev1DataPdbAdj, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')
georep.addHexBins(data=ev1DataPdbLap, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')

georep.addComment('2) Lower resolution distributes differently, which has no meaning in reality')
georep.addComment("Distances are not dependent on experimental measurement, so we can't trust lower res")
georep.addComment('So reduce to <0.85A Resolution')
georep.addScatter(data=ev2DataPdbCut, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Original coords',categorical=False,palette='viridis_r', sort='DESC')
georep.addScatter(data=ev2DataPdbAdj, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Den Adj coords',categorical=False,palette='viridis_r', sort='DESC')
georep.addScatter(data=ev2DataPdbLap, geoX='C:O', geoY='Dist_Den_Lap', hue='RES',title='Lap Adj coords',categorical=False,palette='viridis_r', sort='DESC')

georep.addHexBins(data=ev2DataPdbCut, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')
georep.addHexBins(data=ev2DataPdbAdj, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')
georep.addHexBins(data=ev2DataPdbLap, geoX='C:O', geoY='Dist_Den_Lap', hue='COUNT',title='',palette='cubehelix_r')

georep.addComment('Now with evidential data, compare distributions')
georep.addComment('')
georep.addComment('')
georep.addHistogram(data=ev2DataPdbCut, geoX='N:CA',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbAdj, geoX='N:CA',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbLap, geoX='N:CA',title='', hue='ID')

georep.addHistogram(data=ev2DataPdbCut, geoX='CA:C',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbAdj, geoX='CA:C',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbLap, geoX='CA:C',title='', hue='ID')

georep.addHistogram(data=ev2DataPdbCut, geoX='C:O',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbAdj, geoX='C:O',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbLap, geoX='C:O',title='', hue='ID')

georep.addHistogram(data=ev2DataPdbCut, geoX='C:N+1',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbAdj, geoX='C:N+1',title='', hue='ID')
georep.addHistogram(data=ev2DataPdbLap, geoX='C:N+1',title='', hue='ID')

georep.addScatter(data=ev2DataPdbCut, geoX='C:O', geoY='C:N+1', hue='TAU+1',title='Original coords',categorical=False,palette='jet', sort='RAND')
georep.addScatter(data=ev2DataPdbAdj, geoX='C:O', geoY='C:N+1', hue='TAU+1',title='Den Adj coords',categorical=False,palette='jet', sort='RAND')
georep.addScatter(data=ev2DataPdbLap, geoX='C:O', geoY='C:N+1', hue='TAU+1',title='Lap Adj coords',categorical=False,palette='jet', sort='RAND')

georep.addScatter(data=ev2DataPdbCut, geoX='C:O', geoY='Dist_Den_Lap', hue='dssp',title='Original coords',categorical=True,palette='tab10', sort='RAND')
georep.addScatter(data=ev2DataPdbAdj, geoX='C:O', geoY='Dist_Den_Lap', hue='dssp',title='Den Adj coords',categorical=True,palette='tab10', sort='RAND')
georep.addScatter(data=ev2DataPdbLap, geoX='C:O', geoY='Dist_Den_Lap', hue='dssp',title='Lap Adj coords',categorical=True,palette='tab10', sort='RAND')

dsspList = dataPdbCut["dssp"].unique()

for dssp in dsspList:
    print('### ',dssp,' ###')

    onepdbUn = ev2DataPdbCut.query("dssp == '" + dssp + "'")
    onepdbRes = ev2DataPdbAdj.query("dssp == '" + dssp + "'")
    onepdbCut = ev2DataPdbLap.query("dssp == '" + dssp + "'")

    georep.addHistogram(data=onepdbUn, geoX='C:O', title='Original Coords ' + str(dssp), hue='ID')
    georep.addHistogram(data=onepdbRes, geoX='C:O', title='Density Adjusted ' + str(dssp), hue='ID')
    georep.addHistogram(data=onepdbCut, geoX='C:O', title='Laplacian Adjusted ' + str(dssp), hue='ID')



georep.printToHtml(title, 3, fileName)
