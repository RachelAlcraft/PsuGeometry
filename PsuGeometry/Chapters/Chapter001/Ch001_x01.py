import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help

######################################
# Difference reports
######################################

### Data set up ###
pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

fileMergedDiff = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'
fileUnrestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csv'
fileRestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csv'
fileReduced = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'
fileAdjusted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csv'

### 1. ###
### First create the report that shows the differences between the Fo adjusted data and the restricted data

dataCsvFoDiff = pd.read_csv(fileMergedDiff)
geos = [['N:CA',1.35,1.55],['CA:C',1.4,1.65] ,['C:O',1.1,1.37] ,['C:N+1',1.18,1.45]]
georep1 = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
for geol in geos:
    geo = geol[0]
    dataCsvFoDiff[geo + '_Diff'] = dataCsvFoDiff[geo + '_Orig'] - dataCsvFoDiff[geo + '_Adj']

for geol in geos:
    geo = geol[0]
    georep1.addHistogram(data=dataCsvFoDiff, geoX=geo + '_Orig', title='Pdb Atoms', count=True, hue='pdbCode')
for geol in geos:
    geo = geol[0]
    georep1.addHistogram(data=dataCsvFoDiff, geoX=geo + '_Adj', title='Adjusted Atoms', count=True, hue='pdbCode')
for geol in geos:
    geo = geol[0]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Diff', geoY='Resolution', hue='Software', palette='jet_r', sort='RANDOM', categorical=True, title='Resolution and atom differences ' + geo)
for geol in geos:
    geo = geol[0]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Diff', geoY='Software', hue='Resolution', palette='viridis_r', sort='DESC', categorical=False, title='Resolution and atom differences ' + geo)
for geol in geos:
    geo = geol[0]
    georep1.addHexBins(data=dataCsvFoDiff, geoX=geo + '_Diff', geoY='Resolution', title='Count ' + geo, hue='count', palette='cubehelix_r')
for geol in geos:
    geo = geol[0]
    range = [geol[1],geol[2]]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', hue='Software', palette='jet_r',sort='RANDOM', categorical=True,title='Software and atom positions ' + geo,range=range)
for geol in geos:
    geo = geol[0]
    range = [geol[1], geol[2]]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', hue='Resolution', palette='viridis_r', sort='DESC', categorical=False, title='Resolution and atom positions ' + geo,range=range)
for geol in geos:
    geo = geol[0]
    range = [geol[1], geol[2]]
    georep1.addHexBins(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', title='Count ' + geo, hue='count', palette='cubehelix_r',range=range)

georep1.printToHtml('Comparing Geometry: PDB vs Adjusted, set=Fo3', 4, 'Compare_Fo3')

### 2. ###
### Comparing the 4 sets ###
georep2 = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

dataCsvAdjusted = pd.read_csv(fileAdjusted)
dataCsvReduced = pd.read_csv(fileReduced)
dataCsvRestricted = pd.read_csv(fileRestricted)
dataCsvUnrestricted = pd.read_csv(fileUnrestricted)

geos = ['N:CA','CA:C','C:O','C:N+1','TAU','PHI','PSI','OMEGA']
for geo in geos:
    georep2.addHistogram(data=dataCsvUnrestricted, geoX=geo, title=geo + ' Unrestricted')
    georep2.addHistogram(data=dataCsvRestricted, geoX=geo, title=geo + ' Restricted')
    georep2.addHistogram(data=dataCsvReduced, geoX=geo, title=geo + ' Reduced')
    georep2.addHistogram(data=dataCsvAdjusted, geoX=geo, title=geo + ' Adjusted')

georep2.printToHtml('Comparing geoemtry sets: PDB vs adjusted, set=Fo3', 4, 'Compare_Sets')

