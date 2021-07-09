import pandas as pd
from PsuGeometry import GeoReport as psu

### Data set up ###
pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'
csvFileDiff = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/DataFoDiffs.csv'
csvFileMaxima = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/Data_DefensibleWithGeosALL_Fo_ADJ.csv'
csvFileDefensible = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/Data_DefensibleWithGeosALL_RESTRICTED_CUTFo.csv'
csvFileRestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/Data_DefensibleWithGeosALL_RESTRICTED.csv'
csvFileUnrestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/Data_DefensibleWithGeosALL_UNRESTRICTED.csv'


### 1. ###
### First create the report that shows the differences between the Fo adjusted data and the restricted data

dataCsvFoDiff = pd.read_csv(csvFileDiff)
geos = [['CA:C',1.4,1.65] ,['C:O',1.1,1.35] ,['C:N+1',1.2,1.45] ,['N:CA',1.35,1.65]]
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
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Diff', geoY='RES', hue='SOFTWARE', palette='jet_r', sort='RANDOM', categorical=True, title='Resolution and atom differences ' + geo)
for geol in geos:
    geo = geol[0]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Diff', geoY='SOFTWARE', hue='RES', palette='viridis_r', sort='DESC', categorical=False, title='Resolution and atom differences ' + geo)
for geol in geos:
    geo = geol[0]
    georep1.addHexBins(data=dataCsvFoDiff, geoX=geo + '_Diff', geoY='RES', title='Count ' + geo, hue='count', palette='cubehelix_r')
for geol in geos:
    geo = geol[0]
    range = [geol[1],geol[2]]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', hue='SOFTWARE', palette='jet_r',sort='RANDOM', categorical=True,title='Software and atom positions ' + geo,range=range)
for geol in geos:
    geo = geol[0]
    range = [geol[1], geol[2]]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', hue='RES', palette='viridis_r', sort='DESC', categorical=False, title='Resolution and atom positions ' + geo,range=range)
for geol in geos:
    geo = geol[0]
    range = [geol[1], geol[2]]
    georep1.addScatter(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', hue='count', palette='viridis_r', sort='COUNT', categorical=False, title='Resolution and atom positions ' + geo,range=range)
for geol in geos:
    geo = geol[0]
    range = [geol[1], geol[2]]
    georep1.addHexBins(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', title='Count ' + geo, hue='count', palette='cubehelix_r',range=range)
for geol in geos:
    geo = geol[0]
    range = [geol[1], geol[2]]
    georep1.addProbability(data=dataCsvFoDiff, geoX=geo + '_Orig', geoY=geo + '_Adj', title='Count ' + geo, hue='count', palette='cubehelix_r',range=range)
georep1.printToHtml('Comparing atom positions: PDB vs maxima, set=Fo3', 4, 'Compare_Fo3')

### 2. ###
### Comparing the 4 sets ###
georep2 = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

dataCsvMaxima = pd.read_csv(csvFileMaxima)
dataCsvDefensible = pd.read_csv(csvFileDefensible)
dataCsvRestricted = pd.read_csv(csvFileRestricted)
dataCsvUnrestricted = pd.read_csv(csvFileUnrestricted)

geos = ['CA:C' ,'C:O' ,'C:N+1' ,'N:CA','TAU','PHI','PSI','OMEGA']
for geo in geos:
    georep2.addHistogram(data=dataCsvUnrestricted, geoX=geo, title=geo + ' Unrestricted')
    georep2.addHistogram(data=dataCsvRestricted, geoX=geo, title=geo + ' Restricted')
    georep2.addHistogram(data=dataCsvDefensible, geoX=geo, title=geo + ' Defensible')
    georep2.addHistogram(data=dataCsvMaxima, geoX=geo, title=geo + ' Maxima')

georep2.printToHtml('Comparing geoemtry sets: PDB vs adjusted, set=Fo3', 4, 'Compare_Sets')

