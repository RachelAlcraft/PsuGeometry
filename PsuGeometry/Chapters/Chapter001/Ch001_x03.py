#####################################################################
## Make correlation reports #########################################
#####################################################################
import pandas as pd
from PsuGeometry import GeoReport as psu

pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

fileReduced = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'
fileAdjusted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csv'

dataCsvAdjusted = pd.read_csv(fileAdjusted)
dataCsvReduced = pd.read_csv(fileReduced)

georep1 = psu.GeoReport([],pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
georep2 = psu.GeoReport([],pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

aas = dataCsvAdjusted['aa'].values
aas = list(set(aas))
aas.sort()
for aa in aas:
    dataCutA = dataCsvReduced.query("aa ==  '" + aa + "'")
    dataCutB = dataCsvAdjusted.query("aa ==  '" + aa + "'")
    georep1.addScatter(data=dataCutA, geoX='PSI', geoY='N:N+1', hue='TAU', title=aa + ' PSI|N:N+1|TAU Reduced', palette='jet',sort='NON')
    georep2.addScatter(data=dataCutB, geoX='PSI', geoY='N:N+1', hue='TAU', title=aa + ' PSI|N:N+1|TAU Adjusted', palette='jet', sort='NON')

georep1.printToHtml('Reduced scatter', 3, 'ScatterTau1_Reduced')
georep2.printToHtml('Adjusted scatter', 3, 'ScatterTau1_Adjusted')