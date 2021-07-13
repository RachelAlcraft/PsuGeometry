#####################################################################
## Summary correlation reports #########################################
#####################################################################
import pandas as pd
from PsuGeometry import GeoReport as psu

pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

fileMergedDiffAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csvmean.csv'
fileUnrestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csvmean.csv'
fileRestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csvmean.csv'
fileReducedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csvmean.csv'
fileAdjustedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csvmean.csv'

dataCsvAdjustedAv = pd.read_csv(fileAdjustedAv)
dataCsvReducedAv = pd.read_csv(fileMergedDiffAv)
dataCsvRestrictedAv = pd.read_csv(fileRestrictedAv)
dataCsvUnrestrictedAv = pd.read_csv(fileUnrestrictedAv)

dataCsvAdjustedAv['ID'] = dataCsvAdjustedAv['pdbCode']
dataCsvReducedAv['ID'] = dataCsvReducedAv['pdbCode']
dataCsvRestrictedAv['ID'] = dataCsvRestrictedAv['pdbCode']
dataCsvUnrestrictedAv['ID'] = dataCsvUnrestrictedAv['pdbCode']

georep = psu.GeoReport([],pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

georep.addScatter(data=dataCsvUnrestrictedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestrictedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReducedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjustedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestrictedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestrictedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReducedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjustedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Adjusted', palette='jet',sort='NON')


georep.printToHtml('Average scatters bonds compare', 4, 'ScatterAverage1_Compare')
