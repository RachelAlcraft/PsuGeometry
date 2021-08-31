#####################################################################
## SuAverage data #########################################
#####################################################################
import pandas as pd
from PsuGeometry import GeoReport as psu

pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

fileMergedDiffAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csvmean.csv'
fileUnrestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csvmean.csv'
fileRestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csvmean.csv'
fileReducedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataReduced.csvmean.csv' #'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csvmean.csv'
fileAdjustedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csvmean.csv'

fileUnrestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csv'
fileRestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csv'
fileReduced = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'
fileAdjusted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csv'

dataCsvAdjustedAv = pd.read_csv(fileAdjustedAv)
#dataCsvReducedAv = pd.read_csv(fileMergedDiffAv)
dataCsvReducedAv = pd.read_csv(fileReducedAv)
dataCsvRestrictedAv = pd.read_csv(fileRestrictedAv)
dataCsvUnrestrictedAv = pd.read_csv(fileUnrestrictedAv)

dataCsvAdjusted = pd.read_csv(fileAdjusted)
dataCsvReduced = pd.read_csv(fileReduced)
dataCsvRestricted = pd.read_csv(fileRestricted)
dataCsvUnrestricted = pd.read_csv(fileUnrestricted)

dataCsvAdjustedAv['ID'] = dataCsvAdjustedAv['pdbCode']
dataCsvReducedAv['ID'] = dataCsvReducedAv['pdbCode']
dataCsvRestrictedAv['ID'] = dataCsvRestrictedAv['pdbCode']
dataCsvUnrestrictedAv['ID'] = dataCsvUnrestrictedAv['pdbCode']

georep = psu.GeoReport([],pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

georep.addScatter(data=dataCsvUnrestricted, geoX='PHI', geoY='PSI', hue='C:O', title='Avg Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='PHI', geoY='PSI', hue='C:O', title='Avg Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='PHI', geoY='PSI', hue='C:O', title='Avg Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='PHI', geoY='PSI', hue='C:O', title='Avg Adjusted', palette='jet',sort='NON')
'''
georep.addScatter(data=dataCsvUnrestrictedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Avg Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestrictedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Avg Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReducedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Avg Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjustedAv, geoX='N:CA', geoY='CA:C', hue='C:O', title='Avg Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='N:CA', geoY='CA:C', hue='C:O', title='All Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='N:CA', geoY='CA:C', hue='C:O', title='All Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='N:CA', geoY='CA:C', hue='C:O', title='All Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='N:CA', geoY='CA:C', hue='C:O', title='All Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestrictedAv, geoX='N:CA', geoY='CA:C', hue='Resolution', title='Avg Unrestricted', palette='jet',sort='DESC')
georep.addScatter(data=dataCsvRestrictedAv, geoX='N:CA', geoY='CA:C', hue='Resolution', title='Avg Restricted', palette='jet',sort='DESC')
georep.addScatter(data=dataCsvReducedAv, geoX='N:CA', geoY='CA:C', hue='Resolution', title='Avg Reduced', palette='jet',sort='DESC')
georep.addScatter(data=dataCsvAdjustedAv, geoX='N:CA', geoY='CA:C', hue='Resolution', title='Avg Adjusted', palette='jet',sort='DESC')

georep.addScatter(data=dataCsvUnrestricted, geoX='N:CA', geoY='CA:C', hue='Resolution', title='All Unrestricted', palette='jet',sort='DESC')
georep.addScatter(data=dataCsvRestricted, geoX='N:CA', geoY='CA:C', hue='Resolution', title='All Restricted', palette='jet',sort='DESC')
georep.addScatter(data=dataCsvReduced, geoX='N:CA', geoY='CA:C', hue='Resolution', title='All Reduced', palette='jet',sort='DESC')
georep.addScatter(data=dataCsvAdjusted, geoX='N:CA', geoY='CA:C', hue='Resolution', title='All Adjusted', palette='jet',sort='DESC')

georep.addScatter(data=dataCsvUnrestricted, geoX='N:CA', geoY='CA:C', hue='Software', title='All Unrestricted', categorical=True,palette='jet_r',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='N:CA', geoY='CA:C', hue='Software', title='All Restricted', categorical=True,palette='jet_r',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='N:CA', geoY='CA:C', hue='Software', title='All Reduced', categorical=True,palette='jet_r',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='N:CA', geoY='CA:C', hue='Software', title='All Adjusted', categorical=True,palette='jet_r',sort='NON')

georep.addScatter(data=dataCsvUnrestrictedAv, geoX='N:CA', geoY='CA:C', hue='TAU', title='Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestrictedAv, geoX='N:CA', geoY='CA:C', hue='TAU', title='Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReducedAv, geoX='N:CA', geoY='CA:C', hue='TAU', title='Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjustedAv, geoX='N:CA', geoY='CA:C', hue='TAU', title='Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestrictedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Avg Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestrictedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Avg Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReducedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Avg Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjustedAv, geoX='C:O', geoY='C:N+1', hue='TAU', title='Avg Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='C:O', geoY='C:N+1', hue='TAU', title='All Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='C:O', geoY='C:N+1', hue='TAU', title='All Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='C:O', geoY='C:N+1', hue='TAU', title='All Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='C:O', geoY='C:N+1', hue='TAU', title='All Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='CA:C:O', geoY='O:C:N+1', hue='CA:C:N+1', title='All Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='CA:C:O', geoY='O:C:N+1', hue='CA:C:N+1', title='All Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='CA:C:O', geoY='O:C:N+1', hue='CA:C:N+1', title='All Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='CA:C:O', geoY='O:C:N+1', hue='CA:C:N+1', title='All Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='CA:C', geoY='C:N+1', hue='C:O', title='All Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='CA:C', geoY='C:N+1', hue='C:O', title='All Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='CA:C', geoY='C:N+1', hue='C:O', title='All Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='CA:C', geoY='C:N+1', hue='C:O', title='All Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='CA:C', geoY='N:CA', hue='C:N+1', title='All Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='CA:C', geoY='N:CA', hue='C:N+1', title='All Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='CA:C', geoY='N:CA', hue='C:N+1', title='All Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='CA:C', geoY='N:CA', hue='C:N+1', title='All Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='N:CA:C', geoY='CA:C:N+1', hue='C:O', title='All Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='N:CA:C', geoY='CA:C:N+1', hue='C:O', title='All Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='N:CA:C', geoY='CA:C:N+1', hue='C:O', title='All Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='N:CA:C', geoY='CA:C:N+1', hue='C:O', title='All Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='C:O', geoY='C:N+1', hue='TAU+1', title='All Unrestricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='C:O', geoY='C:N+1', hue='TAU+1', title='All Restricted', palette='jet',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='C:O', geoY='C:N+1', hue='TAU+1', title='All Reduced', palette='jet',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='C:O', geoY='C:N+1', hue='TAU+1', title='All Adjusted', palette='jet',sort='NON')

georep.addScatter(data=dataCsvUnrestricted, geoX='C:O', geoY='C:N+1', hue='dssp', title='All Unrestricted', categorical=True, palette='jet_r',sort='NON')
georep.addScatter(data=dataCsvRestricted, geoX='C:O', geoY='C:N+1', hue='dssp', title='All Restricted', categorical=True, palette='jet_r',sort='NON')
georep.addScatter(data=dataCsvReduced, geoX='C:O', geoY='C:N+1', hue='dssp', title='All Reduced', categorical=True, palette='jet_r',sort='NON')
georep.addScatter(data=dataCsvAdjusted, geoX='C:O', geoY='C:N+1', hue='dssp', title='All Adjusted', categorical=True, palette='jet_r',sort='NON')

'''
georep.printToHtml('Average scatters bonds compare', 4, 'ScatterAverage1_Compare')
