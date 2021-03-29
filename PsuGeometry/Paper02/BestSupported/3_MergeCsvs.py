# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
from PsuGeometry import Globals as glob
import time
import pandas as pd
'''
TAU 
I have 3 data sets
Original - the entire glycine dataset form tghe liost of structures I have chosen that are 1.3* avergae with no multiple occupancy
Good - of all the residues above, those with a nearby density peak (good local density)
Best - of all those above, where the tau values are calculated as identical
'''

myWindowsLaptop = True

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
dsspHue='dssp'
includeDSSP = True
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
    includeDSSP = False  # on my windows computer

dataBest = pd.read_csv(printPath + "Set3_BestSupported.csv")
#ID,pdbCode,chain,rid,aa,bfactor,bfactorRatio,disordered,TAU,TAU2,TAU_DIFF

geoLists = []
geoLists.append(['0', ['dssp']])
geoLists.append(['1BOND', ['N:CA','CA:C','C:O','C-1:N','C:N+1']])
geoLists.append(['2ANGS', ['TAU','C-1:N:CA','CA:C:N+1','CA:C:O','O:C:N+1','CA:C:N+1']])
geoLists.append(['3DIHS', ['PHI','PSI','OMEGA','CA-1:C-1:N:CA']])
geoLists.append(['4DIST', ['N:N+1','N:C']])
geoLists.append(['5HB', ['N:O-2','C:O-2','N:CA:C:O-2','N:CA:N+1:O-2']])
geoLists.append(['6HBO', ['N:{O}','C:{O}','N:CA:C:{O}','N:CA:N+1:{O}']])
#geoLists.append(['7WAT', ['N:HOH','C:HOH','N:CA:C:HOH','N:CA:N+1:HOH']])
#geoLists.append(['8XTRA', ['N:HETATM']])

#files.append(['1GLY', ['N:N+1', 'TAU', 'PSI', 'PHI', 'N:C', 'CA:C', 'C:O', 'N:CA', 'C-1:N', 'C:N+1', 'OMEGA']])
#files.append(['2GLY', ['CA:C:O:N+1', 'O:N+1', 'CA:O', 'CA:N+1', 'CA:C:N+1', 'C-1:N:CA', 'N:O-2', 'N:CA:C:O-2']])
#files.append(['3GLY', ['N:O-2:CA', 'N-1:CA:C', 'CA:HOH', 'CA:HETATM', 'N:HETATM:C', 'N:HOH:C', 'N:CA:C:HETATM', 'N:CA:C:HOH']])
#files.append(['4GLY',['O-2:C','O-2:N:CA','O-2:N:CA:N+1']])
#files.append(['6GLY',['N:{O}','C:{O}','CA:N:{O}','N+1:CA:N:{O}']])
#files.append(['7GLY',['N:O-2','C:O-2','CA:N:O-2','N+1:CA:N:O-2']])
#files.append(['8GLY',['N:{OD1,OG1}','C:{OD1,OG1}','CA:N:{OD1,OG1}','N+1:CA:N:{OD1,OG1}']])

for file in geoLists:
    print('Merging with file',file[0])
    fileName = 'CsvGeos_Set' + file[0] + 'ALL.csv'
    fileList = file[1]
    fileList.append('ID')
    print(fileList)
    dataRow = pd.read_csv(printPath + fileName)
    dataRow['rid'] = dataRow['rid'].astype(str)
    dataRow['ID'] = dataRow['pdbCode'] + dataRow['chain'] + dataRow['rid'] + dataRow['aa']
    dataRow = dataRow[fileList]
    dataBest = pd.merge(dataBest,dataRow,left_on='ID',right_on='ID')
    dataBest = dataBest.dropna()

dataBest.to_csv(printPath + "Set4_BestWithGeosALL.csv", index=False)
