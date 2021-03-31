# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
from PsuGeometry import Globals as glob
import time
import pandas as pd
'''
TAU 

'''
pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/NCACO_001_05/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/EvidencedSet/Data/'

dataFirst = pd.read_csv(printPath + "CsvGeos_BEST_Set1BONDALL.csv")
dataFirst['rid'] = dataFirst['rid'].astype(str)
dataFirst['ID'] = dataFirst['pdbCode'] + dataFirst['chain'] + dataFirst['rid'] + dataFirst['aa']

#ID,pdbCode,chain,rid,aa,bfactor,bfactorRatio,disordered,TAU,TAU2,TAU_DIFF

geoLists = []
geoLists.append(['0DSSP', ['dssp']])
#geoLists.append(['1BOND', ['N:CA','CA:C','C:O','C-1:N','C:N+1']])
geoLists.append(['2ANGS', ['TAU','C-1:N:CA','CA:C:N+1','CA:C:O','O:C:N+1','CA:C:N+1']])
#geoLists.append(['3DIHS', ['PHI','PSI','OMEGA','CA-1:C-1:N:CA']])
#geoLists.append(['4DIST', ['N:N+1','N:C']])
#geoLists.append(['5HB', ['N:O-2','C:O-2','N:CA:C:O-2','N:CA:N+1:O-2']])
#geoLists.append(['6HBO', ['N:{O}','C:{O}','N:CA:C:{O}','N:CA:N+1:{O}']])
#geoLists.append(['7WAT', ['N:HOH','C:HOH','N:CA:C:HOH','N:CA:N+1:HOH']])
#geoLists.append(['8XTRA', ['N:HETATM']])

for file in geoLists:
    print('Merging with file',file[0])
    fileName = 'CsvGeos_BEST_Set' + file[0] + 'ALL.csv'
    fileList = file[1]
    fileList.append('ID')
    print(fileList)
    dataRow = pd.read_csv(printPath + fileName)
    dataRow['rid'] = dataRow['rid'].astype(str)
    dataRow['ID'] = dataRow['pdbCode'] + dataRow['chain'] + dataRow['rid'] + dataRow['aa']
    dataRow = dataRow[fileList]
    dataFirst = pd.merge(dataFirst,dataRow,left_on='ID',right_on='ID')
    dataFirst = dataFirst.dropna()

dataFirst.to_csv(printPath + "Data_DefensibleWithGeosALL.csv", index=False)
