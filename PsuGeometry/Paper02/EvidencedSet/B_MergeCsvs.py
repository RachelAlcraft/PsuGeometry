# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
from PsuGeometry import Globals as glob
import time
import pandas as pd
'''
TAU 

'''

def mergeCsvs(pdbSet):
    print('Running MergeCsvs for', pdbSet)

    loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataA/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'

    dataFirst = pd.read_csv(loadPath + 'CsvGeos_BEST_Set1BONDALL_' + pdbSet + '.csv')
    dataFirst['rid'] = dataFirst['rid'].astype(str)
    dataFirst['ID'] = dataFirst['pdbCode'] + dataFirst['chain'] + dataFirst['rid'] + dataFirst['aa']

    #ID,pdbCode,chain,rid,aa,bfactor,bfactorRatio,disordered,TAU,TAU2,TAU_DIFF

    geoLists = []
    geoLists.append(['0DSSPALL', ['dssp']])
    #geoLists.append(['1BONDALL_', ['N:CA','CA:C','C:O','C-1:N','C:N+1']]) # WE use this as out base don't need to load it again
    geoLists.append(['2ANGSALL_' + pdbSet, ['TAU','C-1:N:CA','CA:C:N+1','CA:C:O','O:C:N+1','CA:C:N+1']])
    geoLists.append(['3DIHSALL_' + pdbSet, ['PHI','PSI','OMEGA','CA-1:C-1:N:CA']])
    geoLists.append(['4DISTALL_' + pdbSet, ['N:N+1','N:C']])
    geoLists.append(['5HBALL_' + pdbSet, ['N:O-2','C:O-2','N:CA:C:O-2','N:CA:N+1:O-2']])
    geoLists.append(['6HBALL_' + pdbSet, ['N:O-3', 'C:O-3', 'N:CA:C:O-3', 'N:CA:N+1:O-3']])
    #geoLists.append(['7HBALL_' + pdbSet, ['N:O-4', 'C:O-4', 'N:CA:C:O-4', 'N:CA:N+1:O-4']])
    #geoLists.append(['8HBALL_' + pdbSet, ['N:O-5', 'C:O-5', 'N:CA:C:O-5', 'N:CA:N+1:O-5']])
    #geoLists.append(['9HBOALL_' + pdbSet, ['N:{O}','C:{O}','N:CA:C:{O}','N:CA:N+1:{O}']])
    #geoLists.append(['10CISALL_' + pdbSet, ['CA-1:C-1:N:CA', 'CA-1:CA']])
    #geoLists.append(['7WAT', ['N:HOH','C:HOH','N:CA:C:HOH','N:CA:N+1:HOH']])
    #geoLists.append(['8XTRA', ['N:HETATM']])

    for file in geoLists:
        print('Merging with file',file[0])
        fileName = 'CsvGeos_BEST_Set' + file[0] + '.csv'
        fileList = file[1]
        fileList.append('ID')
        print(fileList)
        dataRow = pd.read_csv(loadPath + fileName)
        dataRow['rid'] = dataRow['rid'].astype(str)
        dataRow['ID'] = dataRow['pdbCode'] + dataRow['chain'] + dataRow['rid'] + dataRow['aa']
        dataRow = dataRow[fileList]
        dataFirst = pd.merge(dataFirst,dataRow,left_on='ID',right_on='ID')
        dataFirst = dataFirst.dropna()

    filePath = printPath + 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'
    print('...printing', filePath)
    dataFirst.to_csv(filePath, index=False)
