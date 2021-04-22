# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
'''
Compare sets and EH and Jaskolski
'''

def compareSets(tag):

    pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataD/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataE/'

    geos = ['N:CA', 'CA:C', 'C:O', 'C:N+1', 'TAU', 'C-1:N:CA', 'CA:C:N+1', 'CA:C:O', 'O:C:N+1', 'CA:C:N+1']
    aas = ['ALL', 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN','ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    fileName = tag + 'Data_SetsSummaryMerged.csv'
    data = pd.read_csv(loadPath + fileName)

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

    for geo in geos:
        dataCut = data.query('geo == "' + geo + '"')
        dataCutCount = dataCut.query('count > 0')
        dataCutCount = dataCutCount.query('aa != "' + 'ALL' + '"')
        dataCutAll = dataCut.query('aa == "' + 'ALL' + '"')
        dataCutPRO = dataCut.query('aa == "' + 'PRO' + '"')
        dataCutGLY = dataCut.query('aa == "' + 'GLY' + '"')

        georep.addScatter(data=dataCutAll, geoX='mean', geoY='set', hue='sd', title=geo + ' ALL (exc gly/pro)', palette='jet', categorical=False,sort='NON')
        georep.addScatter(data=dataCutGLY, geoX='mean', geoY='set', hue='sd', title=geo + ' GLY', palette='jet', categorical=False,sort='NON')
        georep.addScatter(data=dataCutPRO, geoX='mean', geoY='set',hue='sd', title=geo + ' PRO', palette='jet', categorical=False,sort='NON')

        georep.addScatter(data=dataCutCount, geoX='count', geoY='aa', hue='set', title=geo + ' Best Supported counts per aa', palette='jet_r', categorical=True, sort='NON')
        georep.addScatter(data=dataCutCount, geoX='mean', geoY='aa', hue='set', title=geo + ' Best Supported means per aa', palette='jet_r',categorical=True, sort='NON')
        georep.addScatter(data=dataCutCount, geoX='sd', geoY='aa', hue='set', title=geo + ' Best Supported sd per aa', palette='jet_r', categorical=True, sort='NON')



    georep.printToHtml('Best Supported and Engh&Huber Compare', 3, tag + 'Compare_EH_Sets')

