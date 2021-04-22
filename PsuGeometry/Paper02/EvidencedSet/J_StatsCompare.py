# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
'''
EH stats report coparison
'''

def statsCompare(pdbSetA,pdbSetB):

    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/' + pdbSetA + '/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataJ/'

    AFileName = 'Data_DefensibleWithGeosALL_' + pdbSetA + '.csv'
    BFileName = 'Data_DefensibleWithGeosALL_' + pdbSetB + '.csv'
    dataBestA = pd.read_csv(loadPath + AFileName)
    dataBestB = pd.read_csv(loadPath + BFileName)

    aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    geos = ['N:CA','CA:C','C:O','C:N+1','TAU','CA:C:N+1','CA:C:O','O:C:N+1','C-1:N:CA']

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

    for geo in geos:

        outliersA = dataBestA.sort_values(by=[geo])
        outliersA = outliersA.iloc[[0, -1]]
        outliersRA = outliersA[geo].values

        outliersB = dataBestB.sort_values(by=[geo])
        outliersB = outliersB.iloc[[0, -1]]
        outliersRB = outliersB[geo].values

        range = [min(outliersRA[0],outliersRB[0]), max(outliersRA[1],outliersRB[1])]

        for aa in aas:
            dataCutA = dataBestA.query("aa ==  '" + aa + "'")
            dataCutB = dataBestB.query("aa ==  '" + aa + "'")

            georep.addHistogram(data=dataCutA, geoX=geo, title=aa + ' ' + pdbSetA,  range=range)
            georep.addStatsCompare(dataA=dataCutA, dataB=dataCutB,  descA=aa + ' ' + pdbSetA, descB=aa + ' ' + pdbSetB,  geoX=geo, title=aa + ' ' + geo + ' compare')
            georep.addHistogram(data=dataCutB, geoX=geo, title=aa + ' ' + pdbSetB, range=range)



    georep.printToHtml('Stats Compare , set=' + pdbSetA + ' vs' + pdbSetB, 3, 'StatsCompare_' + pdbSetA + '_' + pdbSetB)
