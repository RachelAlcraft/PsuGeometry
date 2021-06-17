# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
EH stats report coparison
'''

def statsCompare(pdbSetA, dataA,pdbSetB,dataB,isAll):


    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSetA + '/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataO/'

    if isAll:
        aas= ['ALL']
    else:
        aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    geos = ['N:CA','CA:C','C:O','C:N+1','TAU','CA:C:N+1','CA:C:O','O:C:N+1','C-1:N:CA']

    for aa in aas:
        georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
        if aa == 'ALL':
            dataCutA = dataA
            dataCutB = dataB
        else:
            dataCutA = dataA.query("aa ==  '" + aa + "'")
            dataCutB = dataB.query("aa ==  '" + aa + "'")

        for geo in geos:
            outliersA = dataA.sort_values(by=[geo])
            outliersA = outliersA.iloc[[0, -1]]
            outliersRA = outliersA[geo].values

            outliersB = dataB.sort_values(by=[geo])
            outliersB = outliersB.iloc[[0, -1]]
            outliersRB = outliersB[geo].values

            range = [min(outliersRA[0],outliersRB[0]), max(outliersRA[1],outliersRB[1])]

            georep.addHistogram(data=dataCutA, geoX=geo, title=aa + ' ' + pdbSetA,  range=range)
            georep.addStatsCompare(dataA=dataCutA, dataB=dataCutB,  descA=aa + ' ' + pdbSetA, descB=aa + ' ' + pdbSetB,  geoX=geo, title=aa + ' ' + geo + ' compare')
            georep.addHistogram(data=dataCutB, geoX=geo, title=aa + ' ' + pdbSetB, range=range)

        georep.printToHtml(aa + ' Stats Compare , set=' + pdbSetA + ' vs ' + pdbSetB, 3, 'StatsCompare_' + pdbSetA + '_' + pdbSetB + '_' + aa)
