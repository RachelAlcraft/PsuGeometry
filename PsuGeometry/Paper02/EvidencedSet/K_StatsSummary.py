# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
'''
EH stats report summary
'''


def applyCis(aa,preomega):
    if aa != 'PRO':
        return aa
    if abs(preomega) > 100:
        return 'PRO'
    else:
        return 'CISPRO'



def statsSummary(pdbSet, data, geos,tag):

    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataK/'

    fileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'

    allAtoms = False
    bFactorFactor = -1
    if pdbSet == 'RESTRICTED':
        allAtoms = True
        bFactorFactor = 1.3

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

    for geo in geos:
        georep.addStatsSummary(data=data, desc=geo + ' ' + pdbSet, geoX=geo, geoY='aa', hue='ID')


    georep.printToHtml('Stats Summary , set=' + pdbSet, 2, 'StatsSummary_' + pdbSet + tag)
