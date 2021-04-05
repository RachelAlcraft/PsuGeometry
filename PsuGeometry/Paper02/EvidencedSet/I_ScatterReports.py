# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
'''
EH stats report coparison
'''

def scatterReports(pdbSet):

    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/EvidencedSet/DataB/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/EvidencedSet/DataI/'

    BestFileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'
    dataBest = pd.read_csv(loadPath + BestFileName)

    aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    geosPairs = [['N:CA','CA:C','C:O'],['C:N+1','TAU','CA:C:N+1'],['CA:C:O','O:C:N+1','C-1:N:CA']]

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

    for geoPair in geosPairs:
        for aa in aas:
            dataCut = dataBest.query("aa ==  '" + aa + "'")
            georep.addScatter(data=dataCut, geoX=geoPair[0], geoY=geoPair[1], hue=geoPair[2], title=aa + ':' + geoPair[0] + ':'+ geoPair[1] , palette='jet', sort='NON')


    georep.printToHtml('Scatters , set=' + pdbSet, 4, 'Defensible_Scatters_' + pdbSet)
