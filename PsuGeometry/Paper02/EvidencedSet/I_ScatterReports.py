# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
EH stats report coparison
'''

def scatterReports(pdbSet, data, trios, perAA=True, tag=''):
    import matplotlib.pyplot as plt
    plt.close('all')
    plt.clf()
    plt.cla()

    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_data/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataI/'

    #BestFileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'
    #dataBest = pd.read_csv(loadPath + BestFileName)

    aas = data['aa'].values
    aas = list(set(aas))
    aas.sort()

    #aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    #geosPairs = [['PHI','PSI','TAU'],['PSI','N:N+1','TAU'],['N:CA','CA:C','C:O'],['C:N+1','TAU','CA:C:N+1'],['CA:C:O','O:C:N+1','C-1:N:CA']]

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

    for trio in trios:
        if perAA:
            for aa in aas:
                dataCut = data.query("aa ==  '" + aa + "'")
                georep.addScatter(data=dataCut, geoX=trio[0], geoY=trio[1], hue=trio[2], title=aa + ':' + trio[0] + ':'+ trio[1] , palette='jet', sort='NON')
        else:
            georep.addScatter(data=data, geoX=trio[0], geoY=trio[1], hue=trio[2],title=trio[0] + ':' + trio[1], palette='jet', sort='NON')


    georep.printToHtml('Scatters , set=' + pdbSet, 4, 'Defensible_Scatters_' + tag)
