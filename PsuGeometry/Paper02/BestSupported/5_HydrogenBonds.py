# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
'''
EH stats report coparison
'''
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/Reports/'

BestFileName = 'Set4_BestWithGeosALL.csv'
dataBest = pd.read_csv(loadPath + BestFileName)

aas = ['GLY']
geos = ['X']

georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
'''
'N:O-2','C:O-2','N:CA:C:O-2','N:CA:N+1:O-2'
'''
for geo in geos:
    for aa in aas:
        #prepare E&H comparison values
        bestCut = dataBest.query("aa ==  '" + aa + "'")

        georep.addScatter(data=bestCut, geoX='PHI', geoY='PSI', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='PHI', geoY='PSI', hue='N:O-2', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='PHI', geoY='PSI', hue='dssp', title='', palette='Dark2_r', categorical=True)

        georep.addScatter(data=bestCut, geoX='PSI', geoY='N:N+1', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='PSI', geoY='N:N+1', hue='N:O-2', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='PSI', geoY='N:N+1', hue='dssp', title='', palette='Dark2_r', categorical=True)

        georep.addScatter(data=bestCut, geoX='N:O-2', geoY='N:CA:C:O-2', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:O-2', geoY='N:CA:C:O-2', hue='PSI', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:O-2', geoY='N:CA:C:O-2', hue='dssp', title='', palette='Dark2_r', categorical=True)

        georep.addScatter(data=bestCut, geoX='N:O-2', geoY='N:CA:N+1:O-2', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:O-2', geoY='N:CA:N+1:O-2', hue='PSI', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:O-2', geoY='N:CA:N+1:O-2', hue='dssp', title='', palette='Dark2_r', categorical=True)

        georep.addScatter(data=bestCut, geoX='N:CA:C:O-2', geoY='N:CA:N+1:O-2', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:CA:C:O-2', geoY='N:CA:N+1:O-2', hue='PSI', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:CA:C:O-2', geoY='N:CA:N+1:O-2', hue='dssp', title='', palette='Dark2_r', categorical=True)

        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:CA:C:{O}', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:CA:C:{O}', hue='PSI', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:CA:C:{O}', hue='dssp', title='', palette='Dark2_r', categorical=True)

        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:CA:N+1:{O}', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:CA:N+1:{O}', hue='PSI', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:CA:N+1:{O}', hue='dssp', title='', palette='Dark2_r',categorical=True)

        georep.addScatter(data=bestCut, geoX='N:CA:C:{O}', geoY='N:CA:N+1:{O}', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:CA:C:{O}', geoY='N:CA:N+1:{O}', hue='PSI', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:CA:C:{O}', geoY='N:CA:N+1:{O}', hue='dssp', title='', palette='Dark2_r',categorical=True)

        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:O-2', hue='TAU_x', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:O-2', hue='PSI', title='', palette='jet')
        georep.addScatter(data=bestCut, geoX='N:{O}', geoY='N:O-2', hue='dssp', title='', palette='Dark2_r',categorical=True)


georep.printToHtml('Hydrogen Bond Correlations', 3, 'BS_HB')