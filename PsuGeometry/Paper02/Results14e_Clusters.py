# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import Categoriser as cluster
import random
import pandas as pd

'''
TAU correlations
'''

fileName = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/DataCsvOfSets/GLY_Categories.csv'

data = pd.read_csv(fileName)
data['MOTIF'] = data['aa-1'] + data['aa+1']
data = data[data['CA:HOH'] != 0]
#psi, phi, NN1,NO2,NCACO2
data['CLUSTER'] = 'X'#cluster.tauCategory(data['PSI'],data['PHI'],data['N:N+1'],data['N:O-2'],data['N:CA:C:O-2'])
#Apply my cluster function
data['CLUSTER'] = data.apply(lambda row: cluster.tauCategory(row['PSI'],row['PHI'],row['N:N+1'],row['N:O-2'],row['N:CA:C:O-2'],row['O-2:C']), axis=1)

#fewer have hetatms so that is a cut down data set
dataHetatm = data.query("hetatm !=  'GLY'")
dataHetatm = dataHetatm[dataHetatm['CA:HETATM'] < 10]
dataHetatm = dataHetatm[dataHetatm['CA:HETATM'] != 0]



pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'

georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

#Available GEOs
'''
pdbCode,chain,rid,aa,id,aa-2,aa-1,aa+1,bfactor,bfactorRatio,disordered,hetatm
N:N+1,TAU,PSI,PHI,OMEGA,CA:C:O:N+1
N:C,CA:C,C:O,N:CA,C-1:N,C:N+1,O:N+1,CA:O,CA:N+1,CA:C:N+1,C-1:N:CA 
N:O-2, N:CA:C:O-2, #hydrogen bonded with 2 previous
N-1:CA:C #previous residue in line with tau?
CA:HOH, N:HOH:C, N:CA:C:HOH # anything interesting with water
CA:HETATM,N:HETATM:C,N:CA:C:HETATM #anything interesting with heavy atoms
'''
#georep.addScatter(data=data, geoX='aa-2', geoY='aa-1', hue='TAU', title='',palette='jet', sort='NON')
#georep.addScatter(data=data, geoX='aa-2', geoY='aa+1', hue='TAU', title='',palette='jet', sort='NON')

histgeos = ['TAU','PSI','PHI','N:N+1','N:O-2','N:CA','CA:C']
clusters = ['X','A1','A2','A3','A4','B1','C1']
for cat in clusters:

    dataCat = data.query("CLUSTER ==  '" + cat +"'")
    dataCatHet = dataHetatm.query("CLUSTER ==  '" + cat + "'")

    dataED = dataCat['ED Request']

    dataED.to_csv(printPath + "Results14e_EDRequest_Category_" + cat + ".csv", index=False,sep=' ')


    georep.addScatter(data=dataCat, geoX='PHI', geoY='PSI', hue='TAU', title='',palette='jet', categorical=False)
    georep.addScatter(data=dataCat, geoX='PHI', geoY='PSI', hue='O-2:N:CA:N+1', title='', palette='jet', categorical=False)
    georep.addScatter(data=dataCat, geoX='PHI', geoY='PSI', hue='N:O-2', title='', palette='jet',categorical=False)
    georep.addScatter(data=dataCat, geoX='PHI', geoY='PSI', hue='dssp', title='', palette='tab10', categorical=True)

    georep.addScatter(data=dataCat, geoX='PSI', geoY='N:N+1', hue='TAU', title='',palette='jet', categorical=False)
    georep.addScatter(data=dataCat, geoX='PSI', geoY='N:N+1', hue='O-2:N:CA:N+1', title='', palette='jet', categorical=False)
    georep.addScatter(data=dataCat, geoX='PSI', geoY='N:N+1', hue='N:O-2', title='', palette='jet', categorical=False)
    georep.addScatter(data=dataCat, geoX='PSI', geoY='N:N+1', hue='dssp', title='',palette='tab10', categorical=True)

    georep.addScatter(data=dataCat, geoX='N:O-2', geoY='O-2:N:CA:N+1', hue='TAU', title='', palette='jet')
    georep.addScatter(data=dataCat, geoX='N:O-2', geoY='O-2:N:CA:N+1', hue='PSI', title='', palette='jet')
    georep.addScatter(data=dataCat, geoX='N:O-2', geoY='O-2:N:CA:N+1', hue='N:N+1', title='', palette='jet')
    georep.addScatter(data=dataCat, geoX='N:O-2', geoY='O-2:N:CA:N+1', hue='dssp', title='', palette='tab10',categorical=True)

    for hgeo in histgeos:
        georep.addHistogram(data=dataCat, geoX=hgeo, title='')

    georep.printToHtml('Results 14e. Gly Best Supported Cluster ' + cat, 4, 'Results14e_GLYBest_' + cat)
