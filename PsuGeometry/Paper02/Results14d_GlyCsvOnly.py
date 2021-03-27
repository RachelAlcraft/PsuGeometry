# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''

fileName = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/DataCsvOfSets/GLY_Categories.csv'

data = pd.read_csv(fileName)
data['MOTIF'] = data['aa-1'] + data['aa+1']
data = data[data['CA:HOH'] != 0]
#data['hetatm'] = data['hetatm'][:2]
#data = data.dropna()

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

georep.addScatter(data=data, geoX='PHI', geoY='PSI', hue='TAU', title='',palette='jet', categorical=False)
georep.addScatter(data=data, geoX='PHI', geoY='PSI', hue='N:O-2', title='',palette='jet', categorical=False)
georep.addScatter(data=data, geoX='PHI', geoY='PSI', hue='N-1:CA:C', title='',palette='jet', categorical=False)

georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='TAU', title='',palette='jet', categorical=False)
georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='N:O-2', title='',palette='jet', categorical=False)
georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='N-1:CA:C', title='',palette='jet')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='N:CA:C:O-2', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='N:O-2', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='CA:HOH', hue='N:CA:C:HOH', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='CA:HOH', geoY='N:CA:C:HOH', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:HOH', geoY='TAU', hue='CA:HOH', title='',palette='jet', sort='NON')

georep.addScatter(data=dataHetatm, geoX='TAU', geoY='CA:HETATM', hue='N:CA:C:HETATM', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='CA:HETATM', geoY='N:CA:C:HETATM', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='N:CA:C:HETATM', geoY='TAU', hue='CA:HETATM', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='N:N+1', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='PSI', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='PSI', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='PSI', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='aa-2', title='',palette='tab20', categorical=True)
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='aa-1', title='',palette='tab20', categorical=True)
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='aa+1', title='',palette='tab20', categorical=True)

georep.addScatter(data=data, geoX='TAU', geoY='CA:HOH', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='CA:HOH', geoY='N:CA:C:HOH', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:HOH', geoY='TAU', hue='N:N+1', title='',palette='jet', sort='NON')

georep.addScatter(data=dataHetatm, geoX='TAU', geoY='CA:HETATM', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='CA:HETATM', geoY='N:CA:C:HETATM', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='N:CA:C:HETATM', geoY='TAU', hue='N:N+1', title='',palette='jet', sort='NON')

####Row 11
georep.addScatter(data=data, geoX='PSI', geoY='N-1:CA:C', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='PHI', geoY='N-1:CA:C', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='TAU', geoY='N-1:CA:C', hue='N:N+1', title='',palette='jet', sort='NON')

####Row 12
georep.addScatter(data=data, geoX='PHI', geoY='N:O-2', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='TAU', geoY='N-1:CA:C', hue='PHI', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='TAU', geoY='N-1:CA:C', hue='N:O-2', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='N:N+1', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='CA:HOH', geoY='N:HOH:C', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:HOH', geoY='CA:HOH', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:HOH:C', geoY='N:CA:C:HOH', hue='TAU', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='CA:HOH', geoY='N:HOH:C', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:HOH', geoY='CA:HOH', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:HOH:C', geoY='N:CA:C:HOH', hue='N:N+1', title='',palette='jet', sort='NON')

georep.addScatter(data=dataHetatm, geoX='CA:HETATM', geoY='N:HETATM:C', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='N:CA:C:HETATM', geoY='CA:HETATM', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='N:HETATM:C', geoY='N:CA:C:HETATM', hue='TAU', title='',palette='jet', sort='NON')

georep.addScatter(data=dataHetatm, geoX='CA:HETATM', geoY='N:HETATM:C', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='N:CA:C:HETATM', geoY='CA:HETATM', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=dataHetatm, geoX='N:HETATM:C', geoY='N:CA:C:HETATM', hue='N:N+1', title='',palette='jet', sort='NON')

georep.printToHtml('Results 14d. Gly Best Supported Plots', 3, 'Results14d_GLYBest')
