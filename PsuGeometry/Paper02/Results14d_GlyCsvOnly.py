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
#data = data.dropna()


pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'

georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

#Available GEOs
'''
pdbCode,chain,rid,aa,id,aa-2,aa-1,aa+1,bfactor,bfactorRatio,disordered,
N:N+1,TAU,PSI,PHI,N:O-2,N:CA:C:O-2,N:C,CA:C,C:O,N:CA,C-1:N,C:N+1,OMEGA,CA:C:O:N+1,O:N+1,CA:O,CA:N+1,CA:C:N+1,C-1:N:CA
'''

#georep.addScatter(data=data, geoX='aa-2', geoY='aa-1', hue='TAU', title='',palette='jet', sort='NON')
#georep.addScatter(data=data, geoX='aa-2', geoY='aa+1', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='TAU', title='',palette='jet', categorical=False)
georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='N:O-2', title='',palette='jet', categorical=False)
georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='N:CA:C:O-2', title='',palette='jet', categorical=False)

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='N:CA:C:O-2', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='TAU', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='N:O-2', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='N:N+1', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='N:N+1', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='PSI', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='PSI', title='',palette='jet', sort='NON')
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='PSI', title='',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:O-2', hue='aa-2', title='',palette='tab20', categorical=True)
georep.addScatter(data=data, geoX='N:O-2', geoY='N:CA:C:O-2', hue='aa-1', title='',palette='tab20', categorical=True)
georep.addScatter(data=data, geoX='N:CA:C:O-2', geoY='TAU', hue='aa+1', title='',palette='tab20', categorical=True)

georep.printToHtml('Results 14d. Gly Best Supported Plots', 3, 'Results14d_GLYBest')
