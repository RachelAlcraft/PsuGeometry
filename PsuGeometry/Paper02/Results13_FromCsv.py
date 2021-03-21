# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''

fileName = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/DataCsvOfSets/GoodBetterTauSets.csv'

data = pd.read_csv(fileName)


pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'

georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

data0 = data.query('TAU_Diff == 0')

georep.addScatter(data=data, geoX='PHI', geoY='PSI', hue='TAU', title='Full',palette='jet', sort='NON')
georep.addScatter(data=data0, geoX='PHI', geoY='PSI', hue='TAU', title='Best Supported',palette='jet', sort='NON')

georep.addHistogram(data=data, geoX='TAU', title='Full')
georep.addHistogram(data=data0, geoX='TAU', title='Best Supported')

georep.addHistogram(data=data, geoX='N:N+1', title='Full')
georep.addHistogram(data=data0, geoX='N:N+1', title='Best Supported')

georep.addHistogram(data=data, geoX='Abs(PSI)', title='Full, Abs Val')
georep.addHistogram(data=data0, geoX='Abs(PSI)', title='Best Supported, Abs Val')

georep.addHistogram(data=data, geoX='Abs(PHI)', title='Full, Abs Val')
georep.addHistogram(data=data0, geoX='Abs(PHI)', title='Best Supported, Abs Val')

georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='TAU', title='Full',palette='jet', sort='NON')
georep.addScatter(data=data0, geoX='PSI', geoY='N:N+1', hue='TAU', title='Best Supported',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='PSI', hue='N:N+1', title='Full',palette='jet', sort='NON')
georep.addScatter(data=data0, geoX='TAU', geoY='PSI', hue='N:N+1', title='Best Supported',palette='jet', sort='NON')

georep.addScatter(data=data, geoX='TAU', geoY='N:N+1', hue='PSI', title='Full',palette='jet', sort='NON')
georep.addScatter(data=data0, geoX='TAU', geoY='N:N+1', hue='PSI', title='Best Supported',palette='jet', sort='NON')

# all replicated with "category"
georep.addScatter(data=data, geoX='PHI', geoY='PSI', hue='Category', title='Full',palette='tab10', sort='NON',categorical=True)
georep.addScatter(data=data0, geoX='PHI', geoY='PSI', hue='Category', title='Best Supported',palette='tab10', sort='NON',categorical=True)

georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='Category', title='Full',palette='tab10', sort='NON',categorical=True)
georep.addScatter(data=data0, geoX='PSI', geoY='N:N+1', hue='Category', title='Best Supported',palette='tab10', sort='NON',categorical=True)

georep.addScatter(data=data, geoX='TAU', geoY='PSI', hue='Category', title='Full',palette='tab10', sort='NON',categorical=True)
georep.addScatter(data=data0, geoX='TAU', geoY='PSI', hue='Category', title='Best Supported',palette='tab10', sort='NON',categorical=True)

georep.addScatter(data=data, geoX='TAU', geoY='N:N+1', hue='Category', title='Full',palette='tab10', sort='NON',categorical=True)
georep.addScatter(data=data0, geoX='TAU', geoY='N:N+1', hue='Category', title='Best Supported',palette='tab10', sort='NON',categorical=True)


georep.printToHtml('Results 13. Tau Plots from 0 differences', 2, 'Results13_0diff')
