
'''
In this file we compare old and adjusted values
With LapDiff - differene between laplacian minimum and density maximum
'''

import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

print('### LOADING csv files ###')
dataMerged = pd.read_csv(help.loadPath + "MergedEvidenced.csv")
dataMerged = dataMerged.dropna()
dataMergedRes1 = dataMerged.query('RES <= 0.9')
dataMergedRes2 = dataMerged.query('RES <= 0.85')
dataMergedRes3 = dataMerged.query('RES <= 0.8')

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

print('### Creating reports ###')
georep.addScatter(data=dataMerged, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff', hue='dssp',categorical=True,sort='RAND',palette='jet_r')
georep.addScatter(data=dataMergedRes1, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.9', hue='dssp',categorical=True,sort='RAND',palette='jet_r')
georep.addScatter(data=dataMergedRes2, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.85', hue='dssp',categorical=True,sort='RAND',palette='jet_r')
georep.addScatter(data=dataMergedRes3, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.8', hue='dssp',categorical=True,sort='RAND',palette='jet_r')

georep.addScatter(data=dataMerged, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep.addScatter(data=dataMergedRes1, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.9', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep.addScatter(data=dataMergedRes2, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.85', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep.addScatter(data=dataMergedRes3, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.8', hue='RES',categorical=False,sort='DESC',palette='viridis_r')

georep.addProbability(data=dataMerged, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff',palette='cubehelix_r')
georep.addProbability(data=dataMergedRes1, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff<=0.9',palette='cubehelix_r')
georep.addProbability(data=dataMergedRes2, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff<=0.85',palette='cubehelix_r')
georep.addProbability(data=dataMergedRes3, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff<=0.8',palette='cubehelix_r')

georep.addScatter(data=dataMerged, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff', hue='bfactor',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=dataMergedRes1, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.9', hue='bfactor',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=dataMergedRes2, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.85', hue='bfactor',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=dataMergedRes3, geoX='C:O_Adj',geoY='LapDiff', title='C:O against maxima diff <=0.8', hue='bfactor',categorical=False,sort='RAND',palette='jet')


georep.printToHtml('Reporting changed values', 4, 'maxadiff')
