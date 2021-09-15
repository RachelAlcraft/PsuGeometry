
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

dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bbden_adjusted.csv")
dataPdbLap = pd.read_csv(help.loadPath + "bblap_adjusted.csv")

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

print('### Creating reports ###')
georep.addScatter(data=dataPdbCut, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='dssp',categorical=True,sort='RAND',palette='jet_r')
georep.addScatter(data=dataPdbAdj, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='dssp',categorical=True,sort='RAND',palette='jet_r')
georep.addScatter(data=dataPdbLap, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='dssp',categorical=True,sort='RAND',palette='jet_r')

georep.addScatter(data=dataPdbCut, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep.addScatter(data=dataPdbAdj, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep.addScatter(data=dataPdbLap, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='RES',categorical=False,sort='DESC',palette='viridis_r')

georep.addProbability(data=dataPdbCut, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff',palette='cubehelix_r')
georep.addProbability(data=dataPdbAdj, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff',palette='cubehelix_r')
georep.addProbability(data=dataPdbLap, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff',palette='cubehelix_r')

georep.addScatter(data=dataPdbCut, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='bfactor',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=dataPdbAdj, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='bfactor',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=dataPdbLap, geoX='C:O',geoY='Dist_Den_Lap', title='C:O against maxima diff', hue='bfactor',categorical=False,sort='RAND',palette='jet')

georep.printToHtml('Reporting changed values', 3, 'maxadiff')
