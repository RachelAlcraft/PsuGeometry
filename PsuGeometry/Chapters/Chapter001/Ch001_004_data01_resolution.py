
'''
In this file we compare old and adjusted values
With resolution and software
'''

import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

print('### LOADING csv files ###')
dataMerged = pd.read_csv(help.loadPath + "MergedEvidenced.csv")
dataMerged = dataMerged.dropna()
dataMergedRes = dataMerged.query('RES <= 0.9')
dataMergedRes.to_csv(help.loadPath + "MergedEvidenced09.csv", index=False)

dataPdbCut = pd.read_csv(help.loadPath + "bb_reduced.csv")
dataPdbAdj = pd.read_csv(help.loadPath + "bbden_adjusted.csv")
dataPdbLap = pd.read_csv(help.loadPath + "bblap_adjusted.csv")

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

print('### Creating reports ###')
georep.addScatter(data=dataPdbCut, geoX='N:CA_Orig',geoY='RES', title='Change in N:CA', hue='SOFTWARE',categorical=False,sort='DESC',palette='viridis_r',range=[1.1,1.6])
georep.addScatter(data=dataPdbCut, geoX='CA:C_Orig',geoY='RES', title='Change in CA:C', hue='SOFTWARE',categorical=False,sort='DESC',palette='viridis_r',range=[1.3,1.75])
georep.addScatter(data=dataPdbCut, geoX='C:O_Orig',geoY='RES', title='Change in C=O', hue='SOFTWARE',categorical=False,sort='DESC',palette='viridis_r',range=[0.9,1.45])
georep.addScatter(data=dataPdbCut, geoX='C:N+1_Orig',geoY='RES', title='Change in C:N+1', hue='SOFTWARE',categorical=False,sort='DESC',palette='viridis_r',range=[1.05,1.55])

georep.addScatter(data=dataPdbAdj, geoX='N:CA_MaxDiff',geoY='RES', title='Change in N:CA', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=dataPdbAdj, geoX='CA:C_MaxDiff',geoY='RES', title='Change in CA:C', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=dataPdbAdj, geoX='C:O_MaxDiff',geoY='RES', title='Change in C=O', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=dataPdbAdj, geoX='C:N+1_MaxDiff',geoY='RES', title='Change in C:N+1', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')

georep.addScatter(data=dataPdbLap, geoX='N:CA_LapDiff',geoY='RES', title='Change in N:CA', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=dataPdbLap, geoX='CA:C_LapDiff',geoY='RES', title='Change in CA:C', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=dataPdbLap, geoX='C:O_LapDiff',geoY='RES', title='Change in C=O', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=dataPdbLap, geoX='C:N+1_LapDiff',geoY='RES', title='Change in C:N+1', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')

georep.addScatter(data=dataPdbCut, geoX='N:CA_Orig',geoY='Dist_Den_Lap', title='Change in N:CA', hue='BFactor',categorical=False,sort='DESC',palette='viridis_r',range=[1.1,1.6])
georep.addScatter(data=dataPdbCut, geoX='CA:C_Orig',geoY='CA:Dist_Den_Lap', title='Change in CA:C', hue='BFactor',categorical=False,sort='DESC',palette='viridis_r',range=[1.3,1.75])
georep.addScatter(data=dataPdbCut, geoX='C:O_Orig',geoY='C:Dist_Den_Lap', title='Change in C=O', hue='BFactor',categorical=False,sort='DESC',palette='viridis_r',range=[0.9,1.45])
georep.addScatter(data=dataPdbCut, geoX='C:N+1_Orig',geoY='C:N+Dist_Den_Lap', title='Change in C:N+1', hue='BFactor',categorical=False,sort='DESC',palette='viridis_r',range=[1.05,1.55])

georep.printToHtml('Reporting changed values', 4, 'changes')



