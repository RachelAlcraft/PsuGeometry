
'''
In this file we compare old and adjusted values
With resolution and software
'''

import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

print('### LOADING csv files ###')
dataMerged = pd.read_csv(help.loadPath + "MergedEvidenced.csv")

#dataMerged['PDB'] = dataMerged['pdbCode']
#pdbdata = pd.read_csv('../../PdbLists/Pdbs_100.csv')
#pdbdata = pdbdata[['PDB', 'SOFTWARE', 'RES']]
#dataMerged = dataMerged.set_index('PDB').join(pdbdata.set_index('PDB'))
#dataMerged['SOFTWARE'] = dataMerged['SOFTWARE'].str[:8]
dataMerged = dataMerged.dropna()
dataMergedRes = dataMerged.query('RES <= 0.9')
dataMergedRes.to_csv(help.loadPath + "MergedEvidenced09.csv", index=False)

georep05 = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)
georep05Res = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

print('### Creating reports ###')
georep05.addScatter(data=dataMerged, geoX='N:CA_Diff',geoY='RES', title='Change in N:CA', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep05.addScatter(data=dataMerged, geoX='CA:C_Diff',geoY='RES', title='Change in CA:C', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep05.addScatter(data=dataMerged, geoX='C:O_Diff',geoY='RES', title='Change in C=O', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep05.addScatter(data=dataMerged, geoX='C:N+1_Diff',geoY='RES', title='Change in C:N+1', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')

georep05.addScatter(data=dataMerged, geoX='N:CA_Diff',geoY='SOFTWARE', title='Change in N:CA', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep05.addScatter(data=dataMerged, geoX='CA:C_Diff',geoY='SOFTWARE', title='Change in CA:C', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep05.addScatter(data=dataMerged, geoX='C:O_Diff',geoY='SOFTWARE', title='Change in C=O', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep05.addScatter(data=dataMerged, geoX='C:N+1_Diff',geoY='SOFTWARE', title='Change in C:N+1', hue='RES',categorical=False,sort='DESC',palette='viridis_r')

georep05.addScatter(data=dataMerged, geoX='N:CA_Orig',geoY='N:CA_Adj', title='Change in N:CA', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[1.2,1.6])
georep05.addScatter(data=dataMerged, geoX='CA:C_Orig',geoY='CA:C_Adj', title='Change in CA:C', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[1.4,1.8])
georep05.addScatter(data=dataMerged, geoX='C:O_Orig',geoY='C:O_Adj', title='Change in C=O', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[0.9,1.45])
georep05.addScatter(data=dataMerged, geoX='C:N+1_Orig',geoY='C:N+1_Adj', title='Change in C:N+1', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[1.05,1.55])

georep05.addScatter(data=dataMerged, geoX='N:CA_Orig',geoY='N:CA_Adj', title='Change in N:CA', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[1.1,1.6])
georep05.addScatter(data=dataMerged, geoX='CA:C_Orig',geoY='CA:C_Adj', title='Change in CA:C', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[1.3,1.75])
georep05.addScatter(data=dataMerged, geoX='C:O_Orig',geoY='C:O_Adj', title='Change in C=O', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[0.9,1.45])
georep05.addScatter(data=dataMerged, geoX='C:N+1_Orig',geoY='C:N+1_Adj', title='Change in C:N+1', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[1.05,1.55])

georep05.printToHtml('Reporting changed values', 4, 'changes')

print('### Creating reports for higher resolution ###')
georep05.addScatter(data=dataMergedRes, geoX='N:CA_Diff',geoY='RES', title='Change in N:CA', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep05.addScatter(data=dataMergedRes, geoX='CA:C_Diff',geoY='RES', title='Change in CA:C', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep05.addScatter(data=dataMergedRes, geoX='C:O_Diff',geoY='RES', title='Change in C=O', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')
georep05.addScatter(data=dataMergedRes, geoX='C:N+1_Diff',geoY='RES', title='Change in C:N+1', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r')

georep05.addScatter(data=dataMergedRes, geoX='N:CA_Diff',geoY='SOFTWARE', title='Change in N:CA', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep05.addScatter(data=dataMergedRes, geoX='CA:C_Diff',geoY='SOFTWARE', title='Change in CA:C', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep05.addScatter(data=dataMergedRes, geoX='C:O_Diff',geoY='SOFTWARE', title='Change in C=O', hue='RES',categorical=False,sort='DESC',palette='viridis_r')
georep05.addScatter(data=dataMergedRes, geoX='C:N+1_Diff',geoY='SOFTWARE', title='Change in C:N+1', hue='RES',categorical=False,sort='DESC',palette='viridis_r')

georep05.addScatter(data=dataMergedRes, geoX='N:CA_Orig',geoY='N:CA_Adj', title='Change in N:CA', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[1.2,1.6])
georep05.addScatter(data=dataMergedRes, geoX='CA:C_Orig',geoY='CA:C_Adj', title='Change in CA:C', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[1.4,1.8])
georep05.addScatter(data=dataMergedRes, geoX='C:O_Orig',geoY='C:O_Adj', title='Change in C=O', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[0.9,1.45])
georep05.addScatter(data=dataMergedRes, geoX='C:N+1_Orig',geoY='C:N+1_Adj', title='Change in C:N+1', hue='SOFTWARE',categorical=True,sort='NON',palette='jet_r',range=[1.05,1.55])

georep05.addScatter(data=dataMergedRes, geoX='N:CA_Orig',geoY='N:CA_Adj', title='Change in N:CA', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[1.2,1.6])
georep05.addScatter(data=dataMergedRes, geoX='CA:C_Orig',geoY='CA:C_Adj', title='Change in CA:C', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[1.4,1.8])
georep05.addScatter(data=dataMergedRes, geoX='C:O_Orig',geoY='C:O_Adj', title='Change in C=O', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[0.9,1.45])
georep05.addScatter(data=dataMergedRes, geoX='C:N+1_Orig',geoY='C:N+1_Adj', title='Change in C:N+1', hue='RES',categorical=False,sort='DESC',palette='viridis_r',range=[1.05,1.55])

georep05.printToHtml('Reporting changed values Resolution <= 0.9', 4, 'changes_res')



