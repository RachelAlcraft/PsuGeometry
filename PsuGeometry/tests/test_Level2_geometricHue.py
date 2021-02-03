# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as geor
import pandas as pd
'''
This script runs a correlation report on 142 ultrahigh-resolution structures <=0.9A
It runs correlation reports on proline, which it colours on the hue of CHI1 and CA-1:CA
These geometric measures are proxies for up/down pucker of the proline run (up-pucker=-ve CHI1)
And cis-trans proline, where pre-omega means proline, which corresponds to short CA-1:CA
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/1000Structures/'

pdbList = []
pdbdata = pd.read_csv('structures09.csv')
pdbList = pdbdata['pdb_code']


georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False)

geoList = ['PHI','PSI','TAU','C-1:C','C-1:N:CA','CHI1','CA-1:CA','OMEGA','CA:CA+1','C-1:N:CA:C']
hueList = ['aa','bfactor','rid','resolution','pdbCode']

data = georep.getGeoemtryCsv(geoList, hueList)

georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='CHI1', title='Ramachandran', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='CA-1:CA', title='Ramachandran', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='CHI1', title='Omega-Phi', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='CA-1:CA', title='Omega-Phi', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='CHI1', title='C-1:C/Phi', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='CA-1:CA', title='C-1:C/Phi', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='CHI1', title='Tau angles', palette='viridis_r',sort='DESC',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='CA-1:CA', title='Tau angles', palette='Spectral',sort='DESC',restrictions={'aa':'PRO'})

# And finally create the reort with a file name of choice
title = 'Correlations with geometric hue<br/> Proline up-pucker by proxy of CHI1<br/>Proline cis/trans by proxy of short Calpha'
georep.printToHtml(title,2,'geohue')


