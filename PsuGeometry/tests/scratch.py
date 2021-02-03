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
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/'

pdbList = []
pdbdata = pd.read_csv('structures09.csv')
pdbList = pdbdata['pdb_code']
pdbList = pdbList[:50]


georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)

geoList = ['CA:O','O:N+1','N+1:CA','PSI','CA:C:O:N+1','TAU','PHI']
hueList = ['aa','bfactor','rid','resolution','pdbCode','dssp','ridx','atomNo']

data = georep.getGeoemtryCsv(geoList, hueList)

georep.addScatter(data=data, geoX='CA:O',geoY='O:N+1',hue='PSI', title='', palette='jet',sort='NON')
georep.addScatter(data=data, geoX='O:N+1',geoY='N+1:CA',hue='PSI', title='', palette='jet',sort='NON')
georep.addScatter(data=data, geoX='N+1:CA',geoY='CA:O',hue='PSI', title='', palette='jet',sort='NON')

georep.addScatter(data=data, geoX='CA:C:O:N+1',geoY='PSI',hue='TAU', title='', palette='jet',sort='NON',operation='ABS')
georep.addScatter(data=data, geoX='CA:C:O:N+1',geoY='TAU',hue='PSI', title='', palette='jet',sort='NON',operation='ABS')
georep.addScatter(data=data, geoX='CA:C:O:N+1',geoY='TAU',hue='PHI', title='', palette='jet',sort='NON',operation='ABS')

georep.addScatter(data=data, geoX='TAU',geoY='dssp',hue='atomNo', title='', palette='Spectral',sort='NON')
georep.addScatter(data=data, geoX='ridx',geoY='TAU',hue='dssp', title='', palette='Spectral',sort='NON')
georep.addScatter(data=data, geoX='atomNo',geoY='TAU',hue='dssp', title='', palette='Spectral',sort='NON')


# And finally create the reort with a file name of choice
title = ''
georep.printToHtml(title,3,'scratch')


