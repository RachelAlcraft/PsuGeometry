# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
This script looks at the rarity effect
By looking at some plots with against resolution and probability density
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList1000 = geol.GeoPdbLists().getList1000()
#pdbList1000 = pdbList1000[:50]

geoList = ['PSI','CA-1:CA:CA+1','PHI','C-1:C','TAU','CA-2:CA-1:CA','CA:CA+1:CA+2']
hueList = ['resolution']

georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)
data = georep.getGeoemtryCsv(geoList,hueList)

data = data[data['C-1:C'] < 10]


georep.addScatter(data=data, geoX='PHI',geoY='C-1:C',hue='resolution', title='Phi-C-1:CB', palette='viridis_r')
georep.addProbability(data=data, geoX='PHI',geoY='C-1:C',title='Phi-C-1:CB', palette='cubehelix_r')
georep.addScatter(data=data, geoX='PSI',geoY='CA-1:CA:CA+1',hue='resolution', title='PSI-Backbone tau', palette='viridis_r')
georep.addProbability(data=data, geoX='PSI',geoY='CA-1:CA:CA+1',title='PSI-Backbone tau', palette='cubehelix_r')
georep.addScatter(data=data, geoX='TAU',geoY='PSI',hue='resolution', title='Tau-Psi', palette='viridis_r')
georep.addProbability(data=data, geoX='TAU',geoY='PSI',title='Tau-Psi', palette='cubehelix_r')
georep.addScatter(data=data, geoX='CA-2:CA-1:CA',geoY='CA:CA+1:CA+2',hue='resolution', title='Backbone taus', palette='viridis_r')
georep.addProbability(data=data, geoX='CA-2:CA-1:CA',geoY='CA:CA+1:CA+2',title='backbone taus', palette='cubehelix_r')

title = 'The Rarity Effect'
georep.printToHtml(title,2,'Results4_rarity')