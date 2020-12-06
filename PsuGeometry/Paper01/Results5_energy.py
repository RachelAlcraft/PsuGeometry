# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
This script looks at the energy transitions that are indicated in the ordered plots
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList1000 = geol.GeoPdbLists().getList1000()
#pdbList1000 = pdbList1000[:100]

geoList = ['PSI','N:O']
hueList = ['resolution']

georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)
data = georep.getGeoemtryCsv(geoList,hueList)

georep.addScatter(data=data, geoX='PSI',geoY='N:O',hue='resolution', title='', palette='viridis_r',restrictions={'aa':'GLY'})
georep.addScatter(data=data, geoX='PSI',geoY='N:O',hue='resolution', title='', palette='viridis_r',restrictions={'aa':'PRO'})
georep.addScatter(data=data, geoX='PSI',geoY='N:O',hue='resolution', title='', palette='viridis_r',restrictions={'aa':'ALA'})

georep.addProbability(data=data, geoX='PSI',geoY='N:O',title='', palette='cubehelix_r',restrictions={'aa':'GLY'})
georep.addProbability(data=data, geoX='PSI',geoY='N:O',title='', palette='cubehelix_r',restrictions={'aa':'PRO'})
georep.addProbability(data=data, geoX='PSI',geoY='N:O',title='', palette='cubehelix_r',restrictions={'aa':'ALA'})

title = 'Energy'
georep.printToHtml(title,3,'Results5_energy')