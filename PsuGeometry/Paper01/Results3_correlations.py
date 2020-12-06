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
#pdbList1000 = pdbList1000[:10]
geoList = ['PSI','N:O','CB:O']
hueList = ['resolution','dssp']

georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)
data = georep.getGeoemtryCsv(geoList,hueList)

georep.addScatter(data=data, geoX='PSI',geoY='N:O',hue='dssp', title='Psi-N:O', palette='Set1',sort='NON')
georep.addScatter(data=data, geoX='PSI',geoY='CB:O',hue='dssp', title='Psi-CB:O', palette='Set1',sort='NON')
georep.addScatter(data=data, geoX='N:O',geoY='CB:O',hue='dssp', title='N:O-CB:O', palette='Set1',sort='NON')
# just 1i1w
georep.pdbCodes = ['1i1w']
geoList = ['PSI','N:O','CB:O','PHI','N:CA','CA:C','C-1:C']
hueList = ['resolution','aa']
data = georep.getGeoemtryCsv(geoList,hueList)
georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='aa', title='Psi-N:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='PSI',geoY='N:O',hue='aa', title='Psi-N:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='PSI',geoY='CB:O',hue='aa', title='Psi-CB:O', palette='gist_rainbow',ghost=True)

georep.addScatter(data=data, geoX='N:CA',geoY='CA:C',hue='aa', title='Psi-N:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='N:O',geoY='CB:O',hue='aa', title='N:O-CB:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='PHI',geoY='C-1:C',hue='aa', title='Phi-C-1:C', palette='gist_rainbow',ghost=True)

title = 'The Rarity Effect'
georep.printToHtml(title,3,'Results3_correlations')