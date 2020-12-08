# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
This script looks at the rarity effect
By looking at some plots with against resolution and probability density
This is 1i1w only
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

geoList = ['PSI','N:O','CB:O','PHI','N:CA','CA:C','C-1:C']
hueList = ['resolution','aa','2FoFc','2FoFc','dssp','bfactor']

georep = psu.GeoReport(['1i1w'], pdbDataPath, edDataPath, printPath, ed=True, dssp=True, includePdbs=False)
data = georep.getGeoemtryCsv(geoList,hueList)

georep.addScatter(data=data, geoX='PSI',geoY='N:O',hue='dssp', title='Psi-N:O', palette='Set1',sort='NON')
georep.addScatter(data=data, geoX='PSI',geoY='CB:O',hue='dssp', title='Psi-CB:O', palette='Set1',sort='NON')
georep.addScatter(data=data, geoX='N:O',geoY='CB:O',hue='dssp', title='N:O-CB:O', palette='Set1',sort='NON')

georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='aa', title='Psi-N:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='PSI',geoY='N:O',hue='aa', title='Psi-N:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='PSI',geoY='CB:O',hue='aa', title='Psi-CB:O', palette='gist_rainbow',ghost=True)

georep.addScatter(data=data, geoX='N:CA',geoY='CA:C',hue='aa', title='Psi-N:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='N:O',geoY='CB:O',hue='aa', title='N:O-CB:O', palette='gist_rainbow',ghost=True)
georep.addScatter(data=data, geoX='PHI',geoY='C-1:C',hue='aa', title='Phi-C-1:C', palette='gist_rainbow',ghost=True)

georep.addScatter(data=data, geoX='PSI',geoY='CB:O',hue='2FoFc', title='Psi-CB:O', palette='inferno_r',ghost=True)
georep.addScatter(data=data, geoX='N:O',geoY='CB:O',hue='bfactor', title='N:O-CB:O', palette='inferno',ghost=True)
georep.addScatter(data=data, geoX='PHI',geoY='C-1:C',hue='2FoFc', title='Phi-C-1:C', palette='inferno_r',ghost=True)

title = 'The Rarity Effect'
georep.printToHtml(title,3,'Results3_correlations_1i1w')