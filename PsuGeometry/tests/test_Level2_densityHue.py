# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as geor
'''
This script runs a correlation report on a few structures to demonstrate the use of 2Fo-Fc as a hue
The data is precalculated into a datafame
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Levels/'

pdbList=['1ejg','1us0','1tt8','1i1w','1ucs','6jvv','5nqo']

georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=True, dssp=True)

geoList = ['PHI','PSI','TAU','C-1:C','C-1:N:CA','CHI1','CA-1:CA','OMEGA','CA:CA+1','C-1:N:CA:C']
hueList = ['aa','bfactor','rid','resolution','pdbCode','2FoFc','dssp']

data = georep.getGeoemtryCsv(geoList, hueList)

georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='2FoFc', title='Ramachandran', palette='cubehelix_r',sort='ASC')
georep.addScatter(data=data, geoX='PHI',geoY='PSI',hue='dssp', title='Ramachandran', palette='tab10',sort='NON')
georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='2FoFc', title='Omega-Phi', palette='cubehelix_r',sort='ASC')
georep.addScatter(data=data, geoX='OMEGA',geoY='PHI',hue='dssp', title='Omega-Phi', palette='tab10',sort='NON')
georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='2FoFc', title='C-1:C/Phi', palette='cubehelix_r',sort='ASC')
georep.addScatter(data=data, geoX='C-1:C',geoY='PHI',hue='dssp', title='C-1:C/Phi', palette='tab10',sort='NON')
georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='2FoFc', title='Tau angles', palette='cubehelix_r',sort='ASC')
georep.addScatter(data=data, geoX='C-1:N:CA',geoY='TAU',hue='dssp', title='Tau angles', palette='tab10',sort='NON')

# And finally create the reort with a file name of choice
title = 'Some correlation reports with 2fo-Fc as the hue'
georep.printToHtml(title,2,'densityhue')


