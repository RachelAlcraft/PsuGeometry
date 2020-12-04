# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
'''
A single predefined report
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Levels/'

pdbList = ['1us0','1ejg','2cnq']
georep = psu.GeoReport(pdbList,pdbDataPath,edDataPath,printPath)
georep.printReport('RachelsChoice', 'rachel')

