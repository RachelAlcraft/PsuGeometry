# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
'''
This script looks electron density correlations 
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList = ['2bw4','5nqo','1ejg','6q53']

for pdb in pdbList:
    georep = psu.GeoReport([pdb], pdbDataPath, edDataPath, printPath)
    georep.printReport('DataPerPdb', 'Results9_data')
    georep.printReport('Slow_DensityPeaksPerPdb', 'Results9_density')









