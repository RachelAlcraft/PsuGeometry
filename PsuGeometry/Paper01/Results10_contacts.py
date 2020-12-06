# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
'''
This script looks electron density correlations 
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList = ['5nqo','1ejg','1tt8']
georep = psu.GeoReport(pdbList,pdbDataPath, edDataPath,printPath)

for pdb in pdbList:
    georep.addCloseContact(pdb, 'N', 'O', 8, 2,hue='dssp',categorical=True,palette='tab20_r',title=pdb)
    georep.addCloseContact(pdb, 'CA', 'CA', 8, 2,hue='bfactor',categorical=False,palette='cubehelix_r',title=pdb)
    georep.addCloseContact(pdb, 'CB', 'CB', 10, 2,hue='2FoFc',categorical=False,palette='viridis_r',title=pdb)
    georep.addCloseContact(pdb, 'SG', 'SG', 14, 2, hue='distance', categorical=False, palette='viridis_r', title=pdb)

georep.printToHtml('Close Contacts',4,'Results10_contacts')

