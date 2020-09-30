from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/contacts/'

pdbList = ['5nqo']#,'1ejg','4rek']

for pdb in pdbList:
    pdb = pdb.lower()
    georep = geor.GeoReport([pdb],pdbDataPath, edDataPath,printPath)
    georep.addCloseContact(pdb, 'N', 'O', 8,2,title='')
    georep.addCloseContact(pdb, 'CA', 'CA', 8,2)
    georep.addCloseContact(pdb, 'CB', 'CB', 6,2)
    georep.addCloseContact(pdb, 'N', 'O', 10, 3,hue='dssp',categorical=True,palette='rainbow')
    georep.addCloseContact(pdb, 'CA', 'CA', 10, 3,hue='bfactor',categorical=False,palette='rainbow_r')
    georep.addCloseContact(pdb, 'CB', 'CB', 10, 3,hue='2FoFc',categorical=False,palette='rainbow')
    georep.addCloseContact(pdb, 'SG', 'SG', 8, 2)
    georep.addCloseContact(pdb, 'CG', 'CB', 8, 2)
    georep.addCloseContact(pdb, 'CB', 'O', 8, 2)
    georep.printToHtml('Close Contacts',3,pdb+'_close')

