from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/contacts/'

pdbList = ['5nqo','1ejg','4rek']

for pdb in pdbList:
    pdb = pdb.lower()
    georep = geor.GeoReport([pdb],pdbDataPath, edDataPath,printPath)
    georep.addCloseContact(pdb, 'N', 'O', 6,2)
    georep.addCloseContact(pdb, 'CA', 'CA', 6,2)
    georep.addCloseContact(pdb, 'CB', 'CB', 6,2)
    georep.addCloseContact(pdb, 'N', 'O', 6, 2,hue='dssp',categorical=True,palette='rainbow')
    georep.addCloseContact(pdb, 'CA', 'CA', 6, 2,hue='aa',categorical=True,palette='rainbow')
    georep.addCloseContact(pdb, 'CB', 'CB', 6, 2,hue='2FoFc',categorical=False,palette='rainbow')
    georep.addCloseContact(pdb, 'SG', 'SG', 6, 2)
    georep.addCloseContact(pdb, 'CG', 'CB', 6, 2)
    georep.addCloseContact(pdb, 'CB', 'O', 6, 2)
    georep.printToHtml('Close Contacts',3,pdb+'_close')

