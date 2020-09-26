from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/covid/'

pdbList = ['6lu7']

for pdb in pdbList:
    pdb = pdb.lower()
    georep = geor.GeoReport([pdb],pdbDataPath, edDataPath,printPath)
    georep.printReport('Sp2Planarity',pdb + '_sp2')
    georep.printReport('BackboneOutliers',pdb + '_bbone')
    georep.printReport('RachelsChoice',pdb + '_rae')
    georep.printReport('DataPerPdb', pdb + '_data')
    georep.printReport('Slow_DensityPeaksPerPdb', pdb + '_pkden')
