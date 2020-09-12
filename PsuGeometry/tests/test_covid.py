from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/covid/'

pdbList = ['6lu7']

for pdb in pdbList:
    pdb = pdb.lower()
    geoPdb = geop.GeoPdb(pdb, pdbDataPath, edDataPath)
    georep = geor.GeoReport([geoPdb])
    georep.printReport('Sp2Planarity',printPath,geoPdb.pdbCode + '_sp2')
    georep.printReport('BackboneOutliers', printPath,geoPdb.pdbCode + '_bbone')
    georep.printReport('RachelsChoice', printPath, geoPdb.pdbCode + '_rae')
    georep.printReport('DataPerPdb', printPath,geoPdb.pdbCode + '_data')
    #if geoPdb.atoms[0].values['resolution'] < 1.7:
    georep.printReport('Slow_DensityPeaksPerPdb', printPath, geoPdb.pdbCode + '_pkden')
    #else:
    #georep.printReport('Slow_DensityPointsPerPdb', printPath, geoPdb.pdbCode + '_poden')
    geoPdb = None
    geoRep = None