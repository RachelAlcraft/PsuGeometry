
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/'

# split list in 2 for memory purposes
#pdbList = ['1ejg','1us0','1tt8','1i1w','1ucs','1yk4','1yk4','1hje','1r6j']
#pdbList = ['2bw4','3nir','3x2m','2VB1','3A39','2b97','2OV0','2WFI']
#pdbList = ['4ZM7','4REK','5D8V','5NW3','5qkw']
pdbList = ['6jvv','6rr2','6E6O','6S2M','6shk','6fgz','6ctd','6fwf','6q53']
# split list in 2 for memory purposes

pdbList = ['3o4p','1pjx']

runIndividualReports = True


if not runIndividualReports:
    geoList = []
    for pdb in pdbList:
        pdb = pdb.lower()
        geoPdb = geop.GeoPdb(pdb, pdbDataPath,edDataPath)
        geoList.append(geoPdb)

    georep = geor.GeoReport(geoList)

    geoName = 'all'
    georep.printReport('Ramachandran',  printPath,geoName + '_rama' )
    georep.printReport('Sp2Planarity',printPath,geoName+ '_sp2')
    georep.printReport('Sp3Tetrahedra',printPath,geoName + '_sp3')
    georep.printReport('BackboneOutliers', printPath,geoName+ '_bbone')
    georep.printReport('OmegaCis', printPath,geoName + '_ocis')
    georep.printReport('RachelsChoice', printPath,geoName + '_rae')
    georep.printReport('DataPerPdb', printPath,geoName + '_data')

else:
    for pdb in pdbList:
        pdb = pdb.lower()
        geoPdb = geop.GeoPdb(pdb, pdbDataPath, edDataPath)
        georep = geor.GeoReport([geoPdb])
        georep.printReport('Ramachandran',  printPath,geoPdb.pdbCode + '_rama' )
        georep.printReport('Sp2Planarity',printPath,geoPdb.pdbCode + '_sp2')
        georep.printReport('Sp3Tetrahedra',printPath,geoPdb.pdbCode + '_sp3')
        georep.printReport('BackboneOutliers', printPath,geoPdb.pdbCode + '_bbone')
        georep.printReport('OmegaCis', printPath,geoPdb.pdbCode + '_ocis')
        georep.printReport('RachelsChoice', printPath,geoPdb.pdbCode + '_rae')
        georep.printReport('DataPerPdb', printPath,geoPdb.pdbCode + '_data')
        geoPdb = None
        geoRep = None



