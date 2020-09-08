from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/density/'

### split list in 2 for memory purposes
#pdbList = ['1ejg','1us0','1tt8','1i1w','1ucs','1yk4','1yk4','1hje','1r6j']
#pdbList = ['2bw4','3nir','3x2m','2VB1','3A39','2b97','2OV0','2WFI']
#pdbList = ['4ZM7','4REK','4ZM7','5D8V','5NW3','5qkw']
#pdbList = ['6jvv','6rr2','6E6O','6S2M','6shk','6fgz','6ctd','6fwf','6q53']
pdbList = ['1ejg','1us0','1tt8','1i1w','1ucs','1yk4','1yk4','1hje','1r6j','2bw4','3nir','3x2m','2VB1','3A39','2b97','2OV0','2WFI']
#pdbList = ['4ZM7','4REK','4ZM7','5D8V','5NW3','5qkw','6jvv','6rr2','6E6O','6S2M','6shk','6fgz','6ctd','6fwf','6q53']
### split list in 2 for memory purposes


for pdb in pdbList:
    pdb = pdb.lower()
    geoPdb = geop.GeoPdb(pdb, pdbDataPath, edDataPath)
    georep = geor.GeoReport([geoPdb])
    if geoPdb.atoms[0].values['resolution'] < 1.7:
        georep.printReport('Slow_DensityPeaksPerPdb', printPath, geoPdb.pdbCode + '_den')
    else:
        georep.printReport('Slow_DensityPointsPerPdb', printPath, geoPdb.pdbCode + '_den')
    geoPdb = None
    georep = None