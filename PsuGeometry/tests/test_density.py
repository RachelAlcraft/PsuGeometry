from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/density/'

### split list in 2 for memory purposes
#pdbList = ['1ejg','1us0','1tt8','1i1w','1ucs','1yk4','1yk4','1hje','1r6j']
#pdbList = ['2bw4','3nir','3x2m','2VB1','3A39','2b97','2OV0','2WFI']
#pdbList = ['4ZM7','4REK','4ZM7','5D8V','5NW3','5qkw']
pdbList = ['6jvv','6rr2','6E6O','6S2M','6shk','6fgz','6ctd','6fwf','6q53']
#pdbList = ['1ejg','1us0','1tt8','1i1w','1ucs','1yk4','1yk4','1hje','1r6j','2bw4','3nir','3x2m','2VB1','3A39','2b97','2OV0','2WFI','3o4p','1pjx']
#pdbList = ['4ZM7','4REK','4ZM7','5D8V','5NW3','5qkw','6jvv','6rr2','6E6O','6S2M','6shk','6fgz','6ctd','6fwf','6q53']
### split list in 2 for memory purposes

pdbList=['6fwf','6q53','6ctd','6fgz','6shk','6rr2','6jvv']
peaks = True

for pdb in pdbList  :
    georep = geor.GeoReport([pdb],pdbDataPath, edDataPath,printPath)
    georep.addDensityView(pdb, geoX='c', geoY='r', peaks=peaks, divisor=10, palette='cubehelix_r')
    georep.addDensityView(pdb, geoX='r', geoY='s', peaks=peaks, divisor=10, palette='cubehelix_r')
    georep.addDensityView(pdb, geoX='s', geoY='c', peaks=peaks, divisor=10, palette='cubehelix_r')
    georep.addDensityView(pdb,geoX='x',geoY='y',peaks=peaks,divisor=10,palette='cubehelix_r')
    georep.addDensityView(pdb, geoX='y', geoY='z', peaks=peaks, divisor=10, palette='cubehelix_r')
    georep.addDensityView(pdb, geoX='z', geoY='x', peaks=peaks, divisor=10, palette='cubehelix_r')
    georep.addDataView(pdb,geoX='x', geoY='y', palette='cubehelix_r')
    georep.addDataView(pdb, geoX='y', geoY='z', palette='cubehelix_r')
    georep.addDataView(pdb, geoX='z', geoY='x', palette='cubehelix_r')
    georep.addScatter(geoX='x',geoY='y', palette='cubehelix_r',hue='2FoFc')
    georep.addScatter(geoX='y', geoY='z', palette='cubehelix_r',hue='2FoFc')
    georep.addScatter(geoX='z', geoY='x', palette='cubehelix_r',hue='2FoFc')
    # And finally create the reort with a file name of choice
    georep.printToHtml('Manually Chosen Density Views', 3, pdb + '_den2')



