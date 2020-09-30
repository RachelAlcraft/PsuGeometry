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

#pdbList = ['6shk']
#pdbList = ['6jvv']

peaksList=['6jvv']
pointsList=['6fwf','6q53','6ctd','6fgz','6shk','6rr2','6jvv']
pointsList=['6jvv']

for pdb in peaksList  :
    georep = geor.GeoReport([pdb],pdbDataPath, edDataPath,printPath)
    georep.printReport('Slow_DensityPeaksPerPdb', pdb + '_denpk')


for pdb in pointsList:
    georep = geor.GeoReport([pdb],pdbDataPath, edDataPath,printPath)
    georep.printReport('Slow_DensityPointsPerPdb', pdb + '_denpt')
