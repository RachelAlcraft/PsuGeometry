from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'


pdb= '1ejg'
if True:
    georep = geor.GeoReport([pdb],pdbDataPath, edDataPath,printPath)

    georep.addDataView(pdb, geoX='x', geoY='y', palette='cubehelix_r')
    georep.addScatter(geoX='x', geoY='y', palette='rainbow',hue='rid')


    central = [3.31, 3.77, 0.97]
    linear = [3.27, 4.91, 1.7]
    planar = [2.97, 3.73, -0.21]
    sfc = georep.addDensitySlice(pdb,70,0.05,central,linear,planar,palette='cubehelix_r')
    # And finally create the reort with a file name of choice
    georep.printToHtml('Density Views and Slices', 3, pdb + '_densl')




