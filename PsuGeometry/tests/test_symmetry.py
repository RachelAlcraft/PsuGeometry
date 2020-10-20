from PsuGeometry import GeoSymmetry as geosym
from PsuGeometry import GeoPdb as geopdb

import time


###### User Choices ######################################################
pdbCode= '6q53'
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, True, False)
apdb = pdbmanager.getPdb(pdbCode)
if apdb.hasDensity:
    peaks = apdb.getStructureDensity(False,-1,pdbDataPath,edDataPath)
    geosym.GeoSymmetry(peaks)
