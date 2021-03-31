# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import time
import pandas as pd
'''
TAU 
I have 3 data sets
Original - the entire glycine dataset form tghe liost of structures I have chosen that are 1.3* avergae with no multiple occupancy
Good - of all the residues above, those with a nearby density peak (good local density)
Best - of all those above, where the tau values are calculated as identical

This can currently only be done from windows becuase the sets oif pdb files are  there.
This file only creates the best supported set looking at the rtau compare.
Other geos are appended in from the subsequent script.
'''

#These are the paths
pdbOriginalPath = 'F:/Code/ProteinDataFiles/pdb_data/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'


#TIMER
print('----------start report ----------')
startx = time.time()

#This gets the list of pdbs
pdbReDataPath = 'F:/Code/ProteinDataFiles/pdb_out/NCACO_001_05/'
pdbdata = pd.read_csv('../../PdbLists/Pdbs_Inc1.csv') # This is a list of pdbs <= 1.1A non homologous to 90%
pdbListIn = pdbdata['PDB'].tolist()[0:]
pdbList = []
for pdb in pdbListIn:
    import os.path
    if os.path.isfile((pdbReDataPath + 'pdb' + pdb + '.ent').lower()):
        pdbList.append(pdb.lower())
    else:
        print('No file:',(pdbReDataPath + 'pdb' + pdb + '.ent').lower())

print('#######################')
print('')
for pdb in pdbList:
    print(pdb)