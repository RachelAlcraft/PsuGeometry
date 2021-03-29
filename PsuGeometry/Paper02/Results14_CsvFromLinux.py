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
'''

myWindowsLaptop = True

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
dsspHue='dssp'
includeDSSP = True
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/DataCsvOfSets/'
    includeDSSP = False  # on my windows computer


#From linux I am only interested in creating the unrestricted data - it can be appended to the other data on windows

#TIMER
print('----------start report 14----------')
startx = time.time()

#This gets the list of pdbs
pdbdata = pd.read_csv('../PdbLists/Pdbs_TauPaper.csv') # This is a list of pdbs <= 1.1A non homologous to 90%
pdbList = pdbdata['PDB'].tolist()[0:2]
pdbList = ['4a7u']

#This is all the data we are going to be looking at
geoList = [
            #'N:N+1','TAU','PSI','PHI','N:C','CA:C','C:O','N:CA','C-1:N','C:N+1','OMEGA',
            #'CA:C:O:N+1','O:N+1','CA:O','CA:N+1','CA:C:N+1','C-1:N:CA','N:O-2','N:CA:C:O-2',
            #'N:O-2:CA','N-1:CA:C','CA:HOH','CA:HETATM','N:HETATM:C','N:HOH:C','N:CA:C:HETATM','N:CA:C:HOH',
            #'O-2:C','O-2:N:CA','O-2:N:CA:N+1',
            #'N:{O,OD1,OG1}','{O,OD1,OG1}:C','{O,OD1,OG1}:N:CA','{O,OD1,OG1}:N:CA:N+1',
            'N:{O}','C:{O}','CA:N:{O}','N+1:CA:N:{O}',
            'N:O-2','C:O-2','CA:N:O-2','N+1:CA:N:O-2',
            #'N:{OD1,OG1}','C:{OD1,OG1}','CA:N:{OD1,OG1}','N+1:CA:N:{OD1,OG1}',
            ]
hueList = ['aa', 'rid', 'bfactor','pdbCode','bfactorRatio','disordered','dssp']


print('Creating CSV files anew')
georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False,keepDisordered=True)
print('Create unrestricted csv')
dataUnrestricted = georep.getGeoemtryCsv(geoList, hueList, -1,allAtoms=True)
#print(dataUnrestricted)
if myWindowsLaptop:
    dataUnrestricted.to_csv(printPath + "Results14_UnrestrictedCsvFromWindowsx.csv", index=False)
else:
    dataUnrestricted.to_csv(printPath + "Results14_UnrestrictedCsvFromLinux.csv", index=False)

print('----------Finished----------')
endx = time.time()
time_diff = endx - startx
timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
print(timestring)

