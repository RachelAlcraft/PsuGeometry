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
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
dsspHue='dssp'
includeDSSP = True
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
    includeDSSP = False  # on my windows computer


#From linux I am only interested in creating the unrestricted data - it can be appended to the other data on windows

#TIMER
print('----------start report 14----------')
startx = time.time()

#This gets the list of pdbs
pdbdata = pd.read_csv('../../PdbLists/Pdbs_TauPaper.csv') # This is a list of pdbs <= 1.1A non homologous to 90%
pdbList = pdbdata['PDB'].tolist()[0:]

#This is all the data we are going to be looking at
geoLists = []
#geoLists.append(['TAU'])#for dssp from linux only
#geoLists.append(['0', ['TAU']])
# Bond lengths
geoLists.append(['1BOND', ['N:CA','CA:C','C:O','C-1:N','C:N+1']])
# Angles
geoLists.append(['2ANGS', ['TAU','C-1:N:CA','CA:C:N+1','CA:C:O','O:C:N+1','CA:C:N+1']])
# Dihedrals
geoLists.append(['3DIHS', ['PHI','PSI','OMEGA','CA-1:C-1:N:CA']])
# Distances
geoLists.append(['4DIST', ['N:N+1','N:C']])
# Hydrogen bond distances and dihedrals O-2
geoLists.append(['5HB', ['N:O-2','C:O-2','N:CA:C:O-2','N:CA:N+1:O-2']])
# Hydrogen bond distances and dihedrals nearest O
geoLists.append(['6HBO', ['N:{O}','C:{O}','N:CA:C:{O}','N:CA:N+1:{O}']])
# Water
geoLists.append(['7WAT', ['N:HOH','C:HOH','N:CA:C:HOH','N:CA:N+1:HOH']])
# Other!
geoLists.append(['8XTRA', ['N:HETATM']])


hueList = ['aa', 'rid', 'bfactor','pdbCode','bfactorRatio','disordered','dssp']
aas = ['ALL']

print('Creating CSV files anew')
for geoListT in geoLists:
    geoList = geoListT[1]
    set = geoListT[0]
    for aa in aas:
        tag = 'Set' + set + aa
        georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False,keepDisordered=True)
        print('Create unrestricted csv', geoList)
        dataUnrestricted = georep.getGeoemtryCsv(geoList, hueList, -1,allAtoms=True,restrictedAa=aa)
        dataUnrestricted.to_csv(printPath + 'CsvGeos_' + tag + '.csv', index=False)


print('----------Finished----------')
endx = time.time()
time_diff = endx - startx
timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
print(timestring)

