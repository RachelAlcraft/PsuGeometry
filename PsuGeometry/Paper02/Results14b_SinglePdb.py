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

#AM I RECREATING THE DATA SET????

#These are the paths
pdbOriginalPath = 'F:/Code/ProteinDataFiles/pdb_data/'
pdbGoodPath = 'F:/Code/ProteinDataFiles/pdb_out/good/'
pdbBetterPath = 'F:/Code/ProteinDataFiles/pdb_out/better/'
printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'


#TIMER
print('----------start report 14 b----------')
startx = time.time()

#This gets the list of pdbs
pdbReDataPath = 'F:/Code/ProteinDataFiles/pdb_out/good/'
pdbdata = pd.read_csv('../PdbLists/Pdbs_Under1.csv') # This is a list of pdbs <= 1.1A non homologous to 90%
pdbListIn = pdbdata['PDB'].tolist()[0:]
pdbList = []
for pdb in pdbListIn:
    import os.path
    if os.path.isfile((pdbReDataPath + 'pdb' + pdb + '.ent').lower()):
        pdbList.append(pdb.lower())
    else:
        print('No file:',(pdbReDataPath + 'pdb' + pdb + '.ent').lower())

print(pdbList)

#Load them all up
dataUnrestricted = pd.read_csv(printPath + "Results14_UnrestrictedStats.csv")
dataGood = pd.read_csv(printPath + "Results14_GoodStats.csv")
#Cut our data into the columns we wanr
headerList = ['pdbCode','chain','rid','aa','C-1:N_motif','C:N+1_motif','bfactor','bfactorRatio','disordered','N:N+1','CA:C','C:O','N:CA','C-1:N','C:N+1','TAU','PSI','PHI','OMEGA','CA:C:O:N+1']
pdbs = ['5a8c']

#aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
for pdb in pdbs:
    print(pdb)
    georep = psu.GeoReport([pdb], pdbOriginalPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=True, keepDisordered=False)
    aas = ['ALL','GLY','PRO','ALA','HIS']
    for aa in aas:
        dataUnrestrictedCut = dataUnrestricted[headerList]
        dataGood= dataGood[headerList]
        if aa != 'ALL':
            dataGood = dataGood.query("aa ==  '" + aa + "'")

        dataOccupantCut = dataGood.query("disordered ==  'N'")
        data_1_3_Cut = dataOccupantCut.query("bfactorRatio <  1.3")

        # Cut on pdb
        dataGoodPdb = data_1_3_Cut.query("pdbCode ==  '" + pdb + "'")

        georep.addHistogram(data=data_1_3_Cut, geoX='TAU', title=aa + ' Good set, restricted on bfactor and occupant')
        georep.addStatsCompare(dataA=data_1_3_Cut, dataB=dataGoodPdb, descA=aa + ' Good set, restricted on bfactor and occupant', descB=pdb + ' Restricted Data ' + aa, geoX='TAU', title=pdb + ' Tau compare')
        georep.addHistogram(data=dataGoodPdb, geoX='TAU', title=pdb + ' Restricted Data ' + aa)

    georep.printToHtml('Results 14b. Single pdb=' + pdb, 3, 'Results14b_Pdb_' + pdb)


print('----------Finished----------')
endx = time.time()
time_diff = endx - startx
timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
print(timestring)

