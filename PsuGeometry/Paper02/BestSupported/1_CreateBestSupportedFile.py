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
pdbGoodPath = 'F:/Code/ProteinDataFiles/pdb_out/good/'
pdbBetterPath = 'F:/Code/ProteinDataFiles/pdb_out/better/'
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'


#TIMER
print('----------start report 14----------')
startx = time.time()

#This gets the list of pdbs
pdbReDataPath = 'F:/Code/ProteinDataFiles/pdb_out/good/'
pdbdata = pd.read_csv('../../PdbLists/Pdbs_Under1.csv') # This is a list of pdbs <= 1.1A non homologous to 90%
pdbListIn = pdbdata['PDB'].tolist()[0:]
pdbList = []
for pdb in pdbListIn:
    import os.path
    if os.path.isfile((pdbReDataPath + 'pdb' + pdb + '.ent').lower()):
        pdbList.append(pdb.lower())
    else:
        print('No file:',(pdbReDataPath + 'pdb' + pdb + '.ent').lower())

print('*****************\n')
for pdd in pdbList:
    print(pdd)
print('*****************\n')
#This is all the data we are going to be looking at
geoList = ['TAU']
hueList = ['rid','aa', 'rid', 'bfactor','pdbCode','bfactorRatio','disordered']

print('Creating CSV files anew')

#Load original csv unrestricted on occupant and bfactor
pdbmanager = geopdb.GeoPdbs(pdbOriginalPath, edDataPath, False, False, True)

georep = psu.GeoReport(pdbList, pdbOriginalPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=True)
print('Create unrestricted csv')
dataUnrestricted = georep.getGeoemtryCsv(geoList, hueList, -1,restrictedAa='ALL')
print(dataUnrestricted)
#Add id column
dataUnrestricted['rid'] = dataUnrestricted['rid'].astype(str)
dataUnrestricted['ID'] = dataUnrestricted['pdbCode'] + dataUnrestricted['chain'] + dataUnrestricted['rid'] + dataUnrestricted['aa']
dataUnrestricted.to_csv(printPath + "Set1_UnrestrictedStats.csv", index=False)

#Load good csv
pdbmanager.clear()
pdbmanager = geopdb.GeoPdbs(pdbGoodPath, edDataPath, False,False,False)
georep = psu.GeoReport(pdbList, pdbGoodPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
print('Create good csv')
dataGood = georep.getGeoemtryCsv(geoList, hueList,1.3,restrictedAa='ALL')
#Add id column
dataGood['rid'] = dataGood['rid'].astype(str)
dataGood['ID'] = dataGood['pdbCode'] + dataGood['chain'] + dataGood['rid'] + dataGood['aa']
dataGood.to_csv(printPath + "Csv_GoodStats.csv", index=False)

#Load better csv
pdbmanager.clear()
pdbmanager = geopdb.GeoPdbs(pdbBetterPath, edDataPath, False,False,False)
georep = psu.GeoReport(pdbList, pdbBetterPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
print('Create better csv')
dataBetter = georep.getGeoemtryCsv(geoList, hueList,1.3,restrictedAa='ALL')
print(dataBetter)
dataBetter['rid'] = dataBetter['rid'].astype(str)
dataBetter['ID'] = dataBetter['pdbCode'] + dataBetter['chain'] + dataBetter['rid'] + dataBetter['aa']
print(dataBetter)
dataBetter.to_csv(printPath + "Csv_BetterStats.csv", index=False)

#Load them all up anyway even if we already saved them, keep things consistent
#dataUnrestricted = pd.read_csv(printPath + "Results14_UnrestrictedStats.csv")
#dataGood = pd.read_csv(printPath + "Results14_GoodStats.csv")
#dataBetter = pd.read_csv(printPath + "Results14_BetterStats.csv")

#To find best supported we need to find a tau with no difference
headerList = ['ID','pdbCode','chain','rid','aa','bfactor','bfactorRatio','disordered','TAU']
dataUnrestrictedCut = dataUnrestricted[headerList]
dataGoodCut = dataGood[headerList]
dataBetterCut = dataBetter[headerList]

#Create the bfactor and occupant restricted sets
dataOccupantCut = dataUnrestrictedCut.query("disordered ==  'N'")
data_1_3_Cut = dataOccupantCut.query("bfactorRatio <  1.3")
data_1_6_Cut = dataOccupantCut.query("bfactorRatio <  1.6")
data_interim_Cut = dataOccupantCut.query("bfactorRatio <  1.6 and bfactorRatio >=  1.3")

#create a data set where the taus are identical
dataBestSupported = dataGoodCut
dataBestSupported['TAU2'] = dataBetterCut['TAU']
dataBestSupported['TAU_DIFF'] = abs(dataBestSupported['TAU2'] - dataBestSupported['TAU'])
dataNotSupported = dataBestSupported.query('TAU_DIFF > 0')
dataBestSupported = dataBestSupported.query('TAU_DIFF == 0')

dataBestSupported.to_csv(printPath + "Set3_BestSupported.csv", index=False)
data_1_3_Cut.to_csv(printPath + "Set2_Restricted.csv", index=False)

print('----------Finished----------')
endx = time.time()
time_diff = endx - startx
timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
print(timestring)

