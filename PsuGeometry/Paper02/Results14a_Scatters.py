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
print('----------start report 14----------')
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
dataBetter = pd.read_csv(printPath + "Results14_BetterStats.csv")

#To find best supported we need to find a tau with no difference
#Cut our data into the columns we wanr
georep = psu.GeoReport([], pdbOriginalPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

#aas = ['GLY','PRO','ALA']
#aas = ['GLY']
aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
for aa in aas:
    headerList = ['pdbCode','chain','rid','aa','C-1:N_motif','C:N+1_motif','bfactor','bfactorRatio','disordered','N:N+1','CA:C','C:O','N:CA','C-1:N','C:N+1','TAU','PSI','PHI','OMEGA','CA:C:O:N+1']
    dataUnrestrictedCut = dataUnrestricted[headerList]
    dataGoodCut = dataGood[headerList]
    dataBetterCut = dataBetter[headerList]

    #Cut on amino acid
    dataUnrestrictedCut = dataUnrestrictedCut.query("aa ==  '" + aa + "'")
    dataGoodCut = dataGoodCut.query("aa ==  '" + aa + "'")
    dataBetterCut = dataBetterCut.query("aa ==  '" + aa + "'")
    #Create the bfactor and occupant restricted sets
    dataOccupantCut = dataUnrestrictedCut.query("disordered ==  'N'")
    data_1_3_Cut = dataOccupantCut.query("bfactorRatio <  1.3")
    #create a data set where the taus are identical
    dataBestSupported = dataGoodCut
    dataBestSupported['TAU2'] = dataBetterCut['TAU']
    dataBestSupported['TAU_DIFF'] = abs(dataBestSupported['TAU2'] - dataBestSupported['TAU'])
    dataNotSupported = dataBestSupported.query('TAU_DIFF > 0')
    dataBestSupported = dataBestSupported.query('TAU_DIFF == 0')

    #Ans val the OMEGAS
    dataGoodCut['ABS_OMEGA'] = abs(dataGoodCut['OMEGA'])
    dataGoodCut['ABS_CA:C:O:N+1'] = abs(dataGoodCut['CA:C:O:N+1'])
    dataGoodCut['ABS_PSI'] = abs(dataGoodCut['PSI'])
    dataOmegaCut = dataGoodCut.query("ABS_OMEGA >  150")

    #Against ALL DATA
    georep.addScatter(data=dataGoodCut, geoX='PSI', geoY='N:N+1', hue='TAU', title=aa + ' Good', palette='jet', sort='NON')
    georep.addScatter(data=dataOmegaCut, geoX='PSI', geoY='ABS_OMEGA', hue='TAU', title=aa + ' Good', palette='jet', sort='NON')
    georep.addScatter(data=dataGoodCut, geoX='PSI', geoY='ABS_CA:C:O:N+1', hue='TAU', title=aa + ' Good', palette='jet', sort='NON')
    georep.addScatter(data=dataOmegaCut, geoX='ABS_CA:C:O:N+1', geoY='ABS_OMEGA', hue='ABS_PSI', title=aa + ' Good', palette='jet', sort='NON')
    georep.addScatter(data=dataBestSupported, geoX='PSI', geoY='N:N+1', hue='TAU', title=aa + ' Best Supported', palette='jet', sort='NON')
    georep.addScatter(data=dataBestSupported, geoX='N:CA', geoY='CA:C', hue='TAU', title=aa + ' Best Supported', palette='jet', sort='NON')
    georep.addScatter(data=dataBestSupported, geoX='N:CA', geoY='C:O', hue='TAU', title=aa + ' Best Supported', palette='jet', sort='NON')
    georep.addScatter(data=dataBestSupported, geoX='CA:C', geoY='C:O', hue='TAU', title=aa + ' Best Supported', palette='jet', sort='NON')


    #georep.addScatter(data=dataBestSupported, geoX='N:N+1', geoY='TAU', hue='PSI', title=aa + ' Best Supported',palette='jet', sort='NON')
    #georep.addScatter(data=dataBestSupported, geoX='TAU', geoY='PSI', hue='N:N+1', title=aa + ' Best Supported', palette='jet', sort='NON')


georep.printToHtml('Results 14a. Tau Planarity Scatters', 4, 'Results14a_PlanarityScattersAll')


print('----------Finished----------')
endx = time.time()
time_diff = endx - startx
timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
print(timestring)

