# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import _Helpers as help
import time
import pandas as pd
'''
TAU 
I have 3 data sets
Original - the entire glycine dataset form tghe liost of structures I have chosen that are 1.3* avergae with no multiple occupancy
Good - of all the residues above, those with a nearby density peak (good local density)
Best - of all those above, where the tau values are calculated as identical
'''

def createGeosFile(pdbSet, cutOff):
    print('----------start create geos csv----------')
    startx = time.time()

    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataA/'
    print('Running CreateGeosFile on',pdbSet)
    allAtoms = False
    bFactorFactor = 1.3
    if pdbSet == 'UNRESTRICTED':
        allAtoms = True
        bFactorFactor = -1

    geos = ['N:CA', 'CA:C', 'C:O', 'C-1:N', 'C:N+1']
    print('Create csv 1', pdbSet)
    csv1 = help.getCsv(pdbSet, geos,badAtoms,True, True,aa='ALL',includeCis=False,allAtoms=allAtoms, bFactorFactor=bFactorFactor,cutoff=cutOff)
    csv1.to_csv(printPath + 'CsvGeos_BEST_' + 'Set1BONDALL_' + pdbSet + '.csv', index=False)

    geos = ['TAU', 'C-1:N:CA', 'CA:C:N+1', 'CA:C:O', 'O:C:N+1', 'CA:C:N+1']
    print('Create csv 2', pdbSet)
    csv2 = help.getCsv(pdbSet, geos, badAtoms,False, True, aa='ALL', includeCis=False, allAtoms=allAtoms, bFactorFactor=bFactorFactor, cutoff=cutOff)
    csv2.to_csv(printPath + 'CsvGeos_BEST_' + 'Set2ANGSALL_' + pdbSet + '.csv', index=False)

    geos = ['PHI', 'PSI', 'OMEGA', 'CA-1:C-1:N:CA']
    print('Create csv 3', pdbSet)
    csv3 = help.getCsv(pdbSet, geos, badAtoms,False, True, aa='ALL', includeCis=False, allAtoms=allAtoms, bFactorFactor=bFactorFactor, cutoff=cutOff)
    csv3.to_csv(printPath + 'CsvGeos_BEST_' + 'Set3DIHSALL_' + pdbSet + '.csv', index=False)

    geos = ['N:N+1', 'N:C']
    print('Create csv 4', pdbSet)
    csv4 = help.getCsv(pdbSet, geos, badAtoms,False, True, aa='ALL', includeCis=False, allAtoms=allAtoms,bFactorFactor=bFactorFactor, cutoff=cutOff)
    csv4.to_csv(printPath + 'CsvGeos_BEST_' + 'Set4DISTALL_' + pdbSet + '.csv', index=False)

    geos = ['N:O-2', 'C:O-2', 'N:CA:C:O-2', 'N:CA:N+1:O-2']
    print('Create csv 5', pdbSet)
    csv5 = help.getCsv(pdbSet, geos, badAtoms,False, True, aa='ALL', includeCis=False, allAtoms=allAtoms,bFactorFactor=bFactorFactor, cutoff=cutOff)
    csv5.to_csv(printPath + 'CsvGeos_BEST_' + 'Set5HBALL_' + pdbSet + '.csv', index=False)

    geos = ['N:O-3', 'C:O-3', 'N:CA:C:O-3', 'N:CA:N+1:O-3']
    print('Create csv 6', pdbSet)
    csv6 = help.getCsv(pdbSet, geos, badAtoms,False, True, aa='ALL', includeCis=False, allAtoms=allAtoms,bFactorFactor=bFactorFactor, cutoff=cutOff)
    csv6.to_csv(printPath + 'CsvGeos_BEST_' + 'Set6HBALL_' + pdbSet + '.csv', index=False)

    print('----------Finished----------')
    endx = time.time()
    time_diff = endx - startx
    timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
    print(timestring)



def createGeosFileOld(pdbSet, cutOff):
    print('Running CreateGeosFile on',pdbSet)

    setName = 'BEST'
    #pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/NCACO_001_05/'
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    keepDisordered = False
    bFactorFactor = -1 #no need to restrict on bfactor and electron density but we do not want occupant verions

    if pdbSet == 'RESTRICTED':
        pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
        keepDisordered = False
        bFactorFactor = 1.3
    elif pdbSet == 'UNRESTRICTED':
        pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
        keepDisordered = True
        bFactorFactor = -1

    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataA/'


    #TIMER
    print('----------start report 14----------')
    startx = time.time()

    #This gets the list of pdbs
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_Evidenced.csv') # This is a list of pdbs <= 1.1A non homologous to 90%
    pdbListIn = pdbdata['PDB'].tolist()[0:]
    if cutOff > 0:
        pdbListIn = pdbdata['PDB'].tolist()[0:cutOff]

    pdbList = []
    for pdb in pdbListIn:
        import os.path
        if os.path.isfile((pdbDataPath + 'pdb' + pdb + '.ent').lower()):
            pdbList.append(pdb.lower())
        else:
            print('No file:',(pdbDataPath + 'pdb' + pdb + '.ent').lower())


    #Clear the pdb manager cache
    pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, False, False, keepDisordered)
    pdbmanager.clear()



    #This is all the data we are going to be looking at
    geoLists = []
    #geoLists.append(['TAU'])#for dssp from linux only
    #geoLists.append(['0', ['TAU']])
    #We currently have these
    geoLists.append(['1BOND', ['N:CA','CA:C','C:O','C-1:N','C:N+1']])    
    geoLists.append(['2ANGS', ['TAU','C-1:N:CA','CA:C:N+1','CA:C:O','O:C:N+1','CA:C:N+1']])    
    geoLists.append(['3DIHS', ['PHI','PSI','OMEGA','CA-1:C-1:N:CA']])    
    geoLists.append(['4DIST', ['N:N+1','N:C']])    
    geoLists.append(['5HB', ['N:O-2','C:O-2','N:CA:C:O-2','N:CA:N+1:O-2']])
    geoLists.append(['6HB', ['N:O-3', 'C:O-3', 'N:CA:C:O-3', 'N:CA:N+1:O-3']])
    #geoLists.append(['9HBO', ['N:{O}','C:{O}','N:CA:C:{O}','N:CA:N+1:{O}']])

    #geoLists.append(['10CIS', ['CA-1:C-1:N:CA', 'CA-1:CA']])
    # geoLists.append(['7HB', ['N:O-4', 'C:O-4', 'N:CA:C:O-4', 'N:CA:N+1:O-4']])
    # geoLists.append(['8HB', ['N:O-5', 'C:O-5', 'N:CA:C:O-5', 'N:CA:N+1:O-5']])
    # Water
    #geoLists.append(['7WAT', ['N:HOH','C:HOH','N:CA:C:HOH','N:CA:N+1:HOH']])
    # Other!
    #geoLists.append(['8XTRA', ['N:HETATM']])

    hueList = ['aa', 'rid', 'bfactor','pdbCode','bfactorRatio','disordered']
    aas = ['ALL']

    print('Creating CSV files anew')
    for geoListT in geoLists:
        geoList = geoListT[1]
        set = geoListT[0]
        for aa in aas:
            tag = 'Set' + set + aa + '_' + pdbSet
            georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=keepDisordered)
            print('Create csv', pdbDataPath,geoList)
            dataUnrestricted = georep.getGeoemtryCsv(geoList, hueList, bFactorFactor,allAtoms=True,restrictedAa=aa)
            dataUnrestricted.to_csv(printPath + 'CsvGeos_' + setName + '_' + tag + '.csv', index=False)


    print('----------Finished----------')
    endx = time.time()
    time_diff = endx - startx
    timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
    print(timestring)

