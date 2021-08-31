'''
This script creates the csv files that we use for the analysis from the pdb files.
There are 4 ultimate sets
1) Unrestricted from the pdb files, any bfactor or occupancy
2) Restricted to single occupancy and a bfactor factor of 1.3 (of average of the CA of the structre)
3) Restricted cut - cut to be the evidenced set only which is also single occupant and within 0.05 of a maximum
4) As above, but the atom centres are adjusted to the maxima
- Note this does limit the extremes by the 0,05 cutoff
'''

###########################################################
import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help
################################################################

runSetCreator = False
runDiffMaker = False
runChangeProToCis = False
runMakeAverages = False
runAppendDssp = True

### The input data, the names of the output sets and the names of the input pdb files ##############
### INPUTS ###
setName = 'UNRESTRICTED' #Options are UNRESTRICTED RESTRICTED REDUCED ADJUSTED
filesPDBRoot ='C:/Dev/Github/ProteinDataFiles/pdb_data/'
filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fo_ADJ/' #adjusted on Fo at 3 degrees thevenaz
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/' # only works on linux
geos = ['N:CA', 'CA:C', 'C:O', 'C-1:N', 'C:N+1','N:N+1','N:O','N:C','TAU', 'TAU-1', 'TAU+1', 'CA:C:O', 'O:C:N+1', 'CA:C:N+1','PHI', 'PSI', 'OMEGA','PREOMEGA']
hueList = ['aa', 'rid', 'bfactor', 'pdbCode', 'bfactorRatio', 'disordered']
### OUTPUTS ###
outUnrestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csv'
outRestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csv'
outReduced = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataReduced.csv'
outAdjusted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csv'
outMerged = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
####################################################################################################

################ HELPER FUNCTIONS ##########################


#######################################################################################

#The script only runs one at a time. Choose which you want to run. Note the final sets for cut and adjusted need to be made consistent
if runSetCreator:
    #We would only run the set creation on UNRESTRICTED or ADJUSTED
    # unrestricted automatically creates the RESTRICTED based on chosen paramaters and reduced based on ATOM list

    filesRoot = filesPDBRoot
    outPath = outUnrestricted
    keepDisordered = False
    badAtoms = []  # Create the badlist from csv files for mxima differences and occupancy

    # Get our chosen pdb list, this is 70% homology <=1A
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_70.csv')
    pdbListA = pdbdata['PDB'].tolist()[0:]
    pdbList = []
    for pdb in pdbListA:
        import os.path
        filePdb = filesADJRoot + 'pdb' + pdb + '.ent'
        if os.path.isfile((filePdb).lower()):
            pdbList.append(pdb.lower())
        else:
            print('No file:', filesADJRoot, pdb)

    print('len=',len(pdbList))

    if setName == 'ADJUSTED':
        filesRoot = filesADJRoot
        outPath = outAdjusted
        badAtoms = help.getBadAtomsList(filesADJRoot, pdbList, 0.05)  # Get the bad atoms list we will use to reduce the list further
    elif setName == 'REDUCED':
        outPath = outReduced
        badAtoms = help.getBadAtomsList(filesADJRoot, pdbList, 0.05)  # Get the bad atoms list we will use to reduce the list further
    elif setName == 'UNRESTRICTED':
        keepDisordered = True

    from PsuGeometry import GeoPdb as geopdb
    pdbmanager = geopdb.GeoPdbs(filesRoot, edDataPath, False, False, False, badAtoms)
    georep = psu.GeoReport(pdbList, filesRoot, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=True)
    dataCsv = georep.getGeoemtryCsv(geos, hueList)
    try:
        dataCsv['rid'] = dataCsv['rid'].astype(str)
        dataCsv['ID'] = dataCsv['pdbCode'] + dataCsv['chain'] + dataCsv['rid'] + dataCsv['aa']
    except:
        print('empty csv')

    #embellish with software and resolution
    pdbdata = pdbdata[['PDB','SOFTWARE','RES']]
    pdbdata.columns = ['PDB','Software','Resolution']
    dataCsv['PDB'] = dataCsv['pdbCode']
    dataCsv = dataCsv.set_index('PDB').join(pdbdata.set_index('PDB'))
    dataCsv['Software'] = dataCsv['Software'].str[:8]

    if setName != 'UNRESTRICTED':
        dataCsv = help.applyRestrictions(dataCsv)
        dataCsv.to_csv(outPath, index=False)
    else:
        dataCsv.to_csv(outPath, index=False)
        dataReduced = help.applyRestrictions(dataCsv)
        dataReduced.to_csv(outRestricted, index=False)

##########################################################################################
# This creates the file that incldues the diffs between reduced and adjusted
##########################################################################################
if runDiffMaker:
    dataA =  pd.read_csv(outReduced)
    dataB = pd.read_csv(outAdjusted)
    locs = []
    locs.append('ID')
    for geo in geos:
        dataA[geo+'_Orig'] = dataA[geo]
        dataB[geo + '_Adj'] = dataB[geo]
        locs.append(geo + '_Adj')
    dataB = dataB[locs]
    mergedDataSet = dataA
    mergedDataSet = mergedDataSet.set_index('ID').join(dataB.set_index('ID'))
    mergedDataSet = mergedDataSet.dropna()
    mergedDataSet['rid'] = mergedDataSet['rid'].astype(str)
    mergedDataSet['ID'] = mergedDataSet['pdbCode'] + mergedDataSet['chain'] + mergedDataSet['rid'] + mergedDataSet['aa']
    mergedDataSet.to_csv(outMerged, index=False)


def applyCis(aa, preomega):
    if aa != 'PRO':
        return aa
    if abs(preomega) > 100:
        return 'PRO'
    else:
        return 'CIS'

if runChangeProToCis:
    dataCsvAdjusted = pd.read_csv(outAdjusted)
    dataCsvReduced = pd.read_csv(outReduced)
    dataCsvRestricted = pd.read_csv(outRestricted)
    dataCsvUnrestricted = pd.read_csv(outUnrestricted)
    dataCsvMerged = pd.read_csv(outMerged)

    dataCsvAdjusted['aa'] = dataCsvAdjusted.apply(lambda row: applyCis(row['aa'], row['PREOMEGA']), axis=1)
    dataCsvReduced['aa'] = dataCsvReduced.apply(lambda row: applyCis(row['aa'], row['PREOMEGA']), axis=1)
    dataCsvRestricted['aa'] = dataCsvRestricted.apply(lambda row: applyCis(row['aa'], row['PREOMEGA']), axis=1)
    dataCsvUnrestricted['aa'] = dataCsvUnrestricted.apply(lambda row: applyCis(row['aa'], row['PREOMEGA']), axis=1)
    dataCsvMerged['aa'] = dataCsvMerged.apply(lambda row: applyCis(row['aa'], row['PREOMEGA']), axis=1)

    dataCsvAdjusted.to_csv(outAdjusted, index=False)
    dataCsvReduced.to_csv(outReduced, index=False)
    dataCsvRestricted.to_csv(outRestricted, index=False)
    dataCsvUnrestricted.to_csv(outUnrestricted, index=False)
    dataCsvMerged.to_csv(outMerged, index=False)

if runMakeAverages:
    dataCsvAdjusted = pd.read_csv(outAdjusted)
    dataCsvReduced = pd.read_csv(outReduced)
    dataCsvRestricted = pd.read_csv(outRestricted)
    dataCsvUnrestricted = pd.read_csv(outUnrestricted)
    dataCsvMerged = pd.read_csv(outMerged)

    dataCsvAdjustedAv = dataCsvAdjusted.groupby('pdbCode').mean()
    dataCsvReducedAv = dataCsvReduced.groupby('pdbCode').mean()
    dataCsvRestrictedAv = dataCsvRestricted.groupby('pdbCode').mean()
    dataCsvUnrestrictedAv = dataCsvUnrestricted.groupby('pdbCode').mean()
    dataCsvMergedAv = dataCsvMerged.groupby('pdbCode').mean()

    dataCsvAdjustedAv['aa'] = 'ALL'
    dataCsvReducedAv['aa'] = 'ALL'
    dataCsvRestrictedAv['aa'] = 'ALL'
    dataCsvUnrestrictedAv['aa'] = 'ALL'
    dataCsvMergedAv['aa'] = 'ALL'

    dataCsvAdjustedAv.to_csv(outAdjusted + 'mean.csv', index=True)
    dataCsvReducedAv.to_csv(outReduced + 'mean.csv', index=True)
    dataCsvRestrictedAv.to_csv(outRestricted + 'mean.csv', index=True)
    dataCsvUnrestrictedAv.to_csv(outUnrestricted + 'mean.csv', index=True)
    dataCsvMergedAv.to_csv(outMerged + 'mean.csv', index=True)

if runAppendDssp:
    pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

    fileUnrestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csv'
    fileRestricted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csv'
    fileReduced = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataMerged.csv'
    fileAdjusted = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csv'

    fileUnrestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataUnrestricted.csvmean.csv'
    fileRestrictedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataRestricted.csvmean.csv'
    fileReducedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataReduced.csvmean.csv'
    fileAdjustedAv = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/DataAdjusted.csvmean.csv'

    dataCsvAdjustedAv = pd.read_csv(fileAdjustedAv)
    dataCsvReducedAv = pd.read_csv(fileReducedAv)
    dataCsvRestrictedAv = pd.read_csv(fileRestrictedAv)
    dataCsvUnrestrictedAv = pd.read_csv(fileUnrestrictedAv)

    dataCsvAdjusted = pd.read_csv(fileAdjusted)
    dataCsvReduced = pd.read_csv(fileReduced)
    dataCsvRestricted = pd.read_csv(fileRestricted)
    dataCsvUnrestricted = pd.read_csv(fileUnrestricted)

    allList = []
    allList.append([dataCsvAdjusted,fileAdjusted])
    allList.append([dataCsvReduced,fileReduced])
    allList.append([dataCsvRestricted,fileRestricted])
    allList.append([dataCsvUnrestricted,fileUnrestricted])

    avgList = []
    avgList.append([dataCsvAdjustedAv,fileAdjustedAv])
    avgList.append([dataCsvReducedAv,fileReducedAv])
    avgList.append([dataCsvRestrictedAv,fileRestrictedAv])
    avgList.append([dataCsvUnrestrictedAv,fileUnrestrictedAv])

    # embellish with software and resolution
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_70.csv')
    pdbdata = pdbdata[['PDB', 'SOFTWARE', 'RES']]
    pdbdata.columns = ['PDB', 'Software', 'Resolution']

    #for dataBlob in avgList:
    #    dataCsv = dataBlob[0]
    #    dataPath = dataBlob[1]
    #    dataCsv['PDB'] = dataCsv['pdbCode']
    #    dataCsv = dataCsv.set_index('PDB').join(pdbdata.set_index('PDB'))
    #    dataCsv['Software'] = dataCsv['Software'].str[:8]
    #    dataCsv.to_csv(dataPath, index=False)

    #embellish with dssp - the dssp file was created ages ago from the linux laptop
    pdbdssp = pd.read_csv('C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/CsvGeos_BEST_Set0DSSPALL.csv')
    try:
        pdbdssp['rid'] = pdbdssp['rid'].astype(str)
        pdbdssp['DSSPID'] = pdbdssp['pdbCode'] + pdbdssp['chain'] + pdbdssp['rid']
        pdbdssp = pdbdssp[['DSSPID', 'dssp']]

    except:
        print('empty csv')

    for dataBlob in allList:
        dataCsv = dataBlob[0]
        dataPath = dataBlob[1]
        try:
            dataCsv['rid'] = dataCsv['rid'].astype(str)
            dataCsv['DSSPID'] = dataCsv['pdbCode'] + dataCsv['chain'] + dataCsv['rid']
        except:
            print('empty csv')

        dataCsv = dataCsv.set_index('DSSPID').join(pdbdssp.set_index('DSSPID'))
        dataCsv.to_csv(dataPath, index=False)