
import pandas as pd
from PsuGeometry import GeoReport as psu

rootPath = 'C:/Dev/Github/'

def applyCis(aa,preomega):
    if aa != 'PRO':
        return aa
    if abs(preomega) > 100:
        return 'PRO'
    else:
        return 'CISPRO'

def applyDifference(id, queryKey, dataMaxi):
    print(id)
    dataMx = dataMaxi.query('ID == "' + id + '"')
    srs = dataMx[queryKey]
    return srs.mean()


def getList(listName,cutoff,pdbSet):
    pdbList = []
    if listName == 'SMALLEST':
        pdbList = ['6eex', '6mw1', '6mw2', '6mvz','2ol9','1akg','5vsg','6mw0']
    elif listName == 'EVIDENCED':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_Evidenced.csv')  # This is a list of pdbs <= 1.1A non homologous to 90%
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == '30':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_30.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == '40':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_40.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == '50':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_50.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == '70':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_70.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == '90':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_90.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == '95':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_95.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == '100':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_100.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]
    elif listName == 'TOP20':
        pdbdata = pd.read_csv('../../PdbLists/Pdbs_Top20.csv')
        pdbList = pdbdata['PDB'].tolist()[0:]


    if cutoff > 0:
        pdbList = pdbList[0:cutoff]

    if pdbSet == '':
        return pdbList
    else:
        pdbDataPath = rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet + '/'
        if pdbSet == 'PDB':
            pdbDataPath = rootPath + '/ProteinDataFiles/pdb_data/'
        pdbListCut = []
        for pdb in pdbList:
            import os.path
            filePdb = pdbDataPath + 'pdb' + pdb + '.ent'
            #print(filePdb)
            if os.path.isfile((filePdb).lower()):
                pdbListCut.append(pdb.lower())
            else:
                print('No file:', pdbDataPath, pdb)
        return pdbListCut

def getPdbEmbellishment(pdb,field):
    pdbList = []
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_100.csv')
    pdbRow = pdbdata.query("PDB == '" + pdb + "'")
    pdbData = pdbRow[field]
    return pdbData


def addMaximaDiffs(dataPdb, dataMaxima):
    try:
        dataMaxima['ResNo'] = dataMaxima['ResNo'].astype(str)
        dataMaxima['ID'] = dataMaxima['pdbCode'] + dataMaxima['Chain'] + dataMaxima['ResNo'] + dataMaxima['AA']
        dataDiffs = dataMaxima.groupby(['ID']).agg({'GridDistance':['mean']}).reset_index()
        dataDiffs.columns = ['ID','AvgMaximaDiffs']
        #print(dataDiffs)
        print('Applying maxima differences...')
        dataPdb['AvgMaximaDiffs'] = dataDiffs['AvgMaximaDiffs']#dataPdb.apply(lambda row: applyDifference(row['ID'],'Difference', dataMaxima), axis=1)
        #print('Applying sample differences...')
        #dataPdb['PdbDiff'] = dataPdb.apply(lambda row: applyDifference(row['ID'], 'GridDistance', dataMaxima), axis=1)
        #print('Applying interp differences...')
        #dataPdb['InterpDiff'] = dataPdb.apply(lambda row: applyDifference(row['ID'], 'BGridDistance', dataMaxima), axis=1)

        return dataPdb
    except:
        print('empty csv')
        return dataPdb

def getMaximaDiffs(pdbSet, pdbList, includeFake):
    pdbDataPath = rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet
    allRealPdbs = []
    allFakePdbs = []
    badRealPdbs = []
    badFakePdbs = []
    occRealPdbs = []
    for pdb in pdbList:
        #The differences that were found
        realFileName = '_ADJ/MaximaDifferences_' + pdb + '.csv'
        print('Reading ',pdbDataPath + realFileName)
        realData = pd.read_csv(pdbDataPath + realFileName)
        allRealPdbs.append(realData)
        #And the bad data that failed
        badRealFileName = '_ADJ/MaximaDifferences_' + pdb + '_BAD_.csv'
        badRealData = pd.read_csv(pdbDataPath + badRealFileName)
        #print(badRealData)
        badRealPdbs.append(badRealData)
        # And the occupancy checker file
        if not includeFake:
            occRealFileName = '_ADJ/OccupancyMaxima_' + pdb + '.csv'
            print(pdbDataPath + occRealFileName)
            occRealData = pd.read_csv(pdbDataPath + occRealFileName)
            occRealData = occRealData.query('Fraction < 1')
            occRealPdbs.append(occRealData)

        if includeFake:
            # The differences that were found
            fakeFileName = '_ADJ/MaximaDifferences_' + pdb + '.csv'
            print('Reading ', pdbDataPath + fakeFileName)
            fakeData = pd.read_csv(pdbDataPath + fakeFileName)
            allFakePdbs.append(fakeData)
            # And the bad data that failed
            badFakeFileName = '_ADJ/MaximaDifferences_' + pdb + '_BAD_.csv'
            badFakeData = pd.read_csv(pdbDataPath + badFakeFileName)
            badFakePdbs.append(badFakeData)

    # append them all
    realCsv = pd.concat(allRealPdbs, axis=0, sort=False)
    badRealCsv = pd.concat(badRealPdbs, axis=0, sort=False)
    if not includeFake:
        occRealCsv = pd.concat(occRealPdbs, axis=0, sort=False)

    if includeFake:
        fakeCsv = pd.concat(allFakePdbs, axis=0, sort=False)
        badFakeCsv = pd.concat(badFakePdbs, axis=0, sort=False)
        return realCsv,badRealCsv,fakeCsv,badFakeCsv
    else:
        return realCsv, badRealCsv, occRealCsv


def getBadList(dataDiffs, dataBad,dataOcc, maxDiff):
    # Make the file into the format 5kwmA222HG23 and filter on maxDiff
    badpdbs = []
    # 1.  deal with the list of bad first
    pdbsBad = dataBad['pdbCode'].values
    chainsBad = dataBad['Chain'].values
    resBad = dataBad['ResNo'].values
    atomsBad = dataBad['AtomType'].values
    # 2.  deal with the list of multiple occupants
    pdbsOcc = dataOcc['pdbCode'].values
    chainsOcc = dataOcc['Chain'].values
    resOcc = dataOcc['ResNo'].values
    atomsOcc = dataOcc['AtomType'].values
    #print(dataBad)
    #print(atomsBad)
    # Create list of unwanted atoms
    for i in range(0,len(pdbsOcc)):
        badpdbs.append((pdbsOcc[i] + chainsOcc[i] + str(resOcc[i]) + str(atomsOcc[i])))
    for i in range(0,len(pdbsBad)):
        badpdbs.append((pdbsBad[i] + chainsBad[i] + str(resBad[i]) + str(atomsBad[i])))

    pdbsBad = dataDiffs['pdbCode'].values
    chainsBad = dataDiffs['Chain'].values
    resBad = dataDiffs['ResNo'].values
    atomsBad = dataDiffs['AtomType'].values
    diffsBad = dataDiffs['Difference'].values
    for i in range(0, len(pdbsBad)):
        diff = diffsBad[i]
        if diff > maxDiff:
            #print('B',pdbsBad[i], chainsBad[i], resBad[i], atomsBad[i])
            badpdbs.append((pdbsBad[i] + chainsBad[i] + str(resBad[i]) + atomsBad[i]))

    return badpdbs


def getCsv(pdbSet, pdbListIn,geos,badAtoms,reloadPdb, reloadCsv,aa='ALL',includeCis=False,allAtoms=False, bFactorFactor=1.3,cutoff=0):
    print('Getting CSV for',pdbSet)
    pdbDataPath = rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    if pdbSet == 'PDB':
        pdbDataPath = rootPath + '/ProteinDataFiles/pdb_data/'

    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataK/'

    fileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'

    if reloadCsv:
        from PsuGeometry import GeoPdb as geopdb
        pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, False, False, False,badAtoms)
        if reloadPdb:
            pdbmanager.clear()
            pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, False, False, False, badAtoms)

        #pdbdata = pd.read_csv('../../PdbLists/Pdbs_Evidenced.csv')  # This is a list of pdbs <= 1.1A non homologous to 90%
        #pdbListIn = pdbdata['PDB'].tolist()[0:]
        #if cutoff > 0:
        #    pdbListIn = pdbdata['PDB'].tolist()[0:cutoff]


        pdbList = []
        for pdb in pdbListIn:
            import os.path
            filePdb = pdbDataPath + 'pdb' + pdb + '.ent'
            #print('- Adding to csv',filePdb)
            if os.path.isfile((filePdb).lower()):
                pdbList.append(pdb.lower())
            else:
                print('No file:', pdbDataPath, pdb)

        pdbList.sort()

        hueList = ['aa', 'rid', 'bfactor', 'pdbCode', 'bfactorRatio', 'disordered']
        georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=allAtoms)

        if includeCis:
            geos.append('CA-1:C-1:N:CA')

        print('geoList',geos)
        dataBest = georep.getGeoemtryCsv(geos, hueList, bFactorFactor, allAtoms=allAtoms,restrictedAa=aa)
        try:
            dataBest['rid'] = dataBest['rid'].astype(str)
            dataBest['ID'] = dataBest['pdbCode'] + dataBest['chain'] + dataBest['rid'] + dataBest['aa']
        except:
            print('empty csv')
    else:
        dataBest = pd.read_csv(loadPath + fileName)

    #aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    if includeCis:
        dataBest['aa'] = dataBest.apply(lambda row: applyCis(row['aa'], row['CA-1:C-1:N:CA']), axis=1)

    if aa !='ALL':
        dataBest = dataBest.query('aa == "' + aa + '"')

    return dataBest

def comparisonDataSet(dataA, dataB, geos,filePath):

    locs = []
    locs.append('ID')
    for geo in geos:
        dataA[geo+'_Orig'] = dataA[geo]
        dataB[geo + '_Adj'] = dataB[geo]
        locs.append(geo + '_Adj')

    dataB = dataB[locs]


    mergedDataSet = dataA
    mergedDataSet = mergedDataSet.set_index('ID').join(dataB.set_index('ID'))

    mergedDataSet['PDB'] =mergedDataSet['pdbCode']

    pdbdata = pd.read_csv('../../PdbLists/Pdbs_100.csv')
    pdbdata = pdbdata[['PDB','SOFTWARE','RES']]

    mergedDataSet = mergedDataSet.set_index('PDB').join(pdbdata.set_index('PDB'))
    mergedDataSet['SOFTWARE'] = mergedDataSet['SOFTWARE'].str[:8]

    mergedDataSet = mergedDataSet.dropna()

    #mergedDataSet['SOFTWARE'] = mergedDataSet.apply(lambda row: getPdbEmbellishment(row['pdbCode'], 'SOFTWARE'), axis=1)

    mergedDataSet.to_csv(filePath, index=False)
    return mergedDataSet

