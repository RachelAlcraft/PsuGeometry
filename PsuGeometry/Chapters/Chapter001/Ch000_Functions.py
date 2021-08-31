'''
Shared functions file
'''
import pandas as pd
from PsuGeometry import GeoReport as psu

##### globals #############################
loadPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'
printPathLx = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
filesPDBRoot ='C:/Dev/Github/ProteinDataFiles/pdb_data/'
pdbDataPathLx =  '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
#filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fo5_ADJ/' #adjusted on Fo at 5 degrees thevenaz (11.8.21)
#filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fo3_ADJ/' #adjusted on Fo at 3 degrees thevenaz (18.8.21) fix to numerical gap
filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/FO5_ADJ/' #adjusted on Fo at 5 degrees thevenaz (19.8.21) fix to numerical gap
filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fo509_ADJ/' #adjusted on Fo at 5 degrees thevenaz (27.8.21) with <=0.9 A and added side chain atoms for further O analysis
############################################

def getBadAtomsListFromFile():
    pdbbaddata = pd.read_csv(loadPath + 'AllBadAtoms_Maxima.csv')
    pdboccdata = pd.read_csv(loadPath + 'AllBadAtoms_Occupant.csv')
    expdbbaddata = pd.read_csv(loadPath + 'ExtraBadAtoms_Maxima.csv')
    pdbbadList = pdbbaddata['BAD'].tolist()[0:]
    pdboccList = pdboccdata['BAD'].tolist()[0:]
    expdbbadList = expdbbaddata['BAD'].tolist()[0:]
    return pdbbadList + expdbbadList + pdboccList

def getPDBList100():
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_100.csv')
    pdbList = pdbdata['PDB'].tolist()[0:]
    return pdbList

def getPDBList():
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_70.csv')
    pdbListA = pdbdata['PDB'].tolist()[0:]
    #pdbListA = ['5xvt', '6fmc', '4u9h', '4y9w']
    pdbListIn = []
    for pdb in pdbListA:
        import os.path
        isFile = os.path.isfile((filesADJRoot + 'pdb' + pdb + '.ent').lower())
        isOcc = os.path.isfile((filesADJRoot + 'OccupancyMaxima_' + pdb + '.csv').lower())
        if isFile and isOcc:
            pdbListIn.append(pdb.lower())
        else:
            print('No file:', (filesADJRoot + 'pdb' + pdb + '.ent').lower())
    #pdbListIn.remove('1ot6')#2 models
    #pdbListIn.remove('2ce2')#2 models
    print(pdbListIn)
    return pdbListIn

'''
def getBadAtomsList(pdbDataPath, pdbList,maxDiff):
    allRealPdbs = []
    badRealPdbs = []
    occRealPdbs = []
    #printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
    for pdb in pdbList:
        #The differences that were found
        realFileName = 'MaximaDifferences_' + pdb + '.csv'
        print('Reading ',pdbDataPath + realFileName)
        realData = pd.read_csv(pdbDataPath + realFileName)
        allRealPdbs.append(realData)
        #And the bad data that failed
        badRealFileName = 'MaximaDifferences_' + pdb + '_BAD_.csv'
        badRealData = pd.read_csv(pdbDataPath + badRealFileName)
        badRealPdbs.append(badRealData)
        # And the occupancy checker file
        occRealFileName = 'OccupancyMaxima_' + pdb + '.csv'
        print(pdbDataPath + occRealFileName)
        occRealData = pd.read_csv(pdbDataPath + occRealFileName)
        occRealData = occRealData.query('Fraction < 1')
        occRealPdbs.append(occRealData)
    # append them all
    realCsv = pd.concat(allRealPdbs, axis=0, sort=False)
    badRealCsv = pd.concat(badRealPdbs, axis=0, sort=False)
    occRealCsv = pd.concat(occRealPdbs, axis=0, sort=False)

    # Make the file into the format 5kwmA222HG23 and filter on maxDiff
    badpdbs = []
    # 1.  deal with the list of bad first
    pdbsBad = badRealCsv['pdbCode'].values
    chainsBad = badRealCsv['Chain'].values
    resBad = badRealCsv['ResNo'].values
    atomsBad = badRealCsv['AtomType'].values
    # 2.  deal with the list of multiple occupants
    pdbsOcc = occRealCsv['pdbCode'].values
    chainsOcc = occRealCsv['Chain'].values
    resOcc = occRealCsv['ResNo'].values
    atomsOcc = occRealCsv['AtomType'].values
    # 3. Create list of unwanted atoms
    for i in range(0,len(pdbsOcc)):
        badpdbs.append((pdbsOcc[i] + chainsOcc[i] + str(resOcc[i]) + str(atomsOcc[i])))
    for i in range(0,len(pdbsBad)):
        badpdbs.append((pdbsBad[i] + chainsBad[i] + str(resBad[i]) + str(atomsBad[i])))

    pdbsBad = realCsv['pdbCode'].values
    chainsBad = realCsv['Chain'].values
    resBad = realCsv['ResNo'].values
    atomsBad = realCsv['AtomType'].values
    diffsBad = realCsv['Difference'].values
    for i in range(0, len(pdbsBad)):
        diff = diffsBad[i]
        if diff > maxDiff:
            badpdbs.append((pdbsBad[i] + chainsBad[i] + str(resBad[i]) + atomsBad[i]))

    # print them
    with open(printPath + 'badatoms.csv', 'w') as f:
        f.write("BAD" + '\n')
        for atm in badpdbs:
            f.write(atm + '\n')

    return badpdbs 
'''
def applyRestrictions(dataCsv,pdbs,bfactor,disorder,occupancy):

    dataReduced = dataCsv
    if disorder:
        dataReduced = dataReduced.query("disordered=='N'")

    if bfactor:
        dataReduced = dataReduced.query("bfactor <= 18")

    if pdbs:
        # these have multiple models and the occupants are not marked as disordered
        dataReduced = dataReduced.query("pdbCode != '1ot6'")
        dataReduced = dataReduced.query("pdbCode != '2ce2'")
        # these have average bond lengths that are clearly non-standard
        dataReduced = dataReduced.query("pdbCode != '2cnq'")
        dataReduced = dataReduced.query("pdbCode != '5o45'")
        dataReduced = dataReduced.query("pdbCode != '3ql9'")
        dataReduced = dataReduced.query("pdbCode != '2fn3'")
        dataReduced = dataReduced.query("pdbCode != '4r5r'")

    if occupancy:
        dataReduced = dataReduced.query("occupancy == 1")  # some occupnats are not marked as disordered


    return dataReduced

def makeCsv(pdbSet, pdbListIn,geos,badAtoms,disordered):
    print('Getting CSV for',pdbSet)
    pdbDataPath = filesPDBRoot
    if pdbSet == 'ADJUSTED':
        pdbDataPath = filesADJRoot

    from PsuGeometry import GeoPdb as geopdb
    pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, False, False, disordered,badAtoms)
    pdbmanager.clear()
    pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, False, False, disordered, badAtoms)

    pdbList = [] #this is so we don't default to getting a pdb file from somewhere we don;t want
    for pdb in pdbListIn:
        import os.path
        filePdb = pdbDataPath + 'pdb' + pdb + '.ent'
        #print('- Adding to csv',filePdb)
        if os.path.isfile((filePdb).lower()):
            pdbList.append(pdb.lower())
        else:
            print('No file:', pdbDataPath, pdb)
    pdbList.sort()
    hueList = ['aa', 'rid', 'bfactor', 'pdbCode', 'bfactorRatio', 'disordered','occupancy']
    georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=disordered)
    print('geoList',geos)
    dataBest = georep.getGeoemtryCsv(geos, hueList, -1, allAtoms=True)
    try:
        dataBest['rid'] = dataBest['rid'].astype(str)
        dataBest['ID'] = dataBest['pdbCode'] + dataBest['chain'] + dataBest['rid'] + dataBest['aa']
    except:
        print('empty csv')
    return dataBest

def embellishCsv(csvData):
    #we assume the data is a csv data file and has fields dataBest['ID'] = dataBest['pdbCode'] + dataBest['chain'] + dataBest['rid'] + dataBest['aa']
    # the dssp file was created ages ago from the linux laptop
    pdbdssp = pd.read_csv('C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/CsvGeos_BEST_Set0DSSPALL.csv')
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_100.csv')
    atomdata = pd.read_csv(loadPath + "AllAtoms_Maxima.csv")

    # embellish with dssp, resolution and software and Maxima Diffeerence
    #try:
    if True:
        # make the atom data into a rid data using the C distance as the MaxaDiff distance we will use for hues (distance between Laplacian and Radiant maxima)
        atomdata['ResNo'] = atomdata['ResNo'].astype(str)
        atomdata = atomdata.query("AtomType=='C'")
        atomdata['ROWID'] = atomdata['pdbCode'] + atomdata['Chain'] + atomdata['ResNo']
        atomdata = atomdata[['ROWID', 'LapDiff']]
        #atomdata.columns = [['ROWID', 'MaxaDiff']]

        # make the dssp data into the right shape
        pdbdssp['rid'] = pdbdssp['rid'].astype(str)
        pdbdssp['PDB'] = pdbdssp['pdbCode']
        pdbdssp['ROWID'] = pdbdssp['pdbCode'] + pdbdssp['chain'] + pdbdssp['rid']
        pdbdssp['DSSPID'] = pdbdssp['pdbCode'] + pdbdssp['chain'] + pdbdssp['rid']
        pdbdssp = pdbdssp[['DSSPID', 'dssp','ROWID','PDB']]

        #Combine atom and dssp
        atomdata.to_csv(loadPath + "TempATMa.csv", index=False)
        pdbdssp.to_csv(loadPath + "TempDSSPa.csv", index=False)
        pdbdssp = pdbdssp.set_index('ROWID').join(atomdata.set_index('ROWID'))

        #make the pdb data in shape
        pdbdata = pdbdata[['PDB', 'SOFTWARE', 'RES']]

        #Join with pdb and atoms
        pdbdssp = pdbdssp.set_index('PDB').join(pdbdata.set_index('PDB'))
        pdbdssp['SOFTWARE'] = pdbdssp['SOFTWARE'].str[:8]
        pdbdssp = pdbdssp.dropna()

        pdbdssp.to_csv(loadPath + "TempDSSPb.csv", index=False)

        #make the input csv data the right shape
        csvData.to_csv(loadPath + "TempDATAa.csv", index=False)
        csvData['rid'] = csvData['rid'].astype(str)
        csvData['DSSPID'] = csvData['pdbCode'] + csvData['chain'] + csvData['rid']
        #DEBUG
        csvData.to_csv(loadPath + "TempDATAb.csv", index=False)
        #DUBUG
        csvData = csvData.set_index('DSSPID').join(pdbdssp.set_index('DSSPID'))
    #except:
    #    print('empty csv')
    return csvData


def trioReports(nameDataA, nameDataB, nameDataC, nameDataD ,trios, title,printPath,fileName, perCategory=''):

    dataA = nameDataA[1]
    dataB = nameDataB[1]
    dataC = nameDataC[1]
    dataD = nameDataD[1]

    nameA = nameDataA[0]
    nameB = nameDataB[0]
    nameC = nameDataC[0]
    nameD = nameDataD[0]


    cats = []
    if perCategory != '':
        dataA[perCategory] = dataA[perCategory].fillna('?')
        cats = dataA[perCategory].values
        cats = list(set(cats))
        print('Splitting on categories', cats)
        try:
            cats.sort()
        except:
            print('unsorted')
        print('Splitting on categories',cats)

    georep = psu.GeoReport([],"", "", printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

    for trio in trios:
        if perCategory!='':
            for cat in cats:
                dataCutA = dataA.query(perCategory + " ==  '" + cat + "'")
                dataCutB = dataB.query(perCategory + " ==  '" + cat + "'")
                dataCutC = dataC.query(perCategory + " ==  '" + cat + "'")
                dataCutD = dataD.query(perCategory + " ==  '" + cat + "'")
                if len(trio) == 4:
                    pal = 'jet'
                    if trio[3]:
                        pal = 'jet_r'

                    georep.addScatter(data=dataCutA, geoX=trio[0], geoY=trio[1], hue=trio[2], title=cat + ':' + trio[0] + ':'+ trio[1] , categorical = trio[3],palette=pal, sort='NON')
                    georep.addScatter(data=dataCutB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=cat + ':' + trio[0] + ':' + trio[1], categorical = trio[3],palette=pal, sort='NON')
                    georep.addScatter(data=dataCutC, geoX=trio[0], geoY=trio[1], hue=trio[2],title=cat + ':' + trio[0] + ':' + trio[1], categorical = trio[3],palette=pal, sort='NON')
                    georep.addScatter(data=dataCutD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=cat + ':' + trio[0] + ':' + trio[1], categorical = trio[3],palette=pal, sort='NON')
                else:
                    georep.addHistogram(data=dataCutA, geoX=trio[0],title=nameA + ' ' + cat + ':' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutB, geoX=trio[0],title=nameB + ' ' + cat + ':' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutC, geoX=trio[0],title=nameC + ' ' + cat + ':' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutD, geoX=trio[0],title=nameD + ' ' + cat + ':' + trio[0], hue='ID')
        else:
            if len(trio) == 4:
                cat = False
                pal = 'jet'
                if trio[3]:
                    cat = True
                    pal = 'jet_r'
                if trio[2] == "COUNT":
                    georep.addHexBins(data=dataA, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ ' ' + nameA, palette='cubehelix_r')
                    georep.addHexBins(data=dataB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1] + ' ' + nameB, palette='cubehelix_r')
                    georep.addHexBins(data=dataC, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1] + ' ' + nameC, palette='cubehelix_r')
                    georep.addHexBins(data=dataD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1] + ' ' + nameD, palette='cubehelix_r')
                else:
                    georep.addScatter(data=dataA, geoX=trio[0], geoY=trio[1], hue=trio[2],title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + nameA, categorical = cat,palette=pal, sort='NON')
                    georep.addScatter(data=dataB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + nameB, categorical = cat, palette=pal, sort='NON')
                    georep.addScatter(data=dataC, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + nameC, categorical = cat, palette=pal, sort='NON')
                    georep.addScatter(data=dataD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=trio[0] + '|' + trio[1]+ '|' + trio[2] + ' ' + nameD, categorical = cat, palette=pal, sort='NON')
            else:
                georep.addHistogram(data=dataA, geoX=trio[0], title=nameA + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataB, geoX=trio[0], title=nameB + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataC, geoX=trio[0], title=nameC + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataD, geoX=trio[0], title=nameD + ' ' + trio[0], hue='ID')



    georep.printToHtml(title, 4, fileName)

def trioHexbins(nameDataA, nameDataB, nameDataC, nameDataD ,trios, title,printPath,fileName, perCategory=''):

    dataA = nameDataA[1]
    dataB = nameDataB[1]
    dataC = nameDataC[1]
    dataD = nameDataD[1]

    nameA = nameDataA[0]
    nameB = nameDataB[0]
    nameC = nameDataC[0]
    nameD = nameDataD[0]


    cats = []
    if perCategory != '':
        dataA[perCategory] = dataA[perCategory].fillna('?')
        cats = dataA[perCategory].values
        cats = list(set(cats))
        print('Splitting on categories', cats)
        try:
            cats.sort()
        except:
            print('unsorted')
        print('Splitting on categories',cats)

    georep = psu.GeoReport([],"", "", printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

    pal = 'jet'
    for trio in trios:
        if perCategory!='':
            for cat in cats:
                dataCutA = dataA.query(perCategory + " ==  '" + cat + "'")
                dataCutB = dataB.query(perCategory + " ==  '" + cat + "'")
                dataCutC = dataC.query(perCategory + " ==  '" + cat + "'")
                dataCutD = dataD.query(perCategory + " ==  '" + cat + "'")
                if len(trio) == 3:
                    georep.addHexBins(data=dataCutA,geoX=trio[0],geoY=trio[1],hue=trio[2],title=cat + ":" + nameA + ':' + trio[0] + ':'+ trio[1],palette=pal)
                    georep.addHexBins(data=dataCutB, geoX=trio[0], geoY=trio[1], hue=trio[2],title=cat + ":" + nameB + ':' + trio[0] + ':' + trio[1], palette=pal)
                    georep.addHexBins(data=dataCutC, geoX=trio[0], geoY=trio[1], hue=trio[2],title=cat + ":" + nameC + ':' + trio[0] + ':' + trio[1], palette=pal)
                    georep.addHexBins(data=dataCutD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=cat + ":" + nameD + ':' + trio[0] + ':' + trio[1], palette=pal)
                else:
                    georep.addHistogram(data=dataCutA, geoX=trio[0],title=nameA + ' ' + cat + ':' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutB, geoX=trio[0],title=nameB + ' ' + cat + ':' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutC, geoX=trio[0],title=nameC + ' ' + cat + ':' + trio[0], hue='ID')
                    georep.addHistogram(data=dataCutD, geoX=trio[0],title=nameD + ' ' + cat + ':' + trio[0], hue='ID')
        else:
            if len(trio) == 3:
                georep.addHexBins(data=dataA, geoX=trio[0], geoY=trio[1], hue=trio[2], title=nameA + ':' + trio[0] + ':' + trio[1], palette=pal)
                georep.addHexBins(data=dataB, geoX=trio[0], geoY=trio[1], hue=trio[2], title=nameB + ':' + trio[0] + ':' + trio[1], palette=pal)
                georep.addHexBins(data=dataC, geoX=trio[0], geoY=trio[1], hue=trio[2],  title=nameC + ':' + trio[0] + ':' + trio[1], palette=pal)
                georep.addHexBins(data=dataD, geoX=trio[0], geoY=trio[1], hue=trio[2], title=nameD + ':' + trio[0] + ':' + trio[1], palette=pal)

            else:
                georep.addHistogram(data=dataA, geoX=trio[0], title=nameA + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataB, geoX=trio[0], title=nameB + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataC, geoX=trio[0], title=nameC + ' ' + trio[0], hue='ID')
                georep.addHistogram(data=dataD, geoX=trio[0], title=nameD + ' ' + trio[0], hue='ID')



    georep.printToHtml(title, 4, fileName)