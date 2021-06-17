# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
import _Helpers as help
'''
TAU correlations
'''
###############################################################################################
def createGoodDensitySlices(pdbSet,geoset):
    import matplotlib.pyplot as plt
    plt.close('all')
    plt.clf()
    plt.cla()

    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_data/'
    if not 'RESTRICTED' in pdbSet:
        pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet + '/'

    edDataPath = help.rootPath + '/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesG/'
    loadPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'

    setFileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'
    print('Loading',loadPath + setFileName)
    dataBest = pd.read_csv(loadPath + setFileName)


    #Loop through the csv file and get the outliers
    aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN','ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    for atoms,geos in geoset:

        atomCe = atoms[0]
        atomLi = atoms[1]
        atomPl = atoms[2]

        ceAdd = 0
        liAdd = 0
        plAdd = 0

        cePlus = atomCe.find('+')
        ceMinus = atomCe.find('-')
        liPlus = atomLi.find('+')
        liMinus = atomLi.find('-')
        plPlus = atomPl.find('+')
        plMinus = atomPl.find('-')

        if cePlus > 0:
            ceAdd = int(atomCe[cePlus + 1:])
            atomCe = atomCe[:cePlus]
        if liPlus > 0:
            liAdd = int(atomLi[liPlus + 1:])
            atomLi = atomLi[:liPlus]
        if plPlus > 0:
            plAdd = int(atomPl[plPlus + 1:])
            atomPl = atomPl[:plPlus]
        if ceMinus > 0:
            ceAdd = -1 * int(atomCe[ceMinus + 1:])
            atomCe = atomCe[:ceMinus]
        if liMinus > 0:
            liAdd = -1 * int(atomLi[liMinus + 1:])
            atomLi = atomLi[:liMinus]
        if plMinus > 0:
            plAdd = -1 * int(atomPl[plMinus + 1:])
            atomPl = atomPl[:plMinus]

        print(atomCe, ceAdd)
        print(atomLi, liAdd)
        print(atomPl, plAdd)


        for geo in geos:
            dics = []
            slicesList = []
            for aa in aas:
                bestCut = dataBest.query("aa ==  '" + aa + "'")
                if len(bestCut) > 0:
                    print(aa,geo)
                    bestCut = bestCut.sort_values(geo)
                    dfHead = bestCut.head(1)
                    #print('HEAD',dfHead)
                    pdb = dfHead['pdbCode'].values[0]
                    rid = dfHead['rid'].values[0]
                    chain = dfHead['chain'].values[0]
                    if not [pdb, chain, rid] in slicesList:
                        slicesList.append([pdb, chain, rid,aa])

                    dicH = {}
                    dicH['pdbCode'] = pdb
                    dicH['chain'] = chain
                    dicH['rid'] = rid
                    dicH['geo'] = geo
                    dicH['aa'] = aa
                    dicH['value'] = dfHead[geo].values[0]

                    dfTail = bestCut.tail(1)
                    #print('TAIL',dfTail)
                    pdb = dfTail['pdbCode'].values[0]
                    rid = dfTail['rid'].values[0]
                    chain = dfTail['chain'].values[0]
                    if not [pdb, chain, rid] in slicesList:
                        slicesList.append([pdb, chain, rid,aa])

                    dicT = {}
                    dicT['pdbCode'] = pdb
                    dicT['chain'] = chain
                    dicT['rid'] = rid
                    dicT['aa'] = aa
                    dicT['geo'] = geo
                    dicT['value'] = dfTail[geo].values[0]

                    dics.append(dicH)
                    dics.append(dicT)
                    #print(dicH)
                    #print(dicT)

            dataFrame = pd.DataFrame.from_dict(dics)
            geoF = geo.replace('-', '')
            geoF = geoF.replace('+', '')
            geoF = geoF.replace(':', '')
            filePath = printPath + 'GoodOutliers_' + pdbSet + '_' + geoF +  '.csv'
            print('...printing', filePath)
            dataFrame.to_csv(filePath, index=False)
            print(dataFrame)




            bigstring = ""

            for sl in slicesList:
                print(sl)
                georep = psu.GeoReport([sl[0]],pdbDataPath,edDataPath,printPath,ed=False,dssp=False)
                pdbmanager = geopdb.GeoPdbs(pdbDataPath,edDataPath,ed=False,dssp=False)
                apdb = pdbmanager.getPdb(sl[0],True)
                pdbcsv = apdb.getDataFrame()

                queryC = 'rid==' + str(sl[2]+ceAdd) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomCe + '"'
                queryL = 'rid==' + str(sl[2]+liAdd) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomLi + '"'
                queryP = 'rid==' + str(sl[2]+plAdd) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomPl + '"'
                dataC = pdbcsv.query(queryC)
                dataL = pdbcsv.query(queryL)
                dataP = pdbcsv.query(queryP)

                if len(dataC) > 0 and len(dataL) > 0 and len(dataP)>0:
                    cx = round(dataC['x'].values[0],3)
                    cy = round(dataC['y'].values[0], 3)
                    cz = round(dataC['z'].values[0], 3)
                    lx = round(dataL['x'].values[0], 3)
                    ly = round(dataL['y'].values[0], 3)
                    lz = round(dataL['z'].values[0], 3)
                    px = round(dataP['x'].values[0], 3)
                    py = round(dataP['y'].values[0], 3)
                    pz = round(dataP['z'].values[0], 3)

                    row = sl[0] + "," + sl[3] + sl[1] + str(sl[2]) + "," + str(cx) + "," + str(cy) + "," + str(cz)
                    row +=  "," + str(lx) + "," + str(ly) + "," + str(lz)
                    row += "," + str(px) + "," + str(py) + "," + str(pz)

                    print(row)

                    bigstring += row + '\n'

            print("########RESULTS#########")
            print("")
            print(bigstring)

            geoG = atomCe + atomLi + atomPl
            geoG = geoG.replace('-', '')
            geoG = geoG.replace('+', '')
            geoG = geoG.replace(':', '')

            f = open(printPath + 'GoodSlice_' + pdbSet + '_' + geoF + '_' + geoG + '.txt', "w")
            f.write(bigstring)
            f.close()


########################################
def getExtremeAdjustedResidues(pdbSet,geo, geoRange):
    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    loadPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesC/'
    setFileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'
    print('Loading', loadPath + setFileName)
    dataBest = pd.read_csv(loadPath + setFileName)
    # we are now going to cut this down to the extremes
    query = "(`C:O` < " + str(geoRange[0]) + ") or (`C:O` > " + str(geoRange[1]) + ")"
    print(query)
    dataCut = dataBest.query(query)
    dataCut = dataCut.sort_values(geo)

    # Loop through the csv file and get the outliers
    aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG',
           'SER', 'THR', 'VAL', 'TRP', 'TYR']





    dics = []
    slicesList = []
    for aa in aas:
        bestCut = dataCut.query("aa ==  '" + aa + "'")
        pdbs = bestCut['pdbCode'].values
        rids = bestCut['rid'].values
        chains = bestCut['chain'].values
        vals = bestCut[geo].values

        for i in range(0, len(pdbs)):
            pdb = pdbs[i]
            rid = rids[i]
            chain = chains[i]
            val = vals[i]
            print(aa, geo)
            if not [pdb, chain, rid] in slicesList:
                slicesList.append([pdb, chain, rid, aa])

            dicH = {}
            dicH['pdbCode'] = pdb
            dicH['chain'] = chain
            dicH['rid'] = rid
            dicH['geo'] = geo
            dicH['aa'] = aa
            dicH['value'] = val

            dics.append(dicH)

    dataFrame = pd.DataFrame.from_dict(dics)
    geoF = geo.replace('-', '')
    geoF = geoF.replace('+', '')
    geoF = geoF.replace(':', '')
    filePath = printPath + 'AdjExtreme_' + pdbSet + '_' + geoF + '.csv'
    print('...printing', filePath)
    dataFrame.to_csv(filePath, index=False)
    print(dataFrame)

    return slicesList, dataFrame


def createExtremeAdjustedSlices(pdbSet,slices,atoms,geo,dataFrame,tag):

    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_data/'
    if not 'RESTRICTED' in pdbSet:
        pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet + '/'

    edDataPath = help.rootPath + '/ProteinDataFiles/ccp4_data/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesC/'

    geoF = geo.replace('-', '')
    geoF = geoF.replace('+', '')
    geoF = geoF.replace(':', '')

    atomCe = atoms[0]
    atomLi = atoms[1]
    atomPl = atoms[2]

    ceAdd = 0
    liAdd = 0
    plAdd = 0

    cePlus = atomCe.find('+')
    ceMinus = atomCe.find('-')
    liPlus = atomLi.find('+')
    liMinus = atomLi.find('-')
    plPlus = atomPl.find('+')
    plMinus = atomPl.find('-')

    if cePlus > 0:
        ceAdd = int(atomCe[cePlus + 1:])
        atomCe = atomCe[:cePlus]
    if liPlus > 0:
        liAdd = int(atomLi[liPlus + 1:])
        atomLi = atomLi[:liPlus]
    if plPlus > 0:
        plAdd = int(atomPl[plPlus + 1:])
        atomPl = atomPl[:plPlus]
    if ceMinus > 0:
        ceAdd = -1 * int(atomCe[ceMinus + 1:])
        atomCe = atomCe[:ceMinus]
    if liMinus > 0:
        liAdd = -1 * int(atomLi[liMinus + 1:])
        atomLi = atomLi[:liMinus]
    if plMinus > 0:
        plAdd = -1 * int(atomPl[plMinus + 1:])
        atomPl = atomPl[:plMinus]

    print(atomCe, ceAdd)
    print(atomLi, liAdd)
    print(atomPl, plAdd)

    pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, ed=False, dssp=False)
    pdbmanager.clear()
    pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, ed=False, dssp=False)

    bigstring = ""
    for sl in slices:
        print(sl)

        apdb = pdbmanager.getPdb(sl[0],True)
        pdbcsv = apdb.getDataFrame()

        queryC = 'rid==' + str(sl[2]+ceAdd) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomCe + '"'
        queryL = 'rid==' + str(sl[2]+liAdd) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomLi + '"'
        queryP = 'rid==' + str(sl[2]+plAdd) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomPl + '"'
        dataC = pdbcsv.query(queryC)
        dataL = pdbcsv.query(queryL)
        dataP = pdbcsv.query(queryP)

        if len(dataC) > 0 and len(dataL) > 0 and len(dataP)>0:
            cx = round(dataC['x'].values[0],3)
            cy = round(dataC['y'].values[0], 3)
            cz = round(dataC['z'].values[0], 3)
            lx = round(dataL['x'].values[0], 3)
            ly = round(dataL['y'].values[0], 3)
            lz = round(dataL['z'].values[0], 3)
            px = round(dataP['x'].values[0], 3)
            py = round(dataP['y'].values[0], 3)
            pz = round(dataP['z'].values[0], 3)

            row = sl[0] + "," + sl[3] + sl[1] + str(sl[2]) + "," + str(cx) + "," + str(cy) + "," + str(cz)
            row +=  "," + str(lx) + "," + str(ly) + "," + str(lz)
            row += "," + str(px) + "," + str(py) + "," + str(pz)

            print(row)

            bigstring += row + '\n'

    print("########RESULTS#########")
    print("")
    print(bigstring)

    geoG = atomCe + atomLi + atomPl
    geoG = geoG.replace('-', '')
    geoG = geoG.replace('+', '')
    geoG = geoG.replace(':', '')

    f = open(printPath + 'AdjExtreme' + pdbSet + '_' + geoF + '_' + geoG + '.txt', "w")
    f.write(bigstring)
    f.close()


########################################


