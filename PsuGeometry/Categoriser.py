# -- ©Rachel Alcraft 2021, PsuGeometry --


def tauCategory(psi, phi, NN1,CO2,COX,NO2,NCACO2,NCANN1O2,NOX,NCACOX,NCANN1OX):
    '''
    These are what we are checking
    1. Is it hydrogen bonded to an oxygen 2 away?
        If so, is that oxygen planar
    2. Is there a heavy atom nearby?
        If so is that heavy atom planar
    3. Is psi 0 despite not getting through the others?
    4. Is phi 0?
    '''

    #Images
    # 1. Big O in the right corner, planar
    if NO2 == NOX:
        if NO2 < 3.2 and abs(NCACO2) < 15:# and abs(psi) < 5:
            return 'A1'
        if NO2 < 3.6 and abs(NCACO2) < 15:# and abs(psi) < 5:
            return 'A2'

    if NOX < 3.2 and abs (NCACOX) < 15:
        return 'A3'
    if NOX < 3.6 and abs (NCACOX) < 15:
        return 'A4'

    # 2. Big O in the corner, not quite planar (planar with N+1)
    if NO2 == NOX:
        if NO2 < 3.2 and abs(NCANN1O2) < 15:  # and abs(psi) < 5:
            return 'B1'
        if NO2 < 3.6 and abs(NCANN1O2) < 15:  # and abs(psi) < 5:
            return 'B2'

    if NOX < 3.2 and abs(NCANN1OX) < 15:
        return 'B3'
    if NOX < 3.6 and abs(NCANN1OX) < 15:
        return 'B4'

    # 3. N-1 is linear with N-C


    # 4. Big O on the left

    # 5. Big O in the bottom left


    return 'X'



def clusterEdTauMaker(pdbCode, rid,chain,aa):
    from PsuGeometry import GeoReport as psu
    from PsuGeometry import GeoPdb as geopdb

    pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/Reports/'

    georep = psu.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath, ed=False, dssp=False)
    pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, ed=False, dssp=False)
    apdb = pdbmanager.getPdb(pdbCode, True)

    pdbcsv = apdb.getDataFrame()
    queryC = 'rid==' + str(rid) + ' and chain=="' + chain + '"' + ' and atom=="CA"'
    queryL = 'rid==' + str(rid) + ' and chain=="' + chain + '"' + ' and atom=="N"'
    queryP = 'rid==' + str(rid) + ' and chain=="' + chain + '"' + ' and atom=="C"'
    dataC = pdbcsv.query(queryC)
    dataL = pdbcsv.query(queryL)
    dataP = pdbcsv.query(queryP)

    if len(dataC) > 0 and len(dataL) > 0 and len(dataP) > 0:
        cx = round(dataC['x'].values[0], 3)
        cy = round(dataC['y'].values[0], 3)
        cz = round(dataC['z'].values[0], 3)
        lx = round(dataL['x'].values[0], 3)
        ly = round(dataL['y'].values[0], 3)
        lz = round(dataL['z'].values[0], 3)
        px = round(dataP['x'].values[0], 3)
        py = round(dataP['y'].values[0], 3)
        pz = round(dataP['z'].values[0], 3)

        row = pdbCode + "," + chain + str(rid) + "," + str(cx) + "," + str(cy) + "," + str(cz)
        row += "," + str(lx) + "," + str(ly) + "," + str(lz)
        row += "," + str(px) + "," + str(py) + "," + str(pz)

        return row