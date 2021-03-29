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

    #Is the )-2 oxygen the closest oxygen there is?
    if NO2 == NOX:
        if NO2 < 3.6:
            return 'A1'

    #So now check if the nearest O is hydrogen bonded
    if NOX < 3.6:
        if abs(NCACOX) < 20:
            if COX < 4.5:  # Then it is right hand image
               return 'C1'
            else:
                return 'C2'
        elif abs(NCANN1OX) < 20:
            return 'D1'

        if NO2 < NOX:
            return 'E1'
        if NOX < 3.6:
            return 'E2'
        if abs(NCACO2) < 10:
            return 'E3'
        if abs(NCANN1O2) < 10:
            return 'E4'
        if abs(NCACOX) < 10:
            return 'E5'
        if abs(NCANN1OX) < 10:
            return 'E6'

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
