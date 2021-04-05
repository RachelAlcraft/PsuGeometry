# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''
###############################################################################################
def createBadDensitySlices(pdbSet,atomCe,atomLi,atomPl):

    pdbOriginalPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/EvidencedSet/SlicesF/'

    # This gets the list of pdbs
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_Evidenced.csv')  # This is a list of pdbs <= 1.1A non homologous to 90%
    pdbListIn = pdbdata['PDB'].tolist()[0:]
    #pdbListIn = ['1p1x']
    for pdb in pdbListIn:
        slicesList = []
        fileName = (pdbDataPath + 'pdb' + pdb + '_' + atomCe + atomLi + atomPl + '.bad').lower()
        import os.path
        if os.path.isfile(fileName):
            text_file = open(fileName, "r")
            lines = text_file.read().split('\n')
            print(len(lines))
            text_file.close()
            for line in lines:
                atom = line[12:14].lstrip().rstrip()
                aa = line[14:20].lstrip().rstrip()
                chain = line[20:22].lstrip().rstrip()
                rid = line[22:27].lstrip().rstrip()
                print(pdb,atom,aa,chain,rid,line)
                if rid != '':
                    if [pdb,chain,rid] not in slicesList:
                        slicesList.append([pdb,chain,rid])

        bigstring = ""

        for sl in slicesList:
            print(sl)
            georep = psu.GeoReport([sl[0]],pdbOriginalPath,edDataPath,printPath,ed=False,dssp=False)
            pdbmanager = geopdb.GeoPdbs(pdbOriginalPath,edDataPath,ed=False,dssp=False)
            apdb = pdbmanager.getPdb(sl[0],True)
            pdbcsv = apdb.getDataFrame()
            queryC = 'rid==' + str(sl[2]) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomCe + '"'
            queryL = 'rid==' + str(sl[2]) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomLi + '"'
            queryP = 'rid==' + str(sl[2]) + ' and chain=="' + sl[1] + '"' + ' and atom=="' + atomPl + '"'
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

                row = sl[0] + "," + sl[1] + str(sl[2]) + "," + str(cx) + "," + str(cy) + "," + str(cz)
                row +=  "," + str(lx) + "," + str(ly) + "," + str(lz)
                row += "," + str(px) + "," + str(py) + "," + str(pz)

                print(row)

                bigstring += row + '\n'

        if len(slicesList) > 0:
            print("########RESULTS#########")
            print("")
            print(bigstring)

            f = open(printPath + 'BadSlice_' + pdbSet + '_' + pdb + '.txt', "w")
            f.write(bigstring)
            f.close()


