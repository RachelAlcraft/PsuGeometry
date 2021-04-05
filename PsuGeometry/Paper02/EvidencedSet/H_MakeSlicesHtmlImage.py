# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''
###############################################################################################

def makeSlicesHtml(setName, title,titleType):
    FileDir = setName
    firstRow = 1
    fileName = titleType + '_' + setName


    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    edSlicePath = 'F:/Code/ProteinDataFiles/ccp4_out/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/EvidencedSet/SlicesH/'

    #We are going to load the data that has been created by density flight
    edSlicePath += FileDir + "/"
    print(edSlicePath + "_Results.csv")
    inputdata = pd.read_csv(edSlicePath + "_Results.csv")
    pdbs = inputdata['PdbCode'].values
    tags = inputdata['Tag'].values
    taus = inputdata['Angle'].values

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False)

    #origs = []
    #radiants = []

    for i in range(0,len(pdbs)):
        pdb = pdbs[i]
        tag = tags[i]
        tau = taus[i]
        sliceOrigVal = georep.loadSlice(edSlicePath + pdb + tag + "value_slice.csv")
        sliceOrigRad = georep.loadSlice(edSlicePath + pdb + tag + "radiant_slice.csv")
        georep.addSlice(sliceOrigVal, palette='cubehelix_r',title=pdb +' ' +  tag + ' value, angle=' + str(round(tau,3)))
        georep.addSlice(sliceOrigRad, palette='bone',title=pdb + ' ' + tag + ' radiant',Contour=False)
        #origs.append(sliceOrigVal)
        #radiants.append(sliceOrigRad)
        firstRow += 1

    #georep.addSlices(origs, palette='cubehelix_r', title='Average values', logged=False, centre=False)
    #georep.addSlices(radiants, palette='bone', title='Average radiant', logged=False, centre=False,Contour=False)
    georep.printToHtml(title,4,fileName)


###########
#makeSlicesHtml('NCACO_B02_C_O','C:O Outliers, Set B02 (CA:C:O)')
#makeSlicesHtml('NCACO_B02_CA_C','CA:C Outliers, Set B02 (CA:C:O)')