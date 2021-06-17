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

def makeSlicesHtml(setName, fileName, title,titleType,tag):
    import matplotlib.pyplot as plt
    plt.close('all')
    plt.clf()
    plt.cla()

    FileDir = setName
    firstRow = 1
    #fileName = titleType + '_' + setName


    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_data/'
    edDataPath = help.rootPath + '/ProteinDataFiles/ccp4_data/'
    edSlicePath = help.rootPath + '/ProteinDataFiles/ccp4_out/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesH/'

    #We are going to load the data that has been created by density flight
    edSlicePath += FileDir + "/"
    print(edSlicePath + "_Results.csv")
    inputdata = pd.read_csv(edSlicePath + "_Results.csv")
    pdbs = inputdata['PdbCode'].values
    tags = inputdata['Tag'].values
    taus = inputdata['Angle'].values

    valsPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesG/'
    print(valsPath)

    inputVals = pd.read_csv(valsPath + fileName)
    #inputVals = pd.read_csv(valsPath + "GoodOutliers_" + tag + ".csv")
    #print(inputVals)

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False)


    #origs = []
    #radiants = []

    for i in range(0,len(pdbs)):
        pdb = pdbs[i]
        tag = tags[i]
        tau = taus[i]

        pdbInputVals = inputVals.query("pdbCode == '" + pdb + "'")

        vals = pdbInputVals['value'].values
        vpdbs = pdbInputVals['pdbCode'].values
        vaas = pdbInputVals['aa'].values
        rids = pdbInputVals['rid'].values
        chains = pdbInputVals['chain'].values

        for j in range(0,len(vpdbs)):

            vpdb = vpdbs[j]
            newTag = vaas[j] + chains[j] + str(rids[j])

            print(vpdb,pdb,tag,newTag)

            if pdb == vpdb and newTag in tag:

                sliceOrigVal = georep.loadSlice(edSlicePath + pdb + tag + "value_slice.csv")
                sliceOrigRad = georep.loadSlice(edSlicePath + pdb + tag + "radiant_slice.csv")
                '''
                https://stackoverflow.com/questions/16400241/how-to-redefine-a-color-for-a-specific-value-in-a-matplotlib-colormap/16401183#16401183
                '''
                sliceOrigPos = georep.loadSlice(edSlicePath + pdb + tag + "poses_slice.csv")


                #This takes the data from the electron density only
                imtitle = pdb +' ' +  tag + ' value, angle=' + str(round(tau,3))
                # This takes the data from seperately created outliers file
                imtitle = pdb + ' ' + tag + ' value=' + str(round(vals[j], 3))
                #if pdb != vpdb:
                #    print(pdb,vpdb)
                #    imtitle = title = pdb +' ' +  tag + ' value, angle=' + str(round(tau,3)) + ' Error loading outlier'

                georep.addSlice(sliceOrigVal, palette='cubehelix_r', title=imtitle, YellowDots=sliceOrigPos,Contour=True)
                georep.addSlice(sliceOrigRad, palette='bone',title=pdb + ' ' + tag + ' radiant',Contour=False,YellowDots=sliceOrigPos)

                #origs.append(sliceOrigVal)
                #radiants.append(sliceOrigRad)
                firstRow += 1

    #georep.addSlices(origs, palette='cubehelix_r', title='Average values', logged=False, centre=False)
    #georep.addSlices(radiants, palette='bone', title='Average radiant', logged=False, centre=False,Contour=False)
    georep.printToHtml(title,4,fileName)


###########
#makeSlicesHtml('NCACO_B02_C_O','C:O Outliers, Set B02 (CA:C:O)')
#makeSlicesHtml('NCACO_B02_CA_C','CA:C Outliers, Set B02 (CA:C:O)')