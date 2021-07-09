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

def makeSlicesHtml(setName, fileNameOrig,fileNameAdj, title,titleType,tag):
    import matplotlib.pyplot as plt
    plt.close('all')
    plt.clf()
    plt.cla()

    FileDir = setName
    firstRow = 1
    outfileName = titleType + '_' + setName + "_" + tag


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

    valsPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesC/'
    print(valsPath)

    inputValsOrig = pd.read_csv(valsPath + fileNameOrig)
    inputValsAdj = pd.read_csv(valsPath + fileNameAdj)
    #inputVals = pd.read_csv(valsPath + "GoodOutliers_" + tag + ".csv")
    #print(inputVals)

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False)


    #origs = []
    #radiants = []

    for i in range(0,len(pdbs)):
        pdb = pdbs[i]
        tag = tags[i]
        tau = taus[i]

        pdbInputValsOrig = inputValsOrig.query("pdbCode == '" + pdb + "'")
        pdbInputValsAdj = inputValsAdj.query("pdbCode == '" + pdb + "'")

        valsOrig = pdbInputValsOrig['value'].values
        valsAdj = pdbInputValsAdj['value'].values
        vpdbs = pdbInputValsAdj['pdbCode'].values
        vaas = pdbInputValsAdj['aa'].values
        rids = pdbInputValsAdj['rid'].values
        chains = pdbInputValsAdj['chain'].values

        for j in range(0,len(vpdbs)):

            vpdb = vpdbs[j]
            newTag = vaas[j] + chains[j] + str(rids[j])

            print(vpdb,pdb,tag,newTag)

            if pdb == vpdb and newTag in tag:

                sliceOrigVal = georep.loadSlice(edSlicePath + pdb + tag + "value_slice.csv")
                sliceOrigRad = georep.loadSlice(edSlicePath + pdb + tag + "radiant_slice.csv")
                sliceOrigMag = georep.loadSlice(edSlicePath + pdb + tag + "magnitude_slice.csv")
                '''
                https://stackoverflow.com/questions/16400241/how-to-redefine-a-color-for-a-specific-value-in-a-matplotlib-colormap/16401183#16401183
                '''
                sliceOrigPos = georep.loadSlice(edSlicePath + pdb + tag + "poses_slice.csv")


                #This takes the data from the electron density only
                imtitle = pdb +' ' +  tag + ' value, angle=' + str(round(tau,3))
                # This takes the data from seperately created outliers file
                if '_O' in tag:
                    imtitle = pdb + ' ' + tag + ' value=' + str(round(valsOrig[j], 3))
                else:
                    imtitle = pdb + ' ' + tag + ' value=' + str(round(valsAdj[j], 3))
                #if pdb != vpdb:
                #    print(pdb,vpdb)
                #    imtitle = title = pdb +' ' +  tag + ' value, angle=' + str(round(tau,3)) + ' Error loading outlier'

                georep.addSlice(sliceOrigVal, palette='cubehelix_r', title=imtitle, YellowDots=sliceOrigPos,Contour=True)
                georep.addSlice(sliceOrigRad, palette='bone',title=pdb + ' ' + tag + ' radiant',Contour=False,YellowDots=sliceOrigPos)
                georep.addSlice(sliceOrigMag, palette='bone', title=pdb + ' ' + tag + ' magnitude', Contour=True,YellowDots=sliceOrigPos)

                #origs.append(sliceOrigVal)
                #radiants.append(sliceOrigRad)
                firstRow += 1

    #georep.addSlices(origs, palette='cubehelix_r', title='Average values', logged=False, centre=False)
    #georep.addSlices(radiants, palette='bone', title='Average radiant', logged=False, centre=False,Contour=False)
    georep.printToHtml(title,6,outfileName)


def makeSlicesHtmlFromValues(mainTitle, dirName, lineRuns,row,tag):
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/SlicesI/'
    georep = psu.GeoReport([], "", "", printPath, ed=False, dssp=False)
    #We are going to load the data that has been created by density flight and pasted into a text file
    for lineRun in lineRuns:
        title = lineRun[0]
        palette = lineRun[1]
        fileName = lineRun[2]
        posName = lineRun[3]
        inputVals = georep.loadSlice(dirName + fileName)
        inputPoses = georep.loadSlice(dirName + posName)
        georep.addSlice(inputVals, palette=palette, title=title, YellowDots=inputPoses,Contour=True)

    georep.printToHtml(mainTitle,row,tag)


