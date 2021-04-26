# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
Compare the maxima for a synthetic structure to the real structure
'''

def maximaCompare(pdbSet,pdbList):


    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataM/'

    allRealPdbs = []
    allFakePdbs = []
    for pdb in pdbList:
        realFileName = '_ADJ/MaximaDifferences_' + pdb + '.csv'
        fakeFileName = '_FAKE_ADJ/MaximaDifferences_' + pdb + '.csv'

        realData = pd.read_csv(pdbDataPath + realFileName)
        #realData = realData.query("AtomType in ('C','N','CA','O')")
        #realData = realData.query("AA != 'HOH'")
        allRealPdbs.append(realData)

        fakeData = pd.read_csv(pdbDataPath + fakeFileName)
        allFakePdbs.append(fakeData)

    #append them all
    realCsv = pd.concat(allRealPdbs,axis=0,sort=False)
    fakeCsv = pd.concat(allFakePdbs, axis=0, sort=False)

    realOcc = realCsv.query("Occupancy == 1")
    fakeOcc = fakeCsv.query("Occupancy == 1")

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

    georep.addHistogram(data=realCsv, geoX='Difference', title='Real PDB structures',count=True,hue='pdbCode')
    georep.addHistogram(data=fakeCsv, geoX='Difference', title='Fake PDB structures',count=True,hue='pdbCode')

    georep.addHistogram(data=realOcc, geoX='Difference', title='Real PDB structures, Occ=1', count=True, hue='pdbCode')
    georep.addHistogram(data=fakeOcc, geoX='Difference', title='Fake PDB structures, Occ=1', count=True, hue='pdbCode')

    georep.addScatter(data=realCsv, geoX='Difference', geoY='Width', hue='AA', palette='tab20', categorical=True,sort='NON',title='Real')
    georep.addScatter(data=fakeCsv, geoX='Difference', geoY='Width', hue='AA', palette='tab20', categorical=True, sort='NON',title='Fake')

    georep.addScatter(data=realCsv, geoX='Difference', geoY='Occupancy', hue='Width', palette='cubehelix_r',sort='NON', title='Real')
    georep.addScatter(data=fakeCsv, geoX='Difference', geoY='Occupancy',  hue='Width', palette='cubehelix_r',sort='NON', title='Fake')

    #georep.addHistogram(data=realCsv, geoX='Difference', title='Real PDB structures', count=True,splitKey='pdbCode')
    for pdb in pdbList:
        print(pdb)
        onerealData = realData.query("pdbCode == '" + pdb + "'")
        onefakeData = fakeData.query("pdbCode == '" + pdb + "'")
        georep.addHistogram(data=realOcc, geoX='Difference', title='Real Density, Occ=1', hue='ResNo',count=True, restrictions={'pdbCode':pdb})
        georep.addHistogram(data=fakeOcc, geoX='Difference', title='Fake Density, Occ=1', hue='ResNo',count=True,restrictions={'pdbCode': pdb},palette='LightSeaGreen')

    georep.printToHtml('Maxima differences, set=' + pdbSet, 4, 'Maxima_' + pdbSet)



##########
# RUN TO TEST
#maximaCompare('CUBIC',help.getList('SMALLEST',0))

maximaCompare('LINEAR',help.getList('EVIDENCED',15))
maximaCompare('CUBIC',help.getList('EVIDENCED',15))
maximaCompare('QUINTIC',help.getList('EVIDENCED',15))
maximaCompare('HEPTIC',help.getList('EVIDENCED',15))

#maximaCompare('CUBIC',['2hs1','2vb1','3k34','5jdt','6jgj','2ce2','1oew'])
#maximaCompare('QUINTIC',['2hs1','2vb1','3k34','5jdt','6jgj','2ce2','1oew'])
#maximaCompare('QUINTIC',help.getList('SMALLEST'))