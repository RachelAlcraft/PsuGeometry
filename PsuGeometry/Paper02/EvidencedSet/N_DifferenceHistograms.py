# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
Compare the maxima for a synthetic structure to the real structure
'''

def diffHistograms(pdbList,tag):


    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = help.rootPath + '/ProteinDataFiles/ccp4_out/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataN/'

    allRealPdbs = []
    for pdb in pdbList:
        realFileName = pdb + '_DiffHistogram.csv'
        realData = pd.read_csv(loadPath + realFileName)
        allRealPdbs.append(realData)

    #append them all
    realCsv = pd.concat(allRealPdbs,axis=0,sort=False)

    georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

    georep.addHistogram(data=realCsv, geoX='Percent', title='Difference Histograms - %',count=True,hue='PdbCode')
    georep.addHistogram(data=realCsv, geoX='Diff', title='Difference Histograms - diff',count=True,hue='PdbCode',palette='LightSeaGreen')

    georep.addScatter(data=realCsv, geoX='Main', geoY='Percent', hue='Diff', palette='jet', categorical=False,sort='RANDOM',title='Differences')
    georep.addScatter(data=realCsv, geoX='Diff', geoY='Percent', hue='PdbCode', palette='jet_r', categorical=True, sort='RANDOM', title='Differences')

    georep.printToHtml('Difference Histograms', 2, 'DiffHist_' + tag)



##########
# RUN TO TEST
diffHistograms(['6eex','1ejg'],'1')

