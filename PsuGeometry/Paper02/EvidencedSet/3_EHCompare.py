# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
'''
EH stats report coparison
'''
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/EvidencedSet/Data/'
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/EvidencedSet/Reports/'

#EH_SET,aa,N:CA,N:CA_SD,CA:C,CA:C_SD,C:O,C:O_SD,C:N+1,C:N+1_SD,N:CA:C,N:CA:C_SD,CA:C:N+1,CA:C:N+1_SD,CA:C:O,CA:C:O_SD,O:C:N+1,O:C:N+1_SD,C-1:N:CA,C-1:N:CA_SD
EHFileName = 'Data_EH.csv'
BestFileName = 'Data_DefensibleWithGeosALL.csv'

dataEH = pd.read_csv(loadPath + EHFileName)
dataBest = pd.read_csv(loadPath + BestFileName)

aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
geos = ['N:CA','CA:C','C:O','C:N+1','N:CA:C','CA:C:N+1','CA:C:O','O:C:N+1','C-1:N:CA']
#geos = ['N:CA','CA:C','C:O','C:N+1','N:CA:C']

#specifically looking at the mean and sd of the parameters in comparison to EH
print(dataEH)

georepAA = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
georepSummary = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

for geo in geos:
    #Cut on ALL PRO and GLY
    compareSets = ['1991','2001']
    listCompares = []
    ehALL = dataEH.query("aa ==  'ALL'")
    ehGLY = dataEH.query("aa ==  'GLY'")
    ehPRO = dataEH.query("aa ==  'PRO'")

    bestALLCut = dataBest.query("aa !=  'GLY'")
    bestALLCut = bestALLCut.query("aa !=  'PRO'")
    bestGLYCut = dataBest.query("aa ==  'GLY'")
    bestPROCut = dataBest.query("aa ==  'PRO'")

    titleALL = ''
    titleGLY = ''
    titlePRO = ''
    for comp in compareSets:
        ehSetALL = ehALL.query("EH_SET ==  '" + comp + "'")
        ehSetGLY = ehGLY.query("EH_SET ==  '" + comp + "'")
        ehSetPRO = ehPRO.query("EH_SET ==  '" + comp + "'")
        meanALL = round(ehSetALL[geo].values[0], 3)
        sdALL = ehSetALL[geo + '_SD'].values[0]
        meanGLY = round(ehSetGLY[geo].values[0], 3)
        sdGLY = ehSetGLY[geo + '_SD'].values[0]
        meanPRO = round(ehSetPRO[geo].values[0], 3)
        sdPRO = ehSetPRO[geo + '_SD'].values[0]
        titleALL += 'ALL ' + geo + ' ' + comp + ' Mean=' + str(meanALL) + ' (' + str(sdALL) + ')\n'
        titleGLY += 'GLY ' + geo + ' ' + comp + ' Mean=' + str(meanGLY) + ' (' + str(sdGLY) + ')\n'
        titlePRO += 'PRO ' + geo + ' ' + comp + ' Mean=' + str(meanPRO) + ' (' + str(sdPRO) + ')\n'
        print(titleGLY,titlePRO,titleALL)
    if geo == 'N:CA:C':
        georepSummary.addHistogram(data=bestALLCut, geoX='TAU', title=titleALL)
        georepSummary.addHistogram(data=bestGLYCut, geoX='TAU', title=titleGLY)
        georepSummary.addHistogram(data=bestPROCut, geoX='TAU', title=titlePRO)
    else:
        georepSummary.addHistogram(data=bestALLCut, geoX=geo, title=titleALL)
        georepSummary.addHistogram(data=bestGLYCut, geoX=geo, title=titleGLY)
        georepSummary.addHistogram(data=bestPROCut, geoX=geo, title=titlePRO)


    for aa in aas:
        #prepare E&H comparison values
        mean1991 = ''
        mean2001 = ''
        sd1991 = ''
        sd2001 = ''
        useaa = 'ALL'
        if aa == 'PRO' or aa == 'GLY':
            useaa = aa
        ehCut = dataEH.query("aa ==  '" + useaa + "'")
        eh1991 = ehCut.query("EH_SET ==  '1991'")
        eh2001 = ehCut.query("EH_SET ==  '2001'")
        ehCutMean = dataEH[geo]
        ehCutSD = dataEH[geo + '_SD']
        mean1991 = round(eh1991[geo].values[0],3)
        mean2001 = round(eh2001[geo].values[0],3)
        sd1991 = eh1991[geo + '_SD'].values[0]
        sd2001 = eh2001[geo + '_SD'].values[0]
        title = geo + ' ' + aa + '\nEH 2001: mean=' + str(mean2001) + ' (' +str(sd2001) + ')\n'
        title = title + 'EH 1991: mean=' + str(mean1991) + ' (' + str(sd1991) + ')'
        print(aa,geo,mean1991,sd1991,mean2001,sd2001)
        bestCut = dataBest.query("aa ==  '" + aa + "'")
        if geo == 'N:CA:C':
            georepAA.addHistogram(data=bestCut, geoX='TAU', title=title)
        else:
            georepAA.addHistogram(data=bestCut, geoX=geo, title=title)

georepAA.printToHtml('Best Supported Engh&Huber Compare', 4, 'Defensible_AA_EH')
georepSummary.printToHtml('Best Supported Engh&Huber Compare', 3, 'Defensible_EH')