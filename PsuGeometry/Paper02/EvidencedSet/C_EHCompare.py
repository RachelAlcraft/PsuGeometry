# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
EH stats report coparison
'''

def EHCompare(pdbSet):
    import matplotlib.pyplot as plt
    plt.close('all')
    plt.clf()
    plt.cla()

    pdbDataPath = help.rootPath + '/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPathEH = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/Data/'
    loadPathCsv = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataC/'

    #EH_SET,aa,N:CA,N:CA_SD,CA:C,CA:C_SD,C:O,C:O_SD,C:N+1,C:N+1_SD,N:CA:C,N:CA:C_SD,CA:C:N+1,CA:C:N+1_SD,CA:C:O,CA:C:O_SD,O:C:N+1,O:C:N+1_SD,C-1:N:CA,C-1:N:CA_SD
    EHFileName = 'Data_EH.csv'
    BestFileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'

    dataEH = pd.read_csv(loadPathEH + EHFileName)
    dataBest = pd.read_csv(loadPathCsv + BestFileName)

    aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    geos = ['N:CA','CA:C','C:O','C:N+1','TAU','CA:C:N+1','CA:C:O','O:C:N+1','C-1:N:CA']
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
        ehCIS = dataEH.query("aa ==  'CIS'")

        bestALLCut = dataBest.query("aa !=  'GLY'")
        bestALLCut = bestALLCut.query("aa !=  'PRO'")
        bestGLYCut = dataBest.query("aa ==  'GLY'")
        bestPROCut = dataBest.query("aa ==  'PRO'")
        print(pdbSet,bestPROCut)
        bestPROCut['ABSOMEGA'] = abs(bestPROCut['CA-1:C-1:N:CA'])
        bestPROCis = bestPROCut.query("ABSOMEGA < 120")
        bestPROTrans = bestPROCut.query("ABSOMEGA >= 120")

        titleALL = ''
        titleGLY = ''
        titlePRO = ''
        titleCIS = ''
        for comp in compareSets:
            ehSetALL = ehALL.query("EH_SET ==  '" + comp + "'")
            ehSetGLY = ehGLY.query("EH_SET ==  '" + comp + "'")
            ehSetPRO = ehPRO.query("EH_SET ==  '" + comp + "'")
            ehSetCIS = ehCIS.query("EH_SET ==  '" + comp + "'")
            meanALL = round(ehSetALL[geo].values[0], 3)
            sdALL = ehSetALL[geo + '_SD'].values[0]
            meanGLY = round(ehSetGLY[geo].values[0], 3)
            sdGLY = ehSetGLY[geo + '_SD'].values[0]
            meanPRO = round(ehSetPRO[geo].values[0], 3)
            sdPRO = ehSetPRO[geo + '_SD'].values[0]
            meanCIS = round(ehSetCIS[geo].values[0], 3)
            sdCIS = ehSetCIS[geo + '_SD'].values[0]
            titleALL += 'ALL ' + geo + ' ' + comp + ' Mean=' + str(meanALL) + ' (' + str(sdALL) + ')\n'
            titleGLY += 'GLY ' + geo + ' ' + comp + ' Mean=' + str(meanGLY) + ' (' + str(sdGLY) + ')\n'
            titlePRO += 'PRO ' + geo + ' ' + comp + ' Mean=' + str(meanPRO) + ' (' + str(sdPRO) + ')\n'
            titleCIS += 'CIS ' + geo + ' ' + comp + ' Mean=' + str(meanCIS) + ' (' + str(sdCIS) + ')\n'
            print(titleGLY,titlePRO,titleALL,titleCIS)
        if geo == 'N:CA:C':
            georepSummary.addHistogram(data=bestALLCut, geoX='TAU', title=titleALL)
            georepSummary.addHistogram(data=bestGLYCut, geoX='TAU', title=titleGLY)
            georepSummary.addHistogram(data=bestPROTrans, geoX='TAU', title=titlePRO)
            georepSummary.addHistogram(data=bestPROCis, geoX='TAU', title=titleCIS)
        else:
            georepSummary.addHistogram(data=bestALLCut, geoX=geo, title=titleALL)
            georepSummary.addHistogram(data=bestGLYCut, geoX=geo, title=titleGLY)
            georepSummary.addHistogram(data=bestPROTrans, geoX=geo, title=titlePRO)
            georepSummary.addHistogram(data=bestPROCis, geoX=geo, title=titleCIS)

        '''
        for aa in aas:
            #prepare E&H comparison values
            useaa = 'ALL'
            if aa == 'PRO' or aa == 'GLY':
                useaa = aa
            ehCut = dataEH.query("aa ==  '" + useaa + "'")
            eh1991 = ehCut.query("EH_SET ==  '1991'")
            eh2001 = ehCut.query("EH_SET ==  '2001'")
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
        '''

    georepSummary.printToHtml('Best Supported Engh&Huber Compare, set=' + pdbSet, 4, 'Defensible_EH_' + pdbSet)
    #georepAA.printToHtml('Best Supported Engh&Huber Compare, set=' + pdbSet, 4, 'Defensible_EH_AA_' + pdbSet)