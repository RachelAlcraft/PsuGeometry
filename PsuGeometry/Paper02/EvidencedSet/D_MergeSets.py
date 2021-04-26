# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
import _Helpers as help
'''
Compare sets and EH and Jaskolski
'''

def mergeSets(sets,tag):
    import matplotlib.pyplot as plt
    plt.close('all')
    plt.clf()
    plt.cla()

    loadPathEH = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/Data/'
    loadPathCsv = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = help.rootPath + '/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataD/'

    dics = []
    geos = ['N:CA', 'CA:C', 'C:O', 'C:N+1', 'TAU', 'C-1:N:CA', 'CA:C:N+1', 'CA:C:O', 'O:C:N+1', 'CA:C:N+1']
    aas = ['ALL', 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN','ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

    # 1. First load the E&H and Jaskolski values for comparison
    ehSets = ['2001','1991','Jaskolski']
    ehAas = ['ALL','GLY','PRO']
    setEHFileName = 'Data_EH.csv'
    dataEH = pd.read_csv(loadPathEH + setEHFileName)
    for eh in ehSets:
        dataEHCut = dataEH.query('EH_SET == "' + eh + '"')
        for geo in geos:
            #if geo == 'TAU':
            #    geo = 'N:CA:C'
            for aa in ehAas:
                dataEHAACut = dataEHCut.query('aa == "' + aa + '"')
                print('Merging EH',eh,geo,aa)
                meanEH = round(dataEHAACut[geo].values[0], 3)
                if meanEH > 0:
                    sdEH = dataEHAACut[geo + '_SD'].values[0]
                    dic = {}
                    dic['aa'] = aa
                    dic['geo'] = geo
                    dic['set'] = 'EH_' + eh
                    dic['mean'] = meanEH
                    dic['sd'] = sdEH
                    dics.append(dic)
                    #print(dic)

    # 2. Then load all my data up as sumary statistics



    for geo in geos:
        for aa in aas:
            sql = 'aa == "' + aa + '"'
            for pdbSet in sets:
                print('Merging',geo,aa,pdbSet)
                dic = {}
                setFileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'
                dataBest = pd.read_csv(loadPathCsv + setFileName)
                if aa == 'ALL':
                    dataBest = dataBest.query('aa != "PRO"')
                    dataBest = dataBest.query('aa != "GLY"')
                else:
                    dataBest = dataBest.query('aa == "' + aa + '"')
                dataBestCut = dataBest[['aa',geo]]
                dataDescribed = dataBestCut.describe().values
                dic = {}
                dic['geo'] = geo
                dic['set'] = pdbSet
                dic['aa'] = aa
                dic['count'] = int(dataDescribed[0])
                mean = round(float(dataDescribed[1]),3)
                sd = round(float(dataDescribed[2]),3)
                # sd wants to be described as as relating to the 4th significant figure of the mean.
                if mean < 10:
                    sd = sd * 1000
                    sd = round(sd,0)
                else:
                    sd = sd * 10
                    sd = round(sd,0)
                print(geo, mean, sd)

                dic['mean'] = mean
                dic['sd'] = sd

                dic['min'] = round(float(dataDescribed[3]),3)
                dic['median'] = round(float(dataDescribed[5]),3)
                dic['max'] = round(float(dataDescribed[7]),3)
                dics.append(dic)
    dataFrame = pd.DataFrame.from_dict(dics)
    filePath = printPath + tag + 'Data_SetsSummaryMerged.csv'
    print('...printing', filePath)
    dataFrame.to_csv(filePath, index=False)
    print(dataFrame)



