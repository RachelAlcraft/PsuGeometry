# -- ©Rachel Alcraft 2021, PsuGeometry --
import time
import pandas as pd
from PsuGeometry import GeoReport as psu
'''
EH stats report coparison
'''
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/'
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/1.TauCorrelations/Data/BestSupportedCSVs/Reports/'

#EH_SET,aa,N:CA,N:CA_SD,CA:C,CA:C_SD,C:O,C:O_SD,C:N+1,C:N+1_SD,N:CA:C,N:CA:C_SD,CA:C:N+1,CA:C:N+1_SD,CA:C:O,CA:C:O_SD,O:C:N+1,O:C:N+1_SD,C-1:N:CA,C-1:N:CA_SD
EHFileName = 'Data1_EH.csv'
BestFileName = 'Set4_BestWithGeosALL.csv'

dataEH = pd.read_csv(loadPath + EHFileName)
dataBest = pd.read_csv(loadPath + BestFileName)

aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
geos = ['N:CA','CA:C','C:O','C:N+1','N:CA:C','CA:C:N+1','CA:C:O','O:C:N+1','C-1:N:CA']
#geos = ['N:CA','CA:C','C:O','C:N+1','N:CA:C']

#specifically looking at the mean and sd of the parameters in comparison to EH
print(dataEH)

georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)

for geo in geos:
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
        eh2001 = ehCut.query("EH_SET ==  '1991'")
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
            georep.addHistogram(data=bestCut, geoX='TAU_x', title=title)
        else:
            georep.addHistogram(data=bestCut, geoX=geo, title=title)

georep.printToHtml('Best Supported Engh&Huber Compare', 4, 'BS_EH')