# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = False
bfactorFactor = 1.3
keepDisordered = False
pdbList = ['6aiq','4m7g'] # disordered list
#pdbList = ['3goe','1mxt','1gvw','4g9s','6rhh','4a7u','6g1i','3wcq'] # ordered list
pdbList = ['4g9s'] # ordered list

geoList = ['N:N+1','TAU','PSI']
hueList = ['aa', 'rid', 'bfactor']
aas = ['GLY']
###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'


pdbdata = pd.read_csv('../PdbLists/Pdbs_1_1.csv') # This is a list of pdbs <= 1.1A non homologous to 90%
pdbList = pdbdata['PDB'].tolist()[0:]

print('Creating ordered report')
georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
print('Creating DISordered report')
georepDIS = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=True)

print('Getting dataframes of geometry -- 1')
dataA = georep.getGeoemtryCsv(geoList, hueList,bfactorFactor)
print('Getting dataframes of geometry -- 2')
dataB = georep.getGeoemtryCsv(geoList, hueList,-1)
print('Getting dataframes of geometry -- 3')
dataC = georepDIS.getGeoemtryCsv(geoList, hueList,bfactorFactor)
print('Getting dataframes of geometry -- 4')
dataD = georepDIS.getGeoemtryCsv(geoList, hueList,-1)

for aa in aas:
    print(aa)
    sql = 'aa == "' + aa + '"'
    dataAaa = dataA.query(sql)
    dataBaa = dataB.query(sql)
    dataCaa = dataC.query(sql)
    dataDaa = dataD.query(sql)

    georep.addScatter(data=dataAaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='Ordered factor=1.3 ' + aa, palette='jet', sort='NON')
    georep.addHistogram(data=dataAaa, geoX='TAU', title='Ordered factor=1.3 ' + aa)
    georep.addScatter(data=dataBaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='Ordered factor=any ' + aa, palette='jet', sort='NON')
    georep.addHistogram(data=dataBaa, geoX='TAU', title='Ordered factor=any ' + aa)

    georep.addScatter(data=dataCaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='Disordered factor=1.3 ' + aa, palette='jet', sort='NON')
    georep.addHistogram(data=dataCaa, geoX='TAU', title='Disordered factor=1.3 ' + aa)
    georep.addScatter(data=dataDaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='Disordered factor=any ' + aa, palette='jet', sort='NON')
    georep.addHistogram(data=dataDaa, geoX='TAU', title='Disordered factor=any ' + aa)

    print('Creating reports')
    #georep.printToHtml('Results 11. Tau Plots, Pdbs=' + str(len(pdbList)) , 2, 'Results11_' + aa + str(len(pdbList)))
    georep.printToHtml('Results 11. Tau Plots\nPdbs=' + str(len(pdbList)) + '\nWith bFactorFactor of ' + str(bfactorFactor), 2, 'Results11_' + aa)



