# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = True
bfactorFactor = -1
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


pdbdata = pd.read_csv('../PdbLists/Pdbs_1_1.csv')
pdbList = pdbdata['PDB'].tolist()[1:50]

georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)
data = georep.getGeoemtryCsv(geoList, hueList)

for aa in aas:
    sql = 'aa == "' + aa + '"'
    dataaa = data.query(sql)
    print(dataaa)
    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU' + aa, palette='jet', sort='NON')
    georep.addHistogram(data=dataaa, geoX='TAU', title='Tau Histogram ' + aa)



    print('Creating reports')
    #georep.printToHtml('Results 11. Tau Plots, Pdbs=' + str(len(pdbList)) , 2, 'Results11_' + aa + str(len(pdbList)))
    georep.printToHtml('Results 11. Tau Plots, Pdbs=' + str(len(pdbList)), 2, 'Results11_' + aa)



