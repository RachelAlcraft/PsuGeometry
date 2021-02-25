# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
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


###########################################################################################
georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=True)
pdbmanager = geopdb.GeoPdbs(georep.pdbDataPath, georep.edDataPath, georep.ed, georep.dssp)

data = georep.getGeoemtryCsv(geoList, hueList,bfactorFactor)
#data = data.query('TAU > 100')
#data = data.query('TAU < 125')
dataPsiRange = data.query('PSI > -50')
dataPsiRange = dataPsiRange.query('PSI < 50')


#for pd in pdbList:
#    apdb = pdbmanager.getPdb(pd, True)
#    atomData = apdb.getDataFrame()
#    fullFileName = printPath + 'Results10_' + pd + '.csv'
#    atomData.to_csv(fullFileName, index=False)

for aa in aas:

    sql = 'aa == "' + aa + '"'
    dataaa = data.query(sql)
    dataPsiRangeaa = dataPsiRange.query(sql)

    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='' + aa, palette='jet_r', sort='NON')
    georep.addScatter(data=dataPsiRangeaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='' + aa, palette='jet_r',sort='NON')

    georep.addScatter(data=dataaa, geoX='TAU', geoY='N:N+1', hue='bfactor', title='' + aa, palette='cubehelix_r', sort='RAND')
    georep.addScatter(data=dataPsiRangeaa, geoX='TAU', geoY='N:N+1', hue='bfactor', title='' + aa, palette='cubehelix_r', sort='RAND')

    georep.addScatter(data=dataaa, geoX='TAU', geoY='N:N+1', hue='rid', title='' + aa, palette='jet_r', sort='RAND',  categorical=False)
    georep.addScatter(data=dataPsiRangeaa, geoX='TAU', geoY='N:N+1', hue='rid', title='' + aa, palette='jet_r',sort='RAND', categorical=True)

    print('Creating reports')
    georep.printToHtml('Tau Plots', 2, 'Results10_' + aa + str(bfactorFactor))

    for pdb in pdbList:
        sql = 'pdbCode == "' + pdb + '"'
        datapdb = dataaa.query(sql)
        dataPsiRangepdb = dataPsiRangeaa.query(sql)

        fullFileName = printPath + 'Results10_' + pdb + aa + str(bfactorFactor) + '.csv'
        dataPsiRangeaa.to_csv(fullFileName, index=False)

        georep.addScatter(data=datapdb, geoX='PSI', geoY='N:N+1', hue='TAU', title='' + aa, palette='jet_r',sort='NON')
        georep.addScatter(data=dataPsiRangepdb, geoX='PSI', geoY='N:N+1', hue='TAU', title='' + aa, palette='jet_r', sort='NON')

        georep.addScatter(data=datapdb, geoX='TAU', geoY='N:N+1', hue='bfactor', title='' + aa, palette='cubehelix_r', sort='RAND')
        georep.addScatter(data=dataPsiRangepdb, geoX='TAU', geoY='N:N+1', hue='bfactor', title='' + aa, palette='cubehelix_r',sort='RAND')

        georep.addScatter(data=datapdb, geoX='TAU', geoY='N:N+1', hue='rid', title='' + aa, palette='jet_r', sort='RAND',categorical=False)
        georep.addScatter(data=dataPsiRangepdb, geoX='TAU', geoY='N:N+1', hue='rid', title='' + aa, palette='jet_r',sort='RAND', categorical=True)

        print('Creating reports')
        georep.printToHtml('Results 10 Tau Plots, Pdb=' + pdb, 2,'Results10_' + aa + pdb + str(bfactorFactor))










