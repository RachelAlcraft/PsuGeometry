# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
import random
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = True
keepDisordered = True
bfactorFactor = -1
pdbList1000 = geol.GeoPdbLists().getListPaper()
#random.shuffle(pdbList1000)
pdbList1000 = pdbList1000[:400]

geoList = ['N:N+1','TAU','PHI','PSI','N:O','OMEGA','C-1:N:CA','CA:C:N+1']
hueList = ['dssp','aa', 'rid', 'bfactor']
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
georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, keepDisordered=keepDisordered,includePdbs=False)
data = georep.getGeoemtryCsv(geoList, hueList,bfactorFactor)
#data = data.query('TAU > 100')
#data = data.query('TAU < 125')
dataPsiRange = data.query('PSI > -50')
dataPsiRange = dataPsiRange.query('PSI < 50')

for aa in aas:
    sql = 'aa == "' + aa + '"'
    dataaa = data.query(sql)
    dataPsiRangeaa = dataPsiRange.query(sql)

    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU' + aa, palette='jet', sort='NON')
    georep.addScatter(data=dataPsiRangeaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU' + aa, palette='jet',sort='NON')

    georep.addHexBins(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU' + aa, palette='jet',bins='log', gridsize=50)
    georep.addHexBins(data=dataPsiRangeaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU' + aa,palette='jet', bins='log', gridsize=50)

    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='rid', title='PSI|N:N+1|TAU' + aa, palette='jet_r',sort='NON')
    georep.addScatter(data=dataPsiRangeaa, geoX='PSI', geoY='N:N+1', hue='rid', title='PSI|N:N+1|TAU' + aa,palette='jet_r', sort='NON')

    georep.addScatter(data=dataaa, geoX='TAU', geoY='N:N+1', hue='bfactor', title='' + aa, palette='cubehelix_r', sort='RAND')
    georep.addScatter(data=dataPsiRangeaa, geoX='TAU', geoY='N:N+1', hue='bfactor', title='' + aa,palette='cubehelix_r', sort='RAND')

    append = 'ordered_' + str(bfactorFactor) + '_'
    if keepDisordered:
        append = 'disordered_' + str(bfactorFactor) + '_'
    print('Creating reports')
    georep.printToHtml('Results 9 Tau Plots, Pdbs=' + str(len(pdbList1000)) , 2, 'Results9_residues_' + aa  + append + str(len(pdbList1000)))

    fullFileName = printPath + 'Results9_Data_' + str(bfactorFactor) + '_'+ aa + '.csv'
    dataPsiRangeaa.to_csv(fullFileName, index=False)





