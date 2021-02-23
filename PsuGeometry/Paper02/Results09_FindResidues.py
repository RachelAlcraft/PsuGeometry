# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
import random
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = False
pdbList1000 = geol.GeoPdbLists().getListPaper()
#random.shuffle(pdbList1000)
pdbList1000 = pdbList1000[:100]

geoList = ['N:N+1','CA-2:CA-1:CA:CA+1','TAU','PHI','PSI','CA-1:CA:CA+1:CA+2','N:O','CA-1:CA:CA+1','C-1:C','O-1:O','CA-2:CA:CA+2']
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
georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)
data = georep.getGeoemtryCsv(geoList, hueList)
data = data.query('TAU > 100')
data = data.query('TAU < 125')
dataPsiRange = data.query('PSI > -5')
dataPsiRange = dataPsiRange.query('PSI < 5')

for aa in aas:
    sql = 'aa == "' + aa + '"'
    dataaa = data.query(sqPSI > -5 and l)
    dataPsiRangeaa = dataPsiRange.query(sql)

    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU' + aa, palette='jet', sort='NON')
    georep.addScatter(data=dataPsiRangeaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU' + aa, palette='jet',sort='NON')

    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='pdbCode', title='PSI|N:N+1|TAU' + aa, palette='jet_r',sort='NON')
    georep.addScatter(data=dataPsiRangeaa, geoX='PSI', geoY='N:N+1', hue='pdbCode', title='PSI|N:N+1|TAU' + aa, palette='jet_r',sort='NON')

    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='rid', title='PSI|N:N+1|TAU' + aa, palette='jet_r',sort='NON')
    georep.addScatter(data=dataPsiRangeaa, geoX='PSI', geoY='N:N+1', hue='rid', title='PSI|N:N+1|TAU' + aa,palette='jet_r', sort='NON')

    print('Creating reports')
    georep.printToHtml('Tau Plots, Pdbs=' + str(len(pdbList1000)) , 2, 'Results9_residues2_' + aa + str(len(pdbList1000)))





