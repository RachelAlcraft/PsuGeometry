# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = False
pdbList1000 = geol.GeoPdbLists().getListPaper()
pdbList1000 = pdbList1000[:200]
geoList = ['N:N+1','CA-2:CA-1:CA:CA+1','TAU','PHI','PSI','CA-1:CA:CA+1:CA+2','N:O','CA-1:CA:CA+1']
hueList = ['aa', 'rid', 'bfactor']
aas = ['GLY','ALA']
includeDSSP = False
###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper01/'
    includeDSSP = False  # on my windows computer

###########################################################################################
georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)
data = georep.getGeoemtryCsv(geoList, hueList)
data = data.query('TAU > 100')
data = data.query('TAU < 125')

for aa in aas:
    sql = 'aa == "' + aa + '"'
    dataaa = data.query(sql)

    georep.addScatter(data=dataaa, geoX='PHI', geoY='PSI', hue='CA-2:CA-1:CA:CA+1', title='PHI|PSI|CA-2:CA-1:CA:CA+1 ',palette='jet', sort='NON')
    georep.addScatter(data=dataaa, geoX='N:N+1', geoY='N:O', hue='TAU', title='N:N+1|N:O|TAU', palette='jet', sort='NON')
    georep.addScatter(data=dataaa, geoX='N:N+1', geoY='CA-2:CA-1:CA:CA+1', hue='TAU', title='N:N+1|CA-2:CA-1:CA:CA+1|TAU', palette='jet', sort='NON')
    georep.addScatter(data=dataaa, geoX='N:N+1', geoY='CA-1:CA:CA+1:CA+2', hue='TAU', title='N:N+1|CA-1:CA:CA+1:CA+2|TAU ', palette='jet', sort='NON')
    georep.addScatter(data=dataaa, geoX='CA-1:CA:CA+1:CA+2', geoY='CA-1:CA:CA+1', hue='TAU',title='CA-1:CA:CA+1:CA+2|CA-1:CA:CA+1|TAU ', palette='jet', sort='NON')

    print('Creating reports')
    georep.printToHtml('Tau Plots, Pdbs=' + str(len(pdbList1000)) , 2, 'Results4_specify_' + aa + str(len(pdbList1000)))





