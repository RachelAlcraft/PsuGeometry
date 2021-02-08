# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = False
pdbList1000 = geol.GeoPdbLists().getListPaper()
pdbList1000 = pdbList1000[:100]
pdbList1000.sort()
dihs = ['PSI','PHI']
angles = ['TAU']
distances = ['N:O']
hueList = ['dssp','aa', 'rid', 'bfactor']
hus = ['dssp']
###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
includeDSSP = True
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper01/'
    includeDSSP = False  # on my windows computer

###########################################################################################
georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)
geoList = []
for geo in dihs:
    geoList.append(geo)
for geo in distances:
    geoList.append(geo)
for geo in angles:
    geoList.append(geo)

for pdb in pdbList1000:
    georepPdb = psu.GeoReport([pdb], pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)
    data = georepPdb.getGeoemtryCsv(geoList, hueList)
    georep.addScatter(data=data, geoX='ridx', geoY='TAU', hue='PSI', title=pdb + ' TAU:PSI', palette='rainbow', sort='NON')
    for hu in hus:
        georep.addScatter(data=data, geoX='ridx', geoY='TAU', hue=hu, title=pdb + ' TAU:' + hu, palette='rainbow', sort='NON')
    georep.addScatter(data=data, geoX='ridx', geoY='PSI', hue='TAU', title=pdb + ' PSI:TAU', palette='rainbow', sort='NON')
    for hu in hus:
        georep.addScatter(data=data, geoX='ridx', geoY='PSI', hue=hu, title=pdb + ' PSI:' + hu, palette='rainbow', sort='NON')


print('Creating reports')
georep.printToHtml('Tau Chain Plots, Pdbs=' + str(len(pdbList1000)) , 4, 'Results3_tauchain_' + str(len(pdbList1000)))





