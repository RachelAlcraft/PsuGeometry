# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
import random
'''
Reference Biochemistry 4th Edition, Donald Voel and Judith G. Voel

'''
myWindowsLaptop = True
###############################################################################################

goodList = ['5xqp','6u66','5lun']

pdbList1000 = geol.GeoPdbLists().getListPaper()
random.shuffle(pdbList1000)
pdbList1000 = pdbList1000[:30]
#pdbList1000 = ['5zj8','5tnv','4wka','5p9v','1n4v','1zjy','6hmc']
pdbList1000 = goodList

dihs = ['PSI','PHI']
angles = ['TAU']
distances = ['N:O-4','N:O-3','N:O-2']
#distances = ['N:O-3']
aas = ['GLY','ALA']
hueList = ['aa', 'ridx', 'bfactor']


###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
dsspHue='dssp'
includeDSSP = False
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'
    includeDSSP = False  # on my windows computer

###########################################################################################

georepData = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=True)

geoList = []
for geo in dihs:
    geoList.append(geo)
for geo in distances:
    geoList.append(geo)
for geo in angles:
    geoList.append(geo)

count  = 0
length = len(pdbList1000)

for pdb in pdbList1000:
    print(pdb,' ', count,'/',length)
    count += 1
    # Create the dataframe
    georepPdb = psu.GeoReport([pdb], pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP,includePdbs=True)
    dataX = georepPdb.getGeoemtryCsv(geoList, hueList)
    # Clean the data
    #dataout = dataX.query('TAU <= 100 or TAU >=125')
    #dataX = dataX.query('TAU > 100')
    #dataX = data.query('TAU < 125')

    for geo in distances:
        print(geo)
        dataRed = dataX.copy()
        dataRed = dataRed[dataRed[geo] < 4]
        #georepData.addScatter(data=dataX, geoX='ridx', geoY='TAU', hue=geo, title=geo + ' < 3.5 ' + pdb, palette='jet',sort='NON')
        georepData.addScatter(data=dataX, geoX='ridx', geoY=geo, hue='TAU', title=geo + ' < 4 ' + pdb, palette='rainbow', sort='NON')
        georepData.addScatter(data=dataRed, geoX='ridx', geoY=geo, hue='TAU', title=geo + ' < 4 ' + pdb,palette='rainbow', sort='NON')
        georepData.addScatter(data=dataRed, geoX='ridx', geoY='TAU', hue=geo, title=geo + ' < 4 ' + pdb, palette='jet',sort='NON')


print('Creating reports')
georepData.printToHtml('Tau In A-Helices Plots, Pdbs=' + str(len(pdbList1000)) , 3, 'Results6_tauHelices_' + str(len(pdbList1000)))






