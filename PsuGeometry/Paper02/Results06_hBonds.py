# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
import random
'''
Reference Biochemistry 4th Edition, Donald Voel and Judith G. Voel

'''
myWindowsLaptop = False
###############################################################################################
pdbList1000 = geol.GeoPdbLists().getListPaper()
random.shuffle(pdbList1000)
#pdbList1000 = pdbList1000[:300]
#pdbList1000 = ['5zj8']

dihs = ['PSI','PHI']
angles = ['TAU']
distances = ['N:O-4','N:O-3','N:O-2']
aas = ['GLY','ALA']
hueList = ['dssp','aa', 'rid', 'bfactor']


###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
dsspHue='dssp'
includeDSSP = True
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'
    includeDSSP = False  # on my windows computer
if not includeDSSP:
    hueList = ['aa', 'rid']
    dsspHue = 'aa'
###########################################################################################

georepData = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)

geoList = []
for geo in dihs:
    geoList.append(geo)
for geo in distances:
    geoList.append(geo)
for geo in angles:
    geoList.append(geo)

# Create the dataframe
dataX = georepData.getGeoemtryCsv(geoList, hueList)

#Clean the data
dataout = dataX.query('TAU <= 100 or TAU >=125')
data = dataX.query('TAU > 100')
data = data.query('TAU < 125')
dataRed = data.copy()

dataRed = dataRed[dataRed['N:O-2'] < 4]
dataRed = dataRed[dataRed['N:O-3'] < 4]
dataRed = dataRed[dataRed['N:O-4'] < 4]


for aa in aas:
    sql = 'aa == "' + aa + '"'
    dataaa = data.query(sql)

    for geo in distances:
        print(aa,geo)
        georepData.addHexBins(data=dataaa, geoX='TAU', geoY=geo, hue='count', title=geo + 'HB count ' + aa, palette='cubehelix_r', bins='log', gridsize=50)
        georepData.addHexBins(data=dataRed, geoX='TAU', geoY=geo, hue='count', title=geo + ' < 4 HB count ' + aa,palette='cubehelix_r', bins='log', gridsize=50)
        georepData.addScatter(data=dataRed, geoX='TAU', geoY=geo, hue=dsspHue, title=geo + ' < 4 HB dssp ' + aa, palette='tab10',sort='NON')

    print('Creating reports')
    georepData.printToHtml('Tau Correlations with HB ' + aa + ' Plots, Pdbs=' + str(len(pdbList1000)) , 3, 'Results6_' + aa + '_tauHB_' + str(len(pdbList1000)))





