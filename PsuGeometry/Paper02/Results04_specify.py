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
random.shuffle(pdbList1000)
#pdbList1000 = pdbList1000[:200]

geoList = ['N:N+1','CA-2:CA-1:CA:CA+1','TAU','PHI','PSI','CA-1:CA:CA+1:CA+2','N:O','CA-1:CA:CA+1','C-1:C','O-1:O','CA-2:CA:CA+2']
hueList = ['dssp','aa', 'rid', 'bfactor']
aas = ['GLY','ALA']
includeDSSP = True
dsspHue = 'dssp'
###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'
    includeDSSP = False  # on my windows computer
    hueList = hueList[1:]
    dsspHue = 'aa'



###########################################################################################
georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)
data = georep.getGeoemtryCsv(geoList, hueList)
data = data.query('TAU > 100')
data = data.query('TAU < 125')

dataPhiPlus = data.query('PHI > 0')
dataPhiMinus = data.query('PHI <= 0')

dataPhiMinus['PHIABS'] = dataPhiMinus['PHI']*-1
dataPhiMinus['PSINEG'] = dataPhiMinus['PSI']*-1

for aa in aas:
    sql = 'aa == "' + aa + '"'
    dataaa = data.query(sql)
    dataPhiMinusaa = dataPhiMinus.query(sql)
    dataPhiPlusaa = dataPhiPlus.query(sql)

    #Ramachandran on average tau and dssp
    georep.addHexBins(data=dataaa, geoX='PHI', geoY='PSI', hue='TAU', title='Ramachandran with Average tau ' + aa, palette='jet',bins='log', gridsize=50)
    georep.addScatter(data=dataaa, geoX='PHI', geoY='PSI', hue=dsspHue, title='Ramachandran with dssp ' + aa, palette='tab10', sort='NON')

    georep.addHexBins(data=dataPhiPlusaa, geoX='PHI', geoY='PSI', hue='TAU', title='Rama PSI/+ve PHI with avg tau ' + aa, palette='jet', bins='log', gridsize=50)
    georep.addHexBins(data=dataPhiMinusaa, geoX='PHIABS', geoY='PSINEG', hue='TAU',title='Rama -ve PSI/abs(-ve PHI) abs with avg tau ' + aa, palette='jet', bins='log',gridsize=50)

    georep.addScatter(data=dataPhiPlusaa, geoX='PHI', geoY='PSI', hue='TAU', title='Rama PSI/+ve PHI with tau ' + aa, palette='jet',vmin=100,vmax=125)
    georep.addScatter(data=dataPhiMinusaa, geoX='PHIABS', geoY='PSINEG', hue='TAU', title='Rama -ve PSI/abs(-ve PHI) abs with tau ' + aa, palette='jet',vmin=100,vmax=125)



    #Tau vs PHI AND PSI
    georep.addHexBins(data=dataaa, geoX='TAU', geoY='PSI', hue='count', title='TAU vs PSI Density Plot ' + aa,palette='cubehelix_r', bins='log', gridsize=50)
    georep.addScatter(data=dataaa, geoX='TAU', geoY='PSI', hue=dsspHue, title='TAU vs PSI with dssp ' + aa,  palette='tab10', sort='NON')

    georep.addHexBins(data=dataaa, geoX='TAU', geoY='PHI', hue='count', title='TAU vs PHI Density Plot ' + aa, palette='cubehelix_r', bins='log', gridsize=50)
    georep.addScatter(data=dataaa, geoX='TAU', geoY='PHI', hue=dsspHue, title='TAU vs PHI with dssp ' + aa, palette='tab10', sort='NON')

    georep.addHexBins(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI vs N:N+1 with Average Tau ' + aa, palette='jet', bins='log', gridsize=50)
    georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue=dsspHue, title='PSI vs N:N+1 with dssp ' + aa, palette='tab10', sort='NON')

    georep.addHexBins(data=dataaa, geoX='PHI', geoY='C-1:C', hue='TAU', title='PHI vs C-1:C with Average Tau ' + aa,palette='jet', bins='log', gridsize=50)
    georep.addScatter(data=dataaa, geoX='PHI', geoY='C-1:C', hue=dsspHue, title='PHI vs C-1:C with dssp ' + aa, palette='tab10', sort='NON')

    georep.addHexBins(data=dataaa, geoX='PHI', geoY='O-1:O', hue='TAU', title='PHI vs O-1:O with Average Tau ' + aa,palette='jet', bins='log', gridsize=50)
    georep.addHexBins(data=dataaa, geoX='PHI', geoY='O-1:O', hue='count', title='PHI vs O-1:O with Density Plot ' + aa,palette='cubehelix_r', bins='log', gridsize=50)

    georep.addHexBins(data=dataaa, geoX='CA-1:CA:CA+1', geoY='CA-1:CA:CA+1:CA+2', hue='TAU', title='Kleywegt/Lyons CAlpha Plot with Average Tau ' + aa,palette='jet', bins='log', gridsize=50)
    georep.addHexBins(data=dataaa, geoX='CA-1:CA:CA+1', geoY='CA-1:CA:CA+1:CA+2', hue='count',title='Kleywegt/Lyons CAlpha Density Plot ' + aa, palette='cubehelix_r', bins='log', gridsize=50)

    georep.addScatter(data=dataaa, geoX='CA-1:CA:CA+1', geoY='CA-1:CA:CA+1:CA+2', hue=dsspHue, title='Kleywegt/Lyons CAlpha Plot with dssp ' + aa,palette='tab10', sort='NON')
    georep.addScatter(data=dataaa, geoX='CA-2:CA:CA+2', geoY='CA-2:CA-1:CA:CA+1', hue=dsspHue, title='CAlpha +- 2 Plot with dssp ' + aa,palette='tab10', sort='NON')

    georep.addHexBins(data=dataaa, geoX='CA-2:CA:CA+2', geoY='CA-2:CA-1:CA:CA+1', hue='TAU', title='CAlpha +- 2 Plot with Average Tau ' + aa, palette='jet', bins='log', gridsize=50)
    georep.addHexBins(data=dataaa, geoX='CA-2:CA:CA+2', geoY='CA-2:CA-1:CA:CA+1', hue='count',title='CAlpha +- 2  Density Plot ' + aa, palette='cubehelix_r', bins='log', gridsize=50)





    #georep.addHexBins(data=dataaa, geoX='TAU', geoY='CA-1:CA:CA+1:CA+2', hue='PSI', title='', palette='jet', bins='log',gridsize=50)
    #georep.addHexBins(data=dataaa, geoX='TAU', geoY='CA-1:CA:CA+1:CA+2', hue='PHI', title='', palette='jet',bins='log', gridsize=50)
    #georep.addHexBins(data=dataaa, geoX='CA-1:CA:CA+1', geoY='CA-1:CA:CA+1:CA+2', hue='TAU', title='', palette='jet',bins='log', gridsize=50)

    #georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU', palette='jet',sort='NON')
    #georep.addScatter(data=dataaa, geoX='PSI', geoY='N:N+1', hue='dssp', title='PSI|N:N+1|dssp', palette='jet',sort='NON')
    #georep.addHexBins(data=dataaa, geoX='PSI', geoY='N:N+1', hue='TAU', title='PSI|N:N+1|TAU', palette='jet',bins='log', gridsize=50)

    #georep.addScatter(data=dataaa, geoX='PHI', geoY='C-1:C', hue='TAU', title='PHI|C-1:C|TAU', palette='jet', sort='NON')
    #georep.addScatter(data=dataaa, geoX='PHI', geoY='C-1:C', hue='dssp', title='PHI|C-1:C|dssp', palette='jet',sort='NON')
    #georep.addHexBins(data=dataaa, geoX='PHI', geoY='C-1:C', hue='TAU', title='PHI|C-1:C|TAU', palette='jet',bins='log', gridsize=50)

    #georep.addScatter(data=dataaa, geoX='O-1:O', geoY='C-1:C', hue='TAU', title='O-1:O|C-1:C|TAU', palette='jet',sort='NON')
    #georep.addScatter(data=dataaa, geoX='N:N+1', geoY='C-1:C', hue='TAU', title='N:N+1|C-1:C|TAU', palette='jet',sort='NON')

    #georep.addScatter(data=dataaa, geoX='PHI', geoY='PSI', hue='CA-2:CA-1:CA:CA+1', title='PHI|PSI|CA-2:CA-1:CA:CA+1 ',palette='jet', sort='NON')
    #georep.addScatter(data=dataaa, geoX='N:N+1', geoY='N:O', hue='TAU', title='N:N+1|N:O|TAU', palette='jet', sort='NON')
    #georep.addScatter(data=dataaa, geoX='N:N+1', geoY='CA-2:CA-1:CA:CA+1', hue='TAU', title='N:N+1|CA-2:CA-1:CA:CA+1|TAU', palette='jet', sort='NON')
    #georep.addScatter(data=dataaa, geoX='N:N+1', geoY='CA-1:CA:CA+1:CA+2', hue='TAU', title='N:N+1|CA-1:CA:CA+1:CA+2|TAU ', palette='jet', sort='NON')
    #georep.addScatter(data=dataaa, geoX='CA-1:CA:CA+1:CA+2', geoY='CA-1:CA:CA+1', hue='TAU',title='CA-1:CA:CA+1:CA+2|CA-1:CA:CA+1|TAU ', palette='jet', sort='NON')

    print('Creating reports')
    georep.printToHtml('Tau Plots, Pdbs=' + str(len(pdbList1000)) , 2, 'Results4_specify_' + aa + str(len(pdbList1000)))





