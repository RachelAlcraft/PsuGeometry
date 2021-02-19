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
#pdbList1000 = pdbList1000[:20]

geoList = ['N:N+1','CA-2:CA-1:CA:CA+1','TAU','PHI','PSI','CA-1:CA:CA+1:CA+2','N:O','CA-1:CA:CA+1','C-1:C','O-1:O','CA-2:CA:CA+2']
hueList = ['dssp','aa', 'rid', 'bfactor']
#aas = ['GLY','ALA']
#aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
includeDSSP = False
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
    dsspHue = 'aa'



###########################################################################################
georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)
data = georep.getGeoemtryCsv(geoList, hueList)
data = data.query('TAU > 100')
data = data.query('TAU < 125')


georep.addScatter(data=data, geoX='PHI', geoY='PSI', hue='aa', title='', palette='jet_r', sort='NON')
georep.addScatter(data=data, geoX='TAU', geoY='PSI', hue='aa', title='', palette='jet_r', sort='NON')
georep.addScatter(data=data, geoX='TAU', geoY='PHI', hue='aa', title='', palette='jet_r', sort='NON')
georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='aa', title='', palette='jet_r', sort='NON')
georep.addScatter(data=data, geoX='TAU', geoY='N:N+1', hue='aa', title='', palette='jet_r', sort='NON')
georep.addScatter(data=data, geoX='TAU', geoY='N:O', hue='aa', title='', palette='jet_r', sort='NON')

georep.addDifference(dataA=data, dataB=data, geoX='TAU', geoY='PSI', restrictionsA={'aa': 'SER'},restrictionsB={'aa': 'VAL'})

georep.addScatter(data=data, geoX='PHI', geoY='PSI', hue='aa', title='', palette='plasma', sort='RAND',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='PSI', hue='aa', title='', palette='plasma', sort='RAND',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='PHI', hue='aa', title='', palette='plasma', sort='RAND',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='PSI', geoY='N:N+1', hue='aa', title='', palette='plasma', sort='RAND',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='N:N+1', hue='aa', title='', palette='plasma', sort='RAND',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='N:O', hue='aa', title='', palette='plasma', sort='RAND',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='N:O', hue='aa', title='', palette='plasma', sort='ASC',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='N:O', hue='aa', title='', palette='plasma', sort='DESC',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='PSI', hue='aa', title='', palette='plasma', sort='ASC',restrictions={'aa':'SER,VAL'})
georep.addScatter(data=data, geoX='TAU', geoY='PSI', hue='aa', title='', palette='plasma', sort='DESC',restrictions={'aa':'SER,VAL'})


print('Creating reports')
georep.printToHtml('Tau Plots, Pdbs=' + str(len(pdbList1000)) , 3, 'Results8_aa_' + str(len(pdbList1000)))





