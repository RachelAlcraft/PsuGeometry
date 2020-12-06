# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
'''
This script looks at the distribution of N-CA for proline and glycine at different resolutions
Due to the large numbe fo structures it runs in modes load and save 
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList1000 = geol.GeoPdbLists().getListPaper()

#pdbListA = pdbList1000[:700]
#pdbListB = pdbList1000[500:]

geoList = ['CA-1:CA', 'CA-1:C-1:N:CA']
hueList = ['aa', 'resolution']

georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False)
dataA = georep.getGeoemtryCsv(geoList,hueList)
#georep.pdbCodes = pdbListB
#dataB = georep.getGeoemtryCsv(geoList,hueList)

#serialsie for posterity
#dataA.to_csv(printPath + 'Results8a.csv', index=False)
dataA.to_csv(printPath + 'Results8a.csv', index=False)

dataA = dataA[dataA['CA-1:CA'] < 10]


georep.addScatter(data=dataA,geoX='CA-1:CA',geoY='CA-1:C-1:N:CA',title='Omega and CA',hue='aa',palette='nipy_spectral',sort='ASC')
georep.addScatter(data=dataA,geoX='CA-1:CA',geoY='CA-1:C-1:N:CA',title='Omega and CAn',hue='aa',palette='nipy_spectral',sort='Desc')
georep.addScatter(data=dataA,geoX='CA-1:CA',geoY='CA-1:C-1:N:CA',title='Omega and CA',hue='resolution',palette='viridis_r')
#georep.addScatter(data=dataB,geoX='CA-1:CA',geoY='CA-1:C-1:N:CA',title='Lower resolution',hue='aa',palette='nipy_spectral',sort='ASC')
#georep.addScatter(data=dataB,geoX='CA-1:CA',geoY='CA-1:C-1:N:CA',title='Lower resolution',hue='resolution',palette='viridis_r')

georep.printToHtml('N-CA Distributions at different resolutions', 3, 'Results8_omega_ca')