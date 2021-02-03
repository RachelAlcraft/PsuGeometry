# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdbLists as geol
'''
This script looks at distributions aournf the carbonyl atom 
'''
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/'

pdbList = geol.GeoPdbLists().getList100()
georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath,ed=False,dssp=False,includePdbs=False)
#geoList = ['CA:C','C:O','C:N+1','CA:N+1','O:N+1','CA:O']
#hueList = ['dssp','aa','bfactor','resolution']
#data = georep.getGeoemtryCsv(geoList,hueList)

#bond lengths
georep.addHistogram(geoX='CA:C', title='CA:C')
georep.addHistogram(geoX='C:O', title='C:O')
georep.addHistogram(geoX='C:N+1', title='C:N+1')
#distances
georep.addHistogram(geoX='CA:N+1', title='CA:N+1')
georep.addHistogram(geoX='O:N+1', title='O:N+1')
georep.addHistogram(geoX='CA:O', title='CA:O')
#angles
georep.addHistogram(geoX='CA:C:N+1', title='CA:C:N+1')
georep.addHistogram(geoX='N+1:C:O', title='N+1:C:O')
georep.addHistogram(geoX='O:C:CA', title='O:C:CA')
#dihedrals
georep.addHistogram(geoX='CA:C:N+1:O', title='CA:C:N+1:O',operation='ABS')
georep.addHistogram(geoX='N+1:C:O:CA', title='N+1:C:O:CA',operation='ABS')
georep.addHistogram(geoX='N:CA', title='N:CA')

#next datom
georep.addHistogram(geoX='N:N+1', title='N:N+1')
georep.addHistogram(geoX='CA:CA+1', title='CA:CA+1')
georep.addHistogram(geoX='C:C+1', title='C:C+1')

#Check my closecontact assumption for C
georep.addCloseContact('4rek', 'C', 'CB', 6,1)
georep.addCloseContact('4rek', 'C', 'O', 6,1)
georep.addCloseContact('4rek', 'N', 'CA', 6,0)

#TEST
georep.addHistogram(geoX='N:CA', title='N:CA')
georep.addHistogram(geoX='CA:C', title='CA:C')
georep.addHistogram(geoX='C:O', title='C:O')
georep.addHistogram(geoX='CA:C:O', title='CA:C:O')
georep.addHistogram(geoX='N:CA:C)', title='N:CA:C')
georep.addHistogram(geoX='CA:O', title='CA:O')
georep.addHistogram(geoX='N:O', title='N:O')
georep.addHistogram(geoX='N:C', title='N:C')
georep.addHistogram(geoX='N:CA:C:O', title='N:CA:C:O')




georep.printToHtml('Carbonyl analysis',3,'carbonyl')
