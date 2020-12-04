# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdbLists as geol
import pandas as pd
'''
This script looks at the distribution of N-CA for proline and glycine at different resolutions 
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList = geol.GeoPdbLists().getListJask()
georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath,ed=False,dssp=False,includePdbs=False)
geoList = ['N:CA','CA:C']
hueList = ['dssp','aa','bfactor','resolution']
data = georep.getGeoemtryCsv(geoList,hueList)

#Jaskolski
georep.addScatter(data=data,geoX='N:CA',geoY='CA:C',title='Jaskolski N-CA distribution',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')

#Compare to low res
pdbList = geol.GeoPdbLists().getList17()
georep.pdbCodes = pdbList
data1 = georep.getGeoemtryCsv(geoList,hueList)
georep.addScatter(data=data1,geoX='N:CA',geoY='CA:C',title='<0.9 N-CA distribution',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')

# Compare the jaskolski results to current top 1000 at different resolutions
'''
pdbList = geol.GeoPdbLists().getList1000()
georep.pdbCodes = pdbList
data = georep.getGeoemtryCsv(geoList,hueList)


#<=0.9
data1 = data[data['resolution'] <= 0.9]
georep.addScatter(data=data1,geoX='N:CA',geoY='CA:C',title='<0.9 N-CA distribution',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')

#=0.9-1.0
data2 = data[data['resolution'] > 0.9]
data2 = data2[data2['resolution'] <= 1.0]
georep.addScatter(data=data2,geoX='N:CA',geoY='CA:C',title='0.9-1.0 N-CA distribution',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')

#=1.0
data3 = data[data['resolution'] > 1.0]
georep.addScatter(data=data3,geoX='N:CA',geoY='CA:C',title='>1.0 N-CA distribution',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')
'''
georep.printToHtml('N-CA Distributions at different resolutions', 2, 'nca')
