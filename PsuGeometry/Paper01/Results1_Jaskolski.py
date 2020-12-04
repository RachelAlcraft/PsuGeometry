# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdbLists as geol
import pandas as pd
'''
This script replicates results from Jaskolski et al 2007 
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList = geol.GeoPdbLists().getListJask()
georep = geor.GeoReport(pdbList, pdbDataPath, edDataPath, printPath,ed=False,dssp=False,includePdbs=False)
geoList = ['N:CA','CA:C','C:O','C-1:N']
hueList = ['dssp','aa','bfactor','resolution']
data = georep.getGeoemtryCsv(geoList,hueList)

for pdb in pdbList:
    datapdb = data[data['pdbCode'] == pdb]
    georep.addHistogram(data=datapdb,geoX='N:CA', title=pdb + ' N:CA', exclusions={'aa': 'PRO,GLY'})
    georep.addHistogram(data=datapdb,geoX='CA:C', title=pdb + ' CA:C', exclusions={'aa': 'GLY'})
    #georep.addHistogram(data=datapdb, geoX='CA:C', title=pdb + ' CA:C', restrictions={'aa': 'GLY'})
    georep.addHistogram(data=datapdb,geoX='C-1:N', title=pdb + ' C-1:N', exclusions={'aa': 'PRO'})
    #georep.addHistogram(data=datapdb, geoX='C-1:N', title=pdb + ' C-1:N', restrictions={'aa': 'PRO'})
    georep.addHistogram(data=datapdb,geoX='C:O', title=pdb + ' C:O')

georep.addHistogram(data=data,geoX='N:CA', title='All N:CA', exclusions={'aa': 'PRO,GLY'})
georep.addHistogram(data=data,geoX='CA:C', title='All CA:C', exclusions={'aa': 'GLY'})
georep.addHistogram(data=data,geoX='C-1:N', title='All C-1:N', exclusions={'aa': 'PRO'})
georep.addHistogram(data=data,geoX='C:O', title='All C:O')

# Compare the jaskolski results to current top 1000 at different resolutions
pdbList = geol.GeoPdbLists().getList1000()
georep.pdbCodes = pdbList
data = georep.getGeoemtryCsv(geoList,hueList)

#<=0.9
data1 = data[data['resolution'] <= 0.9]
georep.addHistogram(data=data1,geoX='N:CA', title='1000<=0.9 N:CA', exclusions={'aa': 'PRO,GLY'})
georep.addHistogram(data=data1,geoX='CA:C', title='1000<=0.9 CA:C', exclusions={'aa': 'GLY'})
georep.addHistogram(data=data1,geoX='C-1:N', title='1000<=0.9 C-1:N', exclusions={'aa': 'PRO'})
georep.addHistogram(data=data1,geoX='C:O', title='1000<=0.9 C:O')

#=0.9-1.0
data2 = data[data['resolution'] > 0.9]
data2 = data2[data2['resolution'] <= 1.0]
georep.addHistogram(data=data2,geoX='N:CA', title='1000 0.9-1 N:CA', exclusions={'aa': 'PRO,GLY'})
georep.addHistogram(data=data2,geoX='CA:C', title='1000 0.9-1 CA:C', exclusions={'aa': 'GLY'})
georep.addHistogram(data=data2,geoX='C-1:N', title='1000 0.9-1 C-1:N', exclusions={'aa': 'PRO'})
georep.addHistogram(data=data2,geoX='C:O', title='1000 0.9-1 C:O')

#=1.0
data3 = data[data['resolution'] > 1.0]
georep.addHistogram(data=data3,geoX='N:CA', title='1000>1.0 N:CA', exclusions={'aa': 'PRO,GLY'})
georep.addHistogram(data=data3,geoX='CA:C', title='1000>1.0 CA:C', exclusions={'aa': 'GLY'})
georep.addHistogram(data=data3,geoX='C-1:N', title='1000>1.0 C-1:N', exclusions={'aa': 'PRO'})
georep.addHistogram(data=data3,geoX='C:O', title='1000>1.0 C:O')

georep.printToHtml('Jaskolski Table 2', 4, 'jask_t2')
