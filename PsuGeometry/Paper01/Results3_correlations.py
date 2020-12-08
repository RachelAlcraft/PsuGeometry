# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
import time
'''
This script shows examples of all the "ordered correlations"
By looking at some plots with against resolution and probability density
This scrupt takes a very long time becauyse I do not create a dataframe to use for all plots
Each plot creates a seperate dataframe.
The advantage is that we do not lose any residues for not having CB or +- 1 or CG etc
But it takes ages
But, we only need to run it once
'''

#start timings
from datetime import datetime
import time
start = time.time()
print("RESULTS-- start time:",datetime.now().strftime("%H:%M:%S"))


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList1000 = geol.GeoPdbLists().getListPaper()
#pdbList1000 = pdbList1000[:50]

#geoList = ['C-1:C']
#hueList = ['resolution','aa','pdbCode']

print("RESULTS--", datetime.now().strftime("%H:%M:%S"),"Make report")
georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=True, includePdbs=False)

print("RESULTS--", datetime.now().strftime("%H:%M:%S"),"Row1")
georep.addScatter(geoX='PSI',geoY='N:O',hue='dssp', title='Psi/N:O', palette='gist_ncar',sort='NON')
georep.addScatter(geoX='PSI',geoY='CB:O',hue='dssp', title='Psi/CB:O', palette='gist_ncar',sort='NON')
georep.addScatter(geoX='N:O',geoY='CB:O',hue='dssp', title='N:O/CB:O', palette='gist_ncar',sort='NON')

print("RESULTS--", datetime.now().strftime("%H:%M:%S"),"Row2")
georep.addScatter(geoX='PHI',geoY='C-1:C',hue='dssp', title='PHI/C-1:C', palette='gist_ncar',sort='NON')
georep.addScatter(geoX='PHI',geoY='C-1:CB',hue='dssp', title='PHI/C-1:C', palette='gist_ncar',sort='NON')
georep.addScatter(geoX='CA-1:CA',geoY='CA-1:C-1:N:CA',hue='dssp', title='Cis/Pro', palette='gist_ncar',sort='NON')

print("RESULTS--", datetime.now().strftime("%H:%M:%S"),"Row3")
georep.addScatter(geoX='CA-2:CA-1:CA',geoY='CA:CA+1:CA+2',hue='dssp', title='The Square Plot', palette='gist_ncar',sort='NON')
georep.addScatter(geoX='PSI',geoY='N:CA:C:O',hue='dssp', title='Psi Lines', palette='gist_ncar',sort='NON')
georep.addScatter(geoX='CHI2',geoY='CHI3',hue='dssp', title='Histidine Lines', palette='gist_ncar',sort='NON',restrictions={'aa':'HIS'})

print("RESULTS--", datetime.now().strftime("%H:%M:%S"),"Row4")
georep.addScatter(geoX='CHI1',geoY='CA:CB:CG',hue='dssp', title='Proline 1', palette='gist_ncar',sort='NON',restrictions={'aa':'PRO'})
georep.addScatter(geoX='CHI1',geoY='CHI2',hue='dssp', title='Proline 2', palette='gist_ncar',sort='NON',restrictions={'aa':'PRO'})
georep.addScatter(geoX='CHI3',geoY='CHI4',hue='dssp', title='Proline 3', palette='gist_ncar',sort='NON',restrictions={'aa':'PRO'})

print("RESULTS--", datetime.now().strftime("%H:%M:%S"),"Row5")
georep.addScatter(geoX='C-1:N:CA',geoY='TAU',hue='dssp', title='Proline 4', palette='gist_ncar',sort='NON',restrictions={'aa':'PRO'})
georep.addScatter(geoX='PSI',geoY='CA-1:CA:CA+1',hue='dssp', title='Psi/CA', palette='gist_ncar',sort='NON')
georep.addScatter(geoX='TAU',geoY='PSI',hue='dssp', title='Bimodal Tau', palette='gist_ncar',sort='NON')

print("RESULTS--", datetime.now().strftime("%H:%M:%S"),"Print to html")
title = 'Correlations Collection'
georep.printToHtml(title,3,'Results3_correlations')

#END timings
end = time.time()
print("RESULTS-- end time:",datetime.now().strftime("%H:%M:%S"))
time_diff = end - start
timestring = str(int(time_diff / 60)) + "m " + str(int(time_diff % 60)) + "s"
print('RESULTS-- time taken',timestring)