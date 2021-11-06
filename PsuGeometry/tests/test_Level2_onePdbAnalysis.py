# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
import pandas as pd
'''
This script analyses a single pdb on a few correlation reports.
It uses the ghost structure so that if there are few residues it can be seen if it is in the correct region.

The structure looks at 6xe9, citation:
Yang, S., Tiwari, P., Lee, K.H. et al. Cryo-EM structure of the inhibited (10S) form of myosin II. Nature (2020). https://doi.org/10.1038/s41586-020-3007-0
'''
#########################################################################################
#pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
#edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
#printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Levels/'
pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/pdb_data/'
edDataPath = 'C:/Dev/Github/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/'
singlePdb = '7a6a'# Analyse a single pdb
#########################################################################################

print('---Set data lists up...')
pdbList = [singlePdb]
geoList = ['CA:C','N:CA','C:O','PHI','PSI','OMEGA','TAU','C-1:N','C:N+1','N:O','C-1:C','N:CA:C:O','TAU+1']
geoListCB = ['N:CA','PHI','PSI','CB:O','N:O','C-1:CB'] # we need a seperate CB as the results are only for measures that exists for ALL and thus we would miss glycine our of our reports otherwise
hueList = ['dssp','aa','bfactor']

print('---Create report object...')
georep = psu.GeoReport(pdbList,pdbDataPath,edDataPath,printPath,ed=False,dssp=False)

print('---Build dataframes of geometric data...')
dataCB = georep.getGeoemtryCsv(geoListCB, hueList)
data = georep.getGeoemtryCsv(geoList, hueList)

print('---Add the data plots...')
printList = []
georep.addHistogram(data=data,geoX='N:CA',title='N-CA',ghost=True,hue='pdbCode')
georep.addHistogram(data=data,geoX='CA:C',title='CA-C',ghost=True,hue='pdbCode')
georep.addHistogram(data=data,geoX='C:O',title='C-O',ghost=True)

georep.addDataView(singlePdb, geoX='x',geoY='y',title='X-Y Coordinates',hue='ridx',palette='gnuplot2_r',sort='DESC')
georep.addDataView(singlePdb, geoX='y',geoY='z',title='Y-Z Coordinates',hue='ridx',palette='gnuplot2_r',sort='ASC')
georep.addDataView(singlePdb, geoX='z',geoY='x',title='Z-X Coordinates',hue='ridx',palette='gnuplot2_r',sort='DESC')

georep.addScatter(data=data,geoX='PHI',geoY='PSI',title='Ramachandran',hue='aa',palette='gist_ncar',ghost=True)
georep.addScatter(data=data,geoX='PHI',geoY='C-1:N',title='Phi-Minus',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(data=data,geoX='PSI',geoY='C:N+1',title='Psi-Plus',hue='dssp',palette='gist_rainbow',ghost=True)

georep.addScatter(data=data,geoX='C:O',geoY='C:N+1',title='',hue='TAU+1',palette='copper_r',ghost=True)
georep.addScatter(data=data,geoX='TAU+1',geoY='C:O',title='',hue='C:N+1',palette='copper_r',ghost=True)
georep.addScatter(data=data,geoX='C:N+1',geoY='TAU+1',title='',hue='C:O',palette='copper_r',ghost=True)

georep.addScatter(data=data,geoX='N:CA',geoY='CA:C',title='Outliers-1',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')
georep.addScatter(data=data,geoX='N:CA',geoY='C:O',title='Outliers',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(data=data,geoX='OMEGA',geoY='TAU',title='Omega-Tau',hue='dssp',palette='gist_rainbow',ghost=True,sort='ASC')

georep.addScatter(data=data,geoX='PSI',geoY='N:O',title='Psi-NO',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')
georep.addScatter(data=dataCB,geoX='PSI',geoY='CB:O',title='Psi-CBO',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(data=dataCB,geoX='N:O',geoY='CB:O',title='Ellipse',hue='dssp',palette='gist_rainbow',ghost=True,sort='ASC')

georep.addScatter(data=dataCB,geoX='PHI',geoY='C-1:CB',title='Phi-CB',hue='aa',palette='gist_ncar',ghost=True,sort='ASC')
georep.addScatter(data=data,geoX='PHI',geoY='C-1:C',title='Phi-C',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(data=data,geoX='PSI',geoY='N:CA:C:O',title='Psi-Line',hue='dssp',palette='gist_rainbow',ghost=True,sort='ASC')

print('---Print final html report to file...')
georep.printToHtml('Analysis of ' + singlePdb,3,'singleppdb_' + singlePdb)
