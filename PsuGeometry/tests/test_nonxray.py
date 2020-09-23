
from PsuGeometry import GeoReport as psu

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/nonxray/'

pdbList = ['5xsg','6j60','6uos','6kj2','6bzm','6m9j','6cf4','6axz'] # high res cryoem
geoName = 'all'

# Create the geoemtric data
geoPsi = ['N:O','CB:O','N:CA:C:N+1']
geoListMain = ['CA:C','N:CA','C:O']
hueList = ['dssp','aa','bfactor','2FoFc','rid'] # note the hues are the sum od the atoms

georep = psu.GeoReport(pdbList,pdbDataPath,edDataPath,printPath)

dataPsi = georep.getGeoemtryCsv(geoPsi, hueList)
dataMain = georep.getGeoemtryCsv(geoListMain, hueList)

#Create the geoplots
printList = []
georep.addHistogram(data=dataMain,geoX='N:CA',title='N-CA',ghost=True,hue='pdbCode')
georep.addHistogram(data=dataMain,geoX='CA:C',title='CA-C',ghost=True,hue='pdbCode')
georep.addHistogram(data=dataMain,geoX='C:O',title='C-O',ghost=True)

georep.addScatter(geoX='C-1:N:CA:C',geoY='N:CA:C:N+1',title='Ramachandran',hue='aa',palette='gist_rainbow',ghost=True)
georep.addScatter(geoX='C-1:N:CA:C',geoY='C-1:N',title='Phi-Minus',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(geoX='N:CA:C:N+1',geoY='C:N+1',title='Psi-Plus',hue='pdbCode',palette='gist_rainbow',ghost=True)

georep.addScatter(data=dataMain,geoX='N:CA',geoY='CA:C',title='Outliers-1',hue='aa',palette='gist_rainbow',ghost=True)
georep.addScatter(data=dataMain,geoX='N:CA',geoY='C:O',title='Outliers',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(geoX='CA:C:N+1:CA+1',geoY='N:CA:C',title='Omega-Tau',hue='pdbCode',palette='gist_rainbow',ghost=True)

georep.addScatter(data=dataPsi,geoX='N:CA:C:N+1',geoY='N:O',title='Psi-NO',hue='aa',palette='gist_rainbow',ghost=True)
georep.addScatter(data=dataPsi,geoX='N:CA:C:N+1',geoY='CB:O',title='Psi-CBO',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(data=dataPsi,geoX='N:O',geoY='CB:O',title='Ellipse',hue='pdbCode',palette='gist_rainbow',ghost=True)

georep.addScatter(geoX='N:CA:C:N+1',geoY='N:CA:C:O',title='Psi-Line',hue='aa',palette='gist_rainbow',ghost=True)
georep.addScatter(geoX='C-1:N:CA:C',geoY='C-1:C',title='Phi-C',hue='bfactor',palette='copper_r',ghost=True)
georep.addScatter(geoX='C-1:N:CA:C',geoY='C-1:CB',title='Phi-CB',hue='pdbCode',palette='gist_rainbow',ghost=True)

georep.printToHtml('Cryo EM Structures',3,'nonxray')
