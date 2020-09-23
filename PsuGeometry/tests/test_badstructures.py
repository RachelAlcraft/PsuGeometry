
from PsuGeometry import GeoReport as psu

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/bad/'


# Create the geoemtric data
geoPsi = ['N:O','CB:O','N:CA:C:N+1']
geoListMain = ['CA:C','N:CA','C:O']
hueList = ['dssp','aa','bfactor','2FoFc','rid'] # note the hues are the sum od the atoms

pdbList = ['2cnq','1i1w'] # structures with errors

for pdbCode in pdbList:

    georep = psu.GeoReport([pdbCode],pdbDataPath,edDataPath,printPath)

    dataPsi = georep.getGeoemtryCsv(geoPsi, hueList)
    dataMain = georep.getGeoemtryCsv(geoListMain, hueList)

    #Create the geoplots
    printList = []
    georep.addHistogram(data=dataMain,geoX='N:CA',title='N-CA',ghost=True,splitKey='pdbCode')
    georep.addHistogram(data=dataMain,geoX='CA:C',title='CA-C',ghost=True)
    georep.addHistogram(data=dataMain,geoX='C:O',title='C-O',ghost=True)

    georep.addScatter(geoX='C-1:N:CA:C',geoY='N:CA:C:N+1',title='Ramachandran',hue='dssp',palette='gist_rainbow',ghost=True)
    georep.addScatter(geoX='C-1:N:CA:C',geoY='C-1:N',title='Phi-Minus',hue='pdbCode',palette='gist_rainbow',ghost=True)
    georep.addScatter(geoX='N:CA:C:N+1',geoY='C:N+1',title='Psi-Plus',hue='pdbCode',palette='gist_rainbow',ghost=True)

    georep.addScatter(data=dataMain,geoX='N:CA',geoY='CA:C',title='Outliers-1',hue='dssp',palette='gist_rainbow',ghost=True)
    georep.addScatter(data=dataMain,geoX='N:CA',geoY='C:O',title='Outliers',hue='dssp',palette='gist_rainbow',ghost=True)
    georep.addScatter(geoX='CA:C:N+1:CA+1',geoY='N:CA:C',title='Omega-Tau',hue='aa',palette='gist_rainbow',ghost=True)

    georep.addScatter(data=dataPsi,geoX='N:CA:C:N+1',geoY='N:O',title='Psi-NO',hue='dssp',palette='gist_rainbow',ghost=True)
    georep.addScatter(data=dataPsi,geoX='N:CA:C:N+1',geoY='CB:O',title='Psi-CBO',hue='bfactor',palette='copper_r',ghost=True)
    georep.addScatter(data=dataPsi,geoX='N:O',geoY='CB:O',title='Ellipse',hue='2FoFc',palette='viridis_r',ghost=True)

    georep.addScatter(geoX='N:CA:C:N+1',geoY='N:CA:C:O',title='Psi-Line',hue='dssp',palette='gist_rainbow',ghost=True)
    georep.addScatter(geoX='C-1:N:CA:C',geoY='C-1:C',title='Phi-C',hue='bfactor',palette='copper_r',ghost=True)
    georep.addScatter(geoX='C-1:N:CA:C',geoY='C-1:CB',title='Phi-CB',hue='2FoFc',palette='viridis_r',ghost=True)

    geoName = pdbCode.lower() + "_bad"
    georep.printToHtml('Structures With Unexpected Geometry',3,geoName)
