
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
from PsuGeometry import GeoPlot as geopl


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/'

pdbList = ['2cnq']

#Create the main report object
geoList = []
for pdb in pdbList:
    pdb = pdb.lower()
    geoPdb = geop.GeoPdb(pdb, pdbDataPath,edDataPath)
    geoList.append(geoPdb)

georep = geor.GeoReport(geoList)
geoName = 'all'

# Create the geoemtric data
geoList = ['N:O','CB:O','N:CA:C:N+1','N:CA:C:O','N:CA','CA:C','CA:C:N+1:CA+1','N:CA:C']
geoListPhi = ['C-1:N:CA:C','C-1:C','C-1:CB']
geoListRama = ['N:CA:C:N+1','C-1:N:CA:C']
hueList = ['dssp','aa','bfactor','2FoFc','rid'] # note the hues are the sum od the atoms
data = georep.getGeoemtryCsv(geoList, hueList)
dataPhi = georep.getGeoemtryCsv(geoListPhi, hueList)
dataRama = georep.getGeoemtryCsv(geoListRama, hueList)

#Create the geoplots
printList = []

printList.append(geopl.GeoOverlay(geopl.GeoPlot(dataRama,'C-1:N:CA:C',geoY='N:CA:C:N+1',title='Ramachandran',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))
printList.append(geopl.GeoOverlay(geopl.GeoPlot(data,'N:CA',geoY='CA:C',title='Outliers',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))
printList.append(geopl.GeoOverlay(geopl.GeoPlot(data,'CA:C:N+1:CA+1',geoY='N:CA:C',title='Omega-Tau',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))

geoList = ['N:O','CB:O','C-1:N:CA:C','N:CA:C:N+1','N:CA:C:O','N:CA','CA:C','CA:C:N+1:CA+1','N:CA:C','C-1:C','C-1:CB']
printList.append(geopl.GeoOverlay(geopl.GeoPlot(data,'N:CA:C:N+1',geoY='N:O',title='Psi-NO',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))
printList.append(geopl.GeoOverlay(geopl.GeoPlot(data,'N:CA:C:N+1',geoY='CB:O',title='Psi-CBO',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))
printList.append(geopl.GeoOverlay(geopl.GeoPlot(data,'N:O',geoY='CB:O',title='Ellipse',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))

geoList = ['N:O','CB:O','C-1:N:CA:C','N:CA:C:N+1','N:CA:C:O','N:CA','CA:C','CA:C:N+1:CA+1','N:CA:C','C-1:C','C-1:CB']
printList.append(geopl.GeoOverlay(geopl.GeoPlot(data,'N:CA:C:N+1',geoY='N:CA:C:O',title='Psi-Line',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))
printList.append(geopl.GeoOverlay(geopl.GeoPlot(dataPhi,'C-1:N:CA:C',geoY='C-1:C',title='Phi-C',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))
printList.append(geopl.GeoOverlay(geopl.GeoPlot(dataPhi,'C-1:N:CA:C',geoY='C-1:CB',title='Phi-CB',hue='pdbCode',palette='gist_rainbow'),'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath))

# And finally create the reort with a file name of choice
georep.printCsvToHtml(printList,georep.pdbs,'Bad Geo Structures',3,printPath,'bad')

