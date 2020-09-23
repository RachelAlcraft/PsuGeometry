
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
from PsuGeometry import GeoPlot as geopl


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/'

pdbList = ['5xsg','6j60','6uos','6kj2','6bzm','6m9j','6cf4','6axz'] # high res cryoem
#pdbList = ['2chh','2cnq','2ggc','3ccd','3rwn']

# Create the dummy report object
geoDummy = geop.GeoPdb('ghost', pdbDataPath,edDataPath)
dummyReport = geor.GeoReport([geoDummy])

'''
Currently deprecated as I have added the isGhost which is the only use of overlay
'''

#Create the main report object
geoList = []
for pdb in pdbList:
    pdb = pdb.lower()
    geoPdb = geop.GeoPdb(pdb, pdbDataPath,edDataPath)
    geoList.append(geoPdb)

georep = geor.GeoReport(geoList)
geoName = 'all'

# Create the geoemtric data
geoList = ['N:O','CB:O']
hueList = ['dssp','aa','bfactor','2FoFc','rid'] # note the hues are the sum od the atoms
data = georep.getGeoemtryCsv(geoList, hueList)
dummydata = dummyReport.getGeoemtryCsv(geoList, hueList)

#Create the geoplots
geoData = geopl.GeoPlot(data,'N:O',geoY='CB:O',title='Main',hue='pdbCode',palette='gist_ncar')
geoDummy = geopl.GeoPlot(dummydata,'N:O',geoY='CB:O',title='ghost',hue='pdbCode',palette='Greys')
geoOverlay = geopl.GeoOverlay(geoDummy,geoData,title='Overlay')
geoOverlay2 = geopl.GeoOverlay(geoData,'',title='ghost',pdbDataPath=pdbDataPath,edDataPath=edDataPath)

printList = []
printList.append(geoData)
printList.append(geoDummy)
printList.append(geoOverlay)
printList.append(geoOverlay2)

# And finally create the reort with a file name of choice
georep.printCsvToHtml(printList,georep.pdbs,'Testing the overlay',2,printPath,'overlay')

