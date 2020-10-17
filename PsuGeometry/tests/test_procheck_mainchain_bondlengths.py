
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/procheck/'

pdbList = ['2cnq','1i1w'] # structures with errors
pdbList = ['2cnq'] # structures with errors

for pdbCode in pdbList:

    georep = psu.GeoReport([pdbCode],pdbDataPath,edDataPath,printPath)

    #Create the geoplots
    printList = []
    georep.addHistogram(geoX='C-1:N', title='C-N', ghost=True, hue='rid',count=True,exclusions={'aa':'PRO'})
    georep.addHistogram(geoX='C-1:N', title='C-N', ghost=True, hue='rid', count=True,restrictions={'aa': 'PRO'})
    georep.addHistogram(geoX='C:O', title='C-O', ghost=True, hue='rid',count=True)
    georep.addHistogram(geoX='CA:C', title='CA-C', ghost=True, hue='rid',count=True, exclusions={'aa': 'GLY'})
    georep.addHistogram(geoX='CA:C', title='CA-C', ghost=True, hue='rid',count=True, restrictions={'aa': 'GLY'})
    georep.addHistogram(geoX='CA:CB', title='CA-CB', ghost=True, hue='rid',count=True, restrictions={'aa': 'ALA'})
    georep.addHistogram(geoX='CA:CB', title='CA-CB', ghost=True, hue='rid',count=True, restrictions={'aa': 'ILE,THR,VAL'})
    georep.addHistogram(geoX='CA:CB', title='CA-CB', ghost=True, hue='rid',count=True, exclusions={'aa': 'ALA,ILE,THR,VAL'})
    georep.addHistogram(geoX='N:CA', title='N-CA', ghost=True, hue='rid',count=True, exclusions={'aa': 'GLY,PRO'})
    georep.addHistogram(geoX='N:CA', title='N-CA', ghost=True, hue='rid',count=True, restrictions={'aa': 'GLY'})
    georep.addHistogram(geoX='N:CA', title='N-CA', ghost=True, hue='rid',count=True, restrictions={'aa': 'PRO'})


    geoName = pdbCode.lower() + "_procheck1"
    georep.printToHtml('Compare to Procheck Main-chain bond lengths',3,geoName)



