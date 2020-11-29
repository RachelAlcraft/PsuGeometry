
from PsuGeometry import GeoReport as psu

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/'

pdbList = ['1i1w']

georep = psu.GeoReport(pdbList,pdbDataPath,edDataPath,printPath,ed=True,dssp=True)

georep.addScatter(geoX='PHI',geoY='PSI',title='Ramachandran Plot',hue='2FoFc',palette='cubehelix_r',ghost=True)
georep.addScatter(geoX='N:O',geoY='CB:O',title='NO-CBO',hue='2FoFc',palette='cubehelix_r',ghost=True)
georep.addScatter(geoX='PSI',geoY='N:O',title='PSI-NO',hue='2FoFc',palette='cubehelix_r',ghost=True)
georep.addScatter(geoX='PSI',geoY='CB:O',title='PSI-CBO',hue='2FoFc',palette='cubehelix_r',ghost=True)


georep.printToHtml('Simple Correlations',2,'SimpleCorr')
