# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
'''
An example of every type of plot that can be produced
Some parameters over the defaults are added
'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Document s/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Levels/'
pdbCode = '6jvv'

pdbList = [pdbCode]
georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath)
georep.addScatter(geoX='PHI', geoY='PSI', title='Ramachandran Plot')
georep.addProbability(geoX='PHI', geoY='PSI', title='Ramachandran Probability',palette='ocean_r')
georep.addHistogram(geoX='N:CA',title='N:CA Differences')
georep.addDataView(pdbCode,'x','y',hue='2FoFc',palette='cubehelix_r',title='XY Coordinates of the solved structure')
georep.addDensityView(pdbCode,'x','y',hue='2FoFc',palette='cubehelix_r',peaks=True,title='XY Coordinates of the density peaks') # warning this is slow as it computes density peaks
georep.addCloseContact(pdbCode,'N','O',6,2,palette='terrain')
georep.addDifference(geoX='TAU',geoY='PSI',restrictionsA={'aa':'PRO'},exclusionsB={'aa':'PRO'})
georep.printToHtml('An Example of Each Plot Type', 3, 'each')