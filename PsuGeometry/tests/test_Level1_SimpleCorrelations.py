# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
'''
Simple correlation report with minimum inputs

'''

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Levels/'

pdbList = ['1i1w']
georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=True, dssp=True)
georep.addScatter(geoX='PHI', geoY='PSI', title='Ramachandran Plot')
georep.addScatter(geoX='N:O', geoY='CB:O', title='NO-CBO')
georep.addScatter(geoX='PSI', geoY='N:O', title='PSI-NO')
georep.addScatter(geoX='PSI', geoY='CB:O', title='PSI-CBO')

georep.printToHtml('Simple Correlations', 2, 'SimpleCorr')
