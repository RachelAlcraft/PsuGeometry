# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
'''
This script looks electron density correlations 
'''


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper01/'

pdbList = ['2bw4','5nqo','1ejg','6q53']

geoList = ['PHI','PSI','TAU']
hueList = ['aa','bfactor','rid','resolution','pdbCode','2FoFc','dssp']


georep2 = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath)

for pdb in pdbList:
    georep = psu.GeoReport([pdb], pdbDataPath, edDataPath, printPath)
    #georep.printReport('DataPerPdb', 'Results9_data')
    #georep.printReport('Slow_DensityPeaksPerPdb', 'Results9_density')


    data = georep.getGeoemtryCsv(geoList, hueList)

    georep2.addScatter(data=data, geoX='ridx', geoY='2FoFc', hue='TAU', title=pdb + ' Backbone tau', palette='jet',sort='NON', vmin=106, vmax=116)
    georep2.addScatter(data=data, geoX='ridx', geoY='2FoFc', hue='PHI', title=pdb + ' Backbone phi', palette='jet',sort='NON', vmin=-170, vmax=170)
    georep2.addScatter(data=data, geoX='ridx', geoY='2FoFc', hue='PSI', title=pdb + ' Backbone psi', palette='jet',sort='NON', vmin=-170, vmax=170)
    georep2.addScatter(data=data, geoX='ridx', geoY='2FoFc', hue='dssp', title=pdb + ' Backbone dssp', palette='jet_r', sort='NON')


georep2.printToHtml('Backbone data', 4, 'Results11_tau_phi_psi_data')










