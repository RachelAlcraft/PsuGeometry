
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/'

from PsuGeometry import GeoReport as geor

# Create the GeoPdb object and the report object with just that single pdb
pdbs = ['2bw4','1ejg','1us0']
pdbs = ['4rek']

georep = geor.GeoReport(pdbs,pdbDataPath,edDataPath,printPath,False,False)
# Choose the geometric calculations desired and the hues we might want to look at
title = 'Protein Warhol'
georep.addProbability(geoX='N:CA:C:CB',geoY='N:CA:C:N+1',title='', palette='Spectral')#,restriction={'aa':'PRO'})
georep.addProbability(geoX='N:CA:C:CB',geoY='N:CA:C:N+1',title='', palette='twilight_shifted')
georep.addProbability(geoX='N:CA:C:CB',geoY='N:CA:C:N+1',title='', palette='inferno')
georep.addProbability(geoX='N:CA:C:CB',geoY='N:CA:C:N+1',title='', palette='viridis_r')
# And finally create the reort with a file name of choice
georep.printToHtml(title,2,'warhol')

# Choose the geometric calculations desired and the hues we might want to look at
title = 'Protein Halo'
georep.addProbability(geoX='N:O',geoY='CB:O',title='', palette='Spectral')
georep.addProbability(geoX='N:O',geoY='CB:O',title='', palette='twilight_shifted')
georep.addProbability(geoX='N:O',geoY='CB:O',title='', palette='inferno')
georep.addProbability(geoX='N:O',geoY='CB:O',title='', palette='nipy_spectral_r')
# And finally create the reort with a file name of choice
georep.printToHtml(title,2,'angel')
