
from PsuGeometry import GeoReport as psu

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/'


# Create the geoemtric data
geoX,geoY = 'PHI','PSI'
pdbList = ['1ucs'] # structures with errors
georep = psu.GeoReport(pdbList,pdbDataPath,edDataPath,printPath,ed=True,dssp=True)

georep.addScatter(geoX=geoX,geoY=geoY,title='Ramachandran Plot',hue='dssp',palette='gist_rainbow',ghost=True)
georep.addProbability(geoX=geoX,geoY=geoY,title='Ramachandran Plot',palette='gist_stern_r',ghost=True)
georep.addScatter(geoX=geoX,geoY=geoY,title='Ramachandran Plot',hue='2FoFc',palette='gnuplot_r',ghost=True)
georep.addScatter(geoX='CHI1',geoY='CHI2',title='Ramachandran Plot',hue='aa',palette='gist_rainbow')
georep.addScatter(geoX='CHI3',geoY='CHI4',title='Ramachandran Plot',hue='aa',palette='gist_rainbow')
georep.addHistogram(geoX='aa',count=True,palette='RosyBrown')

georep.printToHtml('Ramachandran Plot',3,'rama')
