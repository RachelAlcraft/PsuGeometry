
from PsuGeometry import GeoReport as psu

pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/PdbFiles/'
printPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/'
edDataPath = ''

pdbA = '1ejg'
pdbB = 'density_1ejg'
pdbC = 'syndensity_1ejg'
geoName = '1ejg'

# Create the geoemtric data
geos = ['PHI','PSI','N:CA','CA:C','C:O','C:N+1','TAU+1','N:N+1','TAU']
hueList = ['dssp','aa','bfactor','rid','atomNo'] # note the hues are the sum od the atoms

georepA = psu.GeoReport([pdbA], pdbDataPath, edDataPath, printPath,ed=False,dssp=False)
georepB = psu.GeoReport([pdbB], pdbDataPath, edDataPath, printPath,ed=False,dssp=False)
georepC = psu.GeoReport([pdbC], pdbDataPath, edDataPath, printPath,ed=False,dssp=False)
dataA = georepA.getGeoemtryCsv(geos, hueList)
dataB = georepB.getGeoemtryCsv(geos, hueList)
dataC = georepC.getGeoemtryCsv(geos, hueList)

#Create the geoplots

georepA.addHistogram(data=dataA,geoX='N:CA',title=pdbA,hue='ridx')
georepA.addHistogram(data=dataB, geoX='N:CA', title=pdbB,hue='ridx')
georepA.addHistogram(data=dataC, geoX='N:CA', title=pdbC,hue='ridx')

georepA.addHistogram(data=dataA, geoX='CA:C', title=pdbA,hue='ridx')
georepA.addHistogram(data=dataB, geoX='CA:C', title=pdbB,hue='ridx')
georepA.addHistogram(data=dataC, geoX='CA:C', title=pdbC,hue='ridx')

georepA.addHistogram(data=dataA, geoX='C:O', title=pdbA,hue='ridx')
georepA.addHistogram(data=dataB, geoX='C:O', title=pdbB,hue='ridx')
georepA.addHistogram(data=dataC, geoX='C:O', title=pdbC,hue='ridx')

georepA.addHistogram(data=dataA, geoX='C:N+1', title=pdbA,hue='ridx')
georepA.addHistogram(data=dataB, geoX='C:N+1', title=pdbB,hue='ridx')
georepA.addHistogram(data=dataC, geoX='C:N+1', title=pdbC,hue='ridx')

georepA.addScatter(data=dataA,geoX='PHI',geoY='PSI',title=pdbA,hue='aa',palette='gist_rainbow',ghost=True)
georepA.addScatter(data=dataB, geoX='PHI', geoY='PSI', title=pdbB, hue='aa', palette='gist_rainbow', ghost=True)
georepA.addScatter(data=dataC, geoX='PHI', geoY='PSI', title=pdbC, hue='aa', palette='gist_rainbow', ghost=True)

georepA.addScatter(data=dataA, geoX='N:CA', geoY='CA:C', title=pdbA, hue='aa', palette='gist_rainbow', ghost=False)
georepA.addScatter(data=dataB, geoX='N:CA', geoY='CA:C', title=pdbB, hue='aa', palette='gist_rainbow', ghost=False)
georepA.addScatter(data=dataC, geoX='N:CA', geoY='CA:C', title=pdbC, hue='aa', palette='gist_rainbow', ghost=False)

georepA.addScatter(data=dataA, geoX='C:O', geoY='C:N+1', title=pdbA, hue='aa', palette='gist_rainbow', ghost=False)
georepA.addScatter(data=dataB, geoX='C:O', geoY='C:N+1', title=pdbB, hue='aa', palette='gist_rainbow', ghost=False)
georepA.addScatter(data=dataC, geoX='C:O', geoY='C:N+1', title=pdbC, hue='aa', palette='gist_rainbow', ghost=False)

georepA.addScatter(data=dataA, geoX='PSI', geoY='N:N+1', title=pdbA, hue='TAU', palette='gist_rainbow', ghost=True)
georepA.addScatter(data=dataB, geoX='PSI', geoY='N:N+1', title=pdbB, hue='TAU', palette='gist_rainbow', ghost=True)
georepA.addScatter(data=dataC, geoX='PSI', geoY='N:N+1', title=pdbC, hue='TAU', palette='gist_rainbow', ghost=True)

georepA.printToHtml('Comparison of 3 - ' + pdbA + ':' + pdbB + ':' + pdbC,3,'Comp3_' + pdbA)
