
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/'


# hand crafted report for 2bw4 non planar omega
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
from PsuGeometry import GeoPlot as geopl
#printPath = '???' #'wherever on your computer you wish to save the file'

# Create the GeoPdb object and the report object with just that single pdb
geoPdb = geop.GeoPdb('2bw4', pdbDataPath,edDataPath)
georep = geor.GeoReport([geoPdb])

# Choose the geometric calculations desired and the hues we might want to look at
geoList = ['CA:C:N+1:CA+1','N:CA:C','C:N+1','CA:CA+1']
hueList = ['dssp','aa','bfactor','2FoFc','rid'] # note the hues are the sum od the atoms
data = georep.getGeoemtryCsv(geoList, hueList)
printList = []
printList.append(geopl.GeoPlot(data,'CA:C:N+1:CA+1',geoY='N:CA:C',title='All residues',hue='rid',palette='gist_ncar'))
printList.append(geopl.GeoPlot(data, 'CA:C:N+1:CA+1',geoY='N:CA:C',title='All residues',hue='aa',palette='tab20'))
printList.append(geopl.GeoPlot(data, 'CA:C:N+1:CA+1',geoY='N:CA:C',title='All residues', palette='gist_ncar',plot='probability'))

# Narrow down the data to a zoomed in area aorund the residue of interest
riddata = data[data['CA:C:N+1:CA+1'] > 100]
riddata = riddata.query('aa =="HIS" or aa=="ASN"')
riddata = riddata.sort_values(by='rid', ascending=True)
printList.append(geopl.GeoPlot(riddata,'CA:C:N+1:CA+1',geoY='N:CA:C',title='Zoomed in',hue='rid',palette='gist_ncar'))
printList.append(geopl.GeoPlot(riddata, 'CA:C:N+1:CA+1',geoY='N:CA:C',title='Zoomed in', hue='aa',palette='tab20'))
printList.append(geopl.GeoPlot(riddata, 'CA:C:N+1:CA+1',geoY='N:CA:C',title='Zoomed in', palette='cubehelix_r'))

# Another validation check of omega against Ca distance
printList.append(geopl.GeoPlot(data, 'CA:C:N+1:CA+1',geoY='CA:CA+1',title='Omega CA', hue='rid',palette='gist_ncar'))
printList.append(geopl.GeoPlot(data, 'CA:C:N+1:CA+1',geoY='CA:CA+1',title='Omega CA', hue='aa',palette='tab20'))
printList.append(geopl.GeoPlot(data, 'CA:C:N+1:CA+1',geoY='C:N+1',title='Omega CA', palette='cubehelix_r'))

# Narrow down the data to a zoomed in area aorund the residue of interest
serdata = data[data['CA:CA+1'] > 5]
printList.append(geopl.GeoPlot(serdata, 'CA:C:N+1:CA+1',geoY='CA:CA+1',title='Zoomed in', hue='rid',palette='Accent'))
printList.append(geopl.GeoPlot(serdata, 'CA:C:N+1:CA+1',geoY='N:CA:C',title='Zoomed in', hue='aa',palette='Accent'))
printList.append(geopl.GeoPlot(serdata, 'CA:C:N+1:CA+1',geoY='C:N+1',title='Zoomed in', palette='Accent'))


# And finally create the reort with a file name of choice
georep.printCsvToHtml(printList,georep.pdbs,'Non Planar Omega in 2bw4',3,printPath,'2bw4_omega')
