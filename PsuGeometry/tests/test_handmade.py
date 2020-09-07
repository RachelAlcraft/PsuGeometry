
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/'


# hand crafted report for 2bw4 non planar omega
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
#printPath = '???' #'wherever on your computer you wish to save the file'

# Create the GeoPdb object and the report object with just that single pdb
geoPdb = geop.GeoPdb('2bw4', pdbDataPath,edDataPath)
georep = geor.GeoReport([geoPdb])

# Choose the geometric calculations desired and the hues we might want to look at
geoList = ['CA:C:N+1:CA+1','N:CA:C','C:N+1','CA:CA+1']
hueList = ['dssp','aa','bfactor','2FoFc','rid'] # note the hues are the sum od the atoms
data = georep.getGeoemtryCsv(geoList, hueList)
printList = []
printList.append(['All residues',data,'CA:C:N+1:CA+1','N:CA:C','','rid','gist_ncar',False,0,0])
printList.append(['All residues', data, 'CA:C:N+1:CA+1', 'N:CA:C', '', 'aa', 'tab20', False, 0, 0])
printList.append(['All residues', data, 'CA:C:N+1:CA+1', 'N:CA:C', '', '2FoFc', 'cubehelix_r', False, 0, 0])

# Narrow down the data to a zoomed in area aorund the residue of interest
riddata = data[data['CA:C:N+1:CA+1'] > 100]
riddata = riddata.query('aa =="HIS" or aa=="ASN"')
riddata = riddata.sort_values(by='rid', ascending=True)
printList.append(['Zoomed in',riddata,'CA:C:N+1:CA+1','N:CA:C','','rid','gist_ncar',False,0,0])
printList.append(['Zoomed in', riddata, 'CA:C:N+1:CA+1', 'N:CA:C', '', 'aa', 'tab20', False, 0, 0])
printList.append(['Zoomed in', riddata, 'CA:C:N+1:CA+1', 'N:CA:C', '', '2FoFc', 'cubehelix_r', False, 0, 0])

# Another validation check of omega against Ca distance
printList.append(['Omega CA', data, 'CA:C:N+1:CA+1', 'CA:CA+1', '', 'rid', 'gist_ncar', False, 0, 0])
printList.append(['Omega CA', data, 'CA:C:N+1:CA+1', 'CA:CA+1', '', 'aa', 'tab20', False, 0, 0])
printList.append(['Omega CA', data, 'CA:C:N+1:CA+1', 'C:N+1', '', '2FoFc', 'cubehelix_r', False, 0, 0])

# Narrow down the data to a zoomed in area aorund the residue of interest
serdata = data[data['CA:CA+1'] > 5]
printList.append(['Zoomed in', serdata, 'CA:C:N+1:CA+1', 'CA:CA+1', '', 'rid', 'Accent', False, 0, 0])
printList.append(['Zoomed in', serdata, 'CA:C:N+1:CA+1', 'N:CA:C', '', 'aa', 'Accent', False, 0, 0])
printList.append(['Zoomed in', serdata, 'CA:C:N+1:CA+1', 'C:N+1', '', '2FoFc', 'Accent', False, 0, 0])

# And finally create the reort with a file name of choice
georep.printCsvToHtml(printList,georep.pdbs,'Non Planar Omega in 2bw4',3,printPath,'2bw4_omega')
