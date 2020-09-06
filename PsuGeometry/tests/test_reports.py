#import sys
#sys.path.append("/home/rachel/Documents/Bioinformatics/BbkProject/Project2/PythonCode/PdbGeometry/")

from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results/'

# split list in 2 for memory purposes
#pdbList = ['1ejg','1us0','1tt8','1i1w','1ucs','1yk4','1yk4','1hje','1r6j']
#pdbList = ['2bw4','3nir','3x2m','2VB1','3A39','2b97','2OV0','2WFI']
#pdbList = ['4ZM7','4REK','4ZM7','5D8V','5NW3','5qkw']
pdbList = ['6jvv','6rr2','6E6O','6S2M','6shk','6fgz','6ctd','6fwf','6q53']

# split list in 2 for memory purposes

#pdbList = ['3blm']
#pdbList = ['2b97','2OV0','2WFI']
#pdbList = ['3nir','1ejg','1us0','3x2m']
#pdbList = ['1ejg','6ctd']


if False:
    geoList = []
    for pdb in pdbList:
        pdb = pdb.lower()
        geoPdb = geop.GeoPdb(pdb, pdbDataPath,edDataPath)
        geoList.append(geoPdb)

    georep = geor.GeoReport(geoList)

    geoName = 'all'
    georep.printReport('Ramachandran',  printPath,geoName + '_rama' )
    georep.printReport('Sp2Planarity',printPath,geoName+ '_sp2')
    georep.printReport('Sp3Tetrahedra',printPath,geoName + '_sp3')
    georep.printReport('BackboneOutliers', printPath,geoName+ '_bbone')
    georep.printReport('OmegaCis', printPath,geoName + '_ocis')
    georep.printReport('RachelsChoice', printPath,geoName + '_rae')
    georep.printReport('DataPerPdb', printPath,geoName + '_data')
    georep.printReport('DensityPerPdb',printPath,geoName +  '_denPoints') # This one is slow
    georep.printReport('DensityPeaksPerPdb', printPath, geoName + '_denPeaks')  # This one is slow

if False:
    for pdb in pdbList:
        geoPdb = geop.GeoPdb(pdb, pdbDataPath, edDataPath)
        georep = geor.GeoReport([geoPdb])
        georep.printReport('Ramachandran',  printPath,geoPdb.pdbCode + '_rama' )
        georep.printReport('Sp2Planarity',printPath,geoPdb.pdbCode + '_sp2')
        georep.printReport('Sp3Tetrahedra',printPath,geoPdb.pdbCode + '_sp3')
        georep.printReport('BackboneOutliers', printPath,geoPdb.pdbCode + '_bbone')
        georep.printReport('OmegaCis', printPath,geoPdb.pdbCode + '_ocis')
        georep.printReport('RachelsChoice', printPath,geoPdb.pdbCode + '_rae')
        georep.printReport('DataPerPdb', printPath,geoPdb.pdbCode + '_data')
        georep.printReport('DensityPerPdb',printPath,geoPdb.pdbCode +  '_denPoints') # This one is slow
        georep.printReport('DensityPeaksPerPdb', printPath, geoPdb.pdbCode + '_denPeaks')  # This one is slow

# hand crafted report for 2bw4 non planar omega
if False:

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


if True: #Jaskolski test My computer cannot cope with the memory requirements!!
    pdbList = ['1ejg', '1ucs', '1us0', '1yk4', '1r6j', '1hje', '3al1', '2b97', '1gci', '1x6z']
    pdbList = ['1ejg', '1ucs', '1us0']
    geopList = []
    # Create a list of GeoPdb objects
    for pdb in pdbList:
        pdb = pdb.lower()
        geoPdb = geop.GeoPdb(pdb, pdbDataPath,edDataPath)
        geopList.append(geoPdb)
    #Init the report object
    georep = geor.GeoReport(geopList)
    #specify the paramaters
    geoList = ['C-1:N','N:CA','CA:C','C:N+1','C:O','N:CA:C']
    hueList = ['dssp','aa','bfactor','2FoFc','rid'] # note the hues are the sum od the atoms
    #Create the dataframe
    data = georep.getGeoemtryCsv(geoList, hueList)
    #Specify the plots required
    printList = []
    printList.append(['C-1:N', data, 'C-1:N', ''])
    printList.append(['N:CA', data, 'N:CA', ''])
    printList.append(['CA:C', data, 'CA:C', ''])
    printList.append(['C:N+1', data, 'C:N+1', ''])
    printList.append(['C:O', data, 'C:O', ''])
    printList.append(['N:CA:C', data, 'N:CA:C', ''])
    #Print the report
    georep.printReport('RachelsChoice', printPath, 'rae')
    georep.printCsvToHtml(printList,georep.pdbs,'Backbone Report, Jaskolski Pdbs',2,printPath,'jaskolski')