from PsuGeometry import GeoReport as psu


#printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'
printPath = 'F:/Code/ProteinDataFiles/ccp4_data/'

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'


georep = psu.GeoReport([],pdbDataPath,edDataPath,printPath,ed=True,dssp=True)


#slice0d = georep.loadSlice(printPath + "1ejg_diff.ccp4_slice.csv")
#slice1 = georep.loadSlice(printPath + "1ejg.ccp4_slice1.csv")
slice00 = georep.loadSlice(printPath + "1ejg.ccp4_slice0.csv")
slice01 = georep.loadSlice(printPath + "1ejg.ccp4_slice1.csv")
slice02 = georep.loadSlice(printPath + "1ejg.ccp4_slice2.csv")

slice00x = georep.loadSlice(printPath + "1p7h.ccp4_slice0.csv")
slice01x = georep.loadSlice(printPath + "1p7h.ccp4_slice1.csv")
slice02x = georep.loadSlice(printPath + "1p7h.ccp4_slice2.csv")
#slice26 = georep.loadSlice(printPath + "61ejg.ccp4_slice2.csv")


georep.addSlice(slice00,palette='cubehelix_r')
georep.addSlice(slice01,palette='seismic')
georep.addSlice(slice02,palette='cubehelix')

georep.addSlice(slice00x,palette='cubehelix_r')
georep.addSlice(slice01x,palette='seismic')
georep.addSlice(slice02x,palette='cubehelix')
#georep.addSlice(slice06,palette='cubehelix_r')
#georep.addSlice(slice26,palette='cubehelix')




georep.printToHtml('Slices from Density Flight',3,'ccp4_slices_diffs')