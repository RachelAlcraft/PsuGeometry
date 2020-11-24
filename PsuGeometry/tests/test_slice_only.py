from PsuGeometry import GeoSlice as gsl
#printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'
printPath = 'F:/Code/ProteinDataFiles/ccp4_data/'

gs = gsl.GeoSlice()
slice = gs.load(printPath + "1ejg_diff.ccp4_slice.csv")

gs.addDensitySlice(slice,palette='cubehelix_r')

gs.printToHtml('Test slice only',4,printPath + '1ejg_diff.ccp4_slice.html')