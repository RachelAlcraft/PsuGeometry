from PsuGeometry import GeoSlice as gsl
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

gs = gsl.GeoSlice()
slice = gs.load(printPath + "/slice.csv")

gs.addDensitySlice(slice,palette='cubehelix_r')

gs.printToHtml('Test slice only',4,printPath + "/sliceonly.html")