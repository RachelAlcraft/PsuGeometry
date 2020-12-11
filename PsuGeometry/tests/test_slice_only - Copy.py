from PsuGeometry import GeoSlice as gsl
#printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'
printPath = 'F:/Code/ProteinDataFiles/ccp4_data/'

gs = gsl.GeoSlice()
slice1 = gs.load(printPath + "1ejg.ccp4_slice.csv")
slice1d = gs.load(printPath + "1ejg_diff.ccp4_slice.csv")
slice2 = gs.load(printPath + "1us0.ccp4_slice.csv")
slice2d = gs.load(printPath + "1us0_diff.ccp4_slice.csv")
slice3 = gs.load(printPath + "6jvv.ccp4_slice.csv")
slice3d = gs.load(printPath + "6jvv_diff.ccp4_slice.csv")

gs.addDensitySlice(slice1,palette='cubehelix_r')
gs.addDensitySlice(slice1,palette='cubehelix_r',logged=True)
gs.addDensitySlice(slice1d,palette='seismic',centre=True)

gs.addDensitySlice(slice2,palette='cubehelix_r')
gs.addDensitySlice(slice2,palette='cubehelix_r',logged=True)
gs.addDensitySlice(slice2d,palette='seismic',centre=True)

gs.addDensitySlice(slice3,palette='cubehelix_r')
gs.addDensitySlice(slice3,palette='cubehelix_r',logged=True)
gs.addDensitySlice(slice3d,palette='seismic',centre=True)

gs.printToHtml('Slices from Density Flight',3,printPath + 'ccp4_slices.html')