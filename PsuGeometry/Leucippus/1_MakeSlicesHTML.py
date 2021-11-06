import FilesAndCSVs as fac
from PsuGeometry import GeoReport as psu

pdbDataPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/PdbFiles/'
edDataPath = 'C:/Dev/Github/ProteinDataFiles/ccp4_data/'
printPath = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Analysis/'

#### For Adjusted 7a6a ######
#leuFileName = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Analysis/7a6a_SLICESFILE_DensityAdjusted.csv'
#tags = ['Chain A','Chain B','Chain E','Chain e','Chain r','Chain G','Chain I','Chain M','Chain O','Chain Q','Chain S','Chain U','Chain W','Chain Y','Chain 2','Chain 4','Chain F','Chain H','Chain P','Chain X','Chain 6']
#header = 'Slices from 7a6a GLU 134 in 21 chains - Density Adjusted'
#htmlName = '7a6a_slices_adj'

#### For Original 7a6a ######
leuFileName = 'C:/Dev/Github/ProteinDataFiles/LeicippusTesting/Analysis/7a6a_SLICESFILE_Original.csv'
tags = ['Chain A','Chain 1','Chain K','Chain a','Chain B','Chain E','Chain e','Chain r','Chain G','Chain I','Chain M','Chain O','Chain Q','Chain S','Chain U','Chain W','Chain Y','Chain 2','Chain 4','Chain F','Chain H','Chain P','Chain X','Chain 6']
header = 'Slices from 7a6a GLU 134 in 21 chains - Original Positions'
htmlName = '7a6a_slices_orig'


numSlices, results = fac.getCsvFromCppResults_Slices(leuFileName)
slicesList = []
radiantsList = []
laplaciansList = []

for i in range(numSlices):
    ID = 'DENSITYSLICE_' + str(i)
    df = results[ID]
    mtx = fac.DataFrameToMatrix(df, 'Density')
    slicesList.append(mtx)

    ID = 'RADIANTSLICE_' + str(i)
    df = results[ID]
    mtx = fac.DataFrameToMatrix(df, 'Radiant')
    radiantsList.append(mtx)

    ID = 'LAPLACIANSLICE_' + str(i)
    df = results[ID]
    mtx = fac.DataFrameToMatrix(df, 'Laplacian')
    laplaciansList.append(mtx)

georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=True, dssp=True)
ch = 0
for s in range(len(slicesList)):
    dens = slicesList[s]
    rads = radiantsList[s]
    laps = laplaciansList[s]
    title = tags[ch]
    ch += 1
    georep.addSlice(dens, palette='magma_r',title=title + ' Density ' + str(dens[0,0]))
    georep.addSlice(rads, palette='bone', title=title + ' Radiant ' + str(rads[0,0]),Contour=False)
    georep.addSlice(laps, palette='magma', title=title + ' Laplacian ' + str(laps[0,0]))


georep.addSlices(slicesList, palette='magma_r', title='Average Density', logged=False, centre=False)
georep.addSlices(radiantsList, palette='bone', title='Average Radiant', logged=False, centre=False,Contour=False)
georep.addSlices(laplaciansList, palette='magma', title='Average Laplacian', logged=False, centre=False)

georep.printToHtml(header, 3, htmlName)