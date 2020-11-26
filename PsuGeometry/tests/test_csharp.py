


import pdb_eda


from PsuGeometry import GeoTransformation as trans
from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoSpace as space

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
matDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/mat_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

pal = "cubehelix_r"

pdbCode = '1us0'
central = [16.748,-5.111,11.617]
linear = [15.845,-5.415,11.313]
planar = [15.865,-6.151,11.989]

georepa = geor.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath)
sfcp = georepa.addDensitySlice(pdbCode,2,-1,2,1,central,linear,planar,interp='linear',palette=pal)

1
georep = geor.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath)
sfc = georep.addDensitySlice(pdbCode,2,-1,8,0.1,central,linear,planar,interp='linear',palette=pal,logged=False)
sfc = georep.addDensitySlice(pdbCode,2,-1,8,0.1,central,linear,planar,interp='linear',palette=pal,logged=True)
#sfc = georep.addDensitySlice(pdbCode,2,-1,4,0.02,central,linear,planar,interp='nearest',palette=pal)
#sfc = georep.addDensitySlice(pdbCode,2,-1,4,0.02,central,linear,planar,interp='linear',palette=pal,logged=True)
#sfc = georep.addDensitySlice(pdbCode,2,-1,4,0.02,central,linear,planar,interp='nearest',palette=pal,logged=True)
#sfc = georep.addDensitySlice(pdbCode,2,-1,2,0.2,central,linear,planar,interp='linear',palette=pal)
#sfc = georep.addDensitySlice(pdbCode,2,-1,2,1,central,linear,planar,interp='linear',palette=pal)
#sfc = georep.addDensitySlice(pdbCode,2,-1,4,0.2,central,linear,planar,interp='nearest',palette=pal)
#sfc = georep.addDensitySlice(pdbCode,2,-1,2,0.2,central,linear,planar,interp='linear',palette=pal)
#sfc = georep.addDensitySlice(pdbCode,2,-1,2,0.2,central,linear,planar,interp='nearest',palette=pal)
#sfc = georep.addDensitySlice(pdbCode,2,-1,1,0.2,central,linear,planar,interp='linear',palette=pal)
#sfc = georep.addDensitySlice(pdbCode,2,-1,1,0.2,central,linear,planar,interp='nearest',palette=pal)

georep.printToHtml(pdbCode.upper(), 4, pdbCode + '_csharp')

analyser = pdb_eda.densityAnalysis.fromPDBid(pdbCode)

width=1
slice = space.GeoSpace()
squares = 3
slice = slice.getSquare(squares, 1, central, linear, planar)


print('The xyz coords')
print('------------------------------')
for i in range(len(slice[0])-1,-1,-1):
    for j in range(0,len(slice[0])):
        print('(',end='')
        print(round(slice[0][i,j],2),end=",")
        print(round(slice[1][i,j],2),end=",")
        print(round(slice[2][i,j],2),end=",")
        if j == len(slice[0]) - 1:
            print(')')
        else:
            print(')', end='')

print('')
print('The CRS coords')
print('------------------------------')
for i in range(len(slice[0])-1,-1,-1):
    for j in range(0,len(slice[0])):
        c, r, s = analyser.densityObj.header.xyz2crsCoord([slice[0][i,j],slice[1][i,j],slice[2][i,j]])
        print('(',end='')
        print(round(c,0),end=",")
        print(round(r,0),end=",")
        print(round(s,0),end=",")
        if j == len(slice[0]) - 1:
            print(')')
        else:
            print(')', end='')

print('')
print('The calculated values coords')
print('------------------------------')
for i in range(len(slice[0])-1,-1,-1):
    for j in range(0,len(slice[0])):

        print(round(sfcp[i,j],4),end=",")
        if j == len(slice[0]) - 1:
            print('')
        else:
            print('', end='')

#print(sfc)

