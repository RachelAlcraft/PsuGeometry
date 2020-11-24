

import subprocess as sub
# Testing some values for interpolation
import pdb_eda

from PsuGeometry import GeoTransformation as trans
from PsuGeometry import GeoReport as geor

pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
matDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/mat_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

pdbCode = '1us0'
x,y,z = 4.98,6.17,-0.14

pdb_eda.densityAnalysis.ccp4folder = edDataPath
pdb_eda.densityAnalysis.pdbfolder = pdbDataPath
analyser = pdb_eda.densityAnalysis.fromPDBid(pdbCode)
transformer = trans.transformation([0,0,0],[1,1,1],[0,1,2])

x,y,z = 33.15, 2.66, 22.43
c,r,s = analyser.densityObj.header.xyz2crsCoord([x,y,z])
den = analyser.densityObj.density[s,r,c]
print("ca",c,r,s,den)
x,y,z = 32.21, 1.83, 21.81
c,r,s = analyser.densityObj.header.xyz2crsCoord([x,y,z])
den = analyser.densityObj.density[s,r,c]
print("li",c,r,s,den)
x,y,z = 32.3, 0.47, 22.08
c,r,s = analyser.densityObj.header.xyz2crsCoord([x,y,z])
den = analyser.densityObj.density[s,r,c]
print("pl",c,r,s,den)

c,r,s = 147,226,137
den = analyser.densityObj.density[s,r,c]
print("ca2",c,r,s,den)
c,r,s = 144,226,133
den = analyser.densityObj.density[s,r,c]
print("li2",c,r,s,den)
c,r,s = 145,222,127
den = analyser.densityObj.density[s,r,c]
print("pl2",c,r,s,den)


matrix_dimensions = str(x) + "_" + str(y) + "_" + str(z)
matrix_list = ""
coord_list = ""


print(len(analyser.densityObj.densityArray))
print(analyser.densityObj.density.shape)
print(analyser.densityObj.density[0,0,0])
print(analyser.densityObj.density[0,1,0])
print(analyser.densityObj.density[1,0,0])
print(analyser.densityObj.density[1,1,1])
x,y,z = analyser.densityObj.header.crs2xyzCoord([1,2,3])
print(x,y,z)
c,r,s = analyser.densityObj.header.xyz2crsCoord([x,y,z])
print(c,r,s)



georep = geor.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath)
sfc = georep.addDensitySlice(pdbCode,2,-1,2,0.5,[33.15, 2.66, 22.43],[32.21, 1.83, 21.81],[32.3, 0.47, 22.08],interp='NEAREST')
print(sfc)

'''
if False:
    detailsName = matDataPath + pdbCode + "_details.csv"
    f = open(detailsName, "w")
    f.write('Details,Values' + '\n')
    dims = 'Dimensions,' + matrix_dimensions
    f.write(dims + '\n')
    f.close()

if False: #we only need to savwe the matrix data once

    matName = matDataPath + pdbCode +val "_vals.csv"
    f = open(matName, "w")
    f.write('CRS,Values' + '\n')

    for a in range(0,x):
        for b in range(0, y):
            for c in range(0, z):
                abc = str(a) + '_' + str(b) + '_' + str(c)
                if b == 0 or c == 0:
                    print(a)

                row = abc + ',' + str(analyser.densityObj.density[a, b, c])
                f.write(row + '\n')


    f.close()

points_dimensions = "2,2,1"
points_list = "1,1,1:2,2,2"

exeFile = "/home/rachel/Documents/Bioinformatics/BbkProject/Project2/CPlusPlusCode/CMakeProjects/ElectronDensity/build/default/MakeElectronDensityx,y,z = analyser.densityObj.density.shape

pig = sub.Popen([exeFile,'SLICE',pdbCode, matDataPath,points_dimensions,points_list],stdout=sub.PIPE)
result = pig.communicate(input=b"this is sample text.\n")
print(result[0])
'''