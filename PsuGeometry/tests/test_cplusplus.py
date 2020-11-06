

import subprocess as sub
# Testing some values for interpolation
import pdb_eda
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
matDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/mat_data/'
pdbCode = '1ejg'
x,y,z = 4.98,6.17,-0.14
pdb_eda.densityAnalysis.ccp4folder = edDataPath
pdb_eda.densityAnalysis.pdbfolder = pdbDataPath
analyser = pdb_eda.densityAnalysis.fromPDBid(pdbCode)
exeFile = "/home/rachel/Documents/Bioinformatics/BbkProject/Project2/CPlusPlusCode/CMakeProjects/ElectronDensity/build/default/MakeElectronDensity"
x,y,z = analyser.densityObj.density.shape
matrix_dimensions = str(x) + "_" + str(y) + "_" + str(z)
matrix_list = ""
coord_list = ""

if False:
    detailsName = matDataPath + pdbCode + "_details.csv"
    f = open(detailsName, "w")
    f.write('Details,Values' + '\n')
    dims = 'Dimensions,' + matrix_dimensions
    f.write(dims + '\n')
    f.close()

if False: #we only need to savwe the matrix data once

    matName = matDataPath + pdbCode + "_vals.csv"
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


pig = sub.Popen([exeFile,'SLICE',pdbCode, matDataPath,points_dimensions,points_list],stdout=sub.PIPE)
result = pig.communicate(input=b"this is sample text.\n")
print(result[0])
