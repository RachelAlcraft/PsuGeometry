# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import random
import pandas as pd
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = True

makeCsv = False

###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/DensityFlight/'
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    edSlicePath = 'F:/Code/ProteinDataFiles/ccp4_out/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'


s1 = [ #The extremes about 0
['1qow','B',9],
['6evg','A',499],
['2car','B',82],
['1yqs','A',238],
['2gj3','A',76],
['3jyo','A',102],
['1gkm','A',54],
['2pvx','B',127],
['4e9s','A',436],
['3bvx','A',771]
]

s2 = [# wide N, +ve psi
['1oai','B',16],
['3nbc','A',43],
['6j93','A',137],
['3nbc','B',43],
['3nbc','A',129],
['3nbc','B',129],
['3su6','A',1012],
['1uwc','A',210],
['3pb6','X',232],
['5b8d','A',1154],
['6m80','D',3],
['3w07','A',102],
['3a1h','F',15],
['2x1p','C',105],
['1qow','B',10],
['3a1h','D',6]
]

s3 = [# wide N, -ve psi
['3g5s','A',134],
['5nq0','A',57],
['5y45','E',19],
['5x9l','A',144],
['5nwp','B',49],
['6jk4','A',107],
['5y0m','B',34],
['3i2z','B',57],
['1rg8','B',6],
['3nbc','A',5],
['1rg8','A',6],
]

s4 = [#short N psi < 5
['4co3','A',84],
['4co3','B',84],
['4mzd','A',503],
['3bvx','A',59],
['4mij','A',189],
['2bf6','A',672],
['1n62','E',77],
['3mvs','A',142],
['5d7w','A',398],
['4w7l','B',254],
['4qi8','B',191],
['6map','B',139],
['1n62','B',77],
['5aq0','B',417],
]

s5 = [#short N psi > 10
['3zuc','A',1],
['3eww','A',403],
['4e9s','A',436],
['3g5s','A',280],
['1l3k','A',56],
['3bvx','A',849],
['4izx','A',50],
['6gke','A',243],
['1n62','C',33],
['2cnq','A',123],
['1l3k','A',20],
['1lqt','A',85],
['1n62','C',230],
]


#Choose which slices list
slicesList = s5


# Create a csv file with the residues we are interested in
bigstring = ""

for sl in slicesList:
    georep = psu.GeoReport([sl[0]],pdbDataPath,edDataPath,printPath,ed=False,dssp=False)
    pdbmanager = geopdb.GeoPdbs(pdbDataPath,edDataPath,ed=False,dssp=False)
    apdb = pdbmanager.getPdb(sl[0],True)
    pdbcsv = apdb.getDataFrame()
    queryC = 'rid==' + str(sl[2]) + ' and chain=="' + sl[1] + '"' + ' and atom=="CA"'
    queryL = 'rid==' + str(sl[2]) + ' and chain=="' + sl[1] + '"' + ' and atom=="N"'
    queryP = 'rid==' + str(sl[2]) + ' and chain=="' + sl[1] + '"' + ' and atom=="C"'
    dataC = pdbcsv.query(queryC)
    dataL = pdbcsv.query(queryL)
    dataP = pdbcsv.query(queryP)

    if len(dataC) > 0 and len(dataL) > 0 and len(dataP)>0:
        cx = round(dataC['x'].values[0],3)
        cy = round(dataC['y'].values[0], 3)
        cz = round(dataC['z'].values[0], 3)
        lx = round(dataL['x'].values[0], 3)
        ly = round(dataL['y'].values[0], 3)
        lz = round(dataL['z'].values[0], 3)
        px = round(dataP['x'].values[0], 3)
        py = round(dataP['y'].values[0], 3)
        pz = round(dataP['z'].values[0], 3)

        row = sl[0] + "," + sl[1] + str(sl[2]) + "," + str(cx) + "," + str(cy) + "," + str(cz)
        row +=  "," + str(lx) + "," + str(ly) + "," + str(lz)
        row += "," + str(px) + "," + str(py) + "," + str(pz)

        print(row)

        bigstring += row + '\n'

print("########RESULTS#########")
print("")
print(bigstring)

