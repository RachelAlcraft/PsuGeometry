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

s6 = [#biggest tau variation for better at 3 degree 0.25 gaps, top 25
['1p1x','B',1193],
['2nrl','A',53],
['4ga2','A',127],
['4r5r','B',59],
['1k5c','A',143],
['4ga2','A',143],
['4rj2','D',81],
['1p1x','A',41],
['4r5r','B',50],
['4rj2','F',44],
['1m1q','A',82],
['4rj2','E',89],
['4r5r','A',24],
['4mtu','A',21],
['4a02','A',46],
['6j93','A',109],
['2g6f','X',41],
['6rk0','A',51],
['5xqv','A',31],
['1f94','A',7],
['2jfr','A',38],
['1k5c','A',305],
['6rk0','B',141],
]

s7 = [#no tau variation for better at 3 degree 0.25 gaps
['6kfn','A',309],
['6kfn','A',329],
['1xmk','A',341],
['4uyr','A',34],
['2y78','A',26],
['2y78','A',35],
['2y78','A',58],
['2y78','A',60],
['2y78','A',77],
['2y78','A',90],
['2y78','A',93],
['4bct','A',82],
['4bct','A',138],
['1vbw','A',17],
['4iau','A',53],
['4iau','A',95],
['4iau','A',112],
['1ucs','A',52],
['6mu9','A',192],
['6mu9','A',253],
['6eio','A',63],
['6eio','A',96],
['6eio','A',207],
['6eio','A',214],
['6eio','A',232],
]

s8 = [#highest res largest diffs
['2jfr','A',38],
['4y9w','A',21],
['2jfr','A',177],
['4y9w','A',205],
['2jfr','A',175],
['2o7a','A',113],
['1pq7','A',193],
['2ixt','A',56],
['2jfr','A',162],
['5kwm','A',130],
['4y9w','A',11],
['5kwm','A',131],
['2pwa','A',9],
['2jfr','A',221],
['4y9w','A',207],
]



Cat1 = [#best supported set Category 1: Long N:N+1 and psi near +/- 180
['4r2x','F',90],
['5otn','A',26],
['5xvt','A',472],
['4u9h','S',119],
['5xvt','A',475],
['4nsv','B',102],
['4a7u','A',108],
['3fil','B',9],
['5emb','A',63],
['4g9s','A',15],
['5f82','B',40],
['6ri6','A',39],
['3x34','A',82],
['6eio','A',96],
['3fsa','A',9],
['3vn3','A',136],
['2xjp','A',142],
['2e4t','A',290],
['6fmc','A',314],
['4nsv','A',143],
]


Cat2 = [#best supported set Category 2: Between values: -100< PSI <100 with N:N+1 between 3-3.4A
['6ri6','A',79], #There are only 12
['1r6j','A',210],
['5o2x','A',118],
['4u9h','L',481],
['2xfr','A',296],
['1pjx','A',93],
['4e3y','A',1121],
['4iau','A',112],
['3zoj','A',222],
['1muw','A',102],
['2gud','B',66],
['2gud','B',26],
]

Cat3 = [#best supported set Category 3: 0 psi and the bottom of the curve, tau<114
['2ddx','A',205],
['2e4t','A',364],
['5xvt','A',424],
['5xvt','A',251],
['1ug6','A',46],
['3qr7','A',123],
['5kwm','A',177],
['2pwa','A',214],
['6eio','A',63],
['5a71','A',41],
['3qr7','A',166],
['4kqp','A',326],
['6q4g','A',259],
['1xmk','A',341],
['4qb3','A',108],
['6q4g','A',114],
['4y9v','A',381],
['4y9v','A',589],
['4rek','A',288],
['2v8t','B',186],
]

Cat4 = [#best supported set Category 4: Shortest N:N+1 < 2.8 at tau > 114
['2ddx','A',205],
['4u9h','L',230],
['3q8j','A',20],
['2wfi','A',77],
['4ua6','A',156],
['2xjp','A',213],
['4y9w','A',166],
['4wpk','A',78],
['6fmc','A',331],
['4iau','A',53],
['4kqp','A',375],
['2e4t','A',352],
['4g9s','B',79],
['3m5q','A',280],
['4r2x','D',162],
['4r2x','B',179],
['3wgx','C',295],
['4nsv','A',195],
['4q4g','X',1059],
['2pwa','A',66],
['4ekf','A',16],
]

Cat5 = [#best supported set Category 5: 2.8 < N:N+1 < 3 at tau > 114
['1muw','A',145],#only 23
['5x9l','A',27],
['1v0l','A',75],
['2zpm','A',303],
['4g78','A',87],
['4ea9','A',15],
['1muw','A',138],
['6kfn','A',180],
['1c75','A',38],
['1n9b','A',281],
['4unu','A',31],
['5xvt','A',161],
['2nrl','A',12],
['5xvt','A',125],
['3zoj','A',230],
['5xvt','A',33],
['3vla','A',376],
['5u3a','A',205],
['6j93','A',25],
['5tif','A',75],
]

Cat6 = [#best supported set Category 6: Psi <25 N:N+1 > 2.8
['5akr','A',268],
['2zq7','A',217],
['4awt','A',242],
['5akr','A',315],
['4r2x','D',82],
['4wka','A',254],
['4o6u','A',176],
['2pve','A',10],
['1rb9','A',43],
['4r2x','C',82],
['2v8t','A',31],
['2v8t','B',228],
['2ddx','A',144],
['5xvt','A',470],
['3wdn','A',94],
['5otn','A',28],
['5xvt','A',209],
['5xvt','A',117],
['2bw4','A',225],
['1bxo','A',168],
]


#Choose which slices list
slicesList = Cat6


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

