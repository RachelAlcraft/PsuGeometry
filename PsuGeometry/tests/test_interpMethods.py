from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

import time

###### User Choices ######################################################
pdbCodes= ['1ejg','2cnq','1us0','6q53','6jvv','4rek']
pdbCodes= ['6jvv','6q53','4rek']
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

# runs are: num of Fo; num of Fc; palette; centre on zero; logg image, differentiation
#numFo,numFc,palette,zero?log?interp,degree,diff
runs = [
[2,-1,'cubehelix_r',False,False,'linear',0,0],
[2,-1,'cubehelix_r',False,False,'spline',1,0],
[2,-1,'cubehelix_r',False,False,'spline',3,0]
]

length = 6
gaps = 0.1
central_atom = 'C'
linear_atom = 'N'
linear_offset = 1
planar_atom = 'O'
planar_offset = 0
restricted_aa = ''
excluded_aa = ''
restrictednumber = 10 # 0 = all of them


#########################################################################


start = time.time()

for pdbCode in pdbCodes:

    georep = geor.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath)
    georep.addDataView(pdbCode, geoX='x', geoY='y', palette='cubehelix_r', hue='2FoFc')
    georep.addDataView(pdbCode, geoX='y', geoY='z', palette='Spectral', hue='bfactor')
    georep.addDataView(pdbCode, geoX='z', geoY='x', palette='rainbow', hue='atomNo')

    pdbmanager = geop.GeoPdbs(pdbDataPath, edDataPath)
    apdb = pdbmanager.getPdb(pdbCode)
    csv = apdb.getDataFrame()

    if restricted_aa != '':
        qry = 'aa =="' + restricted_aa + '"'
        csv = csv.query(qry)
    if excluded_aa != '':
        qry = 'aa =="' + excluded_aa + '"'
        csv = csv.query(qry)

    chains = csv['chain'].unique()
    rids = csv['rid'].unique()
    num = len(chains) * len(rids)
    count = 0
    for ch in chains:
        for rid in rids:
            print(count, '/', num)
            count = count + 1
            #if True:
            if restrictednumber == 0 or count <= restrictednumber:
                cqry = 'rid=="' + str(rid) + '" and atom=="' + central_atom + '"' + ' and chain=="' + ch + '"'
                lqry = 'rid=="' + str(
                    rid + linear_offset) + '" and atom=="' + linear_atom + '"' + ' and chain=="' + ch + '"'
                pqry = 'rid=="' + str(
                    rid + planar_offset) + '" and atom=="' + planar_atom + '"' + ' and chain=="' + ch + '"'

                centralxyz = csv.query(cqry)
                linearxyz = csv.query(lqry)
                planarxyz = csv.query(pqry)

                if len(centralxyz) > 0 and len(linearxyz) and len(planarxyz) > 0:
                    cx = round(centralxyz['x'].values[0], 3)
                    cy = round(centralxyz['y'].values[0], 3)
                    cz = round(centralxyz['z'].values[0], 3)
                    caa = centralxyz['aa'].values[0]

                    lx = round(linearxyz['x'].values[0], 3)
                    ly = round(linearxyz['y'].values[0], 3)
                    lz = round(linearxyz['z'].values[0], 3)
                    laa = linearxyz['aa'].values[0]

                    px = round(planarxyz['x'].values[0], 3)
                    py = round(planarxyz['y'].values[0], 3)
                    pz = round(planarxyz['z'].values[0], 3)
                    paa = planarxyz['aa'].values[0]

                    central, linear, planar = [cx, cy, cz], [lx, ly, lz], [px, py, pz]

                    c = 0
                    for run in runs:
                        Fos = run[0]
                        Fcs = run[1]
                        palette = run[2]
                        zero = run[3]
                        logged = run[4]
                        interp=run[5]
                        degree = run[6]
                        differ = run[7]

                        FoFc = str(Fos) + 'Fo'  + str(Fcs) + 'Fc'
                        title = interp + ' Degree=' + str(degree) + '\nDerivative ' + str(differ) + '\nSize=' + str(length) + 'Å Gaps=' + str(gaps) + 'Å\n' + ch + str(rid) + ' ' + FoFc + '\n' + 'c=' + caa + ':' + str(central) + '\n' + 'l=' + laa + ':' + str(linear) + '\n' + 'p=' + paa + ':' + str(planar)

                        sfc = georep.addDensitySlice(pdbCode,Fos,Fcs,length,gaps,central,linear,planar,palette=palette,title=title,logged=logged,centre=zero,interp=interp,differ=differ,degree=degree)
                        c = c+1

    georep.printToHtml(pdbCode.upper() + ' Density Interp Methods', 3,  pdbCode + '_interps')

    end = time.time()
    time_diff = end - start
    timestring = str(int(time_diff/60)) + "m " + str(int(time_diff%60)) + "s"
    print('Time taken',timestring)




