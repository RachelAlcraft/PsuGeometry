from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

import time

###### User Choices ######################################################
pdbCodes= ['1ejg','2cnq','1us0','6q53','6jvv','4rek']
pdbCodes= ['6jvv','6q53','4rek']
pdbCodes= ['1ejg','2cnq','5nqo','6jvv','6q53','4rek','1us0']
pdbCodes= ['1ejg']
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

# runs are: num of Fo; num of Fc; palette; centre on zero; logg image, differentiation
#numFo,numFc,palette,zero?log?interp,degree,diff
runs = [
[2,-1,'PuBuGn',False,True,'nearest',0,0],
[2,-1,'PuBuGn',False,True,'linear',1,0],
[2,-1,'PuBuGn',False,True,'spline',3,0],
[2,-1,'PuBuGn',False,True,'splinexyz',3,0]
]

length = 6
gaps = 0.05
central_atom = 'C'
linear_atom = 'N'
linear_offset = 1
planar_atom = 'O'
planar_offset = 0
restricted_aa = ''
excluded_aa = ''
restrictednumber = 0 # 0 = all of them


#########################################################################


start = time.time()

for pdbCode in pdbCodes:

    georep = geor.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath)
    if len(runs) > 0:
        georep.addDataView(pdbCode, geoX='x', geoY='y', palette='cubehelix_r', hue='2FoFc')
    if len(runs) > 1:
        georep.addDataView(pdbCode, geoX='y', geoY='z', palette='Spectral', hue='bfactor')
    if len(runs) > 2:
        georep.addDataView(pdbCode, geoX='z', geoY='x', palette='rainbow', hue='atomNo')
    if len(runs) > 3:
        georep.addDataView(pdbCode, geoX='x', geoY='z', palette='seismic', hue='FoFc', centre=True)
    if len(runs) > 4:
        georep.addDataView(pdbCode, geoX='y', geoY='z', palette='rainbow', hue='electrons')

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
    averages = {}
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
                        if c in averages:
                            averages[c].append(sfc)
                        else:
                            averages[c] = [sfc]
                        c = c+1

    for c in averages:
        avg = averages[c]
        palette = runs[c][2]
        zero = runs[c][3]
        interp = runs[c][5]
        den1 = georep.addDensitySlices(avg, palette=palette, title='Average ' + interp, logged=True, centre=zero)

    georep.printToHtml(pdbCode.upper() + ' Density Interp Methods', len(runs),  pdbCode + '_interps')

    end = time.time()
    time_diff = end - start
    timestring = str(int(time_diff/60)) + "m " + str(int(time_diff%60)) + "s"
    print('Time taken',timestring)




