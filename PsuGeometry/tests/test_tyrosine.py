from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

import time

###### User Choices ######################################################
pdbCodes= ['1ejg','2cnq','1us0','6q53','6jvv','4rek']
pdbCodes= ['1us0']
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

# runs are: num of Fo; num of Fc; palette; centre on zero; logg image, differentiation
runs = [
[2,-1,'cubehelix_r',False,True,0],
[2,-1,'seismic',True,False,1],
[2,-1,'seismic',True,False,2],
[2,-1,'seismic',True,False,3],
[2,-1,'seismic',True,False,4]
]

length = 8.5
gaps = 0.1
interpmethod = 'spline' #linear or nearest or spline or sphere
degree = 5
#tyr ring
central_atom = 'CE1'
linear_atom = 'CD2'
linear_offset = 0
planar_atom = 'OH'
planar_offset = 0
restricted_aa = 'TYR'

##
excluded_aa = ''


#########################################################################


start = time.time()

for pdbCode in pdbCodes:

    dens1 = []
    dens2 = []
    dens3 = []
    dens4 = []
    dens5 = []

    georep = geor.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath)
    georep.addDataView(pdbCode, geoX='x', geoY='y', palette='cubehelix_r', hue='2FoFc')
    georep.addDataView(pdbCode, geoX='y', geoY='z', palette='Spectral', hue='bfactor')
    georep.addDataView(pdbCode, geoX='z', geoY='x', palette='rainbow', hue='atomNo')
    if len(runs) > 3:
        georep.addDataView(pdbCode, geoX='x', geoY='z', palette='Spectral', hue='FoFc',centre=True)
    if len(runs) > 4:
        georep.addDataView(pdbCode, geoX='y', geoY='x', palette='rainbow', hue='aa',categorical=True)

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
    count = 1
    for ch in chains:
        for rid in rids:
            print(count, '/', num)
            count = count + 1
            if True:
            # if count < 10:
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
                        differ = run[5]

                        FoFc = str(Fos) + 'Fo'  + str(Fcs) + 'Fc'
                        title = 'Spline Degree=' + str(degree) + '\nDerivative ' + str(differ) + '\nSize=' + str(length) + 'Å Gaps=' + str(gaps) + 'Å\n' + ch + str(rid) + ' ' + FoFc + '\n' + 'c=' + caa + ':' + str(central) + '\n' + 'l=' + laa + ':' + str(linear) + '\n' + 'p=' + paa + ':' + str(planar)

                        sfc = georep.addDensitySlice(pdbCode,Fos,Fcs,length,gaps,central,linear,planar,palette=palette,title=title,logged=logged,centre=zero,interp=interpmethod,differ=differ,degree=degree)
                        if c==0:
                            dens1.append(sfc)
                        elif c==1:
                            dens2.append(sfc)
                        elif c==2:
                            dens3.append(sfc)
                        elif c==3:
                            dens4.append(sfc)
                        elif c==4:
                            dens5.append(sfc)
                        c = c+1

    den1 = georep.addDensitySlices(dens1, palette=runs[0][2], title='Average Derivative 0', logged=False, centre=runs[0][3])
    den2 = georep.addDensitySlices(dens2, palette=runs[1][2], title='Average Derivative 1', logged=False, centre=runs[1][3])
    den3 = georep.addDensitySlices(dens3, palette=runs[2][2], title='Average Derivative 2', logged=False, centre=runs[2][3])
    if len(dens4)>0:
        den4 = georep.addDensitySlices(dens4, palette=runs[3][2], title='Average Derivative 3', logged=False,centre=runs[3][3])
    if len(dens5) > 0:
        den5 = georep.addDensitySlices(dens5, palette=runs[4][2], title='Average Derivative 4', logged=False,centre=runs[4][3])


    georep.printToHtml(pdbCode.upper() + ' Peptide Bond Density Derivatives ' + interpmethod, len(runs),  pdbCode + '_' +interpmethod + '_' + str(degree) + 'degree')

    end = time.time()
    time_diff = end - start
    timestring = str(int(time_diff/60)) + "m " + str(int(time_diff%60)) + "s"
    print('Time taken',timestring)




