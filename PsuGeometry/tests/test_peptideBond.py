from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop

###### User Choices ######################################################
pdbCodes= ['1ejg','2cnq','1us0','6q53','6jvv','4rek']
#pdbCodes= ['4rek']
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'

# runs are num of Fo num of Fc and palette
runs = [
        [3,-2,'YlGnBu',False]
        #,[1,-1,'RdGy',True]
        ]

logged = True
length = 6
gaps = 0.1
central_atom = 'C'
linear_atom = 'N'
linear_offset = 1
planar_atom = 'O'
planar_offset = 0
restricted_aa = ''
excluded_aa = ''

#########################################################################

for pdbCode in pdbCodes:

    for run in runs:
        Fos = run[0]
        Fcs = run[1]
        palette = run[2]
        zero = run[3]

        surfaces = []
        surfacespro = []

        georep = geor.GeoReport([pdbCode],pdbDataPath, edDataPath,printPath)
        georep.addDataView(pdbCode, geoX='x', geoY='y', palette='cubehelix_r',hue='2FoFc')
        georep.addDataView(pdbCode, geoX='y', geoY='z', palette='Spectral',hue='bfactor')
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
        slices = []
        proslices = []
        num = len(chains)*len(rids)
        count = 1
        for ch in chains:
            for rid in rids:
                print(count, '/', num)
                count = count + 1
                cqry = 'rid=="' + str(rid) + '" and atom=="' + central_atom + '"'+ ' and chain=="' + ch + '"'
                lqry = 'rid=="' + str(rid+linear_offset) + '" and atom=="' + linear_atom + '"' + ' and chain=="' + ch + '"'
                pqry = 'rid=="' + str(rid+planar_offset) + '" and atom=="' + planar_atom + '"' + ' and chain=="' + ch + '"'

                centralxyz = csv.query(cqry)
                linearxyz = csv.query(lqry)
                planarxyz = csv.query(pqry)

                if len(centralxyz) > 0 and len(linearxyz) and len(planarxyz) > 0:

                    cx = round(centralxyz['x'].values[0],3)
                    cy = round(centralxyz['y'].values[0],3)
                    cz = round(centralxyz['z'].values[0],3)
                    caa = centralxyz['aa'].values[0]

                    lx = round(linearxyz['x'].values[0],3)
                    ly = round(linearxyz['y'].values[0],3)
                    lz = round(linearxyz['z'].values[0],3)
                    laa = linearxyz['aa'].values[0]

                    px = round(planarxyz['x'].values[0],3)
                    py = round(planarxyz['y'].values[0],3)
                    pz = round(planarxyz['z'].values[0],3)
                    paa = planarxyz['aa'].values[0]

                    central,linear,planar = [cx,cy,cz],[lx,ly,lz],[px,py,pz]
                    title = 'Size=' + str(length) + 'Å Gaps=' + str(gaps) + 'Å\n' + ch + str(rid) + '\n' + 'c='+ caa + ':' + str(central) + '\n' + 'l=' + laa + ':'+ str(linear) + '\n' + 'p=' + paa + ':'+ str(planar)
                    sfc = georep.addDensitySlice(pdbCode,Fos,Fcs,length,gaps,central,linear,planar,palette=palette,title=title,logged=logged)

                    if laa == 'PRO':
                        proslices.append(sfc)
                    else:
                        slices.append(sfc)
                # And finally create the reort with a file name of choice
        sfca = georep.addDensitySlices(slices, palette=palette, title='No Pro - Averaged',logged=False,centre=zero)
        sfcb = georep.addDensitySlices(slices, palette=palette, title='No Pro - Log Averaged',logged=True,centre=zero)

        sfcc = georep.addDensitySlices(proslices, palette=palette, title='Proline Next - Averaged', logged=False,centre=zero)
        sfcd = georep.addDensitySlices(proslices, palette=palette, title='Proline Next - Log Averaged', logged=True,centre=zero)
        surfaces.append([sfca,zero])
        surfacespro.append([sfcc,zero])
        matstring = str(Fos) + 'Fo' + str(Fcs) + 'Fc'
        georep.printToHtml(pdbCode.upper() + ' ' + matstring + ' Peptide Bond Density Slices', 3, pdbCode + matstring + '_dens')

    georep = geor.GeoReport([pdbCode], pdbDataPath, edDataPath, printPath)




