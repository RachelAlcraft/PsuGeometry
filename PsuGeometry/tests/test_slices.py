from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop
import math


pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/slices/'


pdbCode= '1ejg'
if True:

    georep = geor.GeoReport([pdbCode],pdbDataPath, edDataPath,printPath)
    georep.addDataView(pdbCode, geoX='x', geoY='y', palette='cubehelix_r')
    georep.addScatter(geoX='x', geoY='y', palette='rainbow', hue='rid')

    pdbmanager = geop.GeoPdbs(pdbDataPath, edDataPath)
    apdb = pdbmanager.getPdb(pdbCode)
    csv = apdb.getDataFrame()
    csv = csv.query('aa !="PRO"')
    chains = csv['chain'].unique()
    rids = csv['rid'].unique()
    slices = []
    for ch in chains:
        for rid in rids:
            centralxyz = csv.query('rid=="' + str(rid) + '" and atom=="C"'+ ' and chain=="' + ch + '"')
            linearxyz = csv.query('rid=="' + str(rid+1) + '" and atom=="N"' + ' and chain=="' + ch + '"')
            planarxyz = csv.query('rid=="' + str(rid) + '" and atom=="O"' + ' and chain=="' + ch + '"')
            if len(centralxyz) > 0 and len(linearxyz) and len(planarxyz) > 0:
                cx = round(centralxyz['x'].values[0],3)
                cy = round(centralxyz['y'].values[0],3)
                cz = round(centralxyz['z'].values[0],3)
                lx = round(linearxyz['x'].values[0],3)
                ly = round(linearxyz['y'].values[0],3)
                lz = round(linearxyz['z'].values[0],3)
                px = round(planarxyz['x'].values[0],3)
                py = round(planarxyz['y'].values[0],3)
                pz = round(planarxyz['z'].values[0],3)
                central,linear,planar = [cx,cy,cz],[lx,ly,lz],[px,py,pz]
                title = ch + str(rid) + '\n' + 'c=' + str(central) + '\n' + 'l=' + str(linear) + '\n' + 'p=' + str(planar)
                sfc = georep.addDensitySlice(pdbCode,60,0.05,central,linear,planar,palette='cubehelix_r',title=title)
                slices.append(sfc)
            # And finally create the reort with a file name of choice
    sfc = georep.addDensitySlices(slices, palette='cubehelix_r', title='Averaged')
    sfc = georep.addDensitySlices([sfc], palette='cubehelix_r', title='Log Averaged',norm=True)
    georep.printToHtml('Density Views and Slices', 3, pdbCode + '_densl')




