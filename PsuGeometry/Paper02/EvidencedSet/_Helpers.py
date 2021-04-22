
import pandas as pd
from PsuGeometry import GeoReport as psu

def applyCis(aa,preomega):
    if aa != 'PRO':
        return aa
    if abs(preomega) > 100:
        return 'PRO'
    else:
        return 'CISPRO'



def getCsv(pdbSet, geos,reloadPdb, reloadCsv,aa='ALL',includeCis=False,allAtoms=False, bFactorFactor=1.3,cutoff=0):

    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_out/' + pdbSet + '/'
    edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
    loadPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataB/'
    printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataK/'

    fileName = 'Data_DefensibleWithGeosALL_' + pdbSet + '.csv'

    if reloadCsv:
        from PsuGeometry import GeoPdb as geopdb
        pdbmanager = geopdb.GeoPdbs(pdbDataPath, edDataPath, False, False, False)
        if reloadPdb:
            pdbmanager.clear()

        pdbdata = pd.read_csv('../../PdbLists/Pdbs_Evidenced.csv')  # This is a list of pdbs <= 1.1A non homologous to 90%
        pdbListIn = pdbdata['PDB'].tolist()[0:]
        if cutoff > 0:
            pdbListIn = pdbdata['PDB'].tolist()[0:cutoff]


        pdbList = []
        for pdb in pdbListIn:
            import os.path
            if os.path.isfile((pdbDataPath + 'pdb' + pdb + '.ent').lower()):
                pdbList.append(pdb.lower())
            else:
                print('No file:', (pdbDataPath + 'pdb' + pdb + '.ent').lower())

        pdbList.sort()

        hueList = ['aa', 'rid', 'bfactor', 'pdbCode', 'bfactorRatio', 'disordered']
        georep = psu.GeoReport(pdbList, pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=allAtoms)

        if includeCis:
            geos.append('CA-1:C-1:N:CA')

        dataBest = georep.getGeoemtryCsv(geos, hueList, bFactorFactor, allAtoms=allAtoms,restrictedAa=aa)
        dataBest['rid'] = dataBest['rid'].astype(str)
        dataBest['ID'] = dataBest['pdbCode'] + dataBest['chain'] + dataBest['rid'] + dataBest['aa']
    else:
        dataBest = pd.read_csv(loadPath + fileName)

    #aas = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG','SER', 'THR', 'VAL', 'TRP', 'TYR']
    if includeCis:
        dataBest['aa'] = dataBest.apply(lambda row: applyCis(row['aa'], row['CA-1:C-1:N:CA']), axis=1)

    if aa !='ALL':
        dataBest = dataBest.query('aa == "' + aa + '"')

    return dataBest