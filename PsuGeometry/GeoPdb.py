
import Bio.PDB as bio
from Bio.PDB.DSSP import DSSP
import pandas as pd
import numpy as np

from PsuGeometry import GeoAtom as atm
from PsuGeometry import GeoDensity as den
from PsuGeometry import GeoCalcs as calcs


class GeoPdb:
    def __init__(self,pdbCode,pdbDataPath,edDataPath):
        pdbCode = pdbCode.lower()
        print('PSU: init',pdbCode)
        self.pdbCode = pdbCode
        self.pdbDataPath= pdbDataPath
        self.hasDensity = False
        self.hasDSSP = False
        self.hasPDB = False
        self.atoms = []
        self.geoDen = den.GeoDensity(pdbCode,'fifty',pdbDataPath,edDataPath)
        self.hasDensity = self.geoDen.valid
        self.dataFrame = None

        if self.__gatherAtoms():
            self.__applyDssp()
            #print("Gathered atoms")
            self.dataFrame = pd.DataFrame(columns=('pdbCode', 'resolution',
                                                  'chain', 'rid', 'dssp', 'aa',
                                                  'atom', 'atomNo', 'electrons','element', 'x', 'y', 'z','bfactor','occupant', 'occupancy',
                                                  '2FoFc', 'FoFc', 'Fo', 'Fc'))
            for atom in self.atoms:
                nextrow = len(self.dataFrame)
                self.dataFrame.loc[nextrow] = (atom.values['pdbCode'], atom.values['resolution'],
                                               atom.values['chain'], atom.values['rid'], atom.values['dssp'], atom.values['aa'],
                                               atom.values['atom'], atom.values['atomNo'], atom.values['electrons'], atom.values['element'],atom.values['x'], atom.values['y'], atom.values['z'], atom.values['bfactor'], atom.values['occupant'],atom.values['occupancy'],
                                               atom.values['2FoFc'], atom.values['FoFc'], atom.values['Fo'], atom.values['Fc'])  # switching ijk to crs

    #########################################################################################################################
    ## PRIVATE FUNCTIONS FOR THE CLASS
    #########################################################################################################################
    def __gatherAtoms(self):
        # try:
        if True:
            self.hasPDB = True
            pdbCode = self.pdbCode.lower()
            parser = bio.PDBParser()
            biodl = bio.PDBList()
            structure = None
            try:
                structure = parser.get_structure(pdbCode, self.pdbDataPath + 'pdb' + pdbCode + '.ent')
            except:
                biodl.download_pdb_files([pdbCode], pdir=self.pdbDataPath, file_format='pdb')
                structure = parser.get_structure(pdbCode, self.pdbDataPath + 'pdb' + pdbCode + '.ent')
            resolution = structure.header['resolution']
            atomNo = 0
            for model in structure:
                for chain in model:
                    for residue in chain:
                        r = residue.get_resname()
                        # print('Residue:', r)
                        rid = residue.get_full_id()[3][1]
                        chain = residue.get_full_id()[2]
                        if r != 'HOH':  # bio.is_aa(residue):
                            for atom in residue:
                                # print('Atom:', atom)
                                oneAtom = atm.GeoAtom()
                                oneAtom.setStructureInfo(pdbCode, resolution)
                                oneAtom.setResidueInfo(chain, rid, r)
                                atomNo += 1
                                name = atom.get_name()
                                occupant = atom.get_full_id()[4][1]
                                if occupant == ' ':
                                    occupant = 'A'
                                x = atom.get_vector()[0]
                                y = atom.get_vector()[1]
                                z = atom.get_vector()[2]
                                bfactor = atom.get_bfactor()
                                occupancy = atom.get_occupancy()
                                oneAtom.setAtomInfo(name, atomNo, x, y, z, bfactor, occupant, occupancy)
                                # add density if we can
                                if self.hasDensity:
                                    tFoFc, FoFc, Fo, Fc = self.geoDen.getDensityXYZ(x, y, z)
                                    oneAtom.setDensityInfo(tFoFc, FoFc, Fo, Fc)

                                # print('Atom:',atomNo)
                                self.atoms.append(oneAtom)

        # except:
        #    self.hasPDB = False
        return (self.hasPDB)

    def __applyDssp(self):
        print('PSU: apply dssp')
        p = bio.PDBParser()
        pdbFile = self.pdbDataPath + 'pdb' + self.pdbCode + '.ent'
        structure = p.get_structure(self.pdbCode, pdbFile)
        model = structure[0]
        dssp = DSSP(model, pdbFile)
        for akey in list(dssp.keys()):
            chain = akey[0]
            res_no = akey[1][1]
            row = dssp[akey]
            ss = row[2]
            for atom in self.atoms:
                if atom.values['rid'] == res_no and atom.values['chain'] == chain:
                    atom.setDsspInfo(ss)

    def getStructureCsv(self):
        return (self.data)

    def getGeoemtryCsv(self,geoList, hues):
        # geo in format C-1, C+1, C
        # remove anything that is in anyway
        if 'rid' in hues:
            hues.remove('rid')
        if 'pdbCode' in hues:
            hues.remove('pdbCode')
        if 'chain' in hues:
            hues.remove('chain')
        allAtomsA = self.__getAtomsOccupant('A',self.atoms)

        chainList = self.__getChainsUnique(allAtomsA)
        ridList = self.__getRidUnique(allAtomsA)
        rows = len(ridList)
        chrows = len(chainList)

        # set up the geoData to which we concatenate first
        geoData = pd.DataFrame(columns=('pdbCode', 'chain', 'rid'))
        for hue in hues:
            geoData[hue] = np.nan
        for geo in geoList:
            geoData[geo] = np.nan

        for ch in range(0, chrows):
            thisChain = chainList[ch]
            allAtomsChain = self.__getAtomsChain(thisChain,allAtomsA)

            for r in range(0, rows):
                thisResid = ridList[r]
                allValid = True
                listCalcs = []

                for geo in geoList:
                    geos = geo.split(':')
                    geoPairs = self.__geosToPairs(geos)

                    datasA = []
                    for a in range(0, len(geoPairs)):
                        geoPair = geoPairs[a]
                        geoAtom = geoPair[0]
                        ridA = thisResid + geoPairs[a][1]  # add the offset
                        allAtomsRid = self.__getAtomsRid(ridA, allAtomsChain)
                        allAtomsAtom = self.__getAtomsAtom(geoAtom, allAtomsRid)
                        # There shold now be ONLY 1 atom
                        if len(allAtomsAtom) == 1:
                            datasA.append(allAtomsAtom[0])
                        else:
                            allValid = False
                    listCalcs.append([datasA,geo])


                if allValid:
                    #add a new row to the dataframe
                    df1 = pd.DataFrame([[np.nan] * len(geoData.columns)], columns=geoData.columns)
                    geoData = df1.append(geoData, ignore_index=True)
                    thisRow = 0#len(geoData)-1
                    geoData.loc[thisRow, 'pdbCode'] = self.pdbCode
                    geoData.loc[thisRow, 'chain'] = thisChain
                    geoData.loc[thisRow, 'rid'] = int(thisResid)


                    # add the main data to the data frame
                    reshues = {}
                    for hue in hues:
                        reshues[hue] = ''
                    for oneGeo in listCalcs:
                        datasA = oneGeo[0]
                        geo = oneGeo[1]
                        geoatoms = geo.split(':')
                        geoPairs = self.__geosToPairs([geoatoms])
                        gpCount = 0
                        for gp in geoPairs:
                            offset = geoPairs[0][1]
                            if offset == 0:
                                for hue in hues:
                                    oneHue = datasA[gpCount].values[hue]
                                    if reshues[hue] == '':
                                        try:
                                            float(oneHue)
                                            reshues[hue] = 0
                                        except:
                                           reshues[hue] = oneHue

                        if len(datasA) == 4:  # dihedral
                                valA = calcs.torsion(datasA[0].values['x'], datasA[0].values['y'], datasA[0].values['z'],
                                                     datasA[1].values['x'], datasA[1].values['y'], datasA[1].values['z'],
                                                     datasA[2].values['x'], datasA[2].values['y'], datasA[2].values['z'],
                                                     datasA[3].values['x'], datasA[3].values['y'], datasA[3].values['z'])
                                for hue in hues:
                                    aHue = datasA[0].values[hue]
                                    bHue = datasA[0].values[hue]
                                    cHue = datasA[0].values[hue]
                                    dHue = datasA[0].values[hue]
                                    try:
                                        float(aHue)
                                        thisHue = (aHue + bHue + cHue + dHue)/4
                                        reshues[hue] += thisHue
                                        if reshues[hue] != thisHue:
                                            reshues[hue] = reshues[hue]/2 # we want the average of all the atoms in the calculation
                                    except:
                                        reshues[hue] =reshues[hue]



                        elif len(datasA) == 3:  # angle
                                valA = calcs.angle(datasA[0].values['x'], datasA[0].values['y'], datasA[0].values['z'],
                                                     datasA[1].values['x'], datasA[1].values['y'], datasA[1].values['z'],
                                                     datasA[2].values['x'], datasA[2].values['y'], datasA[2].values['z'])
                                for hue in hues:
                                    aHue = datasA[0].values[hue]
                                    bHue = datasA[0].values[hue]
                                    cHue = datasA[0].values[hue]
                                    try:
                                        float(aHue)
                                        thisHue = (aHue + bHue + cHue)/3
                                        reshues[hue] += thisHue
                                        if reshues[hue] != thisHue:
                                            reshues[hue] = reshues[hue]/2 # we want the average of all the atoms in the calculation
                                    except:
                                        reshues[hue] =reshues[hue]

                        else:
                                valA = calcs.distance(datasA[0].values['x'], datasA[0].values['y'], datasA[0].values['z'],
                                                     datasA[1].values['x'], datasA[1].values['y'], datasA[1].values['z'])
                                for hue in hues:
                                    aHue = datasA[0].values[hue]
                                    bHue = datasA[0].values[hue]
                                    try:
                                        float(aHue)
                                        thisHue = (aHue + bHue)/2
                                        reshues[hue] += thisHue
                                        if reshues[hue] != thisHue:
                                            reshues[hue] = reshues[hue]/2 # we want the average of all the atoms in the calculation
                                    except:
                                        reshues[hue] =reshues[hue]

                        geoData.loc[thisRow, geo] = valA
                        # hue could be an average or an
                        for hue in hues:
                            geoData.loc[thisRow, hue] = reshues[hue]

        return geoData


    def __getAtomsRid(self,rid,atoms):
        newAtoms = []
        for atm in atoms:
            if atm.values['rid'] == rid:
                newAtoms.append(atm)
        return(newAtoms)

    def __getAtomsChain(self, chain,atoms):
        newAtoms = []
        for atm in atoms:
            if atm.values['chain'] == chain:
                newAtoms.append(atm)
        return (newAtoms)

    def __getAtomsOccupant(self, occ,atoms):
        newAtoms = []
        for atm in atoms:
            if atm.values['occupant'] == occ:
                newAtoms.append(atm)
        return (newAtoms)

    def __getAtomsAtom(self, atom,atoms):
        newAtoms = []
        for atm in atoms:
            if atm.values['atom'] == atom:
                newAtoms.append(atm)
        return (newAtoms)


    def __getChainsUnique(self, atoms):
        chains = []
        for atm in atoms:
            if atm.values['chain'] not in chains:
                chains.append(atm.values['chain'])
        return (chains)

    def __getRidUnique(self, atoms):
        vals = []
        for atm in atoms:
            if atm.values['rid'] not in vals:
                vals.append(atm.values['rid'])
        return (vals)

    def __geosToPairs(self,geos):
        # geoX in format C-1, C+1, C
        pairs = []
        for geo in geos:
            atomX = ''
            offX = ''
            pm = 0
            for alpha in geo:
                if alpha == '-':
                    pm = -1
                elif alpha == '+':
                    pm = 1
                elif pm == 0:
                    atomX += alpha
                else:  # it is a number offset
                    offX += alpha
            if pm != 0:
                offX = pm * int(offX)
            else:
                offX = 0
            pairs.append([atomX, offX])

        return (pairs)
