
#import Bio.PDB as bio
import pandas as pd
import numpy as np

from PsuGeometry import GeoAtom as atm
from PsuGeometry import GeoDensity as den
from PsuGeometry import GeoCalcs as calcs


'''
singleton object to manage only loading pdbs once
https://python-3-patterns-idioms-test.readthedocs.io/en/latest/Singleton.html
'''
class GeoPdbs:
    class __GeoPdbs:
        def __init__(self,pdbDirectory,edDirectory,ed=True,dssp=True):
            self.pdbs = {}
            self.pdbDirectory = pdbDirectory
            self.edDirectory = edDirectory
            self.ed = ed
            self.dssp=dssp
        def __getPdb__(self,pdbCode):
            return self.pdbs[pdbCode]
        def __existsPdb__(self,pdbCode):
            return pdbCode in self.pdbs
        def __addPdb__(self,pdbCode,pdb):
            self.pdbs[pdbCode] = pdb
    instance=None
    def __init__(self,pdbDirectory,edDirectory,ed=True,dssp=True):
        if not GeoPdbs.instance:
            GeoPdbs.instance = GeoPdbs.__GeoPdbs(pdbDirectory,edDirectory,ed,dssp)
        #else:
        #    GeoPdbs.instance.pdbDirectory = pdbDirectory
        #    GeoPdbs.instance.edDirectory = edDirectory
        #    GeoPdbs.instance.ed = ed
        #    GeoPdbs.instance.dssp = dssp

    def existsPdb(self, pdbCode):
        pdbCode = pdbCode.lower()
        return self.instance.__existsPdb__(pdbCode)

    def getPdb(self, pdbCode,useAll):
        pdbCode = pdbCode.lower()
        if self.instance.__existsPdb__(pdbCode):
            return self.instance.__getPdb__(pdbCode)
        else:
            gp = GeoPdb(pdbCode,self.instance.pdbDirectory,self.instance.edDirectory,self.instance.ed,self.instance.dssp,useAll)
            self.instance.__addPdb__(pdbCode,gp)
            return gp


class GeoPdb:
    def __init__(self,pdbCode,pdbDataPath,edDataPath,ed,dssp,useAll):
        pdbCode = pdbCode.lower()
        self.pdbCode = pdbCode
        self.pdbDataPath= pdbDataPath
        self.hasDensity = False
        self.hasPDB = False
        self.atoms = []
        self.densCSV = pd.DataFrame()
        self.hasDssp = dssp
        self.dataFrame = pd.DataFrame()
        self.ghost = False
        self.useAll = useAll
        if self.pdbCode == 'ghost':
            self.ghost = True
            #self.pdbCode =  '2q1j'
            self.pdbCode = '4rek'
            self.hasDensity = False
            self.hasDssp = False
        else:
            if ed:
                self.geoDen = den.GeoDensity(pdbCode, 'fifty', pdbDataPath, edDataPath)
                self.hasDensity = self.geoDen.valid
            else:
                self.hasDensity = False


        if self.__gatherAtoms():
            if self.hasDssp:
                self.__applyDssp()
            #self.createDataStructure()

        if self.ghost == True:
            self.pdbCode = 'ghost'

    def createDataStructure(self):
        print('PSU: create data structure',self.pdbCode)
        dicdfs = []
        for atom in self.atoms:
            dic={   'pdbCode':atom.values['pdbCode'],'resolution':atom.values['resolution'],
                    'chain':atom.values['chain'], 'rid':atom.values['rid'],'ridx':atom.values['ridx'],
                    'dssp':atom.values['dssp'], 'aa':atom.values['aa'],
                    'atom':atom.values['atom'], 'atomNo':atom.values['atomNo'],
                    'electrons':atom.values['electrons'], 'element':atom.values['element'],
                    'x':atom.values['x'], 'y':atom.values['y'], 'z':atom.values['z'],
                    'bfactor':atom.values['bfactor'], 'occupant':atom.values['occupant'],
                    'occupancy':atom.values['occupancy'],
                    '2FoFc':atom.values['2FoFc'], 'FoFc':atom.values['FoFc'],
                    'Fo':atom.values['Fo'], 'Fc':atom.values['Fc']}
            dicdfs.append(dic)
        self.dataFrame = pd.DataFrame.from_dict(dicdfs)

    def getDataFrame(self):
        #if self.dataFrame == None:
        if self.dataFrame.empty:
            self.createDataStructure()
        return self.dataFrame

    def getDensitySquare(self,squares,Fos,Fcs,interp,differ,degree):
        xsq = squares[0]
        ysq = squares[1]
        zsq = squares[2]
        x,y = xsq.shape
        squ = np.zeros((x,y))
        for i in range(0,x):
            for j in range(0, y):
                a,b,c = xsq[i,j],ysq[i,j],zsq[i,j]
                den = self.geoDen.getInterpolatedDensity(a,b,c,Fos,Fcs,interp,differ,degree)
                squ[i,j] = den
        return squ






    #########################################################################################################################
    ## PRIVATE FUNCTIONS FOR THE CLASS
    #########################################################################################################################
    def __gatherAtoms(self):
        # try:
        if True:
            import Bio.PDB as bio
            self.hasPDB = True
            pdbCode = self.pdbCode.lower()
            print('PSU: load from BioPython', self.pdbCode)
            parser = bio.PDBParser()
            biodl = bio.PDBList()
            structure = None
            try:
                structure = parser.get_structure(pdbCode, self.pdbDataPath + 'pdb' + pdbCode + '.ent')
            except:
                import time
                biodl.download_pdb_files([pdbCode], pdir=self.pdbDataPath, file_format='pdb')
                time.sleep(1)
                try:
                    structure = parser.get_structure(pdbCode, self.pdbDataPath + 'pdb' + pdbCode + '.ent')
                except:
                    import time
                    time.sleep(10)
                    structure = parser.get_structure(pdbCode, self.pdbDataPath + 'pdb' + pdbCode + '.ent')

            resolution = structure.header['resolution']
            atomNo = 0
            resnum = 1
            for model in structure:
                for chain in model:
                    for residue in chain:
                        r = residue.get_resname()
                        # print('Residue:', r)
                        rid = residue.get_full_id()[3][1]
                        chain = residue.get_full_id()[2]
                        hetatm = residue.get_full_id()[3][0]
                        ridx = resnum
                        resnum = resnum+1
                        #decision as to whether r is to be used. for density maps yes, for geoemtry no
                        if (r in self.getAAList() and 'H' not in hetatm) or (self.useAll and r!='HOH'):# != 'HOH':  # bio.is_aa(residue):
                            for atom in residue:
                                if atom.is_disordered():
                                    if atom.disordered_has_id("A"):
                                        atom.disordered_select("A")
                                # print('Atom:', atom)
                                oneAtom = atm.GeoAtom()
                                oneAtom.setStructureInfo(pdbCode, resolution)
                                oneAtom.setResidueInfo(chain, rid, ridx,r)
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
                                #if rid < 3:
                                #    print(oneAtom)
                                # add density if we can
                                if self.hasDensity:
                                    tFoFc, FoFc, Fo, Fc = self.geoDen.getDensityXYZ(x, y, z)
                                    oneAtom.setDensityInfo(tFoFc, FoFc, Fo, Fc)

                                # print('Atom:',atomNo)
                                self.atoms.append(oneAtom)
            print('PSU: loaded successfully from BioPython', self.pdbCode)


        # except:
        #    self.hasPDB = False
        return (self.hasPDB)

    def __applyDssp(self):
        import Bio.PDB as bio
        print('PSU: applying dssp')
        from Bio.PDB.DSSP import DSSP
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
        print('PSU: applied dssp successfully')


    def getStructureDensity(self,allPoints,divisor,pdbDataPath,edDataPath):
        if self.hasDensity:
            if self.densCSV.empty:
                self.geoDen = den.GeoDensity(self.pdbCode,'fifty',pdbDataPath,edDataPath)
                self.densCSV = self.geoDen.getPeaks(allPoints,divisor)
        return self.densCSV


    def getGeoemtryCsv(self,geoListEntered, hues):
        # geo in format C-1, C+1, C
        #print('PSU: creating geometry dataframe')
        dics = []
        usingAliases = False
        # remove anything that is in anyway
        if 'rid' in hues:
            hues.remove('rid')
        if 'pdbCode' in hues:
            hues.remove('pdbCode')
        if 'chain' in hues:
            hues.remove('chain')

        geoList = []
        geoListIn = []
        for geoa in geoListEntered:
            for aa in self.getAAList():
                geo = self.aliasToGeo(geoa,aa)
                if geo != geoa:
                    usingAliases = True
                #print(geoa,geo,aa,usingAliases)
                if ':' not in geo:
                    if geo not in hues:
                        hues.append(geo)
                else:
                    if geo not in geoList:
                        geoList.append(geo)
                    if geoa not in geoListIn:
                        geoListIn.append(geoa)

        if len(geoList)<2:
            geoList.append('N:CA')
            geoList.append('CA:C')
        if len(geoListIn)<2:
            geoListIn.append('N:CA')
            geoListIn.append('CA:C')

        if 'aa' not in hues:
            hues.append('aa')
        if 'ridx' not in hues:
            hues.append('ridx')
        if 'atomNo' not in hues:
            hues.append('atomNo')
        if 'bfactor' not in hues:
            hues.append('bfactor')

        occList = ['A']#self.__getOccList()
        ridList = self.__getRidList()
        chainList = self.__getChainList()

        rows = len(ridList)
        chrows = len(chainList)
        occs = len(occList)

        #an atom will be uniquely defined by rid, chain, occupant

        # set up the geoData to which we concatenate first
        #geoData = pd.DataFrame(columns=('pdbCode', 'chain', 'rid'))
        #for hue in hues:
        #    geoData[hue] = np.nan

        #for geo in geoListIn: #the column names will be the alias names or whatever we passed in AND the aliases
        #    geoData[geo] = np.nan
        #for geo in geoList: #the column names will be the alias names or whatever we passed in AND the aliases
        #    geoData[geo] = np.nan

        for ch in range(0, chrows):
            thisChain = chainList[ch]
            for occ in range(0,occs):
                thisOcc = occList[occ]
                for rid in range(0, rows):
                    thisResid = ridList[rid]
                    thisResidue = self.__getResidue(thisChain, thisResid,thisOcc)# not really a residue but it does for getting aa
                    if thisResidue == None:
                        #print(thisChain,thisResid,thisOcc)
                        a=2
                    else:
                        allValid = True
                        aa = thisResidue.values['aa']

                        listCalcs = []

                        for geoa in geoListIn:
                            geo = self.aliasToGeo(geoa,aa)
                            geos = geo.split(':')
                            geoPairs = self.__geosToPairs(geos)

                            datasA = []
                            for a in range(0, len(geoPairs)):
                                geoPair = geoPairs[a]
                                geoAtom = geoPair[0]
                                ridA = thisResid + geoPairs[a][1]  # add the offset
                                atomA = self.__getAtom(thisChain, ridA,thisOcc,geoAtom)
                                # There should be 1 atom
                                if atomA != None:
                                    datasA.append(atomA)
                                else:
                                    #print(thisChain, ridA, thisOcc,geoAtom)
                                    allValid = False
                            listCalcs.append([datasA,geo])


                        if allValid:
                            #add a new row to the dataframe
                            #df1 = pd.DataFrame([[np.nan] * len(geoData.columns)], columns=geoData.columns)
                            #geoData = df1.append(geoData, ignore_index=True)
                            thisRow = 0#len(geoData)-1
                            #geoData.loc[thisRow, 'pdbCode'] = self.pdbCode
                            #geoData.loc[thisRow, 'chain'] = thisChain
                            #geoData.loc[thisRow, 'rid'] = int(thisResid)
                            dic = {}
                            dic['pdbCode'] = self.pdbCode
                            dic['chain'] = thisChain
                            dic['rid'] = int(thisResid)


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

                                elif len(datasA) == 2:  # distance
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
                                else: # just some data
                                    print('??',datasA)

                                #geoData.loc[thisRow, geo] = valA
                                dic[geo] = valA
                                # hue could be an average or an
                                for hue in hues:
                                    #geoData.loc[thisRow, hue] = reshues[hue]
                                    dic[hue] = reshues[hue]
                                dic['aa'] = aa

                                #print(usingAliases)
                                if usingAliases:
                                    #aa = geoData['aa'][0]
                                    geoa = self.geoToAlias(geo,aa)
                                    #print(geoa,geo,aa)
                                    if geoa != geo:
                                        #geoData.loc[thisRow, geoa] = valA # we have alias and geo column
                                        dic[geoa] = valA


                            dics.append(dic)
        dataFrame = pd.DataFrame.from_dict(dics)

        return dataFrame


    def __getAtomsRid(self,rid,atoms):
        newAtoms = []
        for atm in atoms:
            if atm.values['rid'] == rid:
                newAtoms.append(atm)
        return(newAtoms)

    def __getAtomsChain(self, chain,atoms):
        newAtoms = []
        for atm in atoms:
            if atm.values['chain'] == chain:# and atm.values['aa'] == aa:
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

    def __getResidue(self, chain, rid, occ):
        for atm in self.atoms:
            if atm.values['chain'] == chain and atm.values['rid'] == rid and atm.values['occupant'] == occ:
                return atm
        return None

    def __getAtom(self, chain, rid, occ,atom):
        # The atom number cannot be less than 1
        if rid < 1:
            return None
        for atm in self.atoms:
            if atm.values['chain'] == chain and atm.values['rid'] == rid and atm.values['occupant'] == occ and atm.values['atom'] == atom:
                return atm
        return None

    def __getChainsUnique(self, atoms):
        chains = []
        for atm in atoms:
            if atm.values['chain'] not in chains:
                chains.append(atm.values['chain'])
        return (chains)

    def __getChainList(self):
        chains = []
        for atm in self.atoms:
            if atm.values['chain'] not in chains:
                chains.append(atm.values['chain'])
        return (chains)

    def __getRidList(self):
        rids = []
        for atm in self.atoms:
            if atm.values['rid'] not in rids:
                rids.append(atm.values['rid'])
        return (rids)

    def __getOccList(self):
        occs = []
        for atm in self.atoms:
            if atm.values['occupant'] not in occs:
                occs.append(atm.values['occupant'])
        return (occs)

    def __getRidUnique(self, atoms):
        vals = []
        for atm in atoms:
            if atm.values['rid'] not in vals:
                vals.append(atm.values['rid'])
        print(vals)
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

    def aliasToGeo(self,alias,aa):
        dic = self.getAliasDictionary()
        if alias + '_' + aa in dic:
            return dic[alias+'_'+aa]
        elif alias in dic:
            return dic[alias]
        else:
            return alias

    def geoToAlias(self,geo,aa):
        dic = self.getAliasDictionary()
        for a,g in dic.items():
            if aa in a and g == geo:
                if '_' in a:
                    return a.split('_')[0]
                else:
                    return a
        for a,g in dic.items():
            if g==geo:
                return a
        return geo

    def getAliasDictionary(self):
        return {
                'PHI':'C-1:N:CA:C',
                'PSI':'N:CA:C:N+1',
                'OMEGA': 'CA:C:N+1:CA+1',
                'TAU':'N:CA:C',
                'CHI1':'N:CA:CB:CG',
                'CHI1_ILE':'N:CA:CB:CG1',
                'CHI1_SER': 'N:CA:CB:OG',
                'CHI1_THR': 'N:CA:CB:OG1',
                'CHI1_VAL': 'N:CA:CB:CG1',
                'CHI1_ALA': 'N:CA:CB:HB1',
                'CHI2': 'CA:CB:CG:CD',
                'CHI2_ASN': 'CA:CB:CG:OD1',
                'CHI2_ASP': 'CA:CB:CG:OD1',
                'CHI2_HIS': 'CA:CB:CG:ND1',
                'CHI2_ILE': 'CA:CB:CG1:CD',
                'CHI2_LEU': 'CA:CB:CG:CD1',
                'CHI2_MET': 'CA:CB:CG:SD',
                'CHI2_PHE': 'CA:CB:CG:CD1',
                'CHI2_TRP': 'CA:CB:CG:CD1',
                'CHI2_TYR': 'CA:CB:CG:CD1',
                'CHI2_VAL': 'CA:CB:CG1:HG11',
                'CHI2_THR': 'CA:CB:CG2:HG21',
                'CHI3':'CB:CG:CD:CE',
                'CHI3_ARG': 'CB:CG:CD:NE',
                'CHI3_GLN': 'CB:CG:CD:OE1',
                'CHI3_GLU': 'CB:CG:CD:OE1',
                'CHI3_HIS': 'CA:CB:CG:CD2',
                'CHI3_MET': 'CB:CG:SD:CE',
                'CHI3_PRO': 'CB:CG:CD:N',
                'CHI3_VAL': 'CA:CB:CG2:HG21',
                'CHI4': 'CG:CD:CE:CZ',
                'CHI4_ARG': 'CG:CD:NE:CZ',
                'CHI4_PRO': 'CG:CD:N:CA',
                'CHI4_LYS': 'CG:CD:CE:NZ',
                'CHI5': 'CD:CE:CZ:NH1',
                'CHI5_PRO': 'CD:N:CA:CB',
                }
    def getAAList(self):
        return ['ALA','CYS','ASP','GLU','PHE',
                'GLY','HIS','ILE','LYS','LEU',
                'MET','ASN','PRO','GLN','ARG',
                'SER','THR','VAL','TRP','TYR']
