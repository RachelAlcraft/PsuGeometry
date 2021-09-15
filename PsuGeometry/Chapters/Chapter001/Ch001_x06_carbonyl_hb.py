'''
This script creates a file for hand chosen geos, so it is slower as they are no serialised.
This particular report is designed to look at hydrogen bonding on the carbonyl oxygen
'''

import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help

filesPDBRoot ='C:/Dev/Github/ProteinDataFiles/pdb_data/'
filesADJRoot ='C:/Dev/Github/ProteinDataFiles/pdb_out/Fov2_ADJ/' #adjusted on Fo at 3 degrees thevenaz
loadPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/'
printPath = 'C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/Data/'

##geos = ['CA:C:N+1','O:N+1', 'CA:O', 'C:O','CA:C', 'C:N+1','CA:O','N+1:C','O:{,N,}','O:{,HOH,}','O:{,ND1,}','O:{,NE2,}','O:{,NZ,}','O:{,ND2,}','O:{,NE2,}','O:{,NE,}','O:{,NH1,}','O:{,NH2,}','O:{,NE1,}']
#geos = ['CA:C:N+1','O:N+1', 'CA:O', 'C:O','CA:C', 'C:N+1','CA:O','N+1:C','O:{,N,HOH,}','O:{,HOH,}','O:{,N,}']
#geos = ['CA:C:N+1','O:N+1', 'CA:O', 'C:O','CA:C', 'C:N+1','CA:O','N+1:C','O:{,HOH,}']
geos = ['TAU','CA:C:N+1','C:O','O:{,HOH,}','O:{,HETATM,}','O:{,N,}','O:{,SG,}','O:{,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}','O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}']


title='Carbonyl Oxygen Report'
fileName = 'carbonyl'

createOrLoad = "LOAD"
if createOrLoad == "CREATE":
    pdbdata = pd.read_csv('../../PdbLists/Pdbs_70.csv')
    pdbListA = pdbdata['PDB'].tolist()[0:]
    pdbListIn = []
    for pdb in pdbListA:
        import os.path
        if os.path.isfile((filesADJRoot + 'pdb' + pdb + '.ent').lower()):
            pdbListIn.append(pdb.lower())
        else:
            print('No file:',(filesADJRoot + 'pdb' + pdb + '.ent').lower())
    print(pdbListIn)
    print("---- Getting bad atom list--------")
    badAtoms = help.getBadAtomsListFromFile(loadPath, "badatoms.csv")  # Get the bad atoms list we will use to reduce the list further
    print("---- Making unrestricted--------")
    dataPdbUn = help.makeCsv('PDB', pdbListIn, geos, [],True)
    print("---- Making unrestricted--------")
    dataPdbRes = help.makeCsv('PDB', pdbListIn, geos, [],False)
    dataPdbRes = help.applyRestrictions(dataPdbRes)
    print("---- Making reduced--------")
    dataPdbCut = help.makeCsv('PDB', pdbListIn, geos, badAtoms,False)
    dataPdbCut = help.applyRestrictions(dataPdbCut)
    print("---- Making adjusted--------")
    dataPdbAdj = help.makeCsv('ADJUSTED', pdbListIn, geos, badAtoms,False)
    dataPdbAdj = help.applyRestrictions(dataPdbAdj)
    #Create the "nearest" key
    dataPdbUn['NEAR'] =dataPdbUn['O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}_motif'].str[3:]
    dataPdbRes['NEAR'] =dataPdbRes['O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}_motif'].str[3:]
    dataPdbCut['NEAR'] =dataPdbCut['O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}_motif'].str[3:]
    dataPdbAdj['NEAR'] =dataPdbAdj['O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}_motif'].str[3:]
    # embellish with dssp - the dssp file was created ages ago from the linux laptop
    pdbdssp = pd.read_csv('C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/CsvGeos_BEST_Set0DSSPALL.csv')
    try:
        pdbdssp['rid'] = pdbdssp['rid'].astype(str)
        pdbdssp['DSSPID'] = pdbdssp['pdbCode'] + pdbdssp['chain'] + pdbdssp['rid']
        pdbdssp = pdbdssp[['DSSPID', 'dssp']]

    except:
        print('empty csv')

    allList = []
    allList.append([dataPdbUn, loadPath + "co_unrestricted.csv"])
    allList.append([dataPdbRes, loadPath + "co_restricted.csv"])
    allList.append([dataPdbCut, loadPath + "co_reduced.csv"])
    allList.append([dataPdbAdj, loadPath + "co_adjusted.csv"])

    for dataBlob in allList:
        dataCsv = dataBlob[0]
        dataPath = dataBlob[1]
        try:
            dataCsv['rid'] = dataCsv['rid'].astype(str)
            dataCsv['DSSPID'] = dataCsv['pdbCode'] + dataCsv['chain'] + dataCsv['rid']
        except:
            print('empty csv')

        dataCsv = dataCsv.set_index('DSSPID').join(pdbdssp.set_index('DSSPID'))
        dataCsv.to_csv(dataPath, index=False)

#Load anyway
dataPdbUn = pd.read_csv(loadPath + "co_unrestricted.csv")
dataPdbRes = pd.read_csv(loadPath + "co_restricted.csv")
dataPdbCut = pd.read_csv(loadPath + "co_reduced.csv")
dataPdbAdj = pd.read_csv(loadPath + "co_adjusted.csv")

dataPdbUn = dataPdbUn.query("bfactor <= 15")
dataPdbRes = dataPdbRes.query("bfactor <= 15")
dataPdbCut = dataPdbCut.query("bfactor <= 15")
dataPdbAdj = dataPdbAdj.query("bfactor <= 15")


geoTrios = [
            ['CA:C:N+1', 'TAU', 'C:O',False],
            ['CA:C:N+1', 'O:{,N,}', 'C:O',False],
            ['CA:C:N+1', 'O:{,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'C:O',False],
            ['CA:C:N+1', 'O:{,SG,}', 'C:O',False],
            ['CA:C:N+1', 'O:{,HOH,}', 'C:O',False],
            ['CA:C:N+1', 'O:{,HETATM,}', 'C:O',False],
            ['CA:C:N+1', 'O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'C:O',False],
            ['O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'C:O','CA:C:N+1',False],
            ['CA:C:N+1', 'O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'C:O',False],
            ['CA:C:N+1', 'O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'NEAR',True],
            ['C:O', 'O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'NEAR',True],
            ['CA:C:N+1', 'O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'aa',True],
            ['C:O', 'O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}', 'dssp',True],
            #['O:{,HOH,HETATM,SG,N,ND1,ND2,NE,NE1,NE2,NZ,NH1,NH2,}','dssp', 'C:O',True],

           ]




help.trioReports(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTrios, title,printPath,fileName)

#help.trioReports(["Unrestricted",dataPdbAdj],                 ["Restricted",dataPdbAdj],                 ["Reduced",dataPdbAdj],                 ["Adjusted",dataPdbAdj],                 geoTrios, title,printPath,fileName)



