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

geos = ['C:O','PHI','PSI']


title='Histograms for DSSP No Chain Ends'
fileName = 'dssp_no_ends'

createOrLoad = "CREATE" # CREATE or LOAD
if createOrLoad == "LOAD":
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
    # embellish with dssp - the dssp file was created ages ago from the linux laptop
    pdbdssp = pd.read_csv('C:/Dev/Github/BbkProject/PhDThesis/5.Chapters/1_Summer/CSV/CsvGeos_BEST_Set0DSSPALL.csv')
    try:
        pdbdssp['rid'] = pdbdssp['rid'].astype(str)
        pdbdssp['DSSPID'] = pdbdssp['pdbCode'] + pdbdssp['chain'] + pdbdssp['rid']
        pdbdssp = pdbdssp[['DSSPID', 'dssp']]

    except:
        print('empty csv')

    allList = []
    allList.append([dataPdbUn, loadPath + "histend_unrestricted.csv"])
    allList.append([dataPdbRes, loadPath + "histend_restricted.csv"])
    allList.append([dataPdbCut, loadPath + "histend_reduced.csv"])
    allList.append([dataPdbAdj, loadPath + "histend_adjusted.csv"])

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

print('Loading csv files')#because of the dssp
dataPdbUn = pd.read_csv(loadPath + "histend_unrestricted.csv")
dataPdbRes = pd.read_csv(loadPath + "histend_restricted.csv")
dataPdbCut = pd.read_csv(loadPath + "histend_reduced.csv")
dataPdbAdj = pd.read_csv(loadPath + "histend_adjusted.csv")

#Choose a bfactor to cut at
dataPdbUn = dataPdbUn.query("`C:O_avbfactor` <= 15")
dataPdbRes = dataPdbRes.query("`C:O_avbfactor` <= 15")
dataPdbCut = dataPdbCut.query("`C:O_avbfactor` <= 15")
dataPdbAdj = dataPdbAdj.query("`C:O_avbfactor` <= 15")

dataPdbUn = dataPdbUn.query("pdbCode != '4r2x'")
dataPdbRes = dataPdbRes.query("pdbCode != '4r2x'")
dataPdbCut = dataPdbCut.query("pdbCode != '4r2x'")
dataPdbAdj = dataPdbAdj.query("pdbCode != '4r2x'")


geoTrios = [['C:O']]

help.trioReports(["Unrestricted",dataPdbUn],
                 ["Restricted",dataPdbRes],
                 ["Reduced",dataPdbCut],
                 ["Adjusted",dataPdbAdj],
                 geoTrios, title,printPath,fileName,'dssp')

groupUn = dataPdbUn[["pdbCode","C:O"]]
groupUn = groupUn.groupby(by=["pdbCode"]).mean()
groupUn.reset_index(inplace=True)
print(groupUn)
groupUn["ID"] = groupUn["pdbCode"]

groupRes = dataPdbRes[["pdbCode","C:O"]]
groupRes = groupRes.groupby(by=["pdbCode"]).mean()
groupRes.reset_index(inplace=True)
groupRes["ID"] = groupUn["pdbCode"]

groupCut = dataPdbCut[["pdbCode","C:O"]]
groupCut = groupCut.groupby(by=["pdbCode"]).mean()
groupCut.reset_index(inplace=True)
groupCut["ID"] = groupUn["pdbCode"]

groupAdj = dataPdbAdj[["pdbCode","C:O"]]
groupAdj = groupAdj.groupby(by=["pdbCode"]).mean()
groupAdj.reset_index(inplace=True)
groupAdj["ID"] = groupUn["pdbCode"]

geoTrios = [['C:O']]

help.trioReports(["Unrestricted",groupUn],
                 ["Restricted",groupRes],
                 ["Reduced",groupCut],
                 ["Adjusted",groupAdj],
                 geoTrios, title,printPath,fileName + "_av")











