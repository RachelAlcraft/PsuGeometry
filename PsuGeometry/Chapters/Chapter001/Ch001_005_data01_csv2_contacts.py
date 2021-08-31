
import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoPdb as geopdb
from PsuGeometry import CloseContact as geocc
from PsuGeometry import GeoReport as psu

pdbListIn = help.getPDBList()
contactslist = []
for pdb in pdbListIn:
    print(pdb)
    try:
        import os.path
        iscc = os.path.isfile((help.loadPath + "CloseContacts/CloseContacts_" + pdb + ".csv").lower())
        if iscc:
            ccdata = pd.read_csv(help.loadPath + "CloseContacts/CloseContacts_" + pdb + ".csv")
            ccdata['ridA'] = ccdata['ridA'].astype(str)
            ccdata['ChRid'] =ccdata['chainA'] +ccdata['ridA']
            cc = ccdata[['pdbCode','ChRid']].groupby('ChRid').agg('count')
            cc['ChRid'] = cc.index
            cc.columns = ['Contacts','ChRid']
            cc['pdbCode'] = pdb
            cc['CID'] = cc['pdbCode'] +cc['ChRid']
            contactslist.append(cc)
    except:
        print('Error with',pdb)



ccall = pd.concat(contactslist)
ccall.to_csv(help.loadPath + "Contacts_List.csv", index=False)
print('Merging')
#Now load some data
mergedDataSet = pd.read_csv(help.loadPath + "MergedEvidenced.csv")
ccContacts = ccall[['CID','Contacts']]
mergedDataSet['rid'] = mergedDataSet['rid'].astype(str)
mergedDataSet['CID'] = mergedDataSet['pdbCode'] + mergedDataSet['chain'] +mergedDataSet['rid']
mergedDataSet = mergedDataSet.set_index('CID').join(ccContacts.set_index('CID'))
mergedDataSet.to_csv(help.loadPath + "MergedWithContacts.csv", index=False)

mergedDataSet = mergedDataSet.dropna()
#qu = "dssp == 'E' or dssp == 'B' or dssp == '-' "
#qu = qu + "or dssp == 'T' or dssp == 'S' or dssp == 'H' or dssp == 'G' or dssp == 'I')"
#mergedDataSet = mergedDataSet.query(qu)

# create a report based on contacts
georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

print('### Creating reports ###')
georep.addScatter(data=mergedDataSet, geoX='N:CA_Orig',geoY='Contacts', title='N:CA Original', hue='dssp',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=mergedDataSet, geoX='N:CA_Diff',geoY='Contacts', title='N:CA Difference', hue='bfactor',categorical=False,sort='NON',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='N:CA_Adj',geoY='Contacts', title='N:CA Adjusted', hue='dssp',categorical=True,sort='NON',palette='jet_r')

georep.addScatter(data=mergedDataSet, geoX='CA:C_Orig',geoY='Contacts', title='CA:C Original', hue='dssp',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=mergedDataSet, geoX='CA:C_Diff',geoY='Contacts', title='CA:C Difference', hue='bfactor',categorical=False,sort='NON',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='CA:C_Adj',geoY='Contacts', title='CA:C Adjusted', hue='dssp',categorical=True,sort='NON',palette='jet_r')

georep.addScatter(data=mergedDataSet, geoX='C:O_Orig',geoY='Contacts', title='C:O Original', hue='dssp',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=mergedDataSet, geoX='C:O_Diff',geoY='Contacts', title='C:O Difference', hue='bfactor',categorical=False,sort='NON',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='C:O_Adj',geoY='Contacts', title='C:O Adjusted', hue='dssp',categorical=True,sort='NON',palette='jet_r')

georep.addScatter(data=mergedDataSet, geoX='C:N+1_Orig',geoY='Contacts', title='C:N+1 Original', hue='dssp',categorical=True,sort='NON',palette='jet_r')
georep.addScatter(data=mergedDataSet, geoX='C:N+1_Diff',geoY='Contacts', title='C:N+1 Difference', hue='bfactor',categorical=False,sort='NON',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='C:N+1_Adj',geoY='Contacts', title='C:N+1 Adjusted', hue='dssp',categorical=True,sort='NON',palette='jet_r')

georep.addScatter(data=mergedDataSet, geoX='C:O_Orig',geoY='Contacts', title='C:O Original', hue='bfactor',categorical=False,sort='DESC',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='C:O_Diff',geoY='Contacts', title='C:O Difference', hue='bfactor',categorical=False,sort='DESC',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='C:O_Adj',geoY='Contacts', title='C:O Adjusted', hue='bfactor',categorical=False,sort='DESC',palette='jet')

georep.addScatter(data=mergedDataSet, geoX='C:O_Orig',geoY='Contacts', title='C:O Original', hue='bfactor',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='C:O_Diff',geoY='Contacts', title='C:O Difference', hue='bfactor',categorical=False,sort='RAND',palette='jet')
georep.addScatter(data=mergedDataSet, geoX='C:O_Adj',geoY='Contacts', title='C:O Adjusted', hue='bfactor',categorical=False,sort='RAND',palette='jet')




georep.printToHtml('Contacts and changed values',3, 'contacts')



