# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
import pandas as pd
'''
This is a helper class to create a python class to provide some predfined lists of pdbs
'''

pdbdata09 = pd.read_csv('structures09.csv')
pdbList09 = pdbdata09['pdb_code'].tolist()

pdbdata1000 = pd.read_csv('structures1000.csv')
pdbList1000 = pdbdata1000['pdb_code'].tolist()

pdbdataJask = pd.read_csv('jask.csv')
pdbListJask = pdbdataJask['pdb_code'].tolist()

pdbdata17 = pd.read_csv('structures17.csv')
pdbList17 = pdbdata17['pdb_code'].tolist()

pdbdata25 = pd.read_csv('structures25.csv')
pdbList25 = pdbdata25['pdb_code'].tolist()

pdbdataPap = pd.read_csv('structuresPaper.csv')
pdbListPap = pdbdataPap['pdb_code'].tolist()

fileText = ''
'''
list1000 = ['xxx','sss','sss']
list100 = ['xxx','sss','aaa']
class GeoPdbLists:
      def getList1000(self):
        return list1000
      def getList100(self):
        return list100
            
'''
fileText +="# -- ©Rachel Alcraft 2020, PsuGeometry --\n\n\n"
fileText +="list100 = ["
isFirst = True
for pdb in pdbList09:
    if not isFirst:
        fileText += ","
    isFirst = False
    fileText +="'" + pdb.lower() + "'"
fileText +="]\n"

fileText +="list1000 = ["
isFirst = True
for pdb in pdbList1000:
    if not isFirst:
        fileText += ","
    isFirst = False
    fileText +="'" + pdb.lower() + "'"
fileText +="]\n"

fileText +="listJaskolski = ["
isFirst = True
for pdb in pdbListJask:
    if not isFirst:
        fileText += ","
    isFirst = False
    fileText +="'" + pdb.lower() + "'"
fileText +="]\n"

fileText +="list17 = ["
isFirst = True
for pdb in pdbList17:
    if not isFirst:
        fileText += ","
    isFirst = False
    fileText +="'" + pdb.lower() + "'"
fileText +="]\n"

fileText +="list25 = ["
isFirst = True
for pdb in pdbList25:
    if not isFirst:
        fileText += ","
    isFirst = False
    fileText +="'" + pdb.lower() + "'"
fileText +="]\n"

fileText +="listPap = ["
isFirst = True
for pdb in pdbListPap:
    if not isFirst:
        fileText += ","
    isFirst = False
    fileText +="'" + pdb.lower() + "'"
fileText +="]\n"

fileText +="class GeoPdbLists:\n"
fileText += "\tdef getList1000(self):\n"
fileText += "\t\treturn list1000\n"
fileText += "\tdef getList100(self):\n"
fileText += "\t\treturn list100\n"
fileText += "\tdef getListJask(self):\n"
fileText += "\t\treturn listJask \n"
fileText += "\tdef getList17(self):\n"
fileText += "\t\treturn list17 \n"
fileText += "\tdef getList25(self):\n"
fileText += "\t\treturn list25 \n"
fileText += "\tdef getListPaper(self):\n"
fileText += "\t\treturn listPap \n"
# create the cpp.html so that all the look and feel is consistent
f= open("../GeoPdbLists.py","w+")
print(fileText)
f.write(fileText)
f.close()