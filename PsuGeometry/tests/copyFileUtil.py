
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
import os
import shutil

fileFrom = 'F:/Code/BbkTransfer/pdbfiles/pdbdata/'
fileTo = 'F:/Code/ProteinDataFiles/pdb_data/'
pdbList1000 = geol.GeoPdbLists().getListPaper()

for pdb in pdbList1000:
    fileNameFrom = fileFrom + pdb.upper() + '.pdb'
    fileNameTo = fileTo + 'pdb' + pdb.lower() + '.ent'

    if not os.path.exists(fileNameTo):
        if os.path.exists(fileNameFrom):
            shutil.copy(fileNameFrom, fileNameTo)
            print('Copying',fileNameTo)
