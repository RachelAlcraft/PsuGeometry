
import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoPdb as geopdb
from PsuGeometry import CloseContact as geocc
from PsuGeometry import GeoReport as psu

pdbListIn = help.getPDBList()
for pdb in pdbListIn:
    georep = psu.GeoReport([pdb], help.filesPDBRoot, help.edDataPath, help.printPath,False, False, False, [])
    df = georep.addCloseContact(pdb,'CA','CA',6,2,palette='terrain',hue='distance')
    df.to_csv(help.loadPath + "CloseContacts/CloseContacts_" + pdb + ".csv", index=False)

#georep.printToHtml('An Example of Each Plot Type', 3, 'each')