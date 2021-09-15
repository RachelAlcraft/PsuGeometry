
import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu
pdblist = help.getPDBList100()
pdblist.sort()
#pdblist = pdblist[:10]
hueList = ['aa', 'rid', 'bfactor', 'pdbCode', 'bfactorRatio', 'disordered','occupancy','dssp']

dsspPrintPath = '../../PdbLists/'

georep = psu.GeoReport(pdblist, help.pdbDataPathLx, help.edDataPath, dsspPrintPath, ed=False, dssp=True, includePdbs=False, keepDisordered=True)
datacsv = georep.getGeoemtryCsv(['N:CA'],hueList)
datacsv = datacsv[['pdbCode','chain','rid','aa','dssp']]
print(datacsv)
if False:#don;t accidentally run this and replace it
    datacsv.to_csv(dsspPrintPath + 'dssp.csv', index=False)
print(datacsv)