# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdb as geopdb
import time
import pandas as pd
'''
Global data like lists
'''

def getGeoLists():
    #Tag then list of geos
    geoLists = []
    geoLists.append(['0ALL',['dssp']])
    geoLists.append(['1GLY',['N:N+1','TAU','PSI','PHI','N:C','CA:C','C:O','N:CA','C-1:N','C:N+1','OMEGA']])
    geoLists.append(['2GLY',['CA:C:O:N+1','O:N+1','CA:O','CA:N+1','CA:C:N+1','C-1:N:CA','N:O-2','N:CA:C:O-2']])
    geoLists.append(['3GLY',['N:O-2:CA','N-1:CA:C','CA:HOH','CA:HETATM','N:HETATM:C','N:HOH:C','N:CA:C:HETATM','N:CA:C:HOH']])
    #geoLists.append(['4GLY',['O-2:C','O-2:N:CA','O-2:N:CA:N+1']])
    #geoLists.append(['5GLY',['N:{O,OD1,OG1}','{O,OD1,OG1}:C','{O,OD1,OG1}:N:CA','{O,OD1,OG1}:N:CA:N+1']])
    #geoLists.append(['6GLY',['N:{O}','C:{O}','CA:N:{O}','N+1:CA:N:{O}']])
    #geoLists.append(['7GLY',['N:O-2','C:O-2','CA:N:O-2','N+1:CA:N:O-2']])
    #geoLists.append(['8GLY',['N:{OD1,OG1}','C:{OD1,OG1}','CA:N:{OD1,OG1}','N+1:CA:N:{OD1,OG1}']])
    return geoLists