'''
This script creates a file for hand chosen geos, so it is slower as they are no serialised.
This particular report is designed to look at hydrogen bonding on the carbonyl oxygen
'''

import pandas as pd
from PsuGeometry import GeoReport as psu
import Ch000_Functions as help
import matplotlib
print(matplotlib.__version__)

geos = ['TAU','TAU+1','TAU-1','CA:C:O','O:C:N+1','CA:C:N+1','CA-1:CA:CA+1',
        'O-1:C-1','C-1:N','N:CA','CA:C','C:O','C:N+1','N+1:CA+1','CA+1:C+1','C+1:O+1',
        'PHI','PSI','OMEGA','CA-1:C-1:N:CA',
        'CA-1:CA','CA:CA+1','C-1:C','C:C+1','N-1:N','N:N+1',
        'CA-1:N','CA-1:O-1','O-1:N','C-1:CA','N:C','CA:O','CA:N+1','O:N+1','C:CA+1','N+1:C+1',
        'O-1:CA','N:O','O:CA+1','N+1:O+1','N-1:O-1']


title='Average Backbone Lengths'
fileName = 'average'


print('### LOADING csv files ###') # bit rubbish but we didn;t change the object references with dssp
descdataPdbUn = pd.read_csv(help.loadPath + "DescribeGeos_Unrestricted.csv")
descdataPdbRes = pd.read_csv(help.loadPath + "DescribeGeos_Restricted.csv")
descdataPdbCut = pd.read_csv(help.loadPath + "DescribeGeos_Cut.csv")
descdataPdbAdj = pd.read_csv(help.loadPath + "DescribeGeos_Adjusted.csv")


print('### Creating scatter files ###')



geoTriosA = [
            ['C:O mean'],['C:O 50%'],
            ['N:CA mean'],['N:CA 50%'],
            ['CA:C mean'],['CA:C 50%'],
            ['C:N+1 mean'],['C:N+1 50%'],
            ['TAU mean'],['TAU 50%'],
            ['C:O mean', 'C:O count', 'C:O 50%', False],
            ['C:O mean', 'N:CA mean', 'CA:C mean',False],
           ]

help.trioReports(["Unrestricted",descdataPdbUn],
                 ["Restricted",descdataPdbRes],
                 ["Reduced05",descdataPdbCut],
                 ["Adjusted05",descdataPdbAdj],
                 geoTriosA, title,help.printPath,fileName + "")



