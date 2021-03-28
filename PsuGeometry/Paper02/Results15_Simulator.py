# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import Categoriser as cluster
import random
import math
import pandas as pd
from PsuGeometry import GeoCalcs as calcs
'''
TAU correlations
'''
pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'

######## DEFINE FUNCTIONS ########

def randomOnSphere(point,radius):
    # 1. generate a random vector
    x = random.randint(-100, 100)
    y = random.randint(-100, 100)
    z = random.randint(-100, 100)

    # 2. scale by radius
    mag = math.sqrt(x**2 + y**2 + z**2)
    x = radius * x/mag
    y = radius * y/mag
    z = radius * z/mag

    # 3. translate by original point
    x = x + point[0]
    y = y + point[1]
    z = z + point[2]

    return x, y, z

'''
Ground state vibrations - look at bafctors
https://www.researchgate.net/post/In_a_typical_diatomic_molecular_vibration_how_much_does_the_bond_length_change_during_natural_vibration
https://arxiv.org/ftp/arxiv/papers/1401/1401.5760.pdf
https://pubs.acs.org/doi/10.1021/jp500586h

ANISOU
http://plato.cgl.ucsf.edu/pipermail/chimera-users/2012-February/007246.html

The six terms actually describe a 3x3 matrix that is symmetric about  
the diagonal.  The eigenvectors and eigenvalues of that matrix are the  
atomic displacement axes and the mean squares of the displacements  
respectively.  
The first thing to know is that the values shown in the  
ANISOU records are scaled by 10**4.  
The other thing to know is that  
the six numbers correspond to these positions in the matrix:  1,1 2,2  
3,3 1,2 1,3 2,3.  Since the eigenvalues are mean squares, one  
typically works with the square roots of the eigenvalues.

ANISOU  720  OE2 GLU A  32     4038   2211   7501    423  -2276    585       O 


https://www3.cmbi.umcn.nl/bdb/theory/

'''
###### PROGRAM BEGINS ########
num = 100000
vals = []
for count in range(0,num):
    if count%100 == 0:
        print(count,num)
    dic = {}
    N = randomOnSphere([-15.359,8.006,1.894],0.017)
    CA = randomOnSphere([-15.97, 7.49, 3.074], 0.015)
    C = randomOnSphere([-16.409, 6.044, 3.006], 0.017)
    tau = calcs.angle(N[0],N[1],N[2],CA[0],CA[1],CA[2],C[0],C[1],C[2])
    n_ca = calcs.distance(N[0],N[1],N[2],CA[0],CA[1],CA[2])
    ca_c = calcs.distance(CA[0],CA[1],CA[2],C[0],C[1],C[2])
    dic['pdbCode'] = 'pdb' + str(count)
    dic['chain'] = 'A'
    dic['rid'] = 1
    dic['TAU'] = tau
    dic['N:CA'] = n_ca
    dic['CA:C'] = ca_c
    vals.append(dic)

dataFrame = pd.DataFrame.from_dict(vals)
georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
georep.addHistogram(data=dataFrame, geoX='TAU', title='')
georep.addHistogram(data=dataFrame, geoX='N:CA', title='')
georep.addHistogram(data=dataFrame, geoX='CA:C', title='')
georep.printToHtml('Results 15. Simulated TAU', 3,'Results15_Sim')