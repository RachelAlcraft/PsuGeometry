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
printPath = 'F:/Code/BbkProject/PhDThesis/0.Papers/3.DefensibleGeometry/EvidencedSet/DataL/'

######## DEFINE FUNCTIONS ########

def randomOnSphere(point,radius):
    # 1. generate a random vector
    x = random.randint(-1000, 1000)
    y = random.randint(-1000, 1000)
    z = random.randint(-1000, 1000)
    v = random.randint(0,1000)
    v = v/1000

    mag = math.sqrt(x**2 + y**2 + z**2)

    # 2. scale by the random scale factor with max = 1
    radius = radius * v

    # 3. scale by radius
    if mag > 0:
        x = (radius * x)/mag
        y = (radius * y)/mag
        z = (radius * z)/mag
    else:
        x,y,z = 0,0,0

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

# Define the atoms and their positions and radii
A1_atom = [-15.359,8.006,1.894]
A2_atom =[-15.97, 7.49, 3.074]
A3_atom = [-16.409, 6.044, 3.006]

A1_atom = [4.028,12.634,2.294]
A2_atom = [4.137,11.461,3.154]
A3_atom = [4.607,11.823,4.564]

num = 1000000
vals = []
for count in range(0,num):
    if count%100 == 0:
        print(count,num)
    dic = {}

    A1 = randomOnSphere(A1_atom, 0.01)#7)
    A2 = randomOnSphere(A2_atom, 0.01)#5)
    A3 = randomOnSphere(A3_atom, 0.01)#7)

    a1a2a3 = calcs.angle(A1[0],A1[1],A1[2],A2[0],A2[1],A2[2],A3[0],A3[1],A3[2])
    a1a2 = calcs.distance(A1[0],A1[1],A1[2],A2[0],A2[1],A2[2])
    a2a3 = calcs.distance(A2[0],A2[1],A2[2],A3[0],A3[1],A3[2])
    dic['pdbCode'] = 'Iter_' + str(count)
    dic['chain'] = 'A'
    dic['rid'] = 1
    dic['ANGLE'] = a1a2a3
    dic['A1:A2'] = a1a2
    dic['A2:A3'] = a2a3
    vals.append(dic)

dataFrame = pd.DataFrame.from_dict(vals)
georep = psu.GeoReport([], pdbDataPath, edDataPath, printPath, ed=False, dssp=False, includePdbs=False,keepDisordered=False)
georep.addHistogram(data=dataFrame, geoX='ANGLE', title='')
georep.addHistogram(data=dataFrame, geoX='A1:A2', title='')
georep.addHistogram(data=dataFrame, geoX='A2:A3', title='')
georep.printToHtml('Simulated Atoms', 3,'SimReport')