# -- ©Rachel Alcraft 2020, PsuGeometry --
from PsuGeometry import GeoReport as psu
from PsuGeometry import GeoPdbLists as geol
from sklearn.model_selection import cross_validate
from sklearn.ensemble import RandomForestRegressor
import random

import seaborn as sns
import matplotlib.pyplot as plt
'''
TAU correlations
'''
###############################################################################################
myWindowsLaptop = False
pdbList1000 = geol.GeoPdbLists().getListPaper()
pdbList1000 = pdbList1000[:400]

#randomise list for regression
random.shuffle(pdbList1000)

geoList = ['N:N+1','TAU','PSI']
hueList = ['aa']
aas = ['GLY','ALA']
includeDSSP = True
###################################################################################
pdbDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/pdb_data/'
edDataPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/ccp4_data/'
printPath = '/home/rachel/Documents/Bioinformatics/ProteinDataFiles/results_psu/Paper02/'
if myWindowsLaptop:
    pdbDataPath = 'F:/Code/ProteinDataFiles/pdb_data/'
    edDataPath = 'F:/Code/ProteinDataFiles/ccp4_data/'
    printPath = 'F:/Code/ProteinDataFiles/results_psu/Paper02/'
    includeDSSP = False  # on my windows computer

###########################################################################################
georep = psu.GeoReport(pdbList1000, pdbDataPath, edDataPath, printPath, ed=False, dssp=includeDSSP, includePdbs=False)
data = georep.getGeoemtryCsv(geoList, hueList)
#datacorr = data.corr()
#sns.heatmap(datacorr, annot=True, cmap="vlag", vmin=-1, vmax=1)
#plt.show()

#Clean the data
data = data.drop('pdbCode', axis=1)
data = data.drop('chain', axis=1)
data = data.drop('rid', axis=1)
data = data.drop('aa', axis=1)
data = data.drop('ridx', axis=1)
data = data.drop('atomNo', axis=1)
data = data.drop('bfactor', axis=1)
data = data.drop('N:CA:C:N+1', axis=1)
data = data.drop('N:CA:C', axis=1)
#randomise the data as it is by pdb res no
data = data.sample(frac=1)

Y =  data['TAU']
X = data.drop(['TAU'],axis=1)
print(X)
print(Y)

test_set_size = 100
X_train = X.iloc[test_set_size:]
Y_train = Y.iloc[test_set_size:]
X_test = X.iloc[:test_set_size]
Y_test = Y.iloc[:test_set_size]

clf = RandomForestRegressor(max_depth=2, random_state=0)

clf = clf.fit(X_train, Y_train)

Y_pred = clf.predict(X_test)
#print('Actual',Y_test)
#print('Predicted',Y_pred)

score = clf.score(X_test, Y_test)

for i in range(0,len(Y_pred)):
    pr = Y_pred[i]
    ac = Y_test.to_numpy()[i]
    percent = abs((pr-ac)*100/ac)
    print(ac,':',pr,' %=',percent)

psi = 176
NN1 = 3.5
print('Prediction for',psi,NN1,clf.predict([[NN1,psi]]))
print('Score on training set=',score)










