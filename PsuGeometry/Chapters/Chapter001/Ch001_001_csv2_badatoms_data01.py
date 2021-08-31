
'''
In this file we load all the maxima differernce files to report on how much
each atom type varies
and to decide on the evidential atoms for each pdb
'''


import pandas as pd
import Ch000_Functions as help
from PsuGeometry import GeoReport as psu

pdbPath = help.loadPath + 'AllAtoms_Maxima.csv'
atomsdata = pd.read_csv(pdbPath)

print('Making numeric Vadj')
pd.to_numeric(atomsdata["Vadj"])
print('Making numeric Radj')
pd.to_numeric(atomsdata["Radj"])
print('Making numeric Vpdb')
pd.to_numeric(atomsdata["Vpdb"])
print('Making numeric Rpdb')
pd.to_numeric(atomsdata["Rpdb"])
print('Making numeric Difference')
pd.to_numeric(atomsdata["Difference"])

pdbListIn = atomsdata["pdbCode"].unique()

georep = psu.GeoReport([], "", "", help.printPath, ed=False, dssp=False, includePdbs=False, keepDisordered=False)

allpdbmaxEv  = atomsdata.query("Radj <= 0.065")
allpdbmaxEv = allpdbmaxEv.query("BFactor < 18")

atomsdataN = atomsdata.query("AtomType == 'N'")
print(atomsdataN)
atomsdataCA = atomsdata.query("AtomType == 'CA'")
atomsdataC = atomsdata.query("AtomType == 'C'")
atomsdataO = atomsdata.query("AtomType == 'O'")

allpdbmaxEvN = allpdbmaxEv.query("AtomType == 'N'")
allpdbmaxEvCA = allpdbmaxEv.query("AtomType == 'CA'")
allpdbmaxEvC = allpdbmaxEv.query("AtomType == 'C'")
allpdbmaxEvO = allpdbmaxEv.query("AtomType == 'O'")

georep.addHistogram(data=atomsdata,geoX='Difference', title='N all',hue='ID')
georep.addHistogram(data=atomsdataCA,geoX='Difference',title='CA all',hue='ID')
georep.addHistogram(data=atomsdataC,geoX='Difference',title='C all',hue='ID')
georep.addHistogram(data=atomsdataO,geoX='Difference',title='O all')

georep.addHistogram(data=allpdbmaxEvN,geoX='Difference',title='N: Radj<=0.065, bfactor<18',hue='ID')
georep.addHistogram(data=allpdbmaxEvCA,geoX='Difference',title='CA: Radj<=0.065, bfactor<18',hue='ID')
georep.addHistogram(data=allpdbmaxEvC,geoX='Difference',title='C: Radj<=0.065, bfactor<18',hue='ID')
georep.addHistogram(data=allpdbmaxEvO,geoX='Difference',title='O: Radj<=0.065, bfactor<18',hue='ID')


pdbextrabaddata = []
pdbtailbaddata = []
for pdb in pdbListIn:
    onepdbmax = atomsdata.query("pdbCode == '" + pdb + "'")

    #we want to choose a cutoff from the density values based on the summary statstics. Current plan, but of everything below the 25% quartile
    statspdb = onepdbmax[['Vadj']].describe()
    #print(statspdb)
    pdb25 = statspdb['Vadj'][4]
    #print('25%',pdb25)

    #Include a cut of 25% for vadj
    onepdbmaxEv  = onepdbmax.query("Radj <= 0.065")
    onepdbmaxEv = onepdbmaxEv.query("BFactor < 18")
    onepdbmaxEv = onepdbmaxEv.query("Vadj >= " + str(pdb25))
    onepdbtail = onepdbmax.query("Vadj<" + str(pdb25))
    onepdbbad = onepdbmax.query("Radj > 0.065 or BFactor >= 18")

    pdbextrabaddata.append(onepdbbad)
    pdbtailbaddata.append(onepdbtail)
    #print(onepdbbad[['Radj','BFactor']])

    rowso = onepdbmax.shape[0]
    rowse = onepdbmaxEv.shape[0]
    print(pdb,'Atoms changed from',rowso, 'to',rowse)

    onepdbmaxN = onepdbmax.query("AtomType == 'N'")
    onepdbmaxCA = onepdbmax.query("AtomType == 'CA'")
    onepdbmaxC = onepdbmax.query("AtomType == 'C'")
    onepdbmaxO = onepdbmax.query("AtomType == 'O'")

    evonepdbmaxN = onepdbmaxEv.query("AtomType == 'N'")
    evonepdbmaxCA = onepdbmaxEv.query("AtomType == 'CA'")
    evonepdbmaxC = onepdbmaxEv.query("AtomType == 'C'")
    evonepdbmaxO = onepdbmaxEv.query("AtomType == 'O'")

    onepdbmaxNdesc = onepdbmaxN[['BFactor','Difference','Vadj','Radj']].describe()
    onepdbmaxCAdesc = onepdbmaxCA[['BFactor','Difference','Vadj','Radj']].describe()
    onepdbmaxCdesc = onepdbmaxC[['BFactor','Difference','Vadj','Radj']].describe()
    onepdbmaxOdesc = onepdbmaxO[['BFactor','Difference','Vadj','Radj']].describe()

    onepdbmaxNdesc.columns = ['Bf','Diff','Vadj','Radj']
    onepdbmaxCAdesc.columns = ['Bf', 'Diff', 'Vadj', 'Radj']
    onepdbmaxCdesc.columns = ['Bf', 'Diff', 'Vadj', 'Radj']
    onepdbmaxOdesc.columns = ['Bf', 'Diff', 'Vadj', 'Radj']

    evonepdbmaxNdesc = evonepdbmaxN[['BFactor','Difference', 'Vadj','Radj']].describe()
    evonepdbmaxCAdesc = evonepdbmaxCA[['BFactor','Difference', 'Vadj','Radj']].describe()
    evonepdbmaxCdesc = evonepdbmaxC[['BFactor','Difference', 'Vadj','Radj']].describe()
    evonepdbmaxOdesc = evonepdbmaxO[['BFactor','Difference', 'Vadj','Radj']].describe()

    evonepdbmaxNdesc.columns = ['Bf', 'Diff', 'Vadj', 'Radj']
    evonepdbmaxCAdesc.columns = ['Bf', 'Diff', 'Vadj', 'Radj']
    evonepdbmaxCdesc.columns = ['Bf', 'Diff', 'Vadj', 'Radj']
    evonepdbmaxOdesc.columns = ['Bf', 'Diff', 'Vadj', 'Radj']

    onepdbmaxNdesc['BfEv'] = evonepdbmaxNdesc['Bf']
    onepdbmaxCAdesc['BfEv'] = evonepdbmaxCAdesc['Bf']
    onepdbmaxCdesc['BfEv'] = evonepdbmaxCdesc['Bf']
    onepdbmaxOdesc['BfEv'] = evonepdbmaxOdesc['Bf']

    onepdbmaxNdesc['DiffEv'] = evonepdbmaxNdesc['Diff']
    onepdbmaxCAdesc['DiffEv'] = evonepdbmaxCAdesc['Diff']
    onepdbmaxCdesc['DiffEv'] = evonepdbmaxCdesc['Diff']
    onepdbmaxOdesc['DiffEv'] = evonepdbmaxOdesc['Diff']

    onepdbmaxNdesc['VEv'] = evonepdbmaxNdesc['Vadj']
    onepdbmaxCAdesc['VEv'] = evonepdbmaxCAdesc['Vadj']
    onepdbmaxCdesc['VEv'] = evonepdbmaxCdesc['Vadj']
    onepdbmaxOdesc['VEv'] = evonepdbmaxOdesc['Vadj']

    onepdbmaxNdesc['REv'] = evonepdbmaxNdesc['Radj']
    onepdbmaxCAdesc['REv'] = evonepdbmaxCAdesc['Radj']
    onepdbmaxCdesc['REv'] = evonepdbmaxCdesc['Radj']
    onepdbmaxOdesc['REv'] = evonepdbmaxOdesc['Radj']

    georep.addCsv(data=onepdbmaxNdesc,title=pdb + ' N')
    georep.addCsv(data=onepdbmaxCAdesc,title=pdb + ' CA')
    georep.addCsv(data=onepdbmaxCdesc, title=pdb + ' C')
    georep.addCsv(data=onepdbmaxOdesc, title=pdb + ' O')


pdbbadcsv = pd.concat(pdbextrabaddata)
pdbtailcsv = pd.concat(pdbtailbaddata)

pdbbadcsv['ResNo'] = pdbbadcsv['ResNo'].astype(str)
pdbbadcsv['BAD'] = pdbbadcsv['pdbCode'] +pdbbadcsv['Chain'] + pdbbadcsv['ResNo'] + pdbbadcsv['AtomType']
pdbtailcsv['ResNo'] = pdbtailcsv['ResNo'].astype(str)
pdbtailcsv['BAD'] = pdbtailcsv['pdbCode'] +pdbtailcsv['Chain'] + pdbtailcsv['ResNo'] + pdbtailcsv['AtomType']

badonly =pdbbadcsv['BAD']
tailonly =pdbtailcsv['BAD']
print('Save BAD to',help.loadPath + "ExtraBadAtoms_Maxima.csv")
badonly.to_csv(help.loadPath + "ExtraBadAtoms_Maxima.csv", index=False)
tailonly.to_csv(help.loadPath + "TailBadAtoms_Maxima.csv", index=False)


georep.printToHtml('How Far Are Atoms Adjusted?', 4, 'atoms_per_pdb')

