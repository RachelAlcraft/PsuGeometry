import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64

from PsuGeometry import GeoAtom as atm
from PsuGeometry import GeoDensity as den
from PsuGeometry import GeoCalcs as calcs
from PsuGeometry import GeoQuery as que


class GeoReport:

    def __init__(self,listPdbs):
        self.pdbs = listPdbs


    def getGeoemtryCsv(self,calcList, hueList):
        dfs = []
        for apdb in self.pdbs:
            data = apdb.getGeoemtryCsv(calcList, hueList)
            dfs.append(data)
        df = pd.concat(dfs, ignore_index=True)
        return (df)

    def getReportCsv(self, reportName):
        hueList = ['2FoFc','FoFc','bfactor','aa','dssp']
        if reportName == 'Ramachandran':
            calcList = ['C-1:N:CA:C', 'N:CA:C:N+1']
        elif reportName == 'Sp2Planarity':
            calcList = ['CA:C:O','O:C:N+1','N+1:C:CA','N+1:O:C:CA']
        elif reportName == 'Sp3Tetrahedra':
            calcList = ['N:CA:C','C:CA:CB','N:CA:CB']
        elif reportName == 'BackboneOutliers':
            calcList = ['C-1:N','N:CA','CA:C','C:N+1','C-1:N:CA','N:CA:C','CA:C:N+1']
        elif reportName == 'OmegaCis':
            calcList = ['CA-1:CA','CA:CA+1','CA:C:N+1:CA+1','CA-1:C-1:N:CA','N:CA:C']
        elif reportName == 'RachelsChoice':
            calcList = ['N:O','CB:O','N:CA:C:N+1','C-1:N:CA:C']
        dfs = []
        for apdb in self.pdbs:
            data = apdb.getGeoemtryCsv(calcList, hueList)
            dfs.append(data)
        df = pd.concat(dfs,ignore_index=True)
        return (df)


    def printReport(self, reportName,printPath,fileName):
        print('PSU: create report',reportName,'for',fileName)
        printList = []
        if reportName == 'BackboneOutliers': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Backbone Outliers Report'
            cols = 3
            printList = []
            #printList.append(GeoQuery(['Bonds', atomData, 'C-1:N', 'N:CA','aa', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(que.GeoQuery(atomData,'C-1:N',geoY='N:CA',title='Bonds', splitKey='aa'))
            printList.append(que.GeoQuery(atomData, 'CA:C', geoY='C:N+1', title='Bonds'))
            printList.append(que.GeoQuery(atomData, 'C-1:N', geoY='C:N+1', title='Bonds'))
            printList.append(que.GeoQuery(atomData, 'C-1:N:CA', geoY='N:CA:C', title='Angles'))
            printList.append(que.GeoQuery(atomData, 'N:CA:C', geoY='CA:C:N+1', title='Angles'))
            printList.append(que.GeoQuery(atomData, 'C-1:N:CA', geoY='CA:C:N+1', title='Angles'))
            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)
        elif reportName == 'RachelsChoice':
            atomData = self.getReportCsv(reportName)
            title = "Rachel's Choice of Correlations"
            cols = 4
            printList = []

            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='dssp',palette='Set1',newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='aa', palette='Set1',newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='pdbCode', palette='Set1',newData=True))

            printList.append(que.GeoQuery(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='aa', palette='Set1',newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='pdbCode', palette='Set1',newData=True))

            printList.append(que.GeoQuery(None, 'N:CA', geoY='CA:C', title='', hue='dssp', palette='Set1',newData=True))
            printList.append(que.GeoQuery(None, 'N:CA', geoY='CA:C', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA', geoY='CA:C', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA', geoY='CA:C', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'N:O', geoY='CB:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:O', geoY='CB:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'N:O', geoY='CB:O', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:O', geoY='CB:O', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='aa', palette='Set1',newData=True))
            printList.append(que.GeoQuery(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='pdbCode', palette='Set1',newData=True))

            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='pdbCode', palette='Set1', newData=True))

            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='aa', palette='Set1',newData=True))
            printList.append(que.GeoQuery(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='pdbCode', palette='Set1',newData=True))

            printList.append(que.GeoQuery(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(que.GeoQuery(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='aa', palette='Set1',newData=True))
            printList.append(que.GeoQuery(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='pdbCode', palette='Set1',newData=True))

            printList.append(que.GeoQuery(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='aa', palette='Set1', newData=True))
            printList.append(que.GeoQuery(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='pdbCode', palette='Set1', newData=True))

            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)


        elif reportName == 'Sp2Planarity': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Sp2 Planarity'
            cols = 4
            printList = []
            printList.append(que.GeoQuery(atomData,'N+1:O:C:CA',title='AbsVal Dihedral',hue='rid',operation='ABS'))
            printList.append(que.GeoQuery(atomData,'CA:C:O',title='CA:C:O',hue='rid'))
            printList.append(que.GeoQuery(atomData,'O:C:N+1',title='O:C:N+1',hue='rid'))
            printList.append(que.GeoQuery(atomData,'N+1:C:CA',title='N+1:C:CA',hue='rid'))

            printList.append(que.GeoQuery(atomData, 'N+1:O:C:CA', geoY='CA:C:O',hue='dssp'))
            printList.append(que.GeoQuery(atomData, 'N+1:O:C:CA', geoY='O:C:N+1'))
            printList.append(que.GeoQuery(atomData, 'N+1:O:C:CA', geoY='N+1:C:CA',hue='bfactor'))
            printList.append(que.GeoQuery(atomData, 'N+1:C:CA', geoY='O:C:N+1',hue='FoFc', palette='Spectral',centre=True))

            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)

        elif reportName == 'DataPerPdb':
            for apdb in self.pdbs:
                print('\tPSU:', reportName, 'for', apdb.pdbCode)
                atomData = apdb.dataFrame
                title = 'General Data Report'
                cols = 3
                printList = []
                printList.append(que.GeoQuery(atomData, 'atomNo', geoY='aa', hue='aa', palette='tab20'))
                printList.append(que.GeoQuery(atomData, 'atomNo', geoY='dssp',hue= 'aa',palette='tab20'))
                printList.append(que.GeoQuery(atomData, '2FoFc', geoY='bfactor',hue= 'element',palette='tab20'))
                printList.append(que.GeoQuery(atomData, 'atomNo', geoY='bfactor',hue= 'element',palette='tab20'))
                printList.append(que.GeoQuery(atomData, 'atomNo', geoY='2FoFc',hue='element',palette='tab20'))
                printList.append(que.GeoQuery(atomData, 'atomNo', geoY='FoFc',hue='element',palette='tab20'))
                printList.append(que.GeoQuery(atomData, 'x', geoY='y',hue='atomNo', palette='plasma_r'))
                printList.append(que.GeoQuery(atomData, 'y', geoY='z',hue='atomNo', palette='plasma_r'))
                printList.append(que.GeoQuery(atomData, 'z', geoY='x',hue='atomNo', palette='plasma_r'))
                self.printCsvToHtml(printList, [apdb], title, cols, printPath, fileName + '_' + apdb.pdbCode)

        elif reportName == 'Slow_DensityPointsPerPdb' or reportName == 'Slow_DensityPeaksPerPdb': # this can only be done per pdb
            for apdb in self.pdbs:
                if apdb.hasDensity:
                    print('\tPSU:', reportName, 'for', apdb.pdbCode)
                    allPoints = True
                    title = 'Density and Atom Points Comparison'
                    if reportName == 'Slow_DensityPeaksPerPdb':
                        allPoints = False
                        title = 'Density and Atom Peaks Comparison'
                    peaksData = apdb.geoDen.getPeaks(allPoints)
                    atomData = apdb.dataFrame
                    atomData['FoFc2'] = atomData['FoFc'] ** 2

                    cols = 3
                    printList = []
                    printList.append(que.GeoQuery(peaksData, 'c', geoY='r', title='Density CR Fo',hue='peakFo',palette='gist_gray_r'))
                    printList.append(que.GeoQuery(peaksData, 'r', geoY='s',title='Density RS Fo',hue='peakFo',palette='gist_gray_r'))
                    printList.append(que.GeoQuery(peaksData, 's', geoY='c',title='Density SC Fo',hue='peakFo',palette='gist_gray_r'))

                    printList.append(que.GeoQuery(peaksData, 'c', geoY='r',title='Density CR Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 'r', geoY='s',title='Density RS Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 's', geoY='c',title='Density SC Fo',hue='peakFo',palette='cubehelix_r'))

                    printList.append(que.GeoQuery(peaksData, 'c', geoY='r',title='Density CR 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 'r', geoY='s',title='Density RS 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 's', geoY='c',title='Density SC 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))

                    printList.append(que.GeoQuery(peaksData, 'c', geoY='r',title='Density CR FC',hue='peakFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 'r', geoY='s',title='Density RS FC',hue='peakFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 's', geoY='c',title='Density SC FC',hue='peakFc',palette='cubehelix_r'))

                    printList.append(que.GeoQuery(peaksData, 'c', geoY='r',title='Density CR FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(que.GeoQuery(peaksData, 'r', geoY='s',title='Density RS FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(que.GeoQuery(peaksData, 's', geoY='c',title='Density SC FoFC',hue='peakFoFc',palette='PiYG',centre=True))

                    printList.append(que.GeoQuery(peaksData, 'x', geoY='y',title='Density XY Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 'y', geoY='z',title='Density YZ Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 'z', geoY='x',title='Density ZX Fo',hue='peakFo',palette='cubehelix_r'))

                    printList.append(que.GeoQuery(peaksData, 'x', geoY='y',title='Density XY FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(que.GeoQuery(peaksData, 'y', geoY='z',title='Density YZ FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(que.GeoQuery(peaksData, 'z', geoY='x', title='Density ZX FoFC',hue='peakFoFc',palette='PiYG',centre=True))

                    printList.append(que.GeoQuery(peaksData, 'x', geoY='y',title='Density XY 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 'y', geoY='z',title='Density YZ 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(peaksData, 'z', geoY='x',title='Density ZX 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))

                    printList.append(que.GeoQuery(atomData, 'x', geoY='y',title='PDB XY 2FoFc',hue='2FoFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'y', geoY='z',title='PDB YZ 2FoFc',hue='2FoFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'z', geoY='x',title='PDB ZX 2FoFc',hue='2FoFc',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'x', geoY='y',title='PDB XY Electrons',hue='electrons',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'y', geoY='z',title='PDB YZ Electrons',hue='electrons',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'z', geoY='x',title='PDB ZX Electrons',hue='electrons',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'x', geoY='y',title='PDB XY FoFc',hue='FoFc',palette='PiYG',centre=True))
                    printList.append(que.GeoQuery(atomData, 'y', geoY='z',title='PDB YZ FoFc',hue='FoFc',palette='PiYG',centre=True))
                    printList.append(que.GeoQuery(atomData, 'z', geoY='x',title='PDB ZX FoFc',hue='FoFc',palette='PiYG',centre=True))
                    printList.append(que.GeoQuery(atomData, 'x', geoY='y',title='PDB XY bfactor',hue='bfactor',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'y', geoY='z',title='PDB YZ bfactor',hue='bfactor',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'z', geoY='x',title='PDB ZX bfactor',hue='bfactor',palette='cubehelix_r'))
                    printList.append(que.GeoQuery(atomData, 'x', geoY='y',title='PDB XY atom nos',hue='atomNo',palette='gist_ncar'))
                    printList.append(que.GeoQuery(atomData, 'y', geoY='z',title='PDB YZ atom nos',hue='atomNo',palette='gist_ncar'))
                    printList.append(que.GeoQuery(atomData, 'z', geoY='x',title='PDB ZX atom nos',hue='atomNo',palette='gist_ncar'))
                    printList.append(que.GeoQuery(atomData, 'x', geoY='y',title='PDB XY amino acids',hue='aa',palette='tab20'))
                    printList.append(que.GeoQuery(atomData, 'y', geoY='z',title='PDB YZ amino acids',hue='aa',palette='tab20'))
                    printList.append(que.GeoQuery(atomData, 'z', geoY='x',title='PDB ZX amino acids',hue='aa',palette='tab20'))
                    printList.append(que.GeoQuery(atomData, 'bfactor', geoY='FoFc2',title='PDB bfactor vs FoFc^2',hue='electrons'))
                    printList.append(que.GeoQuery(atomData, 'Fc', geoY='Fo',title='PDB Fc vs Fc',hue='electrons'))
                    printList.append(que.GeoQuery(atomData, 'electrons', geoY='2FoFc',title='PDB electrons vs 2FoFc',hue='element',palette='tab10'))
                    printList.append(que.GeoQuery(atomData, 'aa',title='Amino Acids'))
                    printList.append(que.GeoQuery(atomData, 'element',title='Atoms'))
                    printList.append(que.GeoQuery(atomData, '2FoFc',title='Peaks in 2FoFc'))
                    self.printCsvToHtml(printList, [apdb], title, cols, printPath, fileName + '_' + apdb.pdbCode)
                else:
                    print('\tPSU:',apdb.pdbCode,'has no density matrix')




    def printCsvToHtml(self, queryList,pdbList,title,cols,printPath,fileName):
        width=str(100/cols)
        html = '<!DOCTYPE html><html lang="en"><head><title>PSU-' + fileName + '-GEO</title>\n'
        #html += '<style> body {background-color:SeaShell;} table {table-layout:fixed;display:table;margin:0 auto;}td {border:1px solid RosyBrown;background-color:SeaShell;}</style>'
        #html += '<style> body {background-color:HoneyDew;} table {background-color:HoneyDew;} .innertable td {border:1px solid MistyRose;background-color:MintCream;}</style>'
        html += '<style> body {text-align:center;background-color:LightSteelBlue ;} img {width:95% }'
        html += 'table {font-size:0.8vw;width:95%;table-layout:fixed;display:table;margin:0 auto;background-color:LightSteelBlue ;}'
        html += ' td {border:1px solid MistyRose;background-color:AliceBlue;}</style>'
        html += '</head>\n'
        html += '<body><h1>' + title + '</h1>\n'
        html += '<h2>PSU: Geometric Correlations</h2>\n'
        html += '<hr/>'

        if len(pdbList) > 0:
            html += '<table><tr><td>PdbCode</td><td>Resolution</td><td>Pdb Link</td><td>PDBe Link</td></tr>\n'
            for apdb in pdbList:
                html += '<tr>\n'
                html += '<td>' + apdb.pdbCode + '</td>\n'
                html += '<td>' + str(apdb.atoms[0].values['resolution']) + '</td>\n'
                html += "<td><a href='https://www.rcsb.org/structure/" + apdb.pdbCode + "' title='PDB Link' target='_blank'>Link to PDB</a></td>\n"
                html += "<td><a href='https://www.ebi.ac.uk/pdbe/entry/pdb/" + apdb.pdbCode + "' title='PDB Link' target='_blank'>Link to PDBe</a></td>\n"
                html += '</tr>\n'
            html += '</table>\n'


        html += '<hr/>\n'

        reportPath = printPath + fileName + ".html"

        count = 0
        #html += '<table style="width:90%">\n'
        html += '<table>\n'
        for geoq in queryList:

            title = geoq.title
            alldata = geoq.data
            geoX = geoq.geoX
            geoY = geoq.geoY
            splitKey = geoq.splitKey
            splitList = ['']
            if splitKey != '':
                splitList = alldata[splitKey].unique()

            for split in splitList:
                geoqSplit = geoq
                if count == 0:
                    html += '<tr>'
                elif count % cols == 0:
                    html += '</tr>\n<tr>'

                count += 1

                if split != '':
                    data = alldata[alldata[splitKey] == split]
                    sptitle = title + ' ' + split
                    geoqSplit.data = data
                    geoqSplit.title = title + ' ' + split

                if geoY == '': # then it is 1d
                    html += self.oneHTMLHistogram(geoqSplit,width)
                else:
                    hue = geoq.hue
                    palette = geoq.palette
                    centre = geoq.centre
                    vmin = geoq.vmin
                    vmax = geoq.vmax
                    html += self.oneHTMLCorrelation(geoqSplit,width)

        html += '</tr></table><hr/><p>Produced by PsuGeometry, written by Rachel Alcraft </p></body>\n'

        # and print
        f = open(reportPath, "w+")
        f.write(html)
        print('PSU: saved file to',reportPath)
        f.close()

    def oneHTMLCorrelation(self,geoQ,width):
        html = '<td width=' + width + '%>' + self.createScatterPlot(geoQ) + '</td>\n'
        return (html)

    def oneHTMLHistogram(self,geoQ,width):
        html = '<td width=' + width + '%>'  + self.createHistogram(geoQ) + '</td>\n'
        return (html)

    def createHistogram(self,geoQ):
        fig, ax = plt.subplots()
        if geoQ.newData:
            calcList = [geoQ.geoX]
            hueList = [geoQ.hue]
            dfs = []
            for apdb in self.pdbs:
                data = apdb.getGeoemtryCsv(calcList, hueList)
                dfs.append(data)
            geoQ.data = pd.concat(dfs, ignore_index=True)

        data = geoQ.data.sort_values(by=geoQ.geoX,ascending=True)
        title = geoQ.title

        if geoQ.operation == 'ABS':
            data = data[data[geoQ.geoX] == abs(data[geoQ.geoX])]
        elif geoQ.operation == 'SQUARE':
            data = data[data[geoQ.geoX] == data[geoQ.geoX]**2]

        firstVal = data.head(1)[geoQ.geoX].values[0]
        lastVal = data.tail(1)[geoQ.geoX].values[0]
        firstHue = data.head(1)[geoQ.hue].values[0]
        lastHue = data.tail(1)[geoQ.hue].values[0]

        try:
            firstVal = round(firstVal, 2)
            lastVal = round(lastVal,2)
        except:
            pass
        try:
            firstHue = round(firstHue, 2)
            lastHue = round(lastHue,2)
        except:
            pass
        title += '\nFirst:' + geoQ.hue + ' ' + str(firstHue) + '='  + str(firstVal)
        title += '\nLast:' + geoQ.hue + ' ' + str(lastHue) + '='  + str(lastVal)


        # sns.distplot(data[xName], norm_hist=True, bins=50, kde=False)
        plt.hist(data[geoQ.geoX], EdgeColor='k', bins=50)

        plt.title(title)
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        encoded = base64.b64encode(img.getvalue())
        #html = '<img width=100% src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        html = '<p><img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        plt.close('all')
        dfdesc = geoQ.data[geoQ.geoX].describe()
        rows = len(dfdesc.index)
        colsNames = list(dfdesc.index)
        html += "<table class='innertable'>\n"
        html += "<tr>\n"
        for r in range(0, rows):
            html += "<td>" + str(colsNames[r]) + "</td>\n"
        html += "</tr>\n"
        html += "<tr>"
        for r in range(0, rows):
            html += "<td>"
            try:
                html += str(round(dfdesc[r],2))
            except:
                html += str(dfdesc[r])
            html += "</td>\n"

        html += "</tr>\n"
        html += "</table></p>\n"

        return html

    def createScatterPlot(self,geoQ):

        if geoQ.newData:
            calcList = [geoQ.geoX,geoQ.geoY]
            hueList = [geoQ.hue]
            dfs = []
            for apdb in self.pdbs:
                datatmp = apdb.getGeoemtryCsv(calcList, hueList)
                dfs.append(datatmp)
            geoQ.data = pd.concat(dfs, ignore_index=True)

        fig, ax = plt.subplots()
        if geoQ.centre:
            geoQ.data[geoQ.hue + '2'] = geoQ.data[geoQ.hue]**2
            data = geoQ.data.sort_values(by=geoQ.hue+'2', ascending=True)
            maxh = max(data[geoQ.hue].max(), -1 * data[geoQ.hue].min())
            minh = maxh * -1
            g = ax.scatter(data[geoQ.geoX], data[geoQ.geoY], c=data[geoQ.hue], cmap=geoQ.palette, vmin=minh, vmax=maxh)
            fig.colorbar(g)
            ax.set_xlabel(geoQ.geoX)
            ax.set_ylabel(geoQ.geoY)
        elif geoQ.vmin < geoQ.vmax:
            data = geoQ.data.sort_values(by=geoQ.hue, ascending=True)
            g = ax.scatter(data[geoQ.geoX], data[geoQ.geoY], c=data[geoQ.hue], cmap=geoQ.palette, vmin=geoQ.vmin, vmax=geoQ.vmax)
            fig.colorbar(g)
            ax.set_xlabel(geoQ.geoX)
            ax.set_ylabel(geoQ.geoY)
        else:

            lw = 0.5
            if geoQ.palette == 'gist_gray_r':
                lw = 0 # this gives a crystollagraphic image look

            if geoQ.hue == 'aa':
                try:
                    geoQ.data = geoQ.data.sort_values(by='2FoFc', ascending=True)
                except:
                    geoQ.data = geoQ.data.sort_values(by=geoQ.hue, ascending=True)
            else:
                geoQ.data = geoQ.data.sort_values(by=geoQ.hue, ascending=True)

            im = sns.scatterplot(x=geoQ.geoX, y=geoQ.geoY, hue=geoQ.hue, data=geoQ.data, alpha=0.65, palette=geoQ.palette,edgecolor='aliceblue',linewidth=lw)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # Put the legend out of the figure

        count = len(geoQ.data.index)
        title = geoQ.title
        if title == '':
            title += 'Count=' + str(count)
        else:
            title += '\nCount=' + str(count)

        plt.title(title)
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        encoded = base64.b64encode(img.getvalue())
        #html = '<img width=100% src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        html = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        plt.close('all')
        return html
