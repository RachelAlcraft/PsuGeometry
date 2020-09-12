import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64

from PsuGeometry import GeoAtom as atm
from PsuGeometry import GeoDensity as den
from PsuGeometry import GeoCalcs as calcs
from PsuGeometry import GeoPlot as geop


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
            printList.append(geop.GeoPlot(atomData,'C-1:N',geoY='N:CA',title='Bonds', splitKey=''))
            printList.append(geop.GeoPlot(atomData, 'CA:C', geoY='C:N+1', title='Bonds'))
            printList.append(geop.GeoPlot(atomData, 'C-1:N', geoY='C:N+1', title='Bonds'))
            printList.append(geop.GeoPlot(atomData, 'C-1:N:CA', geoY='N:CA:C', title='Angles'))
            printList.append(geop.GeoPlot(atomData, 'N:CA:C', geoY='CA:C:N+1', title='Angles'))
            printList.append(geop.GeoPlot(atomData, 'C-1:N:CA', geoY='CA:C:N+1', title='Angles'))
            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)
        elif reportName == 'RachelsChoice':
            atomData = self.getReportCsv(reportName)
            title = "Rachel's Choice of Correlations"
            cols = 4
            printList = []

            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='dssp',palette='Set1',newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue='dssp', palette='Set1',newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='2FoFc', palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='2FoFc', palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)


        elif reportName == 'Sp2Planarity': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Sp2 Planarity'
            cols = 4
            printList = []
            printList.append(geop.GeoPlot(atomData,'N+1:O:C:CA',title='AbsVal Dihedral',hue='rid',operation='ABS'))
            printList.append(geop.GeoPlot(atomData,'CA:C:O',title='CA:C:O',hue='rid'))
            printList.append(geop.GeoPlot(atomData,'O:C:N+1',title='O:C:N+1',hue='rid'))
            printList.append(geop.GeoPlot(atomData,'N+1:C:CA',title='N+1:C:CA',hue='rid'))

            printList.append(geop.GeoPlot(atomData, 'N+1:O:C:CA', geoY='CA:C:O',hue='dssp'))
            printList.append(geop.GeoPlot(atomData, 'N+1:O:C:CA', geoY='O:C:N+1'))
            printList.append(geop.GeoPlot(atomData, 'N+1:O:C:CA', geoY='N+1:C:CA',hue='bfactor'))
            printList.append(geop.GeoPlot(atomData, 'N+1:C:CA', geoY='O:C:N+1',hue='FoFc', palette='Spectral',centre=True))

            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)

        elif reportName == 'DataPerPdb':
            for apdb in self.pdbs:
                print('\tPSU:', reportName, 'for', apdb.pdbCode)
                atomData = apdb.dataFrame
                title = 'General Data Report'
                cols = 3
                printList = []
                printList.append(geop.GeoPlot(atomData, 'atomNo', geoY='aa', hue='aa', categorical=True,palette='nipy_spectral'))
                printList.append(geop.GeoPlot(atomData, 'atomNo', geoY='dssp',hue= 'aa',categorical=True,palette='nipy_spectral'))
                printList.append(geop.GeoPlot(atomData, '2FoFc', geoY='bfactor',hue= 'element',palette='jet',categorical=True))
                printList.append(geop.GeoPlot(atomData, 'atomNo', geoY='bfactor',hue= 'element',palette='jet',categorical=True))
                printList.append(geop.GeoPlot(atomData, 'atomNo', geoY='2FoFc',hue='element',palette='jet',categorical=True))
                printList.append(geop.GeoPlot(atomData, 'atomNo', geoY='FoFc',hue='element',palette='jet',categorical=True))
                printList.append(geop.GeoPlot(atomData, 'x', geoY='y',hue='atomNo', palette='plasma_r'))
                printList.append(geop.GeoPlot(atomData, 'y', geoY='z',hue='atomNo', palette='plasma_r'))
                printList.append(geop.GeoPlot(atomData, 'z', geoY='x',hue='atomNo', palette='plasma_r'))
                self.printCsvToHtml(printList, [apdb], title, cols, printPath, fileName + '_' + apdb.pdbCode)

        elif reportName == 'Slow_DensityPointsPerPdb' or reportName == 'Slow_DensityPeaksPerPdb': # this can only be done per pdb
            for apdb in self.pdbs:
                if apdb.hasDensity:
                    print('\tPSU:', reportName, 'for', apdb.pdbCode)
                    allPoints = True
                    title = 'Density Points and Atoms Comparison'
                    if reportName == 'Slow_DensityPeaksPerPdb':
                        allPoints = False
                        title = 'Density Peaks and Atoms Comparison'
                    peaksData = apdb.geoDen.getPeaks(allPoints)
                    atomData = apdb.dataFrame
                    atomData['FoFc2'] = atomData['FoFc'] ** 2

                    cols = 3
                    printList = []
                    printList.append(geop.GeoPlot(peaksData, 'c', geoY='r', title='Density CR Fo',hue='peakFo',palette='gist_gray_r'))
                    printList.append(geop.GeoPlot(peaksData, 'r', geoY='s',title='Density RS Fo',hue='peakFo',palette='gist_gray_r'))
                    printList.append(geop.GeoPlot(peaksData, 's', geoY='c',title='Density SC Fo',hue='peakFo',palette='gist_gray_r'))

                    printList.append(geop.GeoPlot(peaksData, 'c', geoY='r',title='Density CR Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 'r', geoY='s',title='Density RS Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 's', geoY='c',title='Density SC Fo',hue='peakFo',palette='cubehelix_r'))

                    printList.append(geop.GeoPlot(peaksData, 'c', geoY='r',title='Density CR 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 'r', geoY='s',title='Density RS 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 's', geoY='c',title='Density SC 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))

                    printList.append(geop.GeoPlot(peaksData, 'c', geoY='r',title='Density CR FC',hue='peakFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 'r', geoY='s',title='Density RS FC',hue='peakFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 's', geoY='c',title='Density SC FC',hue='peakFc',palette='cubehelix_r'))

                    printList.append(geop.GeoPlot(peaksData, 'c', geoY='r',title='Density CR FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(geop.GeoPlot(peaksData, 'r', geoY='s',title='Density RS FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(geop.GeoPlot(peaksData, 's', geoY='c',title='Density SC FoFC',hue='peakFoFc',palette='PiYG',centre=True))

                    printList.append(geop.GeoPlot(peaksData, 'x', geoY='y',title='Density XY Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 'y', geoY='z',title='Density YZ Fo',hue='peakFo',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 'z', geoY='x',title='Density ZX Fo',hue='peakFo',palette='cubehelix_r'))

                    printList.append(geop.GeoPlot(peaksData, 'x', geoY='y',title='Density XY FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(geop.GeoPlot(peaksData, 'y', geoY='z',title='Density YZ FoFC',hue='peakFoFc',palette='PiYG',centre=True))
                    printList.append(geop.GeoPlot(peaksData, 'z', geoY='x', title='Density ZX FoFC',hue='peakFoFc',palette='PiYG',centre=True))

                    printList.append(geop.GeoPlot(peaksData, 'x', geoY='y',title='Density XY 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 'y', geoY='z',title='Density YZ 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(peaksData, 'z', geoY='x',title='Density ZX 2FoFC',hue='peak2FoFc',palette='cubehelix_r'))

                    printList.append(geop.GeoPlot(atomData, 'x', geoY='y',title='PDB XY 2FoFc',hue='2FoFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(atomData, 'y', geoY='z',title='PDB YZ 2FoFc',hue='2FoFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(atomData, 'z', geoY='x',title='PDB ZX 2FoFc',hue='2FoFc',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(atomData, 'x', geoY='y',title='PDB XY Electrons',hue='electrons',palette='cubehelix_r',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'y', geoY='z',title='PDB YZ Electrons',hue='electrons',palette='cubehelix_r',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'z', geoY='x',title='PDB ZX Electrons',hue='electrons',palette='cubehelix_r',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'x', geoY='y',title='PDB XY FoFc',hue='FoFc',palette='PiYG',centre=True))
                    printList.append(geop.GeoPlot(atomData, 'y', geoY='z',title='PDB YZ FoFc',hue='FoFc',palette='PiYG',centre=True))
                    printList.append(geop.GeoPlot(atomData, 'z', geoY='x',title='PDB ZX FoFc',hue='FoFc',palette='PiYG',centre=True))
                    printList.append(geop.GeoPlot(atomData, 'x', geoY='y',title='PDB XY bfactor',hue='bfactor',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(atomData, 'y', geoY='z',title='PDB YZ bfactor',hue='bfactor',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(atomData, 'z', geoY='x',title='PDB ZX bfactor',hue='bfactor',palette='cubehelix_r'))
                    printList.append(geop.GeoPlot(atomData, 'x', geoY='y',title='PDB XY atom nos',hue='atomNo',palette='gist_ncar'))
                    printList.append(geop.GeoPlot(atomData, 'y', geoY='z',title='PDB YZ atom nos',hue='atomNo',palette='gist_ncar'))
                    printList.append(geop.GeoPlot(atomData, 'z', geoY='x',title='PDB ZX atom nos',hue='atomNo',palette='gist_ncar'))
                    printList.append(geop.GeoPlot(atomData, 'x', geoY='y',title='PDB XY amino acids',hue='aa',palette='nipy_spectral',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'y', geoY='z',title='PDB YZ amino acids',hue='aa',palette='nipy_spectral',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'z', geoY='x',title='PDB ZX amino acids',hue='aa',palette='nipy_spectral',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'bfactor', geoY='FoFc2',title='PDB bfactor vs FoFc^2',hue='electrons',palette='viridis_r',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'Fc', geoY='Fo',title='PDB Fc vs Fc',hue='electrons',palette='viridis_r',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'electrons', geoY='2FoFc',title='PDB electrons vs 2FoFc',hue='element',palette='viridis_r',categorical=True))
                    printList.append(geop.GeoPlot(atomData, 'aa',title='Amino Acids'))
                    printList.append(geop.GeoPlot(atomData, 'element',title='Atoms'))
                    printList.append(geop.GeoPlot(atomData, '2FoFc',title='Peaks in 2FoFc'))
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
        row = 1
        for geoPl in queryList:

            title = geoPl.title
            alldata = geoPl.data
            geoX = geoPl.geoX
            geoY = geoPl.geoY
            splitKey = geoPl.splitKey
            splitList = ['']
            if splitKey != '':
                splitList = alldata[splitKey].unique()

            for split in splitList:
                geoqSplit = geoPl
                if count == 0:
                    html += '<tr><td colspan=' + '"' + str(cols) + '"' + '>Row ' + str(row) + '</td></tr>\n'
                    row += 1
                    html += '<tr>'
                elif count % cols == 0:
                    html += '</tr>\n'
                    html += '<tr><td colspan=' + '"' + str(cols) + '"' + '>Row ' + str(row) + '</td></tr>\n'
                    row += 1
                    html += '<tr>'


                count += 1

                if split != '':
                    data = alldata[alldata[splitKey] == split]
                    geoqSplit.data = data
                    geoqSplit.title = title + ' ' + split

                html += self.onePlot(geoqSplit,width)


        html += '</tr></table><hr/><p>Produced by PsuGeometry, written by Rachel Alcraft </p></body>\n'

        # and print
        f = open(reportPath, "w+")
        f.write(html)
        print('PSU: saved file to',reportPath)
        f.close()

    def onePlot(self,geoPl,width):
        if geoPl.newData:
            if geoPl.plot == 'histogram':

                calcList = [geoPl.geoX]
                hueList = [geoPl.hue]
                dfs = []
                for apdb in self.pdbs:
                    data = apdb.getGeoemtryCsv(calcList, hueList)
                    dfs.append(data)
                geoPl.data = pd.concat(dfs, ignore_index=True)
            else:
                calcList = [geoPl.geoX, geoPl.geoY]
                hueList = [geoPl.hue]
                dfs = []
                for apdb in self.pdbs:
                    datatmp = apdb.getGeoemtryCsv(calcList, hueList)
                    dfs.append(datatmp)
                geoPl.data = pd.concat(dfs, ignore_index=True)

        html = '<td width=' + width + '%>' + geoPl.getPlot() + '</td>\n'
        return (html)




