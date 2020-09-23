import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64

from PsuGeometry import GeoAtom as atm
from PsuGeometry import GeoDensity as den
from PsuGeometry import GeoCalcs as calcs
from PsuGeometry import GeoPlot as geop
from PsuGeometry import GeoPdb as geopdb
from PIL import Image
import numpy as np


class GeoReport:

    def __init__(self,listPdbs,pdbDataPath,edDataPath,outDataPath):
        self.pdbs = []
        self.pdbDataPath = pdbDataPath
        self.edDataPath = edDataPath
        self.outDataPath = outDataPath
        pdbmanager = geopdb.GeoPdbs(pdbDataPath,edDataPath)
        for pdbCode in listPdbs:
            self.pdbs.append(pdbmanager.getPdb(pdbCode))
        self.plots = []

    def addHistogram(self,geoX='',data=None,title='',ghost=False,operation='',splitKey='',hue=''):
        isNew = False
        if data is None:
            isNew=True
        if hue=='':
            hue='pdbCode'
        gp = geop.GeoPlot(data,geoX,geoY='',title=title,newData=isNew,operation=operation,splitKey=splitKey,plot='histogram',hue=hue)
        if not ghost:
            self.plots.append(gp)
        else:
            self.plots.append(geop.GeoOverlay(gp,'',title='ghost', pdbDataPath=self.pdbDataPath, edDataPath=self.edDataPath))

    def addScatter(self,geoX='',geoY='',data=None,title='',ghost=False,operation='',splitKey='',hue='bfactor',palette='viridis_r',centre=False,vmin=0,vmax=0,categorical=False):
        isNew = False
        if data is None:
            isNew = True
        gp = geop.GeoPlot(data, geoX, geoY=geoY, title=title, newData=isNew, operation=operation,splitKey=splitKey,hue=hue,palette=palette,centre=centre,vmin=vmin,vmax=vmax,categorical=categorical,plot='scatter')
        if not ghost:
            self.plots.append(gp)
        else:
            self.plots.append(
                geop.GeoOverlay(gp, '', title='ghost', pdbDataPath=self.pdbDataPath, edDataPath=self.edDataPath))

    def addProbability(self,geoX='',geoY='',data=None,title='',ghost=False,operation='',splitKey='',hue='bfactor',palette='viridis_r',centre=False,vmin=0,vmax=0,categorical=False):
        isNew = False
        if data is None:
            isNew = True
        gp = geop.GeoPlot(data, geoX, geoY=geoY, title=title, newData=isNew, operation=operation,splitKey=splitKey,hue=hue,palette=palette,centre=centre,vmin=vmin,vmax=vmax,categorical=categorical,plot='probability')
        if not ghost:
            self.plots.append(gp)
        else:
            self.plots.append(
                geop.GeoOverlay(gp, '', title='ghost', pdbDataPath=self.pdbDataPath, edDataPath=self.edDataPath))

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
        elif reportName == 'MainChainHistograms':
            calcList = ['C-1:N','N:CA','CA:C','C:N+1','C-1:N:CA','N:CA:C','CA:C:N+1','C:O','CA-1:CA','CA:CA+1','CA:C:N+1:CA+1','C-1:N:CA:C','N:CA:C:N+1']
        elif reportName == 'OmegaCis':
            calcList = ['CA-1:CA','CA:CA+1','CA:C:N+1:CA+1','CA-1:C-1:N:CA','N:CA:C']
        elif reportName == 'RachelsChoice' or reportName == 'RachelsChoiceNonXRay':
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
            self.printToHtml(printList, self.pdbs, title, cols, printPath, fileName)
        elif reportName == 'RachelsChoice' or reportName == 'RachelsChoiceNonXRay' :
            atomData = self.getReportCsv(reportName)
            # We want the dummy trace correlation plot so we can see if there are areas of interest
            title = "Rachel's Choice of Correlations"
            cols = 4
            printList = []

            densityHue = '2FoFc'
            if reportName == 'RachelsChoiceNonXRay':
                densityHue = 'bfactor'

            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='dssp',palette='Set1',newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue=densityHue, palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='N:CA:C:N+1', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue=densityHue, palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:CB:CG', geoY='CA:CB:CG:CD', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue='dssp', palette='Set1',newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA', geoY='CA:C', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA:CA+1', geoY='CA-1:CA', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='N:CA:C', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:O', geoY='CB:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CB:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='N:CA:C:O', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue=densityHue, palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'N:CA:C:N+1', geoY='CA-1:CA:CA+1', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:C', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA:C', geoY='C-1:CB', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue=densityHue, palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA:C:N+1:CA+1', geoY='CA-1:C-1:N:CA', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue=densityHue, palette='cubehelix_r',newData=True))
            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='aa', palette='nipy_spectral',newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'CA-2:CA-1:CA', geoY='CA:CA+1:CA+2', title='', hue='pdbCode', palette='Set1',newData=True,categorical=True))

            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='dssp', palette='Set1', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue=densityHue, palette='cubehelix_r', newData=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='aa', palette='nipy_spectral', newData=True,categorical=True))
            printList.append(geop.GeoPlot(None, 'C-1:N:CA', geoY='CA:C:N+1', title='', hue='pdbCode', palette='Set1', newData=True,categorical=True))

            self.printToHtml(printList, self.pdbs, title, cols, printPath, fileName)

        elif reportName == 'MainChainHistograms': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Sp2 Planarity'
            cols = 3
            printList = []

            printList.append(geop.GeoPlot(atomData,'C-1:N',title='C-1:N',hue='rid'))
            printList.append(geop.GeoPlot(atomData,'N:CA',title='O:C:N+1',hue='rid'))
            printList.append(geop.GeoPlot(atomData,'CA:C',title='N+1:C:CA',hue='rid'))

            printList.append(geop.GeoPlot(atomData, 'C:O', title='C:0', hue='rid'))
            printList.append(geop.GeoPlot(atomData, 'CA-1:CA', title='CA-1:CA', hue='rid'))
            printList.append(geop.GeoPlot(atomData, 'CA:CA+1', title='CA:CA+1', hue='rid'))

            printList.append(geop.GeoPlot(atomData, 'C-1:N:CA', title='Tau-1', hue='rid'))
            printList.append(geop.GeoPlot(atomData, 'N:CA:C', title='Tau', hue='rid'))
            printList.append(geop.GeoPlot(atomData, 'CA:C:N+1', title='Tau+1', hue='rid'))


            printList.append(geop.GeoPlot(atomData, 'C-1:N:CA:C', title='PHI', hue='rid'))
            printList.append(geop.GeoPlot(atomData, 'N:CA:C:N+1', title='PSI', hue='rid'))
            printList.append(geop.GeoPlot(atomData, 'CA:C:N+1:CA+1', title='AbsVal OMEGA', hue='rid', operation='ABS'))

            self.printToHtml(printList, self.pdbs, title, cols, printPath, fileName)

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

            self.printToHtml(printList, self.pdbs, title, cols, printPath, fileName)

        elif reportName == 'DataPerPdb':
            for apdb in self.pdbs:
                print('\tPSU:', reportName, 'for', apdb.pdbCode)
                atomData = apdb.dataFrame
                title = 'General Data Report'
                cols = 3
                self.addScatter(data=atomData, geoX='atomNo', geoY='aa', hue='aa', categorical=True,palette='nipy_spectral')
                self.addScatter(data=atomData, geoX='atomNo', geoY='dssp',hue= 'aa',categorical=True,palette='nipy_spectral')
                self.addScatter(data=atomData, geoX='2FoFc', geoY='bfactor',hue= 'element',palette='jet',categorical=True)
                self.addScatter(data=atomData, geoX='atomNo', geoY='bfactor',hue= 'element',palette='jet',categorical=True)
                self.addScatter(data=atomData, geoX='atomNo', geoY='2FoFc',hue='element',palette='jet',categorical=True)
                self.addScatter(data=atomData, geoX='atomNo', geoY='FoFc',hue='element',palette='jet',categorical=True)
                self.addScatter(data=atomData, geoX='x', geoY='y',hue='atomNo', palette='plasma_r')
                self.addScatter(data=atomData, geoX='y', geoY='z',hue='atomNo', palette='plasma_r')
                self.addScatter(data=atomData, geoX='z', geoY='x',hue='atomNo', palette='plasma_r')
                self.printToHtml(title, cols, fileName + '_' + apdb.pdbCode)

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
                    self.printToHtml(printList, [apdb], title, cols, printPath, fileName + '_' + apdb.pdbCode)
                else:
                    print('\tPSU:',apdb.pdbCode,'has no density matrix')




    #def printCsvToHtml(self, queryList,pdbList,title,cols,printPath,fileName):
    def printToHtml(self, title, cols, fileName):
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

        if len(self.pdbs) > 0:
            html += '<table><tr><td>PdbCode</td><td>Resolution</td><td>Pdb Link</td><td>PDBe Link</td></tr>\n'
            for apdb in self.pdbs:
                html += '<tr>\n'
                html += '<td>' + apdb.pdbCode + '</td>\n'
                html += '<td>' + str(apdb.atoms[0].values['resolution']) + '</td>\n'
                html += "<td><a href='https://www.rcsb.org/structure/" + apdb.pdbCode + "' title='PDB Link' target='_blank'>Link to PDB</a></td>\n"
                html += "<td><a href='https://www.ebi.ac.uk/pdbe/entry/pdb/" + apdb.pdbCode + "' title='PDB Link' target='_blank'>Link to PDBe</a></td>\n"
                html += '</tr>\n'
            html += '</table>\n'


        html += '<hr/>\n'

        reportPath = self.outDataPath + fileName + ".html"

        count = 0
        #html += '<table style="width:90%">\n'
        html += '<table>\n'
        row = 1
        for geoPl in self.plots:
            if not type(geoPl) is geop.GeoOverlay:
                #html += self.twoPlots(geoPl.plotA,geoPl.plotB,width)

                title = geoPl.title
                alldata = geoPl.data
                geoX = geoPl.geoX
                geoY = geoPl.geoY
                splitKey = geoPl.splitKey
                splitList = ['']
                if splitKey != '':
                    splitList = alldata[splitKey].unique()
            else:
                splitList = ['']


            if True:

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

                    if type(geoqSplit) is geop.GeoOverlay:
                        html += self.twoPlotsOverlay(geoPl.plotA,geoPl.plotB,width)
                    else:
                        html += self.onePlot(geoqSplit, width)

        html += '</tr></table><hr/><p>Produced by PsuGeometry, written by Rachel Alcraft </p></body>\n'

        # and print
        f = open(reportPath, "w+")
        f.write(html)
        print('PSU: saved file to',reportPath)
        f.close()

    def onePlot(self,geoPl,width):
        #try:
        if True:
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

            if geoPl.data.empty:
                html = '<td width=' + width + '%>' + 'No Data:' + geoPl.geoX + ' ' + geoPl.geoY  + '</td>\n'
            else:
                fig, ax = plt.subplots()
                ret = geoPl.plotToAxes(fig, ax)
                encoded = geoPl.getPlotImage(fig, ax)
                htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
                htmlstring += ret

                html = '<td width=' + width + '%>' + htmlstring + '</td>\n'
        #except:
        #    html = '<td width=' + width + '%>' + 'Error:' + geoPl.geoX + ' ' + geoPl.geoY + '</td>\n'

        return (html)


    def twoPlotsOverlay(self,geoPlA,geoPlB,width):
        '''
        https://stackoverflow.com/questions/6871201/plot-two-histograms-on-single-chart-with-matplotlib
        '''
        #try:
        if True:
            fig, ax = plt.subplots()
            if geoPlA.newData:
                geoPlA.getNewData(self.pdbs)
            if geoPlB.newData:
                geoPlB.getNewData(self.pdbs)

            retA = geoPlA.plotToAxes(fig, ax)
            retB = geoPlB.plotToAxes(fig, ax)
            encoded = geoPlA.getPlotImage(fig, ax)

            htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
            htmlstring += retA + retB

            if geoPlA.data.empty:
                html = '<td width=' + width + '%>' + 'No Data:' + geoPlA.geoX + ' ' + geoPlA.geoY  + '</td>\n'
            else:
                html = '<td width=' + width + '%>' + htmlstring + '</td>\n'
        #except:
        #    html = '<td width=' + width + '%>' + 'Error:' + geoPlA.geoX + ' ' + geoPlA.geoY + '</td>\n'

        return (html)




