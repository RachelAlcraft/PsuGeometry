import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64

from PsuGeometry import GeoAtom as atm
from PsuGeometry import GeoDensity as den
from PsuGeometry import GeoCalcs as calcs


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
            calcList = ['CA:C:O','O:C:N+1','N+1:C:CA']
        elif reportName == 'Sp3Tetrahedra':
            calcList = ['N:CA:C','C:CA:CB','N:CA:CB']
        elif reportName == 'BackboneOutliers':
            calcList = ['C-1:N','N:CA','CA:C','C:N+1','C-1:N:CA','N:CA:C','CA:C:N+1']
        elif reportName == 'OmegaCis':
            calcList = ['CA-1:CA','CA:CA+1','CA:C:N+1:CA+1','CA-1:C-1:N:CA','N:CA:C']
        elif reportName == 'RachelsChoice':
            calcList = ['N:O','CB:O','N:CA:C:N+1']
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
            printList.append(['Bonds', atomData, 'C-1:N', 'N:CA','aa', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['Bonds', atomData, 'CA:C', 'C:N+1','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['Bonds', atomData, 'C-1:N', 'C:N+1','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['Angles', atomData, 'C-1:N:CA', 'N:CA:C','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['Angles', atomData, 'N:CA:C', 'CA:C:N+1','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['Angles', atomData, 'C-1:N:CA', 'CA:C:N+1','', '2FoFc', 'viridis_r', False, 0, 0])
            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)
        elif reportName == 'RachelsChoice': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = "Rachel's Choice of Correlations"
            cols = 5
            printList = []
            printList.append(['2FoFc', atomData, 'N:CA:C:N+1', 'N:O','pdbCode', '2FoFc', 'plasma_r', False, 0, 0])
            printList.append(['BFactor', atomData, 'N:CA:C:N+1', 'N:O','', 'bfactor', 'plasma', False, 0, 0])
            printList.append(['Amino Acid', atomData, 'N:CA:C:N+1', 'N:O','', 'aa', 'gist_ncar', False, 0, 0])
            printList.append(['rid', atomData, 'N:CA:C:N+1', 'N:O', '','rid', 'gist_rainbow', False, 0, 0])
            printList.append(['FoFc', atomData, 'N:CA:C:N+1', 'N:O', '','FoFc', 'Spectral', True, 0, 0])

            printList.append(['2FoFc', atomData, 'N:CA:C:N+1', 'CB:O','', '2FoFc', 'plasma_r', False, 0, 0])
            printList.append(['BFactor', atomData, 'N:CA:C:N+1', 'CB:O','', 'bfactor', 'plasma', False, 0, 0])
            printList.append(['Amino Acid', atomData, 'N:CA:C:N+1', 'CB:O','', 'aa', 'gist_ncar', False, 0, 0])
            printList.append(['rid', atomData, 'N:CA:C:N+1', 'CB:O','', 'rid', 'gist_rainbow', False, 0, 0])
            printList.append(['FoFc', atomData, 'N:CA:C:N+1', 'CB:O','', 'FoFc', 'Spectral', True, 0, 0])

            printList.append(['2FoFc', atomData, 'N:O', 'CB:O','', '2FoFc', 'plasma_r', False, 0, 0])
            printList.append(['BFactor', atomData, 'N:O', 'CB:O','', 'bfactor', 'plasma', False, 0, 0])
            printList.append(['Amino Acid', atomData, 'N:O', 'CB:O','', 'aa', 'gist_ncar', False, 0, 0])
            printList.append(['rid', atomData, 'N:O', 'CB:O','', 'rid', 'gist_rainbow','', False, 0, 0])
            printList.append(['FoFc', atomData, 'N:O', 'CB:O','', 'FoFc', 'Spectral', True, 0, 0])
            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)

        elif reportName == 'OmegaCis': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Omega Cis Report'
            cols = 4
            printList = []

            printList.append(['Pre Omega Tau 2FoFc', atomData, 'CA-1:C-1:N:CA', 'N:CA:C', '', 'rid', 'gist_rainbow', False, 0, 0])
            printList.append(['Pre Omega Tau Amino Acid', atomData, 'CA-1:C-1:N:CA', 'N:CA:C', '', 'aa', 'tab20', False, 0, 0])
            printList.append(['Post Omega Tau 2FoFc', atomData, 'CA:C:N+1:CA+1', 'N:CA:C', '', 'rid', 'gist_rainbow', False, 0, 0])
            printList.append(['Post Omega Tau Amino Acid', atomData, 'CA:C:N+1:CA+1', 'N:CA:C', '', 'aa', 'tab20', False, 0, 0])

            printList.append(['CA Report 2FoFC', atomData, 'CA-1:CA', 'CA:CA+1','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['CA Report Amino Acid', atomData, 'CA-1:CA', 'CA:CA+1','', 'aa', 'tab20', False, 0, 0])
            printList.append(['CA Report bfactor', atomData, 'CA-1:CA', 'CA:CA+1','', 'bfactor', 'viridis_r', False, 0, 0])
            printList.append(['CA Report FoFc', atomData, 'CA-1:CA', 'CA:CA+1','', 'FoFc', 'PiYG', True, 0, 0])

            printList.append(['Pre Omega 2FoFc', atomData, 'CA-1:C-1:N:CA', 'CA-1:CA','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['Pre Omega Amino Acid', atomData, 'CA-1:C-1:N:CA', 'CA-1:CA','', 'aa', 'tab20', False, 0, 0])
            printList.append(['Pre Omega bfactor', atomData, 'CA-1:C-1:N:CA', 'CA-1:CA','', 'bfactor', 'viridis_r', False, 0, 0])
            printList.append(['Pre Omega FoFc', atomData, 'CA-1:C-1:N:CA', 'CA-1:CA','', 'FoFc', 'PiYG', True, 0, 0])

            printList.append(['Post Omega 2FoFc', atomData, 'CA:C:N+1:CA+1', 'CA:CA+1','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['Post Omega Amino Acid', atomData, 'CA:C:N+1:CA+1', 'CA:CA+1','', 'aa', 'tab20', False, 0, 0])
            printList.append(['Post Omega bfactor', atomData, 'CA:C:N+1:CA+1', 'CA:CA+1','', 'bfactor', 'viridis_r', False, 0, 0])
            printList.append(['Post Omega FoFc', atomData, 'CA:C:N+1:CA+1', 'CA:CA+1','', 'FoFc', 'PiYG', True, 0, 0])

            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)

        elif reportName == 'Ramachandran': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Ramachandran Report'
            cols = 2
            printList = []
            printList.append(['2FoFc', atomData, 'C-1:N:CA:C', 'N:CA:C:N+1','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['FoFc', atomData, 'C-1:N:CA:C', 'N:CA:C:N+1','', 'FoFc', 'PiYG', True, 0, 0])
            printList.append(['BFactor', atomData, 'C-1:N:CA:C', 'N:CA:C:N+1','', 'bfactor', 'viridis_r', False, 0, 0])
            printList.append(['Residue No', atomData, 'C-1:N:CA:C', 'N:CA:C:N+1','', 'rid', 'viridis_r', False, 0, 0])
            printList.append(['Amino Acid', atomData, 'C-1:N:CA:C', 'N:CA:C:N+1','', 'aa', 'tab20', False, 0, 0])
            printList.append(['DSSP', atomData, 'C-1:N:CA:C', 'N:CA:C:N+1','dssp', 'aa', 'tab20', False, 0, 0])
            printList.append(['PHI', atomData, 'C-1:N:CA:C','', 'dssp'])
            printList.append(['PSI', atomData, 'N:CA:C:N+1', ''])
            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)
        elif reportName == 'Sp2Planarity': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Sp2 Planarity'
            cols = 3
            printList = []
            printList.append(['CA:C:O', atomData, 'CA:C:O','', ''])
            printList.append(['O:C:N+1', atomData, 'O:C:N+1','', ''])
            printList.append(['N+1:C:CA', atomData, 'N+1:C:CA','', 'aa'])
            printList.append(['', atomData, 'CA:C:O', 'O:C:N+1','', '2FoFc', 'viridis_r', False, 0, 0])
            printList.append(['', atomData, 'N+1:C:CA', 'CA:C:O','', 'dssp', 'tab20', False, 0, 0])
            printList.append(['', atomData, 'O:C:N+1', 'N+1:C:CA','', 'bfactor', 'cubehelix_r', False, 0, 0])
            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)
        elif reportName == 'Sp3Tetrahedra': # Sp2Planarity, DensityAtomCompare, OmegaCis
            atomData = self.getReportCsv(reportName)
            title = 'Sp3 Tetrahedra'
            cols = 3
            printList = []
            printList.append(['N:CA:C', atomData, 'N:CA:C','', ''])
            printList.append(['N:CA:CB', atomData, 'N:CA:CB','', ''])
            printList.append(['C:CA:CB', atomData, 'C:CA:CB','', ''])
            self.printCsvToHtml(printList, self.pdbs, title, cols, printPath, fileName)
        elif reportName == 'DataPerPdb':
            for apdb in self.pdbs:
                print('\tPSU:', reportName, 'for', apdb.pdbCode)
                atomData = apdb.dataFrame
                atomData = atomData.sort_values(by='2FoFc', ascending=True)
                title = 'General Data Report'
                cols = 2
                printList = []
                printList.append(['', atomData, 'atomNo', 'aa','', 'aa', 'tab20', False, 0, 0])
                printList.append(['', atomData, '2FoFc', 'bfactor','', 'element', 'tab20', False, 0, 0])
                printList.append(['', atomData, 'atomNo', 'bfactor','', 'element', 'tab20', False, 0, 0])
                printList.append(['', atomData, 'atomNo', '2FoFc','', 'element', 'tab20', False, 0, 0])
                printList.append(['', atomData, 'atomNo', 'FoFc','', 'element', 'tab20', False, 0, 0])
                printList.append(['', atomData, 'x', 'y','', 'atomNo', 'plasma_r', False, 0, 0])
                printList.append(['', atomData, 'y', 'z','', 'atomNo', 'plasma_r', False, 0, 0])
                printList.append(['', atomData, 'z', 'x','', 'atomNo', 'plasma_r', False, 0, 0])
                self.printCsvToHtml(printList, [apdb], title, cols, printPath, fileName + '_' + apdb.pdbCode)

        elif reportName == 'DensityPerPdb' or reportName == 'DensityPeaksPerPdb': # this can only be done per pdb
            for apdb in self.pdbs:
                if apdb.hasDensity:
                    print('\tPSU:', reportName, 'for', apdb.pdbCode)
                    allPoints = True
                    if reportName == 'DensityPeaksPerPdb':
                        allPoints = False
                    peaksData = apdb.geoDen.getPeaks(allPoints)
                    atomData = apdb.dataFrame
                    atomData['FoFc2'] = atomData['FoFc'] ** 2
                    title = 'Density and Atom Comparison'
                    cols = 3
                    printList = []
                    printList.append(['Density CR Fo', peaksData, 'c', 'r', '', 'peakFo', 'gist_gray_r', False, 0, 0])
                    printList.append(['Density RS Fo', peaksData, 'r', 's', '', 'peakFo', 'gist_gray_r', False, 0, 0])
                    printList.append(['Density SC Fo', peaksData, 's', 'c', '', 'peakFo', 'gist_gray_r', False, 0, 0])

                    printList.append(['Density CR Fo', peaksData, 'c', 'r', '', 'peakFo', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density RS Fo', peaksData, 'r', 's', '', 'peakFo', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density SC Fo', peaksData, 's', 'c', '', 'peakFo', 'cubehelix_r', False, 0, 0])

                    printList.append(['Density CR 2FoFC', peaksData, 'c', 'r','', 'peak2FoFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density RS 2FoFC', peaksData, 'r', 's','', 'peak2FoFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density SC 2FoFC', peaksData, 's', 'c','', 'peak2FoFc', 'cubehelix_r', False, 0, 0])

                    printList.append(['Density CR FC', peaksData, 'c', 'r','', 'peakFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density RS FC', peaksData, 'r', 's','', 'peakFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density SC FC', peaksData, 's', 'c','', 'peakFc', 'cubehelix_r', False, 0, 0])

                    printList.append(['Density CR FoFC', peaksData, 'c', 'r', '', 'peakFoFc', 'PiYG', True, 0, 0])
                    printList.append(['Density RS FoFC', peaksData, 'r', 's', '', 'peakFoFc', 'PiYG', True, 0, 0])
                    printList.append(['Density SC FoFC', peaksData, 's', 'c', '', 'peakFoFc', 'PiYG', True, 0, 0])

                    printList.append(['Density XY Fo', peaksData, 'x', 'y', '', 'peakFo', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density YZ Fo', peaksData, 'y', 'z', '', 'peakFo', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density ZX Fo', peaksData, 'z', 'x', '', 'peakFo', 'cubehelix_r', False, 0, 0])

                    printList.append(['Density XY FoFC', peaksData, 'x', 'y', '', 'peakFoFc', 'PiYG', True, 0, 0])
                    printList.append(['Density YZ FoFC', peaksData, 'y', 'z', '', 'peakFoFc', 'PiYG', True, 0, 0])
                    printList.append(['Density ZX FoFC', peaksData, 'z', 'x', '', 'peakFoFc', 'PiYG', True, 0, 0])

                    printList.append(['Density XY 2FoFC', peaksData, 'x', 'y','', 'peak2FoFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density YZ 2FoFC', peaksData, 'y', 'z','', 'peak2FoFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['Density ZX 2FoFC', peaksData, 'z', 'x','', 'peak2FoFc', 'cubehelix_r', False, 0, 0])

                    printList.append(['PDB XY 2FoFc', atomData, 'x', 'y','', '2FoFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB YZ 2FoFc', atomData, 'y', 'z','', '2FoFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB ZX 2FoFc', atomData, 'z', 'x','', '2FoFc', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB XY Electrons', atomData, 'x', 'y','', 'electrons', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB YZ Electrons', atomData, 'y', 'z','', 'electrons', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB ZX Electrons', atomData, 'z', 'x','', 'electrons', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB XY FoFc', atomData, 'x', 'y','', 'FoFc', 'PiYG', True, 0, 0])
                    printList.append(['PDB YZ FoFc', atomData, 'y', 'z','', 'FoFc', 'PiYG', True, 0, 0])
                    printList.append(['PDB ZX FoFc', atomData, 'z', 'x','', 'FoFc', 'PiYG', True, 0, 0])
                    printList.append(['PDB XY bfactor', atomData, 'x', 'y','', 'bfactor', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB YZ bfactor', atomData, 'y', 'z','', 'bfactor', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB ZX bfactor', atomData, 'z', 'x','', 'bfactor', 'cubehelix_r', False, 0, 0])
                    printList.append(['PDB bfactor vs FoFc^2', atomData, 'bfactor', 'FoFc2','', 'electrons', 'viridis_r', False, 0, 0])
                    printList.append(['PDB Fc vs Fc', atomData, 'Fc', 'Fo','', 'electrons', 'viridis_r', False, 0, 0])
                    printList.append(['PDB electrons vs 2FoFc', atomData, 'electrons', '2FoFc','', 'element', 'tab10', False, 0, 0])
                    printList.append(['Amino Acids', atomData, 'aa', '',''])
                    printList.append(['Atoms', atomData, 'element', '',''])
                    printList.append(['Peaks in 2FoFc', atomData, '2FoFc', '',''])
                    self.printCsvToHtml(printList, [apdb], title, cols, printPath, fileName + '_' + apdb.pdbCode)
                else:
                    print('\tPSU:',apdb.pdbCode,'has no density matrix')




    def printCsvToHtml(self, reportsList,pdbList,title,cols,printPath,fileName):
        html = '<!DOCTYPE html><html lang="en"><head><title>PSU-' + fileName + '-GEO</title>'
        #html += '<style> body {background-color:WhiteSmoke;} table,th,td {border:1px solid LightSteelBlue;background-color:WhiteSmoke;}</style></head>'
        html += '<style> body {background-color:SeaShell;} td {border:1px solid RosyBrown;background-color:SeaShell;}</style></head>'
        html += '<body><h1>' + title + '</h1>'
        html += '<h2>PSU: Geometric Correlations</h2>'
        html += '<hr/>'

        if len(pdbList) > 0:
            html += '<table><tr><td>PdbCode</td><td>Resolution</td><td>Pdb Link</td><td>Pdbe Link</td></tr>'
            for apdb in pdbList:
                html += '<tr>'
                html += '<td>' + apdb.pdbCode + '</td>'
                html += '<td>' + str(apdb.atoms[0].values['resolution']) + '</td>'
                html += "<td><a href='https://www.rcsb.org/structure/" + apdb.pdbCode + "' title='PDB Link' target='_blank'>Link to PDB</a></td>"
                html += "<td><a href='https://www.ebi.ac.uk/pdbe/entry/pdb/" + apdb.pdbCode + "' title='PDB Link' target='_blank'>Link to PDBe</a></td>"
                html += '</tr>'
            html += '</table>'


        html += '<hr/>'

        reportPath = printPath + fileName + ".html"

        count = 0
        html += '<table style="width:90%">'
        for report in reportsList:

            title = report[0]
            alldata = report[1]
            geoX = report[2]
            geoY = report[3]
            splitKey = ''
            splitList = ['']
            if len(report) > 4:
                splitKey = report[4]

            if splitKey != '':
                splitList = alldata[splitKey].unique()

            for split in splitList:
                if count == 0:
                    html += '<tr>'
                elif count % cols == 0:
                    html += '</tr><tr>'

                count += 1

                if split != '':
                    data = alldata[alldata[splitKey] == split]
                    sptitle = title + ' ' + split
                else:
                    data = alldata
                    sptitle = title

                if geoY == '': # then it is 1d
                    html += self.oneHTMLHistogram(sptitle,data,geoX)
                else:
                    hue = report[5]
                    palette = report[6]
                    centre = report[7]
                    vmin = report[8]
                    vmax = report[9]
                    html += self.oneHTMLCorrelation(sptitle,data,geoX,geoY,hue,palette,centre,vmin,vmax)

        html += '</tr></table><hr/><p>Produced by PsuGeometry, written by Rachel Alcraft </p><body>'

        # and print
        f = open(reportPath, "w+")
        f.write(html)
        print('PSU: saved file to',reportPath)
        f.close()

    def oneHTMLCorrelation(self,title,data,geoX,geoY,hue,palette,centre,vmin,vmax):
        html = '<td>' + self.createScatterPlot(title, geoX, geoY, hue, data, palette, centre,vmin,vmax) + '</td>'
        return (html)

    def oneHTMLHistogram(self,title,data,geoX):
        html = '<td>' + self.createHistogram(geoX, data, title) + '</td>'
        return (html)

    def createHistogram(self,xName, data, title):
        fig, ax = plt.subplots()

        # sns.distplot(data[xName], norm_hist=True, bins=50, kde=False)
        plt.hist(data[xName], EdgeColor='k', bins=50)
        plt.title(title)
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        encoded = base64.b64encode(img.getvalue())
        html = '<img width = 95% src="data:image/png;base64, {}">'.format(encoded.decode('utf-8'))
        plt.close('all')
        dfdesc = data[xName].describe()
        rows = len(dfdesc.index)
        colsNames = list(dfdesc.index)
        html += "<table>"
        html += "<tr>"
        for r in range(0, rows):
            html += "<td>" + str(colsNames[r]) + "</td>"
        html += "</tr>"
        html += "<tr>"
        for r in range(0, rows):
            html += "<td>" + str(round(dfdesc[r],4)) + "</td>"
        html += "</tr>"
        html += "</table>"

        return html

    def createScatterPlot(self,title, geoX, geoY, hue, data, palette, centre,vmin,vmax):
        fig, ax = plt.subplots()
        if centre:
            data[hue + '2'] = data[hue]**2
            data = data.sort_values(by=hue+'2', ascending=True)
            maxh = max(data[hue].max(), -1 * data[hue].min())
            minh = maxh * -1
            g = ax.scatter(data[geoX], data[geoY], c=data[hue], cmap=palette, vmin=minh, vmax=maxh)
            fig.colorbar(g)
            ax.set_xlabel(geoX)
            ax.set_ylabel(geoY)
        elif vmin < vmax:
            data = data.sort_values(by=hue, ascending=True)
            g = ax.scatter(data[geoX], data[geoY], c=data[hue], cmap=palette, vmin=vmin, vmax=vmax)
            fig.colorbar(g)
            ax.set_xlabel(geoX)
            ax.set_ylabel(geoY)
        else:
            lw = 0.5
            if palette == 'gist_gray_r':
                lw = 0 # this gives a crystollagraphic image look

            data = data.sort_values(by=hue, ascending=True)
            im = sns.scatterplot(x=geoX, y=geoY, hue=hue, data=data, alpha=0.65, palette=palette,edgecolor='aliceblue',linewidth=lw)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # Put the legend out of the figure

        count = len(data.index)
        if title == '':
            title += 'Count=' + str(count)
        else:
            title += '\nCount=' + str(count)

        plt.title(title)
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        encoded = base64.b64encode(img.getvalue())
        html = '<img width = 95% src="data:image/png;base64, {}">'.format(encoded.decode('utf-8'))
        plt.close('all')
        return html
