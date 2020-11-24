import numpy as np
import pandas as pd
from PsuGeometry import GeoPlot as geop
import matplotlib.pyplot as plt

class GeoSlice:
    '''
    A minimal cut down of GeoReport for slices only with no density, pdb or dssp code (ie for my windows machine!)
    '''

    def __init__(self):
        self.plots = []

    def save(self,dataarray, filepath):
        with open(filepath, 'w') as outfile:
            x, y = dataarray.shape
            for i in range(0,x):
                for j in range(0, y):
                    val = dataarray[i,j]
                    if j > 0:
                        outfile.write(str(','))
                    elif i > 0:
                        outfile.write(str('\n'))
                    outfile.write(str(val))

    def load(self,filepath):

        with open(filepath,'r') as f:
            ed_data = f.read().splitlines()

        rows = len(ed_data)
        ed_slice = np.zeros((rows,rows))
        for i in range(0,rows):
            row = ed_data[i].split(',')
            for j in range(0, rows):
                val = float(row[j])
                ed_slice[i,j] = val

        return ed_slice


    def addDensitySlice(self, slice, palette='viridis', title='',logged=False,centre=False):
        gp = geop.GeoPlot(data=None,geoX='',title=title, palette=palette, plot='surface', report=self,centre=centre)
        gp.surface = slice
        gp.logged=logged
        gp.differ=0
        self.plots.append(gp)


    def printToHtml(self, maintitle, cols, fileName):
        print('PSU: formatting to html...')
        width=str(100/cols)
        reportPath = fileName
        count = 0
        #html += '<table style="width:90%">\n'
        html = '<table>\n'
        row = 1
        for geoPl in self.plots:
            if type(geoPl) is geop.GeoPlot:
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

                    print('PSU: plot',count,'/',len(self.plots))
                    if type(geoqSplit) is geop.GeoOverlay:
                        html += self.twoPlotsOverlay(geoPl.plotA,geoPl.plotB,width)
                    elif type(geoqSplit) is geop.GeoDifference:
                        html += self.onePlot(geoqSplit, width)
                    else:
                        html += self.onePlot(geoqSplit, width)

        html += '</tr></table><hr/><p>Produced by PsuGeometry, written by Rachel Alcraft<br/>Please cite the application note...</p></body>\n'
        hhtml = self.getHeaderString(fileName, maintitle)
        # and print
        f = open(reportPath, "w+")
        f.write(hhtml + html)
        print('PSU: saved file to',reportPath)
        self.flush()
        f.close()

    def onePlot(self, geoPl, width):
        # try:
        if True:
            if geoPl.newData:
                geoPl.getNewData()

            geoPl.applyRestrictions()
            geoPl.applyExclusions()

            if geoPl.plot == 'surface' or geoPl.plot == 'surfaces':
                fig, ax = plt.subplots()
                ret = geoPl.plotToAxes(fig, ax)
                encoded = geoPl.getPlotImage(fig, ax)
                htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
                htmlstring += ret
                html = '<td width=' + width + '%>' + htmlstring + '</td>\n'
            elif geoPl.hasMatrix:
                fig, ax = plt.subplots()
                ret = geoPl.plotToAxes(fig, ax)
                encoded = geoPl.getPlotImage(fig, ax)
                htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
                htmlstring += ret
                html = '<td width=' + width + '%>' + htmlstring + '</td>\n'
            elif geoPl.data.empty:
                html = '<td width=' + width + '%>' + 'No Data:' + geoPl.geoX + ' ' + geoPl.geoY + '</td>\n'
            else:
                fig, ax = plt.subplots()
                ret = geoPl.plotToAxes(fig, ax)
                encoded = geoPl.getPlotImage(fig, ax)
                htmlstring = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
                htmlstring += ret
                html = '<td width=' + width + '%>' + htmlstring + '</td>\n'
        # except:
        #    html = '<td width=' + width + '%>' + 'Error:' + geoPl.geoX + ' ' + geoPl.geoY + '</td>\n'

        return (html)

    def getHeaderString(self,fileName,title):
        html = '<!DOCTYPE html><html lang="en"><head><title>PSU-' + fileName + '-GEO</title>\n'
        # html += '<style> body {background-color:SeaShell;} table {table-layout:fixed;display:table;margin:0 auto;}td {border:1px solid RosyBrown;background-color:SeaShell;}</style>'
        # html += '<style> body {background-color:HoneyDew;} table {background-color:HoneyDew;} .innertable td {border:1px solid MistyRose;background-color:MintCream;}</style>'
        html += '<style> body {text-align:center;background-color:LightSteelBlue ;} img {width:95% }'
        html += 'table {font-size:0.8vw;width:95%;table-layout:fixed;display:table;margin:0 auto;background-color:LightSteelBlue ;}'
        html += ' td {border:1px solid MistyRose;background-color:AliceBlue;}</style>'
        html += '</head>\n'
        html += '<body><h1>' + title + '</h1>\n'
        html += '<h2>PSU: Geometric Correlations</h2>\n'
        html += '<hr/>\n'
        return html

    def flush(self):
        self.plots = []
        html = ''

