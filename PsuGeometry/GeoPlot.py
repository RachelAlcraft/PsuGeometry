
import pandas as pd
import numpy as np
from matplotlib import cm
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64
import PIL

from PsuGeometry import GeoReport as geor
from PsuGeometry import GeoPdb as geop


class GeoPlot:
    def __init__(self,data,geoX,geoY='',title='',hue='bfactor',splitKey='',palette='viridis_r',
                 centre=False,vmin=0,vmax=0,operation='',newData=False,plot='scatter',categorical=False):
        self.plot = plot
        self.data = data
        self.geoX = geoX
        self.geoY = geoY
        self.title=title
        self.hue = hue
        self.splitKey=splitKey
        self.palette=palette
        self.centre = centre
        self.vmin=vmin
        self.vmax=vmax
        self.operation=operation
        self.newData = newData
        self.categorical = categorical
        self.numpy = []
        self.hasMatrix = False
        self.restrictions = {}
        self.axX = 0,0
        self.axY = 0, 0
        if self.geoY == '':
            self.plot = 'histogram'
            #if self.hue=='bfactor':
            #    self.hue = 'pdbCode'

    #def getPlot(self,fig, ax):
    #    if self.plot == 'histogram':
    #        return self.plotHistogram(True,fig, ax)
    #    elif self.plot == 'scatter':
    #        return self.plotScatter(True,fig, ax)
    #    elif self.plot == 'probability':
    #        return self.plotProbability(True,fig, ax)

    def plotToAxes(self,fig, ax):
        if self.plot == 'histogram':
            return self.plotHistogram(fig, ax)
        elif self.plot == 'scatter':
            return self.plotScatter(fig, ax)
        elif self.plot == 'probability':
            return self.plotProbability(fig, ax)

    def plotHistogram(self,fig, ax):
        data = self.data.sort_values(by=self.geoX, ascending=True)
        title = self.title

        if self.operation == 'ABS':
            data = data[data[self.geoX] == abs(data[self.geoX])]
        elif self.operation == 'SQUARE':
            data = data[data[self.geoX] == data[self.geoX] ** 2]

        firstVal = data.head(1)[self.geoX].values[0]
        lastVal = data.tail(1)[self.geoX].values[0]
        firstHue = data.head(1)[self.hue].values[0]
        lastHue = data.tail(1)[self.hue].values[0]

        try:
            firstVal = round(firstVal, 2)
            lastVal = round(lastVal, 2)
        except:
            pass
        try:
            firstHue = round(firstHue, 2)
            lastHue = round(lastHue, 2)
        except:
            pass
        title += '\nFirst:' + self.hue + ' ' + str(firstHue) + '=' + str(firstVal)
        title += '\nLast:' + self.hue + ' ' + str(lastHue) + '=' + str(lastVal)

        # sns.distplot(data[xName], norm_hist=True, bins=50, kde=False)
        histCol = 'crimson'
        alpha=1
        bins = min(max(int(len(data[self.geoX])/6),10),50)
        if self.title == 'ghost':
            histCol = 'gainsboro'
            alpha=0.5
            plt.hist(data[self.geoX], EdgeColor='k', bins=bins,color=histCol,alpha=alpha,density=True,label='ghost')
            #sns.distplot(data[self.geoX], label='x', norm_hist=True, bins=50, kde=False,color='gainsboro')
        else:
            #if self.hue != '':
            #    splitList = data[self.hue].unique()
            #    for split in splitList:
            #        dfx = data[data[self.hue] == split]
            #        bins = max(int(len(dfx[self.geoX]) / 6),10)
            #        #plt.hist(dfx[self.geoX], EdgeColor='k', bins=bins, alpha=alpha, density=True)
            #        sns.distplot(dfx[self.geoX], label=split, norm_hist=True, bins=bins, kde=False,hist_kws=dict(alpha=0.5,EdgeColor='silver'))
            #    plt.legend()
            #else:
            #sns.distplot(data[self.geoX], label='', norm_hist=True, bins=bins, kde=False,hist_kws=dict(alpha=0.8,EdgeColor='silver'))
            plt.hist(data[self.geoX], EdgeColor='k', bins=bins, color=histCol, alpha=alpha, density=True)


        plt.title(title)
        plt.xlabel(self.geoX)
        #if returnData:
        #    img = io.BytesIO()
        #    fig.savefig(img, format='png', bbox_inches='tight')
        #    img.seek(0)
        #    encoded = base64.b64encode(img.getvalue())
        #    # html = '<img width=100% src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        #    html = '<p><img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        #    plt.close('all')
        if self.title != 'ghost':
            dfdesc = self.data[self.geoX].describe()
            rows = len(dfdesc.index)
            colsNames = list(dfdesc.index)
            html = "<table class='innertable'>\n"
            html += "<tr>\n"
            for r in range(0, rows):
                html += "<td>" + str(colsNames[r]) + "</td>\n"
            html += "</tr>\n"
            html += "<tr>"
            for r in range(0, rows):
                html += "<td>"
                try:
                    html += str(round(dfdesc[r], 2))
                except:
                    html += str(dfdesc[r])
                html += "</td>\n"

            html += "</tr>\n"
            html += "</table></p>\n"

            return html
        else:
            return ''

    def plotScatter(self,fig, ax):

        #fig, ax = plt.subplots()
        if self.categorical or self.hue == 'dssp':
            #blanksdata = self.data[self.data[self.hue] == '']
            #print(blanksdata)
            # it is possible for errors in dssp assignment in which case we call them X, but it will cover any errors not just dssp
            self.data.loc[self.data[self.hue] == '', self.hue] = 'X'

            gradients = {}
            gradsorig = self.data.sort_values(by=self.hue, ascending=True)[self.hue].unique()
            grads = self.getHueLists(self.hue,gradsorig)
            evenly_spaced_interval = np.linspace(0, 1, len(grads))
            if self.palette is str:
                try:
                    sns.set_palette(sns.color_palette(self.palette, len(grads)))
                    colors = [cm.get_cmap(self.palette)(x) for x in evenly_spaced_interval]
                except:
                    colors = [cm.get_cmap(self.palette)(x) for x in evenly_spaced_interval]
                i = 0
                for g in grads:
                    gradients[g] = colors[i]
                    i = i+1
                self.palette = gradients

        if self.centre:
            self.data[self.hue + '2'] = self.data[self.hue] ** 2
            data = self.data.sort_values(by=self.hue + '2', ascending=True)
            maxh = max(data[self.hue].max(), -1 * data[self.hue].min())
            minh = maxh * -1
            g = ax.scatter(data[self.geoX], data[self.geoY], c=data[self.hue], cmap=self.palette, vmin=minh,vmax=maxh)
            fig.colorbar(g)
            ax.set_xlabel(self.geoX)
            ax.set_ylabel(self.geoY)
        elif self.vmin < self.vmax:
            data = self.data.sort_values(by=self.hue, ascending=True)
            g = ax.scatter(data[self.geoX], data[self.geoY], c=data[self.hue], cmap=self.palette, vmin=self.vmin,
                           vmax=self.vmax)
            fig.colorbar(g)
            ax.set_xlabel(self.geoX)
            ax.set_ylabel(self.geoY)
        else:

            lw = 0.5
            if self.palette == 'gist_gray_r':
                lw = 0  # this gives a crystollagraphic image look

            if self.hue == 'aa':
                try:
                    self.data = self.data.sort_values(by='2FoFc', ascending=True)
                except:
                    self.data = self.data.sort_values(by=self.hue, ascending=True)
            else:
                self.data = self.data.sort_values(by=self.hue, ascending=True)

            alpha=0.8
            if self.title=='ghost':
                alpha = 0.4

            im = sns.scatterplot(x=self.geoX, y=self.geoY, hue=self.hue, data=self.data, alpha=alpha,
                                 palette=self.palette, edgecolor='aliceblue', linewidth=lw)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  # Put the legend out of the figure

        count = len(self.data.index)
        title = self.title
        if title == '':
            title += 'Count=' + str(count)
        else:
            title += '\nCount=' + str(count)

        plt.title(title)
        return ''
        #if returnData:
        #    img = io.BytesIO()
        #    fig.savefig(img, format='png', bbox_inches='tight')
        #    img.seek(0)
        #    encoded = base64.b64encode(img.getvalue())
        #    # html = '<img width=100% src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        #    html = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        #    plt.close('all')
        #    return html

    def getPlotImage(self,fig, ax):
        #fig, ax = plt.subplots()
        img = io.BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        encoded = base64.b64encode(img.getvalue())
        plt.close('all')
        return encoded

    def getAxes(self):
        xMin = min(self.data[self.geoX])
        xMax = max(self.data[self.geoX])
        yMin = min(self.data[self.geoY])
        yMax = max(self.data[self.geoY])
        return ([xMin,xMax,yMin,yMax])

    def plotProbability(self,fig, ax):

        # These shold be settings
        contours = 12
        kde = 0.10
        bins=50
        minX,maxX = self.axX[0],self.axX[1]
        minY,maxY = self.axY[0],self.axY[1]

        if minX == maxX:
            minX = min(self.data[self.geoX])
            maxX = max(self.data[self.geoX])
            minY = min(self.data[self.geoY])
            maxY = max(self.data[self.geoY])


        #fig, ax = plt.subplots()
        plt.axis([minX, maxX, minY, maxY])
        if self.hasMatrix:
            xgrid, ygrid, zgrid = self.numpy
            print(zgrid)
        else:
            xgrid, ygrid, zgrid = self.kde2D_scipy(kde,[minX,maxX,minY,maxY], bins)

        xgrid = np.linspace(minX, maxX, bins)
        ygrid = np.linspace(minY,maxY, bins)

        ax.grid(True, which='major', axis='both', linestyle='-', color=(0.8, 0.8, 0.8), alpha=0.3)

        if self.centre:
            self.data[self.hue + '2'] = self.data[self.hue] ** 2
            self.data = self.data.sort_values(by=self.hue + '2', ascending=True)
            self.vmax = max(self.data[self.hue].max(), -1 * self.data[self.hue].min())
            self.vmin = self.vmax * -1

        alpha=1
        if self.title=='ghost':
            alpha = 0.4
            #self.palette = 'seismic'

        if self.vmin == self.vmax:
            im = plt.pcolormesh(xgrid, ygrid, zgrid, shading='gouraud', cmap=self.palette,alpha=alpha)
            cs = plt.contour(xgrid, ygrid, zgrid, contours, colors='0.7', linewidths=0.4,alpha=alpha)
        else:
            im = plt.pcolormesh(xgrid, ygrid, zgrid, shading='gouraud', cmap=self.palette, vmin=self.vmin, vmax=self.vmax,alpha=alpha)
            cs = plt.contour(xgrid, ygrid, zgrid, contours, colors='tab:purple', linewidths=0.05,alpha=alpha)

        cbar = fig.colorbar(im, ax=ax)
        cbar.remove()


        ax.set_xlabel(self.geoX)
        ax.set_ylabel(self.geoY)

        if self.hasMatrix:
            if self.title == '':
                title = 'Difference Image'
            else:
                title = self.title + '\nDifference Image'
        else:
            count = len(self.data.index)
            title = self.title
            if title == '':
                title += 'Count=' + str(count)
            else:
                title += '\nCount=' + str(count)

        plt.title(title)
        return ''

        #if returnData:
        #    img = io.BytesIO()
        #    fig.savefig(img, format='png', bbox_inches='tight')
        #    img.seek(0)
        #    encoded = base64.b64encode(img.getvalue())

        #    html = '<img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        #    plt.close('all')
        #    return html


    def kde2D_scipy(self,bandwidth, axes, bins):

        xdata = self.data[self.geoX]
        ydata = self.data[self.geoY]
        data = np.vstack([xdata, ydata])
        xgrid = np.linspace(axes[0], axes[1], bins)
        ygrid = np.linspace(axes[2], axes[3], bins)
        Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
        grid_sized = np.vstack([Xgrid.ravel(), Ygrid.ravel()])
        # fit an array of size [Ndim, Nsamples]

        kde = gaussian_kde(data, bw_method=bandwidth)
        # evaluate on a regular grid
        Z = kde.evaluate(grid_sized)
        zgrid = Z.reshape(Xgrid.shape)
        return xgrid, ygrid, zgrid

    def getDifference(self,comparePlot):
        return ''

    def getOverlay(self,compareList):
        return ''

    def getHueLists(self,hue,huelist):
        if hue == 'dssp':
            return ['-','B','E','G','H','I','S','T','X']
        else:
            return huelist

    def getNewData(self,pdbs,hues=None):
        if self.plot == 'histogram':
            calcList = [self.geoX]
        else:
            calcList = [self.geoX, self.geoY]
        hueList = hues
        if hues == None:
            hueList = [self.hue]
        dfs = []
        print('calcList', calcList)
        for apdb in pdbs:
            data = apdb.getGeoemtryCsv(calcList, hueList)
            dfs.append(data)
        self.data = pd.concat(dfs, ignore_index=True)
        # now the data can be restricted as per the restrictions, which is a dictionary of restrictions, eg aa:'THR,PRO'
        if len(self.restrictions)>0:
            dfs = []
            for hue in self.restrictions:
                allowed = self.restrictions[hue]
                data = self.data[self.data[hue] == allowed]
                dfs.append(data)
                if self.title != '':
                    self.title += '/n'
                self.title += hue + ':' + allowed
            self.data = pd.concat(dfs, ignore_index=True)

    def getMatrix(self):

        kde = 0.10
        bins = 50
        minX,maxX,minY,maxY = self.axX[0],self.axX[1],self.axY[0],self.axY[1]
        xgrid, ygrid, zgrid = self.kde2D_scipy(kde, [minX, maxX, minY, maxY], bins)
        return xgrid, ygrid, zgrid



class GeoOverlay:
    def __init__(self,plotA, plotB, title,pdbDataPath='',edDataPath='',outDataPath=''):
        self.title = title
        if title!='ghost':
            self.plotA = plotA
            self.plotB = plotB
        else:#In this case we have only the main plot, so we create the dummy plot
            self.plotB = plotA
            ghostReport = geor.GeoReport(['ghost'],pdbDataPath,edDataPath,outDataPath)
            geoList = []
            geoList.append(self.plotB.geoX)
            if self.plotB.geoY != '':
                geoList.append(self.plotB.geoY)
            ghostdata = ghostReport.getGeoemtryCsv(geoList, ['pdbCode'])
            self.plotA = GeoPlot(ghostdata, self.plotB.geoX, geoY=self.plotB.geoY, title='ghost', hue='pdbCode', palette='Greys',plot=self.plotB.plot,operation=self.plotB.operation)

class GeoDifference:
    def __init__(self,pdbs,dataA,dataB,geoX,geoY='',title='',restrictionsA={},restrictionsB={},newData=False,palette='seismic'):

        #huelist is all of the restrictions
        print('rest2', restrictionsA)
        hues = []
        for hue in restrictionsA:
            if hue not in hues:
                hues.append(hue)


        self.plotA = GeoPlot(dataA,geoX,geoY=geoY,newData=newData,palette=palette,plot='probability')
        self.plotB = GeoPlot(dataB, geoX, geoY=geoY, newData=newData, palette=palette + '_r',plot='probability')

        self.plotA.restrictions = restrictionsA
        self.plotB.restrictions = restrictionsB

        self.plotA.getNewData(pdbs,hues)
        self.plotB.getNewData(pdbs,hues)

        self.plotA.newData = False
        self.plotB.newData = False

        axesA = self.plotA.getAxes()
        axesB = self.plotB.getAxes()
        axesAll = min(axesA[0],axesB[0]),max(axesA[1],axesB[1]),min(axesA[2],axesB[2]),max(axesA[3],axesB[3])
        self.plotA.axX = axesAll[0], axesAll[1]
        self.plotA.axY = axesAll[2], axesAll[3]
        self.plotB.axX = axesAll[0], axesAll[1]
        self.plotB.axY = axesAll[2], axesAll[3]

        arAA = self.plotA.getMatrix()
        arBB = self.plotB.getMatrix()

        print('same',arAA[2] == arBB[2])

        arA = arAA[2]
        arB = arBB[2]
        arDiff = arA - arB
        minVal,maxVal=0,0
        for i in range(arA.shape[0]):
            for j in range(arA.shape[1]):
                maxVal = max(maxVal, arA[i, j])
                minVal = min(minVal, arA[i, j])

        for i in range(arB.shape[0]):
            for j in range(arB.shape[1]):
                maxVal = max(maxVal, arB[i, j])
                minVal = min(minVal, arB[i, j])

        maxVal = max(abs(maxVal), abs(minVal))
        minVal = 0 - maxVal

        self.plotA.vmin = minVal
        self.plotB.vmin = minVal
        self.plotA.vmax = maxVal
        self.plotB.vmax = maxVal


        self.plotDiff = GeoPlot(dataA, geoX, geoY=geoY, newData=False, palette=palette,vmin=minVal,vmax=maxVal,plot='probability')
        self.plotDiff.hasMatrix = True
        self.plotDiff.numpy = arB[0],arB[1],arDiff
        self.plotDiff.vmax = maxVal
        self.plotDiff.vmax = maxVal
        self.plotDiff.axX = axesAll[0], axesAll[1]
        self.plotDiff.axY = axesAll[2], axesAll[3]
