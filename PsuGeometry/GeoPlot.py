
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
    def __init__(self,data,geoX,geoY='',title='',hue='2FoFc',splitKey='',palette='viridis_r',
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
        if self.geoY == '':
            self.plot = 'histogram'

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
        plt.hist(data[self.geoX], EdgeColor='k', bins=50,color='tomato')

        plt.title(title)
        #if returnData:
        #    img = io.BytesIO()
        #    fig.savefig(img, format='png', bbox_inches='tight')
        #    img.seek(0)
        #    encoded = base64.b64encode(img.getvalue())
        #    # html = '<img width=100% src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        #    html = '<p><img src="data:image/png;base64, {}">'.format(encoded.decode('utf-8')) + '\n'
        #    plt.close('all')
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

    def plotScatter(self,fig, ax):

        #fig, ax = plt.subplots()
        if self.categorical or self.hue == 'dssp':
            gradients = {}
            grads = self.data.sort_values(by=self.hue, ascending=True)[self.hue].unique()
            grads = self.getHueLists(self.hue,grads)
            evenly_spaced_interval = np.linspace(0, 1, len(grads))
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
            if self.title=='Dummy':
                alpha = 0.3

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
        minX = 0
        maxX = 0
        minY = 0
        maxY = 0

        if minX == maxX:
            minX = min(self.data[self.geoX])
            maxX = max(self.data[self.geoX])
            minY = min(self.data[self.geoY])
            maxY = max(self.data[self.geoY])


        #fig, ax = plt.subplots()
        plt.axis([minX, maxX, minY, maxY])
        xgrid, ygrid, zgrid = self.kde2D_scipy(kde,[minX,maxX,minY,maxY], bins)
        ax.grid(True, which='major', axis='both', linestyle='-', color=(0.8, 0.8, 0.8), alpha=0.3)

        if self.centre:
            self.data[self.hue + '2'] = self.data[self.hue] ** 2
            self.data = self.data.sort_values(by=self.hue + '2', ascending=True)
            self.vmax = max(self.data[self.hue].max(), -1 * self.data[self.hue].min())
            self.vmin = self.vmax * -1

        if self.vmin == self.vmax:
            im = plt.pcolormesh(xgrid, ygrid, zgrid, shading='gouraud', cmap=self.palette)
            cs = plt.contour(xgrid, ygrid, zgrid, contours, colors='0.7', linewidths=0.4)
        else:
            im = plt.pcolormesh(xgrid, ygrid, zgrid, shading='gouraud', cmap=self.palette, vmin=self.vmin, vmax=self.vmax)
            cs = plt.contour(xgrid, ygrid, zgrid, contours, colors='tab:purple', linewidths=0.05)

        cbar = fig.colorbar(im, ax=ax)
        cbar.remove()


        ax.set_xlabel(self.geoX)
        ax.set_ylabel(self.geoY)

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
            return ['-','B','E','G','H','I','S','T']
        else:
            return huelist

    def getNewData(self,pdbs):
        if self.plot == 'histogram':
            calcList = [self.geoX]
        else:
            calcList = [self.geoX, self.geoY]
        hueList = [self.hue]
        dfs = []
        for apdb in pdbs:
            data = apdb.getGeoemtryCsv(calcList, hueList)
            dfs.append(data)
        self.data = pd.concat(dfs, ignore_index=True)



class GeoOverlay:
    def __init__(self,plotA, plotB, title,pdbDataPath='',edDataPath=''):
        self.title = title
        if title!='ghost':
            self.plotA = plotA
            self.plotB = plotB
        else:#In this case we have only the main plot, so we create the dummy plot
            self.plotB = plotA
            geoDummy = geop.GeoPdb('ghost', pdbDataPath, edDataPath)
            dummyReport = geor.GeoReport([geoDummy])
            geoList = []
            geoList.append(self.plotB.geoX)
            if self.plotB.geoY != '':
                geoList.append(self.plotB.geoY)
            dummydata = dummyReport.getGeoemtryCsv(geoList, ['pdbCode'])
            self.plotA = GeoPlot(dummydata, self.plotB.geoX, geoY=self.plotB.geoY, title='ghost', hue='pdbCode', palette='Greys')




