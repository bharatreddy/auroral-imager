import pandas
import os
import numpy
from scipy.interpolate import interp1d
import datetime
import matplotlib.pyplot as plt
import ssusi_utils
import aacgmv2
from davitpy import utils
import matplotlib.pyplot as plt


if __name__ == "__main__":
    inpDir = "../data/processed/"#/home/bharat/Documents/code/data/ssusi-prcsd/"
    fileDate = datetime.datetime( 2014, 12, 16 )
    inpTime = datetime.datetime( 2014, 12, 16, 18, 30 )
    ssObj = ssusi_utils.UtilsSsusi( inpDir, fileDate )

    fDict = ssObj.filter_data_by_time(inpTime, timeDelta=30)
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1,1,1)
    m = utils.plotUtils.mapObj(boundinglat=40., coords="mag")
    ssObj.overlay_sat_data( fDict, m, ax, satList=["F18"],\
             inpTime=inpTime, vmin=0., vmax=1000., autoScale=False )
    figName = "../figs/ssusi-sats.pdf" 
    fig.savefig(figName,bbox_inches='tight')

class UtilsSsusi(object):
    """
    A class to Download SSUSI data
    given a date and datatype!
    """
    def __init__(self, inpDir, fileDate, satList=[ "F18", "F17", "F16" ]):
        """
        Given a input dir read data from the satellites.
        """
        # Read data from the files and store them in a
        # dict of dataframes
        self.frames = {}
        dirList = []
        self.fileDate = fileDate
        currDtStr = self.fileDate.strftime("%Y%m%d")
        # Loop through the satList and read data 
        # from each satellite.
        # NOTE we expect the sub directory in the
        # parent directory to be the name of the
        # satellite listed in satList
        for sat in satList:
            currFname = inpDir + sat + "/" + currDtStr + ".txt"
            # check if file exists
            if os.path.exists( currFname ):
                print "reading data from--->", currFname
                self.frames[ "ssusi" + sat ] = pandas.read_csv(\
                                 currFname, delim_whitespace=True,\
                                infer_datetime_format=True,\
                                parse_dates=["date"] )
            else:
                print "file not found-->", currFname

    def filter_data_by_time(self, inpTime, hemi="north",\
             timeDelta=15, filterLat=40.):
        """
        Filter the processed data for 
        the desired time (+/- timeDel) of
        a given time and hemisphere.
        NOTE - timeDel is in minutes
        """
        # We'll output the results in a dict
        filteredDict = {}
        for key in self.frames.keys():
            ssusiDF = self.frames[key]
            ssusiDF = ssusiDF.fillna(0.)
            # get the time ranges that confine
            # inptime and timedelta
            timeMin = inpTime - datetime.timedelta(minutes=timeDelta)
            timeMax = inpTime + datetime.timedelta(minutes=timeDelta)
            # Choose DF rows which lie between timeMin and timeMax
            ssusiDF = ssusiDF[ (ssusiDF["date"] >= timeMin) &\
                                (ssusiDF["date"] <= timeMax) ]
            # select data based on hemi
            if hemi == "north":
                evalStr = "(ssusiDF['{0}'] >" + str( int(filterLat) ) + ".)" #
                # select all rows where lats are positive
                # we'll use the eval func for this purpose
                filterCol = [col for col in ssusiDF if col.startswith('mlat')]
                ssusiDF = ssusiDF[eval(" & ".join([\
                        evalStr.format(col) 
                        for col in filterCol]))].reset_index(drop=True)
            else:
                evalStr = "(ssusiDF['{0}'] <" + str( int(-1*filterLat) ) + ".)" #
                filterCol = [col for col in df if col.startswith('mlat')]
                ssusiDF = ssusiDF[eval(" & ".join([\
                        evalStr.format(col) 
                        for col in filterCol]))].reset_index(drop=True)
            if ssusiDF.shape[0] == 0:
                print "********NO DATA FOUND, CHECK FOR A " +\
                         "DIFFERENT TIME OR INCREASE TIMEDEL********"
            filteredDict[key] = ssusiDF
        return filteredDict


    def filter_data_by_orbit(self, inpTime, hemi="north", filterLat=0.):
        """
        Filter the processed data for 
        the desired closest orbit in time
        and hemisphere
        """
        # We'll output the results in a dict
        filteredDict = {}
        for key in self.frames.keys():
            ssusiDF = self.frames[key]
            ssusiDF = ssusiDF.fillna(0.)
            # get min and max times in each orbit
            orbitMin = ssusiDF[ ["date", "sat", "orbitNum"] \
                        ].groupby(["orbitNum", "sat"]).min().reset_index()
            orbitMin.columns = [ "orbitNum", "sat", "date_min" ]
            orbitMax = ssusiDF[ ["date", "sat", "orbitNum"] \
                        ].groupby(["orbitNum", "sat"]).max().reset_index()
            orbitMax.columns = [ "orbitNum", "sat", "date_max" ]
            orbitDF = pandas.merge( orbitMin, orbitMax,\
                         on=["orbitNum", "sat"] )
            selOrbit = orbitDF[ (orbitDF["date_min"] <= inpTime) &\
                 (orbitDF["date_max"] >= inpTime)\
                  ].reset_index(drop=True)
            # Only select the required orbit
            ssusiDF = ssusiDF.merge( selOrbit, on=[ "orbitNum", "sat" ] )
            # select data based on hemi
            if hemi == "north":
                evalStr = "(ssusiDF['{0}'] >" + str( int(filterLat) ) + ".)" #
                # select all rows where lats are positive
                # we'll use the eval func for this purpose
                filterCol = [col for col in ssusiDF if col.startswith('mlat')]
                ssusiDF = ssusiDF[eval(" & ".join([\
                        evalStr.format(col) 
                        for col in filterCol]))].reset_index(drop=True)
            else:
                evalStr = "(ssusiDF['{0}'] <" + str( int(-1*filterLat) ) + ".)" #
                filterCol = [col for col in df if col.startswith('mlat')]
                ssusiDF = ssusiDF[eval(" & ".join([\
                        evalStr.format(col) 
                        for col in filterCol]))].reset_index(drop=True)
            filteredDict[key] = ssusiDF
        return filteredDict

    def cart2pol(self, x, y):
        """
        convert from cartesian to polar coords
        """
        colat = numpy.sqrt(x**2 + y**2)
        lat = 90. - colat
        lon = numpy.rad2deg( numpy.arctan2(y, x) )
        return (lat, lon)

    def pol2cart(self, lat, lon):
        """
        convert from polar to cartesian coords
        """
        colat = 90. - lat
        x = colat * numpy.cos(numpy.deg2rad(lon))
        y = colat * numpy.sin(numpy.deg2rad(lon))
        return (x, y)
        
    def overlay_sat_data(self, filteredDict, mapHandle, ax,\
                        satList=["F18", "F17", "F16"], plotType='d135',\
                        overlayTime=True, overlayTimeInterval=5, timeMarker='o',\
                        timeMarkerSize=2., timeColor="grey", timeFontSize=8.,\
                         plotCBar=True, autoScale=True, vmin=0., vmax=1000.,\
                         plotTitle=True, titleString=None, inpTime=None, markSatName=True):
        """
        Plot SSUSI data on a map
        # overlayTimeInterval is in minutes
        """
        # Loop through and read data
        for key in filteredDict.keys():
            ssusiDF = filteredDict[key]
            satNameKey= key[-3:]
            if satNameKey not in satList:
                continue
            ssusiMlats = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlat')\
                            ]].values
            ssusiMlons = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlon')\
                            ]].values
            ssusiDisk = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith(plotType)\
                            ]].values
            # Need to get max and min for the plot
            # we'll round off to the nearest 500
            # and keep a cap of 1000.
            if autoScale:
                vmin = 0.
                vmax = numpy.round( numpy.max( ssusiDisk )/500. )*500.
            xVecs, yVecs = mapHandle(ssusiMlons, ssusiMlats, coords="mag")
            ssusiPlot = mapHandle.scatter(xVecs, yVecs, c=ssusiDisk, s=10.,\
                       cmap="Greens", alpha=0.7, zorder=5., \
                                 edgecolor='none', marker="s",\
                                  vmin=vmin, vmax=vmax)
            # p = mapHandle.pcolormesh(ssusiMlats, ssusiMlons,\
            #                 ssusiDisk,\
            #                 latlon=True, zorder=1.9,
            #                 vmin=0, vmax=vmax,
            #                 ax=ax, alpha=1, cmap='Greens')
            ssusiPlot.set_rasterized(True)
            # overlay time
            if overlayTime:
                uniqueTimeList = ssusiDF["date"].unique()
                timeDiff = ( uniqueTimeList.max() -\
                             uniqueTimeList.min()\
                              ).astype('timedelta64[m]')
                delRange = overlayTimeInterval*\
                            uniqueTimeList.shape[0]/timeDiff.astype('int')
                # get a list of times every timeinterval
                # for the given day
                nextDayTime = self.fileDate + datetime.timedelta(days=1)
                currDt = self.fileDate
                allDayDatesList = []
                allDayTSList = []
                minDate = datetime.datetime.utcfromtimestamp(\
                                uniqueTimeList.min().tolist()/1e9)
                maxDate = datetime.datetime.utcfromtimestamp(\
                                uniqueTimeList.max().tolist()/1e9)
                while currDt <= nextDayTime:
                    if ( currDt >= minDate ) & ( currDt <= maxDate ):
                        ts = (numpy.datetime64(currDt) - \
                                numpy.datetime64('1970-01-01T00:00:00Z')\
                                ) / numpy.timedelta64(1, 's')
                        allDayTSList.append( ts )
                        allDayDatesList.append( currDt )
                    currDt += datetime.timedelta(minutes=overlayTimeInterval)
                
                timeSSusiMlats = ssusiDF\
                        [ssusiDF.columns[pandas.Series(\
                        ssusiDF.columns).str.startswith('mlat')\
                        ]].values
                timeSSusiMlons = ssusiDF\
                        [ssusiDF.columns[pandas.Series(\
                        ssusiDF.columns).str.startswith('mlon')\
                        ]].values
                timeSSusiTimes = ssusiDF["date"].values
                satTSArr = (timeSSusiTimes - \
                            numpy.datetime64('1970-01-01T00:00:00Z')\
                            ) / numpy.timedelta64(1, 's')
                # Interpolate the values to get times
                for dd in range( len(allDayDatesList)-1 ):
                    for pixel in range(timeSSusiMlons.shape[1]):
                        currPixelMlons = timeSSusiMlons[:,pixel]
                        currPixelMlats = timeSSusiMlats[:,pixel]
                        (x,y) = self.pol2cart( currPixelMlats, currPixelMlons )
                        xArr = numpy.interp(allDayTSList[dd], satTSArr, x)
                        yArr = numpy.interp(allDayTSList[dd], satTSArr, y)
                        (timePlotLatArr, timePlotLonArr) = self.cart2pol( xArr, yArr )
                        xTVecs, yTVecs = mapHandle(timePlotLonArr,\
                                         timePlotLatArr, coords="mag")
                        mapHandle.plot(xTVecs, yTVecs,\
                             marker=timeMarker,color=timeColor,\
                              markersize=timeMarkerSize, zorder=7.)
                        timeStr = allDayDatesList[dd].strftime("%H:%M")
                        # Write Sat names used in plotting
                        if pixel == 0:
                            if markSatName:
                                timeStr = timeStr + " (" + satNameKey + ")"
                            timeXVecs, timeYVecs = mapHandle(timePlotLonArr,\
                                 timePlotLatArr, coords="mag")
                            ax.text(timeXVecs, timeYVecs, timeStr,\
                                fontsize=timeFontSize,fontweight='bold',
                                ha='left',va='center',color='k',\
                                 clip_on=True, zorder=7.)
            # plot colorbar
            if plotCBar:
                cbar = plt.colorbar(ssusiPlot, orientation='vertical')
                cbar.set_label('Rayleighs', size=14)
            # Title
            if plotTitle:
                if titleString is not None:
                    plt.title(titleString)
                else:
                    if inpTime is not None:
                        inpTimeStr = inpTime.strftime("%Y-%m-%d  %H:%M")
                        plt.title( inpTimeStr + " UT" )
                    else:
                        print "***********NEED INPTIME FOR TITLE***********"
            