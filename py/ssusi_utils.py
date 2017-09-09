import pandas
import os
import numpy
import datetime
import matplotlib.pyplot as plt
import ssusi_utils
import aacgmv2
from davitpy import utils
import matplotlib.pyplot as plt


if __name__ == "__main__":
    inpDir = "../data/processed/"
    fileDate = datetime.datetime( 2014, 12, 16 )
    inpTime = datetime.datetime( 2014, 12, 16, 20, 30 )
    ssObj = ssusi_utils.UtilsSsusi( inpDir, fileDate )
    fDict = ssObj.filter_data(inpTime)
    # fig = plt.figure(figsize=(12, 8))
    # ax = fig.add_subplot(1,1,1)
    # m = utils.plotUtils.mapObj(boundinglat=40., coords="mag")
    # ssObj.overlay_sat_data( fDict, m, ax, satList=["F16"] )
    # figName = "../figs/ssusi-sats.pdf" 
    # fig.savefig(figName,bbox_inches='tight')
    ssObj.combine_sat_data(fDict)

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

    def filter_data(self, inpTime, hemi="north", filterLat=0.):
        """
        Filter the processed data for 
        the desired time and hemisphere
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

    def combine_sat_data(self, filteredDict, diskType='d135'):
        """
        Combine data from all the satellites
        and get an integrated image.
        """
        # We'll round-off each latitude and longitude
        # and try to identify the mean values from
        # mlat, mlon pairs that fall into the category!
        # combine data from all the frames
        frames = []
        for key in filteredDict.keys():
            currDF = filteredDict[key]
            frames.append( currDF )
        ssusiDF = pandas.concat( frames )
        # round all values in DF so that 
        # we can aggregate MLATs and MLONs
        ssusiDF = ssusiDF.round()
        mlonList = []
        mlatList = []
        mltList = []
        d121List = []
        d130List = []
        d135List = []
        dLBHSList = []
        dLBHLList = []
        for currMlon in range(-180,181):
            for currMlat in range(-90, 91):
                a = ssusiDF[ ssusiDF ]

        

    def overlay_sat_data(self, filteredDict, mapHandle, ax,\
                        satList=["F18", "F17", "F16"], plotType='d135',\
                        overlayTime=True, overlayTimeInterval=10, timeMarker='D',\
                        timeMarkerSize=2., timeColor="grey", timeFontSize=8.):
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
            vmax = numpy.round( numpy.max( ssusiDisk )/500. )*500.
            if vmax > 1000.:
                vmax = 1000.
            xVecs, yVecs = mapHandle(ssusiMlons, ssusiMlats, coords="mag")
            p = mapHandle.scatter(xVecs, yVecs, c=ssusiDisk, s=10.,\
                       cmap="Greens", alpha=0.7, zorder=5., \
                                 edgecolor='none', marker="s",\
                                  vmin=0., vmax=1000.)
            # p = mapHandle.pcolormesh(ssusiMlats, ssusiMlons,\
            #                 ssusiDisk,\
            #                 latlon=True, zorder=1.9,
            #                 vmin=0, vmax=vmax,
            #                 ax=ax, alpha=1, cmap='Greens')
            p.set_rasterized(True)
            # overlay time
            if overlayTime:
                uniqueTimeList = ssusiDF["date"].unique()
                timeDiff = ( uniqueTimeList.max() -\
                             uniqueTimeList.min()\
                              ).astype('timedelta64[m]')
                delRange = overlayTimeInterval*uniqueTimeList.\
                                shape[0]/timeDiff.astype('int')
                # loop through and overlay times
                for tt in range(0,uniqueTimeList.shape[0],delRange):
                    timeSSusiMlats = ssusiDF[ ssusiDF["date"] ==\
                             uniqueTimeList[tt] ]\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlat')\
                            ]].values
                    timeSSusiMlons = ssusiDF[ ssusiDF["date"] ==\
                             uniqueTimeList[tt] ]\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlon')\
                            ]].values
                    xTVecs, yTVecs = mapHandle(timeSSusiMlons,\
                                     timeSSusiMlats, coords="mag")
                    mapHandle.plot(xTVecs, yTVecs,\
                         marker=timeMarker,color=timeColor,\
                          markersize=timeMarkerSize, zorder=7.)
                    timeStr = pandas.to_datetime(uniqueTimeList[tt]).strftime("%H:%M")
                    timeXVecs, timeYVecs = mapHandle(timeSSusiMlons[-1][-1],\
                         timeSSusiMlats[-1][-1], coords="mag")
                    ax.text(timeXVecs, timeYVecs, timeStr,\
                        fontsize=timeFontSize,fontweight='bold',
                        ha='left',va='center',color='k',\
                         clip_on=True, zorder=7.)