import netCDF4
from cdf.internal import EPOCHbreakdown
import pandas
import os
import numpy
import datetime
import matplotlib.pyplot as plt
import read_ssusi
import aacgmv2


if __name__ == "__main__":
    inpDirs = [ "../data/sdr/f18/20141216/" ]
    outDir = "/home/bharat/Documents/code/auroral-imager/data/processed/"
    ssRdObj = read_ssusi.ProcessData( inpDirs, outDir )
    ssRdObj.process_to_file()
    # ssRdObj.plot_ssusi_data( ssusiDF )

class ProcessData(object):
    """
    A class to Download SSUSI data
    given a date and datatype!
    """
    def __init__(self, inpDirs, outDir):
        """
        Given a list of dirs (SSUSI has multiple files per date
        for the same satellite). Read all files from it.
        """
        # loop through the dir and get a list of the files
        self.fileList = []
        for currDir in inpDirs:
            for root, dirs, files in os.walk(currDir):
                for fNum, fName in enumerate(files):
                    self.fileList.append( root + fName )
        self.outDir = outDir

    def process_to_file(self):
        """
        read the required data into a dataframe
        select only required columns, convert to
        AACGM coords and save data to file!
        """
        # selFname = "PS.APL_V0116S024CB0005_SC.U_DI.A_GP.F18-SSUSI_PA.APL-SDR-DISK_DD.20141216_SN.26612-00_DF.NC"
        for fileInd, currFile in enumerate(self.fileList):
            # if selFname not in currFile:
            #     continue
            # Get Sat name
            print "currently working with file-->", currFile
            print "processing--->", fileInd+1, "/", len(self.fileList), "files"
            satName = "F18"
            if "F17" in currFile:
                satName = "F17"
            if "F16" in currFile:
                satName = "F16"
            currDataSet = netCDF4.Dataset(currFile)
            # get datetime from epoch list
            dtList = numpy.array( [ datetime.datetime( *EPOCHbreakdown( e ) ) \
                        for e in currDataSet.variables["TIME_EPOCH_DAY"][:] ] )
            # get date for filename
            currDate = numpy.unique( numpy.array( [ \
                        x.strftime("%Y%m%d") for x in dtList ] ) )[0]
            # get peircepoints
            prpntLats = currDataSet.variables['PIERCEPOINT_DAY_LATITUDE'][:]
            prpntLons = currDataSet.variables['PIERCEPOINT_DAY_LONGITUDE'][:]
            prpntAlts = currDataSet.variables['PIERCEPOINT_DAY_ALTITUDE'][:]
            # GET DISK intensity data - waveband/color radiance data
            # 5 colors are - 121.6, 130.4, 135.6 nm and LBH short and LBH long
            dskInt121 = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 0]
            dskInt130 = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 1]
            dskInt135 = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 2]
            dskIntLBHS = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 3]
            dskIntLBHL = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 4]
            # We'll store the data in a DF. Now we need to be a little cautious
            # when storing the data in a DF. SSUSI measures flux as swaths, so
            # at each time instance we have multiple lats and lons and disk data.
            # I'm taking a simple approach where I take each lat (lon and other
            #  data) at a time instance as a column and time as rows. So if 
            # the array ishaving a dimention of 42x1632, each of the 42 
            # elements becomes a column and the 1632 time instances become rows.
            latColList = [ "glat." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            lonColList = [ "glon." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            d121ColList = [ "d121." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            d130ColList = [ "d130." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            d135ColList = [ "d135." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            dLBHSColList = [ "dlbhs." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            dLBHLColList = [ "dlbhl." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            # create dataframes with
            dfLat = pandas.DataFrame(prpntLats.T,columns=latColList, index=dtList)
            dfLon = pandas.DataFrame(prpntLons.T,columns=lonColList, index=dtList)
            dfD121 = pandas.DataFrame(dskInt121.T,columns=d121ColList, index=dtList)
            dfD130 = pandas.DataFrame(dskInt130.T,columns=d130ColList, index=dtList)
            dfD135 = pandas.DataFrame(dskInt135.T,columns=d135ColList, index=dtList)
            dfDLBHS = pandas.DataFrame(dskIntLBHS.T,columns=dLBHSColList, index=dtList)
            dfDLBHL = pandas.DataFrame(dskIntLBHL.T,columns=dLBHLColList, index=dtList)
            # Merge the dataframes
            ssusiDF = reduce(lambda left,right: pandas.merge(left,right,\
                         left_index=True, right_index=True), [ dfLat, \
                        dfLon, dfD121, dfD130, dfD135, dfDLBHL, dfDLBHS ])
            ssusiDF["orbitNum"] = currDataSet.variables['ORBIT_DAY'][:]
            # Lets also keep track of the sat name and shape of arrays
            ssusiDF["sat"] = satName
            ssusiDF["shapeArr"] = prpntLats.shape[0]
            # # reset index, we need datetime as a col
            ssusiDF = ssusiDF.reset_index()
            ssusiDF = ssusiDF.rename(columns = {'index':'date'})
            # Now we need to convert the GLAT, GLON into MLAT, MLON and MLT
            ssusiDF = ssusiDF.apply(self.convert_to_aacgm, axis=1)
            ssusiDF = ssusiDF.round(2)
            # We'll only need aacgm coords, discard all geog coords
            mlatColList = [ "mlat." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            mlonColList = [ "mlon." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            mltColList = [ "mlt." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            outCols = ["date", "sat", "orbitNum"] + mlatColList + mlonColList + mltColList + d121ColList + \
                        d130ColList + d135ColList + dLBHSColList + dLBHLColList
            ssusiDF = ssusiDF[ outCols ]
            # We now need to write the processed data to a file
            if not os.path.exists(self.outDir + "/" +satName):
                os.makedirs(self.outDir + "/" + satName)
            # if the file for the date exists append data
            # else create the file and write data!!!
            outFileName = self.outDir + "/" + satName + "/" + currDate + ".txt"
            if not os.path.exists( outFileName ):
                # NOTE we only need header when writing data for the first time!
                with open(outFileName, 'w') as ftB:
                    ssusiDF.to_csv(ftB, header=True,\
                                      index=False, sep=' ' )
            else:
                with open(outFileName, 'a') as ftB:
                    ssusiDF.to_csv(ftB, header=False,\
                                      index=False, sep=' ' )


    def convert_to_aacgm(self, row):
        """
        For the SSUSI DF convert all the 42
        Given glat, glon and date return
        mlat, mlon and mlt
        """
        for i in range( row["shapeArr"] ):
            indStr = str(i+1)
            mlat, mlon = aacgmv2.convert(row["glat." + indStr], row["glon." + indStr],\
                               300, row["date"])
            mlt = aacgmv2.convert_mlt(mlon, row["date"], m2a=False)
            row["mlat." + indStr] = numpy.round( mlat, 2)
            row["mlon." + indStr] = numpy.round( mlon, 2)
            row["mlt." + indStr] = numpy.round( mlt, 2)
        return row

    def plot_ssusi_data(self, ssusiDF, plotType='di135', coords="geo",\
                     figName="../figs/ssusi-test.pdf"):
        """
        Plot SSUSI data on a map
        """
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(1,1,1)
        m = plotUtils.mapObj(boundinglat=40., coords=coords)
        # Need to get max and min for the plot
        # we'll round off to the nearest 500
        # and keep a cap of 1000.
        vmax = numpy.round( numpy.max( ssusiDF[plotType] )/500. )*500.
        if vmax > 1000.:
            vmax = 1000.
        p = m.pcolormesh(ssusiDF["glon"].values.reshape(42,1632), \
                        ssusiDF["glat"].values.reshape(42,1632),\
                        ssusiDF[plotType].values.reshape(42,1632),\
                        latlon=True, zorder=1.9,
                        vmin=0, vmax=vmax,
                        ax=ax, alpha=1, cmap='Greens')
        p.set_rasterized(True)
        fig.savefig(figName,bbox_inches='tight')


