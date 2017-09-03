import netCDF4
from cdf.internal import EPOCHbreakdown
import pandas
import os
import numpy
import datetime
import matplotlib.pyplot as plt
import read_ssusi

if __name__ == "__main__":
    inpDirs = [ "../data/sdr/f18/20141216/" ]
    ssRdObj = read_ssusi.ReadData( inpDirs )
    ssRdObj.read_data()

class ReadData(object):
    """
    A class to Download SSUSI data
    given a date and datatype!
    """
    def __init__(self, inpDirs):
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

    def read_data(self):
        """
        read the required data into a dataframe
        """
        for currFile in self.fileList:
            currDataSet = netCDF4.Dataset(currFile)
            # get datetime from epoch list
            dtList = numpy.array( [ datetime.datetime( *EPOCHbreakdown( e ) ) \
                        for e in currDataSet.variables["TIME_EPOCH_DAY"][:] ] )
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
            # We'll store the data in a DF. We need to convert the dates into a 
            # similar shape as disk radiance.
            dateArr = numpy.asarray( [dtList]*dskInt121.shape[0] )
            # store data in a DF
            ssusiDF = pandas.DataFrame({ 'glat': prpntLats.ravel(),\
                         'glon': prpntLons.ravel(),'di121': dskInt121.ravel(),\
                          'di130': dskInt130.ravel(),'di135': dskInt135.ravel(),\
                          'diLBHS': dskIntLBHS.ravel(),'diLBHL': dskIntLBHL.ravel(),\
                          'date':dateArr.ravel()})
            break
        return ssusiDF


