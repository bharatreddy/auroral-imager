import netCDF4
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
            print currFile
            print currDataSet.variables["TIME_DAY"].shape

            epochList = currDataSet.variables["TIME_EPOCH_DAY"][:]
            # account for difference in seconds between
            # CDF epoch and python's epoch, leap year in there
            # (datetime(1971,1,2) - 
            #      datetime(1,1,1)).total_seconds()*1000
            epochList = epochList - 62167219200000
            a = epochList[0]
            print datetime.datetime.fromtimestamp(a/1e3)
            break

