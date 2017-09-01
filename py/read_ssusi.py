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
            break

