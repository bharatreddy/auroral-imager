import os
import datetime
import dwnld_ssusi
import read_ssusi
import shutil

if __name__ == "__main__":
    tempFileDir = "/home/bharat/Documents/code/data/ssusi-temp/"
    prcsdFileDir = "/home/bharat/Documents/code/data/ssusi-prcsd/"
    satList = [ "f16", "f17", "f18" ]
    ssDwnldObj = dwnld_ssusi.SSUSIDownload(\
                    outBaseDir = tempFileDir)
    dataTypeList = [ "sdr" ]#, "l1b", "edr-aur" ]
    currDate = datetime.datetime( 2011, 4, 9 )
    endDate = datetime.datetime( 2011, 4, 9 )
    tDelta = datetime.timedelta(days=1)
    while currDate <= endDate:
        print "currently downloading files for --> ",\
            currDate.strftime("%Y-%m-%d")
        ssDwnldObj.download_files(currDate, dataTypeList)
        for currSat in satList:
            currDir = tempFileDir + "sdr/" + currSat + "/"
            for root, dirs, files in os.walk(currDir):
                for nd, dd in enumerate(dirs):
                    print "processing data --> ",\
                             currDate.strftime("%Y-%m-%d"), " sat-->", currSat
                    ssRdObj = read_ssusi.ProcessData( [root + dd + "/"],\
                                 prcsdFileDir, currDate )
                    ssRdObj.processed_data_to_file()
                    shutil.rmtree(root + dd + "/")
        currDate += tDelta