import os
import datetime
import dwnld_timed_guvi
import read_timed_guvi
import shutil

if __name__ == "__main__":
    tempFileDir = "/home/bharat/Documents/code/data/timed-temp/"
    prcsdFileDir = "/home/bharat/Documents/code/data/timed-prcsd/"
    satList = [ "f16", "f17", "f18" ]
    tgDwnldObj = dwnld_timed_guvi.TimedGuviDownload(\
                    outBaseDir = tempFileDir)
    currDate = datetime.datetime( 2002, 3, 18 )
    endDate = datetime.datetime( 2002, 3, 18 )
    tDelta = datetime.timedelta(days=1)
    while currDate <= endDate:
        print "currently downloading files for --> ",\
            currDate.strftime("%Y-%m-%d")
        tgDwnldObj.download_files(currDate)
        for root, dirs, files in os.walk(tempFileDir):
            for nd, dd in enumerate(dirs):
                print "processing data --> ",\
                         currDate.strftime("%Y-%m-%d")
                tgRdObj = read_timed_guvi.ProcessTGData( [root + dd + "/"],\
                             prcsdFileDir, currDate )
                tgRdObj.processed_data_to_file()
                shutil.rmtree(root + dd + "/")
        currDate += tDelta