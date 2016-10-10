#!/usr/bin/python3

import sys
import os
import re
from subprocess import check_output

def downloader():
    if len(sys.argv) < 5 :
        raise RuntimeError("Need exactly four arguments!\n" + 
                            "USAGE: tleDownloader <ObjectID> <StartDate> <EndDate> <LoginFileName>\n"
                            + "Date Format: YYYY-MM-DD-HH:MM:SS (Do not use fractional seconds)\n"
                            + "Login file must contain the username on the first line and the password on the second line.")
    objectID = int(sys.argv[1])
    if objectID <= 0:
        raise RuntimeError("Object ID must be positive!")
    
    beginDate = sys.argv[2]
    if len(beginDate) != 19:
        raise RuntimeError("Begin date specification must contain 19 characters. Check the format")        
    
    endDate = sys.argv[3]
    if len(endDate) != 19:
        raise RuntimeError("End date specification must contain 19 characters. Check the format")

    loginFileLocation = os.path.join(os.getcwd(), sys.argv[4])
    if not os.path.exists(loginFileLocation) or not os.path.exists(sys.argv[4]):
        raise RuntimeError("Login file could not be found")
    loginFile = open(sys.argv[4])
    username = loginFile.readline()
    password = loginFile.readline()
    loginFile.close()

    outputName = re.sub('[:]','',"TLE_"+str(objectID).zfill(5)+"_"+beginDate+"--"+endDate+".txt")

    downloadCommand = ("wget --limit-rate=50K "+
    "--post-data='identity="+username+"&password="+password +
    "&query=https://www.space-track.org/basicspacedata/query/"+
    "class/tle/NORAD_CAT_ID/"+str(objectID).zfill(5)+
    "/EPOCH/"+beginDate[0:10]+"%20"+beginDate[11:20]+"--"+
    endDate[0:10]+"%20"+endDate[11:20]+"/orderby/EPOCH asc/format/tle' "+ 
    "--cookies=on --no-check-certificate --keep-session-cookies --save-cookies=cookies.txt "+
    "'https://www.space-track.org/ajaxauth/login' -O "+outputName)
    
    downloadCommand = re.sub('[\n]','',downloadCommand)    
    
    print(downloadCommand)  
    dwnCmd = open("dwnCmd.txt", "w")
    dwnCmd.write(downloadCommand)
    dwnCmd.close()
    
    commandOutput = check_output(downloadCommand, shell=True)
    
    #print(commandOutput)
    
downloader()
