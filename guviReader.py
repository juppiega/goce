#!/usr/bin/python

import sys
import os
import numpy
from scipy.io import readsav, savemat

def guviReader():
    if len(sys.argv) < 2:
        raise RuntimeError("Give the GUVI location as argument. \n" + 
        "This is the folder where all the years reside.")
    
    rootdir = sys.argv[1]
    if not os.path.isdir(rootdir):
        raise RuntimeError(rootdir + " is not a directory!")
    #print "GUVI folder: ", rootdir
    
    outputLocation = os.path.join(os.getcwd(), 'guviData')
    if not os.path.exists(outputLocation):
        os.makedirs(outputLocation)
    
    fileCount = 0
    for root, subdirs, files in os.walk(rootdir):
        for name in files:
            if ".sav" in name:
                fileCount += 1    
    
    i = 1
    for root, subdirs, files in os.walk(rootdir):
        for name in files:
            if ".sav" in name:
                processFile(root, name, outputLocation)
                print str(i) + " / " + str(fileCount)
                i += 1

def processFile(rootName, guviFileName, outputLocation):
    guviFileFullPath = os.path.join(rootName, guviFileName)
    guviFile = readsav(guviFileFullPath, python_dict=True)
    dataArray = guviFile['ndpsorbit']
    varNames = dataArray.dtype.names
    matlabDic = {}
    
    for varName in varNames:
        dataField = dataArray[varName]
        if not isinstance(dataField, numpy.core.records.recarray):
            numRecords = len(dataField)
            recordLen = len(dataField[0])
            varArray = numpy.concatenate(dataField)
            varArray = numpy.reshape(varArray, (recordLen, numRecords),
                                     order='F')
        else:
            varArray = numpy.array(dataField)
        matlabDic[varName] = varArray
    
    matfileName = guviFileName.split(".")[0] + ".mat"
    matfilePath = os.path.join(outputLocation, matfileName)
    savemat(matfilePath, matlabDic)
    
    
    #raise RuntimeError("Stop here!")
    
guviReader()
