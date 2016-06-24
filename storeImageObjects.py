# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 13:31:25 2016
@author: casadei_c
GENERATE diffractionImage OBJECTS AND STORE THEM IN A DICTIONARY IN './Output/ExtractExperimentalInfo/rXXXX_imagesDictionary.pkl'
"""
import getopt
import sys
import os
import pickle
import diffractionImageClass

def storeImageObjectsFunction(myArguments):
    runNumber = ''
    tiltAngle = 0
    imageListDirectory = ''

    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "tiltAngle=", "imageListDirectory="])
    except getopt.GetoptError:
        print 'Usage: python storeImageObjects.py --runNumber <runNumber> --tiltAngle <tiltAngle> --imageListDirectory <imageListDirectory>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python storeImageObjects.py --runNumber <runNumber> --tiltAngle <tiltAngle> --imageListDirectory <imageListDirectory>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value
        elif option == "--tiltAngle":
            tiltAngle = float(value)
        elif option == "--imageListDirectory":
            imageListDirectory = value
         
    # PRODUCE diffractionImage OBJECTS 
    imageList = '%s/r%s_ImageNumbers_Filenames.txt'%(imageListDirectory, runNumber)
    print 'Selected image list file is: %s '%imageList
    if not os.path.exists(imageList):
        print 'Directory does not exist.'
        return        
    else:
        imagesDictionary = {}
        myImageFiles = open(imageList, 'r')
        for line in myImageFiles:
            mySplitLine = line.split()
            imageNumber = mySplitLine[0]
            fileName = mySplitLine[1]
            imageSelectionFlag = 0
            imageObject = diffractionImageClass.diffractionImage(fileName, runNumber, imageNumber, tiltAngle, imageSelectionFlag)
            imagesDictionary['%s'%fileName] = imageObject
        myImageFiles.close()
        
        outputFolder = './Output_r%s/ExtractExperimentalInfo'%runNumber
        if not os.path.exists(outputFolder):   
            os.mkdir(outputFolder)          
        print 'STORING diffractionImage OBJECTS IN FILE %s/r%s_imagesDictionary.pkl\n'%(outputFolder, runNumber)
        f = open('%s/r%s_imagesDictionary.pkl'%(outputFolder, runNumber), 'wb') 
        pickle.dump(imagesDictionary, f)
        f.close()
        
if __name__ == "__main__":
    print "\n**** CALLING storeImageObjects ****"
    storeImageObjectsFunction(sys.argv[1:])