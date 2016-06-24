# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 12:01:02 2015
@author: casadei_c
PRODUCE LIST OF ALL IMAGES OF ONE USER-SELECTED RUN IN ./Output/ImageLists/rXXXX_ImageNumbers_Filenames.txt.
"""

import sys
import getopt
import os

def makeImageListFunction(myArguments):  
    runNumber = ''
    imagesDirectoryName = ''

    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "imagesDirectoryName="])
    except getopt.GetoptError:
        print 'Usage: python makeImageList.py --runNumber <runNumber> --imagesDirectoryName <imagesDirectoryName>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python makeImageList.py --runNumber <runNumber> --imagesDirectoryName <imagesDirectoryName>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value
        elif option == "--imagesDirectoryName":
            imagesDirectoryName = value
    print 'Selected run is: %s'%runNumber
    print 'Selected folder is: %s '%imagesDirectoryName
       
    if not os.path.exists(imagesDirectoryName):
        print 'Directory does not exist.'
        return
        
    # PRODUCE IMAGE LIST IN ALPHABETICAL ORDER
    else:
        outputFolder = './Output_r%s/ImageLists'%runNumber
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
   
        myFiles = os.listdir(imagesDirectoryName)
        myFiles = sorted(myFiles)
        temp = open('%s/r%s_temporaryList.txt'%(outputFolder, runNumber), 'w')
        for f in myFiles:
            temp.write(f)
            temp.write('\n')
        temp.close()
   
        temp = open('%s/r%s_temporaryList.txt'%(outputFolder, runNumber), 'r')
        myFile = open('%s/r%s_imageFilenames.txt'%(outputFolder, runNumber), 'w')
        for l in temp:
            if 'LCLS' in l:
                if 'h5' in l:
                    myFile.write(l)
        temp.close()
        myFile.close()
   
        myFile = open('%s/r%s_imageFilenames.txt'%(outputFolder, runNumber), 'r')
        myFileOut = open('%s/r%s_ImageNumbers_Filenames.txt'%(outputFolder, runNumber), 'w')
        print 'Writing %s/r%s_ImageNumbers_Filenames.txt\n'%(outputFolder, runNumber)
        i = 0
        for l in myFile:
            i = i + 1
            myFileOut.write(str(i) + '\t' + l)
        myFile.close()
        myFileOut.close()
   
        os.remove('%s/r%s_temporaryList.txt'%(outputFolder, runNumber))
        os.remove('%s/r%s_imageFilenames.txt'%(outputFolder, runNumber))
 
if __name__ == "__main__":
    print "\n**** CALLING makeImageList ****"
    makeImageListFunction(sys.argv[1:])