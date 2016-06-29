# -*- coding: utf-8 -*-
import getopt
import sys
import os
import pickle

def verifyAssembledImagesFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "imageNumber=", "imageFolderName="])
    except getopt.GetoptError:
        print 'Usage: python verify_assembledImages.py --runNumber <runNumber> --imageNumber <imageNumber> --imageFolderName <imageFolderName>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python verify_assembledImages.py --runNumber <runNumber> --imageNumber <imageNumber> --imageFolderName <imageFolderName>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value
        elif option == "--imageNumber":
            imageNumber = value
        elif option == "--imageFolderName":
            imageFolderName = value
            
    # DISPLAY ASSEMBLED IMAGE
    pklFile = './Output_r%s/ExtractExperimentalInfo/r%s_imagesDictionary.pkl'%(runNumber, runNumber)
    if not os.path.exists(pklFile):
        print 'File %s not found'%pklFile
    else:
        openPkl = open(pklFile, 'rb')
        imageObjects = pickle.load(openPkl)
        openPkl.close()
        if not os.path.exists(imageFolderName):
            print 'File %s not found'%imageFolderName
        else:
            for imageKey, imageObject in imageObjects.items():
                if imageObject.runNumber == runNumber and imageObject.imageNumber == imageNumber:
                    imageObject.displayImage(imageFolderName)
                    


if __name__ == "__main__":
    print "\n**** SHOW ASSEMBLED IMAGE ****"
    verifyAssembledImagesFunction(sys.argv[1:])