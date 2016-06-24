# *- coding: utf-8 -*-
"""
Created on Tue Jan  5 18:18:28 2016
@author: casadei_c
SHOW INFO STORED IN diffractionImage OBJECTS AND DISPLAY THE ASSOCIATED ASSEMBLED IMAGE.
"""
import getopt
import sys
import pickle
import os

            
def verifyPklContentFunction(myArguments):
    myPklFile = ''
    myImageFileName = ''
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["myPklFile=", "myImageFileName="])
    except getopt.GetoptError:
        print 'Usage: python verify_pklContent.py --myPklFile <myPklFile> --myImageFileName <myImageFileName>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python verify_pklContent.py --myPklFile <myPklFile> --myImageFileName <myImageFileName>'
            sys.exit()
        elif option == "--myPklFile":
            myPklFile = value
        elif option == "--myImageFileName":
            myImageFileName = value
    
    # VERIFY SELECTED IMAGE OBJECT
    if not os.path.exists(myPklFile):
        print 'File %s not found'%myPklFile
    else:
        openPkl = open(myPklFile, 'rb')
        myData = pickle.load(openPkl)
        openPkl.close()
        
        k = 0
        for i, j in myData.items():
            k = k + 1
        print 'Number of dictionary items: %d'%k
                
        print '\n* Image %s *'%myImageFileName
        print '\nFile name: %s'%myData['%s'%myImageFileName].fileName
        print 'Run number: %s'%myData['%s'%myImageFileName].runNumber
        print 'Tilt angle: %s degrees'%myData['%s'%myImageFileName].tiltAngle
        print 'Image number: %s'%myData['%s'%myImageFileName].imageNumber
        print 'Selection flag: %s'%myData['%s'%myImageFileName].selectionFlag

        try:
            print 'Wavelength: %f A'%myData['%s'%myImageFileName].wavelength
        except:
            print 'Wavelength not available'
        try:
            print 'Photon energy: %f eV'%myData['%s'%myImageFileName].photonEnergy
        except:
            print 'Photon energy not available'
        try:
            print 'Number of peaks: %d'%myData['%s'%myImageFileName].nPeaks
        except:
            print 'Number of peaks not available'
        try:
            print 'Peaks matrix:'
            print '              x              y              I'
            for peaksMatrixLine in myData['%s'%myImageFileName].peaksMatrix:
                print '%15.3f%15.3f%15.3f'%(peaksMatrixLine[0], peaksMatrixLine[1], peaksMatrixLine[2])
            print '\nx = column index in unassembled matrix'
            print 'y = row index in unassembled matrix\n'
        except:
            print 'Peaks matrix not available'

            
if __name__ == "__main__":
    print "\n**** VERIFY diffractionImage OBJECTS DICTIONARY CONTENT ****"
    verifyPklContentFunction(sys.argv[1:])