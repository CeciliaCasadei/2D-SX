# -*- coding: utf-8 -*-
import getopt
import sys
import os
import pickle
import curses

def verify_browseAssembledImagesFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "imageFolderName="])
    except getopt.GetoptError:
        print 'Usage: python verify_browseAssembledImages.py --runNumber <runNumber> --imageFolderName <imageFolderName>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python verify_browseAssembledImages.py --runNumber <runNumber> --imageFolderName <imageFolderName>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value
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
            nImages = 0
            for imageKey, imageObject in imageObjects.items():
                nImages = nImages + 1
                
            
            n = 0
            while 1:
                screen = curses.initscr()
                try:
                    curses.noecho()
                    curses.curs_set(0)
                    screen.keypad(1)
                    screen.addstr("\nPress the right arrow to go forth, left arrow to go back or any key to quit.")
                    event = screen.getch()
                finally:
                    curses.endwin()
                
                if event == curses.KEY_LEFT or event == curses.KEY_UP:
                    n = n-1
                    if n < 1:
                        n = 1
                    for imageKey, imageObject in imageObjects.items():
                        if imageObject.runNumber == runNumber and imageObject.imageNumber == str(n):                        
                            imageObject.displayImage(imageFolderName)
                elif event == curses.KEY_RIGHT or event == curses.KEY_DOWN:
                    n = n+1
                    if n > nImages:
                        n = nImages
                    for imageKey, imageObject in imageObjects.items():
                        if imageObject.runNumber == runNumber and imageObject.imageNumber == str(n):                        
                            imageObject.displayImage(imageFolderName)
                            
                else:
                    break
                


                    


if __name__ == "__main__":
    print "\n**** SHOW ASSEMBLED IMAGE ****"
    verify_browseAssembledImagesFunction(sys.argv[1:])