# -*- coding: utf-8 -*-
import sys, getopt

def myFunction(myArguments):
    firstNumber = 0
    secondNumber = 0
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["1stNumber=","2ndNumber="])
    except getopt.GetoptError:
        print 'Usage: python Test-CommandLineArgs.py --1stNumber <1stNumber> --2ndNumber <2ndNumber>'
        sys.exit(2)
   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python Test-CommandLineArgs.py --1stNumber <1stNumber> --2ndNumber <2ndNumber>'
            sys.exit()
        elif option == "--1stNumber":
            firstNumber = float(value)
        elif option == "--2ndNumber":
            secondNumber = float(value)
    result = firstNumber + secondNumber
    print result

  
if __name__ == "__main__":
   myFunction(sys.argv[1:])
   
   
### RUN e.g.: python Test_CommandLineArgs.py --1stNumber 5.8 --2ndNumber 7.2
### RUN e.g.: python Test_CommandLineArgs.py -h