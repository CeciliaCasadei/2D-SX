# -*- coding: utf-8 -*-
# INSTRUCTIONS: python Test-runSetup.py
# SETUP FOR CURRENT RUN

import os

# LOG CURRENT SETTINGS
os.system('cp ./Test-runSetup.py ./Test-runSetup.log')

# PARAMETERS
firstNumber = 6.3
secondNumber = 2.7

# RUN SCRIPT
os.system('python Test-CommandLineArgs.py --1stNumber %f --2ndNumber %f'%(firstNumber, secondNumber))