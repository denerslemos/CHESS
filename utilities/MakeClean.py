#!/usr/bin/env python


import os.path



os.chdir('../codes/trento/')
os.system('rm -r build')

os.chdir('../vhlle/')
os.system('make clean')

os.chdir('../tools/')
os.system('rm -r idw2D')


os.chdir('../therminator2/build/src')
os.system('rm therm2_hbtfit.cxx')
os.chdir('../../')
os.system('rm events.ini hbtfit.ini femto.ini')
os.system('make clean')

