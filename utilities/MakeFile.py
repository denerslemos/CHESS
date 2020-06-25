#!/usr/bin/env python


import os.path


#make TRENTO
os.chdir('../codes/trento/')
os.system('mkdir build')
os.chdir('build/')
os.system('cmake ../')
os.system('make')
#make vHLLE
os.chdir('../../vhlle/')
os.system('make')
#IDW
os.chdir('../tools/')
os.system('g++ -O3 idw2D.cpp fast.cpp -o idw2D')
#choose the HBT fit
HBT = 0
if HBT == 0:
	os.chdir('../therminator2/build/src')
	os.system('cp therm2_hbtfitgauss.cxx therm2_hbtfit.cxx')
elif HBT == 1:
	os.chdir('../therminator2/build/src')
	os.system('cp therm2_hbtfitexp.cxx therm2_hbtfit.cxx')
#elif HBT == 2:
#	os.chdir('../therminator2/build/src')
#	os.system('cp therm2_hbtfitlevy.cxx therm2_hbtfit.cxx')    !FUTURE
#make THERMINATOR2
os.chdir('../../')
os.system('make')

print "Execute the code with: python run.py -o output"
