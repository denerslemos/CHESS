#!/usr/bin/env python


import os.path

#os.chdir('~/')
#ROOT
os.system('sudo apt-get install git dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev \
libxft-dev libxext-dev')

#TRENTO
os.system('sudo apt-get install libboost-dev libboost-{filesystem,program-options}-dev libhdf5-dev')

#Python
os.system('sudo apt-get install python python-numpy python-matplotlib python-scipy python-tk')

