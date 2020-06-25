#!/usr/bin/env python


from math import *
import numpy as np
import matplotlib as mp
import os.path
import glob
import time
import optparse

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-o', '--outDir', dest='outDir', help='output directory', default='Job_CHESS', type='string')
(opt, args) = parser.parse_args()
outdir = opt.outDir

start = time.time()
os.chdir('codes/tools/')
os.system('python aux.py -o '+str(outdir))
end =  time.time()
sec = end-start
a = sec//60//60//24
b = (sec//60//60)%24
c = (sec//60)%60
d = sec%60
print "Execution Time:   "+str(int(a))+" days   "+str(int(b))+" hours   "+str(int(c))+" minutes   "+str(d)+" seconds\n"
os.chdir('../../')
fp=open("runningtime.dat","w")
fp.write("Execution Time:   "+str(int(a))+" days   "+str(int(b))+" hours   "+str(int(c))+" minutes   "+str(d)+" seconds\n")
fp.close()
os.system('mv runningtime.dat Results/'+str(outdir)+'/')

	
