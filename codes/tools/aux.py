#!/usr/bin/env python
import sys, traceback
from math import *
import numpy as np
import matplotlib as mp
import os.path
import glob
import matplotlib.pyplot as plt
import scipy.interpolate
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from random import randint
import multiprocessing as mp
import optparse

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-o', '--outDir', dest='outDir', help='output directory', default='Job_CHESS', type='string')
(opt, args) = parser.parse_args()
outdir = opt.outDir

os.system('python main.py -o '+str(outdir))
IC = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[0,1]
NE = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[4,1]
Npartmax = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[12,1]
Npartmin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[13,1]
Multmax = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[14,1]
Multmin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[15,1]
etas = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[16,1]
zetas = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[17,1]
MC = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[21,1]
HBT = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[23,1]
nth = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[38,1]

#run TRENTO
gridstep = 0.2
gridmax = 10.0
N=np.ceil(2.0*gridmax/gridstep)
xymax= (N-1)*gridstep/2.0
xymin = -xymax
x = np.linspace(xymin,xymax,N)
y = np.linspace(xymin,xymax,N)
os.chdir('../trento/build/src')
os.system('./trento -c param.conf')
os.system('mv param.conf ../../../../Results/'+str(outdir)+'/IC/')
os.chdir('saida/')

a=0
zz=0
ZZ=0.
A=0.
B=0.
C=0.
D=0.
E=0.
F=0.
G=0.
H=0.
if IC == 0:
	for file in glob.glob('*.txt'):
		b = np.genfromtxt(file,comments='!', dtype=float)[1,1]
		mult = np.genfromtxt(file,comments='!', dtype=float)[3,1]
		npart = np.genfromtxt(file,comments='!', dtype=int)[2,1]
		e2 = np.genfromtxt(file,comments='!', dtype=float)[4,1]
		e3 = np.genfromtxt(file,comments='!', dtype=float)[5,1]
		e4 = np.genfromtxt(file,comments='!', dtype=float)[6,1]
		e5 = np.genfromtxt(file,comments='!', dtype=float)[7,1]
		B += b
		C += mult
		D += npart
		E += e2
		F += e3
		G += e4
		H += e5
		zz = zz + 1
		if(zz == NE):
                    break
	meanb = B/NE
	meanc = C/NE
	meand = D/NE
	meane = E/NE
	meanf = F/NE
	meang = G/NE
	meanh = H/NE
	for file in glob.glob('*.dat'):
		s = np.loadtxt(file)[:,:].transpose()
		splot = np.loadtxt(file)[:,:]
		smax = np.max(splot)
		nrows, ncols = int(N), int(N)
#		grid = splot.reshape((nrows, ncols))
#		plt.imshow(grid, extent=(x.min(), x.max(), y.max(), y.min()),
#        interpolation='bessel', cmap=cm.gnuplot)
#		plt.xlabel('x [fm]')
#		plt.ylabel('y [fm]')
#		cb = plt.colorbar(orientation='vertical', label='$s$ [1/fm$^{3}$]')
#		plt.savefig('ic_'+str(a)+'.png')
#		cb.remove()
#		plt.close()
		fp=open("out_"+str(a)+".dat","w")
		fp.write(str(int(N))+"	"+str(int(N))+"\n")
		for i in range(len(x)):
			for j in range(len(y)):
				fp.write(str(x[i])+"  "+str(y[j])+"  "+str(s[i][j])+"\n")
		fp.close()
		os.system('mv out_'+str(a)+'.dat ../../../../vhlle/ic/')
#		os.system('mv ic_'+str(a)+'.png ../../../../../Results/'+str(outdir)+'/IC')
		a=a+1
		if(a == NE):
                    break
        if (a < NE):
            print ("Change the Npart or Mult range!")
            os.system('rm -r ../../../../../Results/'+str(outdir))
            os.chdir('../../../../vhlle/params/')
            os.system('rm *')
            sys.exit(0)
else:
	for file in glob.glob('*.txt'):
		b = np.genfromtxt(file,comments='!', dtype=float)[1,1]
		mult = np.genfromtxt(file,comments='!', dtype=float)[3,1]
		npart = np.genfromtxt(file,comments='!', dtype=int)[2,1]
		e2 = np.genfromtxt(file,comments='!', dtype=float)[4,1]
		e3 = np.genfromtxt(file,comments='!', dtype=float)[5,1]
		e4 = np.genfromtxt(file,comments='!', dtype=float)[6,1]
		e5 = np.genfromtxt(file,comments='!', dtype=float)[7,1]
		B += b
		C += mult
		D += npart
		E += e2
		F += e3
		G += e4
		H += e5
		zz = zz + 1
		if(zz == NE):
                    break
	meanb = B/NE
	meanc = C/NE
	meand = D/NE
	meane = E/NE
	meanf = F/NE
	meang = G/NE
	meanh = H/NE
	for file in glob.glob('*.dat'):
		s = np.loadtxt(file)[:,:].transpose()
		splot = np.loadtxt(file)[:,:]
		A +=  s
		ZZ += splot
		a = a + 1
		if(a == NE):
                    break
        if (a < NE):
            print ("Change the Npart or Mult range!")
            os.system('rm -r ../../../../../Results/'+str(outdir))
            os.chdir('../../../../vhlle/params/')
            os.system('rm *')
            sys.exit(0)
	mean = A / NE
	meanplot = ZZ / NE
	nrows, ncols = int(N), int(N)
#	grid = meanplot.reshape((nrows, ncols))
#	plt.imshow(grid, extent=(x.min(), x.max(), y.max(), y.min()),
 #   interpolation='bessel', cmap=cm.gnuplot)
#	plt.xlabel('x [fm]')
#	plt.ylabel('y [fm]')
#	cb = plt.colorbar(orientation='vertical', label='$s$ [1/fm$^{3}$]')
#	plt.savefig('ic_smooth.png')
	fp=open("out_smooth.dat","w")
	fp.write(str(int(N))+"	"+str(int(N))+"\n")
	for i in range(len(x)):
	   	for j in range(len(y)):
	   		 fp.write(str(x[i])+"  "+str(y[j])+"  "+str(mean[i][j])+"\n")
	fp.close()
	os.system('mv out_smooth.dat ../../../../vhlle/ic/')
#	os.system('mv ic_smooth.png ../../../../../Results/'+str(outdir)+'/IC')
os.chdir('../')
os.system('rm -r saida/')
os.chdir('../../../../Results/')
fpp=open(str(outdir)+"/glauber_results.dat","w")
fpp.write("b                "+str(meanb)+"\n")
fpp.write("Multiplicity     "+str(meanc)+"\n")
fpp.write("Npart            "+str(int(meand))+"\n")
fpp.write("e2               "+str(meane)+ "\n")
fpp.write("e3               "+str(meanf)+ "\n")
fpp.write("e4               "+str(meang)+ "\n")
fpp.write("e5               "+str(meanh)+ "\n")
fpp.close()


#run vHLLE
os.chdir('../codes/vhlle/')
print ("\n")
print ("------------ vHLLE ------------")
print ("-------- Hydro Evolution ------------")

def f(el):
	return os.system('./hlle_visc params/event_'+str(el))
p = mp.Pool(nth)
q = range(0,NE)
if IC == 0:
    resultados = p.map(f, q)
    for i in range(0,NE):
        print ("Event_"+str(i))
        os.system('cp ../tools/idw2D  hydro_output/out_'+str(i))
        os.system('cp ../../Results/'+str(outdir)+'/thermo_params.dat  hydro_output/out_'+str(i))
        os.chdir('hydro_output/out_'+str(i))
        os.system('./idw2D')
        os.system('mv therm.xml therm_'+str(i)+'.xml')
        os.system('mv therm_'+str(i)+'.xml ../../../../Results/'+str(outdir)+'/THERM')
        tt = np.loadtxt('freezeout.dat')[:,0]
        xx = np.loadtxt('freezeout.dat')[:,1]
        yy = np.loadtxt('freezeout.dat')[:,2]
#        fig = plt.figure()
#        axx = fig.gca(projection='3d')
#        axx.scatter(xx, yy, tt,c='b', marker = "o")
#        axx.set_xlabel('x [fm]')
#        axx.set_ylabel('y [fm]')
#        axx.set_zlabel('$\\tau$ [fm]')
#        plt.savefig('freezeout_'+str(i)+'.png')
#        os.system('cp freezeout_'+str(i)+'.png ../../../../Results/'+str(outdir)+'/HYDRO')
        os.chdir('../../')
else:
    os.system('./hlle_visc params/event_smooth')
    os.system('cp ../tools/idw2D  hydro_output/out_smooth')
    os.system('cp ../../Results/'+str(outdir)+'/thermo_params.dat  hydro_output/out_smooth')
    os.chdir('hydro_output/out_smooth')
    os.system('./idw2D')
    os.system('mv therm.xml therm_smooth.xml')
    os.system('mv therm_smooth.xml ../../../../Results/'+str(outdir)+'/THERM/')
    tt = np.loadtxt('freezeout.dat')[:,0]
    xx = np.loadtxt('freezeout.dat')[:,1]
    yy = np.loadtxt('freezeout.dat')[:,2]
#    fig = plt.figure()
#    axx = fig.gca(projection='3d')
#    axx.scatter(xx, yy, tt,c='b', marker = "o")
#    axx.set_xlabel('x [fm]')
#    axx.set_ylabel('y [fm]')
#    axx.set_zlabel('$\\tau$ [fm]')
#    plt.savefig('freezeout_smooth.png')
#    plt.close()
#    os.system('cp freezeout_smooth.png ../../../../Results/'+str(outdir)+'/HYDRO')
    os.chdir('../../')

os.chdir('ic')
os.system('mv *.dat ../../../Results/'+str(outdir)+'/IC')
os.chdir('../params')
os.system('mv event* ../../../Results/'+str(outdir)+'/HYDRO')
os.chdir('../hydro_output')
os.system('rm -r *')
os.chdir('../../../Results/'+str(outdir)+'/THERM')
os.system('cp *.xml ../../../codes/therminator2/')
os.chdir('../../../codes/therminator2')
#run THERMINATOR 2
print ("\n")
print ("------------ THERMINATOR 2 ------------")
print ("-------------- Freezeout ------------")
print ("\n")

def g(ell):
	return os.system('./therm2_events therm_'+str(ell)+'.xml')
def k(elll):
	if etas!=0. or zetas!=0.:
		return os.system('./runfigure.sh events/lhyquid2vbi-therm_'+str(elll)+'/ '+str(int(MC/100)))
	else:
		return os.system('./runfigure.sh events/lhyquid2dbi-therm_'+str(elll)+'/ '+str(int(MC/100)))
def h(ellll):
	if etas!=0. or zetas!=0.:
		return os.system('./runhbt.sh events/lhyquid2vbi-therm_'+str(ellll)+'/ '+str(int(MC/100)))    
	else:
		return os.system('./runhbt.sh events/lhyquid2dbi-therm_'+str(ellll)+'/ '+str(int(MC/100)))    
	
pp = mp.Pool(nth)
qq = range(0,NE)
if IC == 0:
	results = pp.map(g, qq)
	figures = pp.map(k, qq)
	if HBT != 0:
		hbt = pp.map(h, qq)
else:
	os.system('./therm2_events therm_smooth.xml')
	if etas!=0. or zetas!=0.:
		os.system('./runfigure.sh events/lhyquid2vbi-therm_smooth/ '+str(int(MC/100)))
	else:
		os.system('./runfigure.sh events/lhyquid2dbi-therm_smooth/ '+str(int(MC/100)))
	if HBT != 0:
		if etas!=0. or zetas!=0.:
			os.system('./runhbt.sh events/lhyquid2vbi-therm_smooth/ '+str(int(MC/100)))
		else:
			os.system('./runhbt.sh events/lhyquid2dbi-therm_smooth/ '+str(int(MC/100)))

os.system('rm *.xml')
os.chdir('events/')
if IC == 0:
	for i in range(0,NE):
		if etas!=0. or zetas!=0.:
			os.chdir('lhyquid2vbi-therm_'+str(i))
		else:
			os.chdir('lhyquid2dbi-therm_'+str(i))
		os.system('mv *.root fmult* *.eps ../../../../Results/'+str(outdir)+'/OBS/Event_'+str(i))
		os.system('rm  *.ini *.xml')
		os.system('mv fig_disteta.dat ../../../../Results/'+str(outdir)+'/OBS/disteta_'+str(i)+'.dat')
		os.system('mv fig_distpt.dat  ../../../../Results/'+str(outdir)+'/OBS/distpt_'+str(i)+'.dat')
		os.system('mv fig_distptLHC.dat  ../../../../Results/'+str(outdir)+'/OBS/LHCdistpt_'+str(i)+'.dat')
		os.system('mv fig_v1chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v1_'+str(i)+'.dat')
		os.system('mv fig_v2chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v2_'+str(i)+'.dat')
		os.system('mv fig_v3chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v3_'+str(i)+'.dat')
		os.system('mv fig_v4chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v4_'+str(i)+'.dat')
		os.system('mv fig_v5chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v5_'+str(i)+'.dat')
		os.system('mv fig_v6chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v6_'+str(i)+'.dat')
		if HBT != 0:
			os.system('mv hbtradii.dat ../../../../Results/'+str(outdir)+'/OBS/hbtradii_'+str(i)+'.dat')
		os.chdir('../')
else:
	if etas!=0. or zetas!=0.:
		os.chdir('lhyquid2vbi-therm_smooth')
	else:
		os.chdir('lhyquid2dbi-therm_smooth')
	os.system('mv *.root fmult* ../../../../Results/'+str(outdir)+'/OBS/')
	os.system('rm *.eps *.ini *.xml')
	os.system('mv fig_disteta.dat ../../../../Results/'+str(outdir)+'/OBS/disteta_smooth.txt')
	os.system('mv fig_distpt.dat  ../../../../Results/'+str(outdir)+'/OBS/distpty_smooth.txt')
	os.system('mv fig_distptLHC.dat  ../../../../Results/'+str(outdir)+'/OBS/distpteta_smooth.txt')
	os.system('mv fig_v1chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v1_smooth.txt')
	os.system('mv fig_v2chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v2_smooth.txt')
	os.system('mv fig_v3chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v3_smooth.txt')
	os.system('mv fig_v4chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v4_smooth.txt')
	os.system('mv fig_v5chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v5_smooth.txt')
	os.system('mv fig_v6chpt2.dat ../../../../Results/'+str(outdir)+'/OBS/v6_smooth.txt')
	if HBT != 0:
		os.system('mv hbtradii.dat ../../../../Results/'+str(outdir)+'/OBS')
	os.chdir('../')
os.system('rm -r *')
os.chdir('../../tools')
if IC == 0:
	os.system('python mean.py -o '+str(outdir))
else:
	if HBT!=0:
		kt1 = 0.2
		kt2 = 0.3
		kt3 = 0.4
		kt4 = 0.525
#kt=0.2
		n1 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[0,1]
		lamb1 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[1,1]
		Rout1 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[2,1]
		Rside1 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[3,1]
		Rlong1 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[4,1]
		n1err = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[0,2]
		lamb1err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[1,2]
		Rout1err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[2,2]
		Rside1err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[3,2]
		Rlong1err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[4,2]
#kt=0.3
		n2 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[5,1]
		lamb2 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[6,1]
		Rout2 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[7,1]
		Rside2 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[8,1]
		Rlong2 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[9,1]
		n2err = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[5,2]
		lamb2err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[6,2]
		Rout2err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[7,2]
		Rside2err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[8,2]
		Rlong2err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[9,2]
#kt=0.4
		n3 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[10,1]
		lamb3 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[11,1]
		Rout3 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[12,1]
		Rside3 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[13,1]
		Rlong3 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[14,1]
		n3err = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[10,2]
		lamb3err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[11,2]
		Rout3err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[12,2]
		Rside3err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[13,2]
		Rlong3err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[14,2]
#kt=0.525
		n4 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[15,1]
		lamb4 = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[16,1]
		Rout4 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[17,1]
		Rside4 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[18,1]
		Rlong4 =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[19,1]
		n4err = np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[15,2]
		lamb4err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[16,2]
		Rout4err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[17,2]
		Rside4err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[18,2]
		Rlong4err =  np.genfromtxt('../../Results/'+str(outdir)+'/OBS/hbtradii.dat', dtype=float)[19,2]
		fp=open("../../Results/"+str(outdir)+"/OBS/hbt_smooth.txt","w")
		fp.write("#kt         Norm      errNorm        Lamb         Lamberr     Rout        Routerr         Rside       Rsideerr        Rlong       Rlongerr\n")
		fp.write(str(kt1)+"  "+str(n1)+"  "+str(n1err)+"  "+str(lamb1)+"  "+str(lamb1err)+"  "+str(Rout1)+"  "+str(Rout1err)+"  "+str(Rside1)+"  "+str(Rside1err)+"  "+str(Rlong1)+"  "+str(Rlong1err)+"\n")
		fp.write(str(kt2)+"  "+str(n2)+"  "+str(n2err)+"  "+str(lamb2)+"  "+str(lamb2err)+"  "+str(Rout2)+"  "+str(Rout2err)+"  "+str(Rside2)+"  "+str(Rside2err)+"  "+str(Rlong2)+"  "+str(Rlong2err)+"\n")
		fp.write(str(kt3)+"  "+str(n3)+"  "+str(n3err)+"  "+str(lamb3)+"  "+str(lamb3err)+"  "+str(Rout3)+"  "+str(Rout3err)+"  "+str(Rside3)+"  "+str(Rside3err)+"  "+str(Rlong3)+"  "+str(Rlong3err)+"\n")
		fp.write(str(kt4)+"  "+str(n4)+"  "+str(n4err)+"  "+str(lamb4)+"  "+str(lamb4err)+"  "+str(Rout4)+"  "+str(Rout4err)+"  "+str(Rside4)+"  "+str(Rside4err)+"  "+str(Rlong4)+"  "+str(Rlong4err)+"\n")
		fp.close()
