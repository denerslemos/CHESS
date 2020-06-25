#!/usr/bin/env python

from math import *
import numpy as np
import matplotlib as mpl
import os.path
import glob
import scipy.interpolate
import time
import optparse

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-o', '--outDir', dest='outDir', help='output directory', default='Job_CHESS', type='string')
(opt, args) = parser.parse_args()
outdir = opt.outDir

#Read Parameter file
IC = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[0,1]
rIC = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[1,1]
p1 = np.genfromtxt('../../Parameters.dat',comments='!', dtype=str)[2,1]
p2 = np.genfromtxt('../../Parameters.dat',comments='!', dtype=str)[3,1]
NE = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[4,1]
p = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[5,1]
k = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[6,1]
w = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[7,1]
SNN = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[8,1]
kappa = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[9,1]
bmin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[10,1]
bmax = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[11,1]
Npartmax= np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[12,1]
Npartmin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[13,1]
Mmax = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[14,1]
Mmin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[15,1]
etas = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[16,1]
zetas = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[17,1]
Tfo = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[18,1]/1000. # -> convert to GeV
tau0 = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[19,1]
eos = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[20,1]
MC = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[21,1]
rMC = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[22,1]
HBT = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[23,1]
Pair = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[24,1]
Mix = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[25,1]
OP = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[26,1]
Lambda = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[27,1]
Rout = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[28,1]
RoutMin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[29,1]
RoutMax = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[30,1]
Rside = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[31,1]
RsideMin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[32,1]
RsideMax = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[33,1]
Rlong = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[34,1]
RlongMin = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[35,1]
RlongMax = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[36,1]
MaxFitRange = np.genfromtxt('../../Parameters.dat',comments='!', dtype=float)[37,1]


print ("CHESS")
print ("Complete Hydrodynamical Evolution SyStem")
print ("by D. S. Lemos and O. Socolowski Jr.")
print ("If you need help, contact us by denerslemos@ift.unesp.br\n")
print ("Codes:")
print ("IC - TRENTO  --  http://qcd.phy.duke.edu/trento/usage.html and https://arxiv.org/pdf/1412.4708.pdf\n")
print ("Hydro - vHLLE  --  https://github.com/yukarpenko/vhlle and https://arxiv.org/pdf/1312.4160.pdf\n")
print ("Freezeout - THERMINATOR2  --  https://therminator2.ifj.edu.pl/ and https://arxiv.org/pdf/1102.0273.pdf")
print ("\n")

if IC == 0:
  print ("IC: Event-by-event")
else:
  print ("IC: Smooth")
print "System              "+str(p1)+"+"+str(p2)
print "Number of events    "+str(NE)

#Equation of State List

#s95p-v1
if eos==7:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("s95p-v1.dat")[:,0]
  pr = np.loadtxt("s95p-v1.dat")[:,1]
  T = np.loadtxt("s95p-v1.dat")[:,2]
  print ("s95p-v1 EoS")

#FO EoS
elif eos==2:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("FO.dat")[:,0]
  pr = np.loadtxt("FO.dat")[:,1]
  T = np.loadtxt("FO.dat")[:,2]
  print ("First Order Phase Transition EoS")

#Laine EoS
elif eos==4:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("Laine_nf3.dat")[:,0]
  pr = np.loadtxt("Laine_nf3.dat")[:,1]
  T = np.loadtxt("Laine_nf3.dat")[:,2]
  print ("Laine EoS")

#WB EoS
elif eos==5:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("WB.dat")[:,0]
  pr = np.loadtxt("WB.dat")[:,1]
  T = np.loadtxt("WB.dat")[:,2]
  print ("WB EoS")

#HotQCD EoS
elif eos==6:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("HOTQCD.dat")[:,0]
  pr = np.loadtxt("HOTQCD.dat")[:,1]
  T = np.loadtxt("HOTQCD.dat")[:,2]
  print ("HotQCD EoS")

#s95n-v1
elif eos==8:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("s95n-v1.dat")[:,0]
  pr = np.loadtxt("s95n-v1.dat")[:,1]
  T = np.loadtxt("s95n-v1.dat")[:,2]
  print ("s95n-v1 EoS")

#s90f-v1
elif eos==9:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("s90f-v1.dat")[:,0]
  pr = np.loadtxt("s90f-v1.dat")[:,1]
  T = np.loadtxt("s90f-v1.dat")[:,2]
  print ("s95n-v1 EoS")

#PCE EoS
elif eos==10:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("PCE-EOS.dat")[:,0]
  pr = np.loadtxt("PCE-EOS.dat")[:,1]
  T = np.loadtxt("PCE-EOS.dat")[:,2]
  print ("PCE EoS")

#PCE EoS 145
elif eos==11:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("PCE-EOS-145.dat")[:,0]
  pr = np.loadtxt("PCE-EOS-145.dat")[:,1]
  T = np.loadtxt("PCE-EOS-145.dat")[:,2]
  print ("PCE EoS 145")

#PCE EoS 150
elif eos==12:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("PCE-EOS-150.dat")[:,0]
  pr = np.loadtxt("PCE-EOS-150.dat")[:,1]
  T = np.loadtxt("PCE-EOS-150.dat")[:,2]
  print ("PCE EoS 150")

#PCE EoS 155
elif eos==13:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("PCE-EOS-155.dat")[:,0]
  pr = np.loadtxt("PCE-EOS-155.dat")[:,1]
  T = np.loadtxt("PCE-EOS-155.dat")[:,2]
  print ("PCE EoS 155")

#PCE EoS 160
elif eos==14:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("PCE-EOS-160.dat")[:,0]
  pr = np.loadtxt("PCE-EOS-160.dat")[:,1]
  T = np.loadtxt("PCE-EOS-160.dat")[:,2]
  print ("PCE EoS 160")

#quasi particle EoS (not public yet)
elif eos==3:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("quasiEOS.dat")[:,0]
  pr = np.loadtxt("quasiEOS.dat")[:,1]
  T = np.loadtxt("quasiEOS.dat")[:,2]
  print ("quasi particle EoS")

#Hadron Gas Resonance
elif eos==1:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("HRG.dat")[:,0]
  pr = np.loadtxt("HRG.dat")[:,1]
  T = np.loadtxt("HRG.dat")[:,2]
  print ("Hadron Gas Resonance")

#Ideal Ultrarelativistic Gas
elif eos==0:
  os.chdir("../vhlle/eos")
  e = np.loadtxt("IURG.dat")[:,0]
  pr = np.loadtxt("IURG.dat")[:,1]
  T = np.loadtxt("IURG.dat")[:,2]
  print ("Ideal Ultrarelativistic Gas")


#Get cs2
h = 0.00000000001
pnew = scipy.interpolate.interp1d(T, pr) 
enew = scipy.interpolate.interp1d(T, e)
pfo = pnew(Tfo)
efo = enew(Tfo)
cs2 = ((pnew(Tfo + h) - pnew(Tfo))/h)/((enew(Tfo + h) - enew(Tfo))/h)
os.chdir("../../tools")


#Generate parameters file for codes individually

#TRENTO
fp=open("../trento/build/src/param.conf","w")
fp.write("projectile = "+str(p1)+"\n")
fp.write("projectile = "+str(p2)+"\n")
fp.write("number-events = "+str(NE)+" \n")
fp.write("output = saida\n")
fp.write("reduced-thickness = "+str(p)+"\n")
fp.write("fluctuation = "+str(k)+"\n")
fp.write("nucleon-width = "+str(w)+" \n")
fp.write("cross-section = "+str(SNN)+"\n")
fp.write("normalization = "+str(kappa)+"\n")
fp.write("npart-max = "+str(Npartmax)+" \n")
fp.write("npart-min = "+str(Npartmin)+" \n")
fp.write("mult-max = "+str(Mmax)+" \n")
fp.write("mult-min = "+str(Mmin)+" \n")
fp.write("nucleon-min-dist = 0.4 \n") #same used in RHIC and LHC
fp.write("b-min = "+str(bmin)+"\n")
fp.write("b-max = "+str(bmax)+"\n")
fp.write("grid-max = 10.\n")
fp.write("grid-step = 0.2 \n")
if rIC == 0:
	fp.write("random-seed = 3\n")
	fp.close()
else:
	fp.close()

if HBT == 0:
	print ("HBT: off")
else:
	print ("HBT: on")

#vHLLE
if IC == 0:
	for i in range(0,NE): 
		fpv=open("../vhlle/params/event_"+str(i),"w")
		fpv.write("outputDir       hydro_output/out_"+str(i)+"	   ! output directory\n")

		fpv.write("eosType       "+str(eos)+"           ! reserved for future use\n")

		fpv.write("etaS      "+str(etas)+"             !  eta/s value\n")
		fpv.write("zetaS     "+str(zetas)+"            !  zeta/s, bulk viscosity\n")
 		fpv.write("e_crit    "+str(efo)+"    !  criterion for surface finding 153\n")

		fpv.write("nx        161             ! number of cells in X direction\n")
		fpv.write("ny        161             ! number of cells in Y direction\n")
		fpv.write("nz        5               ! number of cells in eta direction\n")
		fpv.write("xmin     -16.0            ! coordinate of the first cell\n")
		fpv.write("xmax      16.0            ! coordinate of the last cell\n")
		fpv.write("ymin     -16.0\n")
		fpv.write("ymax      16.0\n")
		fpv.write("etamin   -4.0\n")
		fpv.write("etamax    4.0\n")

		fpv.write("! icModel: 1=Glauber, 2=Glauber_table+parametrized rapidity for energy conserv. check, 4=Gubser flow\n")
		fpv.write("icModel    2\n")
		fpv.write("glauberVar 1       ! reserved for future use\n")
		fpv.write("icInputFile ic/out_"+str(i)+".dat\n")
		fpv.write("s0ScaleFactor   1.0\n")
		fpv.write("epsilon0   30.0\n")
		fpv.write("impactPar  0.0\n")
		fpv.write("Rgt      0.0\n")
		fpv.write("Rgz      0.0\n")
		fpv.write("tau0      "+str(tau0)+"     ! starting proper time\n")
		fpv.write("tauMax     50.0    ! proper time to stop hydro\n")
		fpv.write("dtau       0.05    ! timestep\n")
		fpv.close()
else:
	fpv=open("../vhlle/params/event_smooth","w")
	fpv.write("outputDir       hydro_output/out_smooth	   ! output directory\n")

	fpv.write("eosType       "+str(eos)+"           ! reserved for future use\n")

	fpv.write("etaS      "+str(etas)+"             !  eta/s value\n")
	fpv.write("zetaS     "+str(zetas)+"            !  zeta/s, bulk viscosity\n")
	fpv.write("e_crit    "+str(efo)+"    !  criterion for surface finding 153\n")

	fpv.write("nx        161             ! number of cells in X direction\n")
	fpv.write("ny        161             ! number of cells in Y direction\n")
	fpv.write("nz        5               ! number of cells in eta direction\n")
	fpv.write("xmin     -16.0            ! coordinate of the first cell\n")
	fpv.write("xmax      16.0            ! coordinate of the last cell\n")
	fpv.write("ymin     -16.0\n")
	fpv.write("ymax      16.0\n")
	fpv.write("etamin   -4.0\n")
	fpv.write("etamax    4.0\n")
	
	fpv.write("! icModel: 1=Glauber, 2=Glauber_table+parametrized rapidity for energy conserv. check, 4=Gubser flow\n")
	fpv.write("icModel    2\n")
	fpv.write("glauberVar 1       ! reserved for future use\n")
	fpv.write("icInputFile ic/out_smooth.dat\n")
	fpv.write("s0ScaleFactor   1.0\n")
	fpv.write("epsilon0   30.0\n")
	fpv.write("impactPar  0.0\n")
	fpv.write("Rgt      0.0\n")
	fpv.write("Rgz      0.0\n")
	fpv.write("tau0       "+str(tau0)+"     ! starting proper time\n")
	fpv.write("tauMax     50.0    ! proper time to stop hydro\n")
	fpv.write("dtau       0.05    ! timestep\n")
	fpv.close()

# THERMINATOR2
#events.ini
fpt=open("../therminator2/events.ini","w")
fpt.write("[FreezeOut]\n")
if etas!=0. or zetas!=0.:
	fpt.write("FreezeOutModel = Lhyquid2VBI\n")
	fpt.write("\n")
else:
	fpt.write("FreezeOutModel = Lhyquid2DBI\n")
	fpt.write("\n")
fpt.write("[Event]\n")
fpt.write("NumberOfEvents = "+str(MC)+"\n")
fpt.write("EventFileType = root")
fpt.write("\n")
fpt.write("[Primordial]\n")
fpt.write("MultiplicityDistribution = Poisson\n")
fpt.write("IntegrateSamples = "+str(1000000)+"\n")
fpt.write("\n")
fpt.write("[Random]\n")
fpt.write("Randomize = "+str(rMC)+"\n")
fpt.write("\n")
fpt.write("[Directories]\n")
fpt.write("ShareDir = share/\n")
fpt.write("FreezeOutDir = fomodel/\n")
fpt.write("MacroDir = macro/\n")
fpt.write("EventDir = events/\n")
fpt.write("\n")
fpt.write("[Logging]\n")
fpt.write("LogFile = therminator.log")
fpt.close()

#femto.ini
fptt=open("../therminator2/femto.ini","w")
fptt.write("[Pair]\n")
fptt.write("\n")
if(Pair!=1):
	fptt.write("PairType = pion-pion\n")
else:
	fptt.write("PairType = kaon-kaon\n")
fptt.write("\n")
fptt.write("[Cuts]\n")
fptt.write("\n")
fptt.write("TimeCut = 500.0\n")
fptt.write("\n")
fptt.write("[Event]\n")
fptt.write("\n")
fptt.write("EventsToMix = "+str(Mix))
fptt.write("\n")
fptt.write("[Switches]\n")
fptt.write("\n")
if(OP!=0):
	fptt.write("EnableOnlyPrimordial = no\n")
else:
	fptt.write("EnableOnlyPrimordial = yes\n")
fptt.write("\n")
fptt.write("EnableSourceHistograms = no\n")
fptt.write("\n")
fptt.write("[Logging]\n")
fptt.write("LogFile = therminator.log")
fptt.close()

#hbt.ini
fph=open("../therminator2/hbtfit.ini","w")
fph.write("[Normalization_Fit]\n")
fph.write("Norm = 1.0\n")
fph.write("NormFitType = fixed\n")
fph.write("\n")
fph.write("[Lambda_Fit]\n")
fph.write("Lambda = "+str(Lambda)+"\n")
fph.write("LambdaFitType = free\n")
fph.write("\n")
fph.write("[Rout_Fit]\n")
fph.write("Rout = "+str(Rout)+"\n")
fph.write("RoutFitType = limit\n")
fph.write("RoutMin = "+str(RoutMin)+"\n")
fph.write("RoutMax = "+str(RoutMax)+"\n")
fph.write("\n")
fph.write("[Rside_Fit]\n")
fph.write("Rside = "+str(Rside)+"\n")
fph.write("RsideFitType = limit\n")
fph.write("RsideMin = "+str(RsideMin)+"\n")
fph.write("RsideMax = "+str(RsideMax)+"\n")
fph.write("\n")
fph.write("[Rlong_Fit]\n")
fph.write("Rlong = "+str(Rlong)+"\n")
fph.write("RlongFitType = limit\n")
fph.write("RlongMin = "+str(RlongMin)+"\n")
fph.write("RlongMax = "+str(RlongMax)+"\n")
fph.write("\n")
fph.write("[QRange]\n")
fph.write("MaxFitRange = "+str(MaxFitRange)+"\n")
fph.write("\n")
fph.write("[Histograms]\n")
fph.write("Numerator = cnuma\n")
fph.write("Denominator = cdena\n")
fph.write("\n")
fph.write("[Logging]\n")
fph.write("LogFile = therminator.log")
fph.close()

print ("\n")
print ("------------ TRENTO ------------")
print ("------- Initial Conditions ------------")
print ("\n")

os.chdir('../../Results/')
os.system("mkdir "+str(outdir))
os.chdir(str(outdir))
os.system("mkdir IC")
os.system("mkdir HYDRO")
os.system("mkdir THERM")
os.system("mkdir OBS")
os.chdir("OBS")
for i in range(0,NE):
    os.system("mkdir Event_"+str(i))
os.chdir('../../')
fpp=open(str(outdir)+"/thermo_params.dat","w")
fpp.write("etas                "+str(etas)+"\n")
fpp.write("zetas               "+str(zetas)+"\n")
fpp.write("Tfo                 "+str(Tfo)+"\n")
fpp.write("Efo                 "+str(efo)+ "\n")
fpp.write("Pfo                 "+str(pfo)+ "\n")
fpp.write("cs2                 "+str(cs2)+ "\n")
fpp.write("tau0                "+str(tau0)+ "\n")
fpp.close()
os.chdir('../')
os.system('cp Parameters.dat Results/'+str(outdir)+'/')

