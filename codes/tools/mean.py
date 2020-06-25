#!/usr/bin/env python


from math import *
import numpy as np
import matplotlib as mp
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

#Mean of event-by-event hydrodynamics
NE = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[4,1]
p1 = np.genfromtxt('../../Parameters.dat',comments='!', dtype=str)[2,1]
p2 = np.genfromtxt('../../Parameters.dat',comments='!', dtype=str)[3,1]
HBT = np.genfromtxt('../../Parameters.dat',comments='!', dtype=int)[34,1]
os.chdir('../../Results/'+str(outdir)+'/OBS')

#All the results are from Charged Hadrons (You can change this in CHESS-version/codes/therminator2/macro/).
#Flows are calculated with Reaction Plane Method {RP}

#Eta
A=0.
B=0.
C=0.
for file in glob.glob('disteta*.dat'):
	eta = np.loadtxt(file)[:,0]
	dndeta = np.loadtxt(file)[:,1]
	erreta = np.loadtxt(file)[:,2]
	A+=eta
	B+=dndeta
	C+=erreta
meaneta    = A/NE
meandndeta = B/NE
meanerreta = C/NE
fp=open("disteta_ebe.txt","w")

for i in range(len(eta)):
	fp.write(str(meaneta[i])+"  "+str(meandndeta[i])+"  "+str(meanerreta[i])+"\n")
fp.close()
	
#dN/ptdptdy
D=0.
E=0.
F=0.
for file in glob.glob('distpt*.dat'):
	pty = np.loadtxt(file)[:,0]
	dndpty = np.loadtxt(file)[:,1]
	errpty = np.loadtxt(file)[:,2]
	D+=pty
	E+=dndpty
	F+=errpty
meanpty    = D/NE
meandndpty = E/NE
meanerrpty = F/NE
fp=open("distpty_ebe.txt","w")
for i in range(len(pty)):
	fp.write(str(meanpty[i])+"  "+str(meandndpty[i])+"  "+str(meanerrpty[i])+"\n")
fp.close()


#dN/ptdptdeta
AA=0.
BB=0.
CC=0.
for file in glob.glob('LHC*.dat'):
	pteta = np.loadtxt(file)[:,0]
	dndpteta = np.loadtxt(file)[:,1]
	errpteta = np.loadtxt(file)[:,2]
	AA+=pteta
	BB+=dndpteta
	CC+=errpteta
meanpteta    = AA/NE
meandndpteta = BB/NE
meanerrpteta = CC/NE
fp=open("distpteta_ebe.txt","w")
for i in range(len(pteta)):
	fp.write(str(meanpteta[i])+"  "+str(meandndpteta[i])+"  "+str(meanerrpteta[i])+"\n")
fp.close()


#Diretic Flow - v1(pt)
aaaa=0.
bbbb=0.
cccc=0.
for file in glob.glob('v1*.dat'):
	ptv1 = np.loadtxt(file)[:,0]
	v1 = np.loadtxt(file)[:,1]
	errv1 = np.loadtxt(file)[:,2]
	aaaa+=ptv1
	bbbb+=v1
	cccc+=errv1
meanptv1    = aaaa/NE
meanv1      = bbbb/NE
meanerrv1   = cccc/NE
fp=open("v1_ebe.txt","w")
for i in range(len(ptv1)):
	fp.write(str(meanptv1[i])+"  "+str(meanv1[i])+"  "+str(meanerrv1[i])+"\n")
fp.close()

#Elliptic Flow - v2(pt)
G=0.
H=0.
I=0.
for file in glob.glob('v2*.dat'):
	ptv2 = np.loadtxt(file)[:,0]
	v2 = np.loadtxt(file)[:,1]
	errv2 = np.loadtxt(file)[:,2]
	G+=ptv2
	H+=v2
	I+=errv2
meanptv2    = G/NE
meanv2      = H/NE
meanerrv2   = I/NE
fp=open("v2_ebe.txt","w")
for i in range(len(ptv2)):
	fp.write(str(meanptv2[i])+"  "+str(meanv2[i])+"  "+str(meanerrv2[i])+"\n")
fp.close()

#Triangular Flow - v3(pt)
J=0.
L=0.
M=0.
for file in glob.glob('v3*.dat'):
	ptv3 = np.loadtxt(file)[:,0]
	v3 = np.loadtxt(file)[:,1]
	errv3 = np.loadtxt(file)[:,2]
	J+=ptv3
	L+=v3
	M+=errv3
meanptv3    = J/NE
meanv3      = L/NE
meanerrv3   = M/NE
fp=open("v3_ebe.txt","w")
for i in range(len(ptv3)):
	fp.write(str(meanptv3[i])+"  "+str(meanv3[i])+"  "+str(meanerrv3[i])+"\n")
fp.close()

#Quadrangular Flow - v4(pt)
N=0.
O=0.
P=0.
for file in glob.glob('v4*.dat'):
	ptv4 = np.loadtxt(file)[:,0]
	v4 = np.loadtxt(file)[:,1]
	errv4 = np.loadtxt(file)[:,2]
	N+=ptv4
	O+=v4
	P+=errv4
meanptv4    = N/NE
meanv4      = O/NE
meanerrv4   = P/NE
fp=open("v4_ebe.txt","w")
for i in range(len(ptv4)):
	fp.write(str(meanptv4[i])+"  "+str(meanv4[i])+"  "+str(meanerrv4[i])+"\n")
fp.close()

# v5(pt)
Q=0.
R=0.
S=0.
for file in glob.glob('v5*.dat'):
	ptv5 = np.loadtxt(file)[:,0]
	v5 = np.loadtxt(file)[:,1]
	errv5 = np.loadtxt(file)[:,2]
	Q+=ptv5
	R+=v5
	S+=errv5
meanptv5    = Q/NE
meanv5      = R/NE
meanerrv5   = S/NE
fp=open("v5_ebe.txt","w")
for i in range(len(ptv5)):
	fp.write(str(meanptv5[i])+"  "+str(meanv5[i])+"  "+str(meanerrv5[i])+"\n")
fp.close()

#v6(pt)
T=0.
U=0.
V=0.
for file in glob.glob('v6*.dat'):
	ptv6 = np.loadtxt(file)[:,0]
	v6 = np.loadtxt(file)[:,1]
	errv6 = np.loadtxt(file)[:,2]
	T+=ptv6
	U+=v6
	V+=errv6
meanptv6    = T/NE
meanv6      = U/NE
meanerrv6   = V/NE
fp=open("v6_ebe.txt","w")
for i in range(len(ptv6)):
	fp.write(str(meanptv6[i])+"  "+str(meanv6[i])+"  "+str(meanerrv6[i])+"\n")
fp.close()

if HBT!=0:
	kt1 = 0.2
	kt2 = 0.3
	kt3 = 0.4
	kt4 = 0.525
	AAA=0.
	BBB=0.
	CCC=0.
	DDD=0.
	EEE=0.
	FFF=0.
	GGG=0.
	HHH=0.
	III=0.
	JJJ=0.
	KKK=0.
	LLL=0.
	MMM=0.
	NNN=0.
	OOO=0.
	PPP=0.
	QQQ=0.
	RRR=0.
	SSS=0.
	TTT=0.
	UUU=0.
	VVV=0.
	XXX=0.
	YYY=0.
	WWW=0.
	ZZZ=0.
	AAAA=0.
	BBBB=0.
	CCCC=0.
	DDDD=0.
	EEEE=0.
	FFFF=0.
	GGGG=0.
	HHHH=0.
	IIII=0.
	JJJJ=0.
	KKKK=0.
	LLLL=0.
	MMMM=0.
	NNNN=0.
	for file in glob.glob('hbt*.dat'):
#kt=0.2
		n1 = np.genfromtxt(file, dtype=float)[0,1]
		lamb1 = np.genfromtxt(file, dtype=float)[1,1]
		Rout1 =  np.genfromtxt(file, dtype=float)[2,1]
		Rside1 =  np.genfromtxt(file, dtype=float)[3,1]
		Rlong1 =  np.genfromtxt(file, dtype=float)[4,1]
		n1err = np.genfromtxt(file, dtype=float)[0,2]
		lamb1err =  np.genfromtxt(file, dtype=float)[1,2]
		Rout1err =  np.genfromtxt(file, dtype=float)[2,2]
		Rside1err =  np.genfromtxt(file, dtype=float)[3,2]
		Rlong1err =  np.genfromtxt(file, dtype=float)[4,2]
#kt=0.3
		n2 = np.genfromtxt(file, dtype=float)[5,1]
		lamb2 = np.genfromtxt(file, dtype=float)[6,1]
		Rout2 =  np.genfromtxt(file, dtype=float)[7,1]
		Rside2 =  np.genfromtxt(file, dtype=float)[8,1]
		Rlong2 =  np.genfromtxt(file, dtype=float)[9,1]
		n2err = np.genfromtxt(file, dtype=float)[5,2]
		lamb2err =  np.genfromtxt(file, dtype=float)[6,2]
		Rout2err =  np.genfromtxt(file, dtype=float)[7,2]
		Rside2err =  np.genfromtxt(file, dtype=float)[8,2]
		Rlong2err =  np.genfromtxt(file, dtype=float)[9,2]
#kt=0.4
		n3 = np.genfromtxt(file, dtype=float)[10,1]
		lamb3 = np.genfromtxt(file, dtype=float)[11,1]
		Rout3 =  np.genfromtxt(file, dtype=float)[12,1]
		Rside3 =  np.genfromtxt(file, dtype=float)[13,1]
		Rlong3 =  np.genfromtxt(file, dtype=float)[14,1]
		n3err = np.genfromtxt(file, dtype=float)[10,2]
		lamb3err =  np.genfromtxt(file, dtype=float)[11,2]
		Rout3err =  np.genfromtxt(file, dtype=float)[12,2]
		Rside3err =  np.genfromtxt(file, dtype=float)[13,2]
		Rlong3err =  np.genfromtxt(file, dtype=float)[14,2]
#kt=0.525
		n4 = np.genfromtxt(file, dtype=float)[15,1]
		lamb4 = np.genfromtxt(file, dtype=float)[16,1]
		Rout4 =  np.genfromtxt(file, dtype=float)[17,1]
		Rside4 =  np.genfromtxt(file, dtype=float)[18,1]
		Rlong4 =  np.genfromtxt(file, dtype=float)[19,1]
		n4err = np.genfromtxt(file, dtype=float)[15,2]
		lamb4err =  np.genfromtxt(file, dtype=float)[16,2]
		Rout4err =  np.genfromtxt(file, dtype=float)[17,2]
		Rside4err =  np.genfromtxt(file, dtype=float)[18,2]
		Rlong4err =  np.genfromtxt(file, dtype=float)[19,2]
#Means
		AAA+=n1
		BBB+=lamb1
		CCC+=Rout1
		DDD+=Rside1
		EEE+=Rlong1
		FFF+=n1err
		GGG+=lamb1err
		HHH+=Rout1err
		III+=Rside1err
		JJJ+=Rlong1err
		KKK+=n2
		LLL+=lamb2
		MMM+=Rout2
		NNN+=Rside2
		OOO+=Rlong2
		PPP+=n2err
		QQQ+=lamb2err
		RRR+=Rout2err
		SSS+=Rside2err
		TTT+=Rlong2err
		UUU+=n3
		VVV+=lamb3
		XXX+=Rout3
		YYY+=Rside3
		WWW+=Rlong3
		ZZZ+=n3err
		AAAA+=lamb3err
		BBBB+=Rout3err
		CCCC+=Rside3err
		DDDD+=Rlong3err
		EEEE+=n4
		FFFF+=lamb4
		GGGG+=Rout4
		HHHH+=Rside4
		IIII+=Rlong4
		JJJJ+=n4err
		KKKK+=lamb4err
		LLLL+=Rout4err
		MMMM+=Rside4err
		NNNN+=Rlong4err
	meanAAA = AAA/NE
	meanBBB = BBB/NE
	meanCCC = CCC/NE
	meanDDD = DDD/NE
	meanEEE = EEE/NE
	meanFFF = FFF/NE
	meanGGG = GGG/NE
	meanHHH = HHH/NE
	meanIII = III/NE
	meanJJJ = JJJ/NE
	meanKKK = KKK/NE
	meanLLL = LLL/NE
	meanMMM = MMM/NE
	meanNNN = NNN/NE
	meanOOO = OOO/NE
	meanPPP = PPP/NE
	meanQQQ = QQQ/NE
	meanRRR = RRR/NE
	meanSSS = SSS/NE
	meanTTT = TTT/NE
	meanUUU = UUU/NE
	meanVVV = VVV/NE
	meanXXX = XXX/NE
	meanYYY = YYY/NE
	meanWWW = WWW/NE
	meanZZZ = ZZZ/NE
	meanAAAA = AAAA/NE
	meanBBBB = BBBB/NE
	meanCCCC = CCCC/NE
	meanDDDD = DDDD/NE
	meanEEEE = EEEE/NE
	meanFFFF = FFFF/NE
	meanGGGG = GGGG/NE
	meanHHHH = HHHH/NE
	meanIIII = IIII/NE
	meanJJJJ = JJJJ/NE
	meanKKKK = KKKK/NE
	meanLLLL = LLLL/NE
	meanMMMM = MMMM/NE
	meanNNNN = NNNN/NE
	fp=open("hbt_ebe.txt","w")
	fp.write("#kt         Norm      errNorm        Lamb         Lamberr     Rout        Routerr         Rside       Rsideerr        Rlong       Rlongerr\n")
	fp.write(str(kt1)+"  "+str(meanAAA)+"  "+str(meanFFF)+"  "+str(meanBBB)+"  "+str(meanGGG)+"  "+str(meanCCC)+"  "+str(meanHHH)+"  "+str(meanDDD)+"  "+str(meanIII)+"  "+str(meanEEE)+"  "+str(meanJJJ)+"\n")
	fp.write(str(kt2)+"  "+str(meanKKK)+"  "+str(meanPPP)+"  "+str(meanLLL)+"  "+str(meanQQQ)+"  "+str(meanMMM)+"  "+str(meanRRR)+"  "+str(meanNNN)+"  "+str(meanSSS)+"  "+str(meanOOO)+"  "+str(meanTTT)+"\n")
	fp.write(str(kt3)+"  "+str(meanUUU)+"  "+str(meanZZZ)+"  "+str(meanVVV)+"  "+str(meanAAAA)+"  "+str(meanXXX)+"  "+str(meanBBBB)+"  "+str(meanYYY)+"  "+str(meanCCCC)+"  "+str(meanWWW)+"  "+str(meanDDDD)+"\n")
	fp.write(str(kt4)+"  "+str(meanEEEE)+"  "+str(meanJJJJ)+"  "+str(meanFFFF)+"  "+str(meanKKKK)+"  "+str(meanGGGG)+"  "+str(meanLLLL)+"  "+str(meanHHHH)+"  "+str(meanMMMM)+"  "+str(meanIIII)+"  "+str(meanNNNN)+"\n")
	fp.close()
#os.system('rm -r *.dat')
