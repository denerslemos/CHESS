import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
import matplotlib as mp
import pylab as plt

e = np.loadtxt("IURG.dat")[:,0]
p = np.loadtxt("IURG.dat")[:,1]
#s = np.loadtxt("FO.dat")[:,2]
T = np.loadtxt("IURG.dat")[:,2]


s = (e+p)/T
#et4 = e/T**(4)



#p_interp = Spline(e, p)
#s_interp = Spline(e, s) 
#T_interp = Spline(e, T)

#e_table = np.arange(7.385+0.001477,794.,0.001477)
#p_table = p_interp(e_table)
#s_table = s_interp(e_table)
#T_table = T_interp(e_table)


#for i in range(1, len(e)):
    #print ("e_table[i] - e_table[i-1]"+str(e_table[i] - e_table[i-1]))
fpp=open("IURGnew.dat","w")
#for i in range(len(e2)):
#    fpp.write(str(e2[i])+"         "+str(p2[i])+"         "+str(s2[i])+"         "+str(T2[i])+"\n")
for j in range(len(e)):
    fpp.write(str(e[j])+"         "+str(p[j])+"         "+str(s[j])+"         "+str(T[j])+"\n")

#fpp.close()
#plt.plot(T, et4,'sr')  

#plt.plot(T2, et42,'-b')  
#plt.show()
