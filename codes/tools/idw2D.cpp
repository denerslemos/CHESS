#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <omp.h>
#include<stdio.h>
#include<stdlib.h>
#include <cstring>
#include <ctime>
#include <sstream>
#include "fast.h"

#define PI ( 4.*atan(1.) )
#define hc ( 0.197326968 )


int  main() {

using namespace std;

//Read Thermo_Params File
int N;
clock_t sec0,sec1;
int a,b,r,c;
const double p = 7.; // interpolation power 
int n = 91; // number of points
double etas, zetas, Tfo, cs2, tau0, efo, pfo;
ifstream file("thermo_params.dat");
if (!file.is_open()) {
    cout << "cannot open parameters file " << "thermo_params.dat" << endl;
    exit(1);
}
char parName[255], parValue[255];
while (file.good()) {
    string line;
    getline(file, line);
    istringstream sline(line);
    sline >> parName >> parValue;
if (strcmp(parName, "etas") == 0)
      etas = atof(parValue);  
else if (strcmp(parName, "zetas") == 0)
      zetas = atof(parValue);  
else if (strcmp(parName, "Tfo") == 0)
      Tfo = atof(parValue);  
else if (strcmp(parName, "Efo") == 0)
      efo = atof(parValue);  
else if (strcmp(parName, "Pfo") == 0)
      pfo = atof(parValue);  
else if (strcmp(parName, "cs2") == 0)
      cs2 = atof(parValue);  
else if (strcmp(parName, "tau0") == 0)
      tau0 = atof(parValue);   
}

// Read Freezeout File
double t, x, y, z, u0, u1, u2, u3;
double pi11, pi12, pi22, Pi;
vector<double> tt, xx, yy, zz,  gamma, ux, uy, uz;
vector<double> pixx, pixy, piyy, pizz, PPi;
tt.clear();
xx.clear();
yy.clear();
zz.clear();
gamma.clear();
ux.clear();
uy.clear();
uz.clear();
PPi.clear();
pixx.clear();
piyy.clear();
pixy.clear();


ifstream fin("freezeout.dat");
sec0 = clock();
while(!fin.eof()){
    fin   >> t 
          >> x  
          >> y  
          >> z 
          >> u0  
          >> u1  
          >> u2  
          >> u3
          >> pi11  
          >> pi12  
          >> pi22   
          >> Pi;
    tt.push_back(t);
    xx.push_back(x);
    yy.push_back(y);
    zz.push_back(z);
    gamma.push_back(u0);
    ux.push_back(u1);
    uy.push_back(u2);
    uz.push_back(u3);
    pixx.push_back(pi11);
    pixy.push_back(pi12);
    piyy.push_back(pi22);
    PPi.push_back(Pi);
};
N=tt.size();
double* tau = new double[N];
double* zeta = new double[N];
double* phif= new double[N];
double* phi = new double[N];
double* phii = new double[N];
double* d = new double[N];
double* dd = new double[N];
double* phivt= new double[N];
double* phiivt= new double[N];
double* vx = new double[N];
double* vy = new double[N];
double* vt = new double[N];

//    change the coordenates
for (int i = 0; i < N; i++){   
//     tau
    tau[i] = tt[i] - tau0;
//        phi
    phii[i] = atan2(yy[i],xx[i]);
	if (phii[i] < 0.0){
      phi[i]= phii[i] + 2 * PI;
    }
    else{
      phi[i]= phii[i];    
    }
//     zeta
    zeta[i] =  atan(tau[i]/sqrt(xx[i]*xx[i]+yy[i]*yy[i]));    
   // velocity
	vx[i] = ux[i] / gamma[i];
	vy[i] = uy[i] / gamma[i];
   //    Fluid  Transverse Velocity
    vt[i] = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
    //        Phivt
    phiivt[i] = atan2(vy[i],vx[i]);
    if (phiivt[i] < 0.0){
      phivt[i] = phiivt[i] + 2 * PI;
    }
    else{
      phivt[i]= phiivt[i];    
    }

    phif[i] = phivt[i] - phi[i];
    dd[i]=sqrt(tau[i]*tau[i]+xx[i]*xx[i]+yy[i]*yy[i]); 
    d[i]=(1/hc)*dd[i]; 
  }

  // Inverse Distance Weighting (IDW) - Interpolation
double W = 0.0; 
double dist = 0.0;
double A = 0.0;
double B = 0.0;
double C = 0.0;
double D = 0.0;
double I = 0.0;
double J = 0.0;
double L = 0.0;
double P = 0.0;

double phinew = 0.0;
double zetanew = 0.0;
double phimin= 0.0;
double phimax= 2.*PI;
double zetamin= 0.0;
double zetamax= PI / 2.;

double *Wnew,*totA,*totB,*totC,*totD,*totI,*totJ,*totL,*totP,*dnew,*vtnew,*phifnew, *pixxnew, *pixynew, *piyynew, *Pinew;
Wnew    = new double [n*n];
totA    = new double [n*n];
totB    = new double [n*n];
totD    = new double [n*n];
totI    = new double [n*n];
totJ    = new double [n*n];
totL    = new double [n*n];
totP    = new double [n*n];
dnew    = new double [n*n];
vtnew   = new double [n*n];
phifnew = new double [n*n];
pixxnew = new double [n*n];
pixynew = new double [n*n];
piyynew = new double [n*n];
Pinew = new double [n*n];
cout << endl; 
cout << "Parametrizing Freezeout Surface" << endl;
cout << endl;  

for ( a = 0; a < n; a++){
  zetanew = zetamin + a*(zetamax - zetamin)/(n-1);
  for ( b = 0; b < n; b++){
     phinew = phimin + b*(phimax - phimin)/(n-1);
     for ( c = 0; c < N; c++){
         dist = pow_fast( (zetanew -zeta[c])*(zetanew - zeta[c])+
                         (phinew -phi[c])*(phinew - phi[c]), 0.5 );
     if (dist != 0.){
            W = pow_fast(dist,-p);
            A =    d[c] * W;
            B =   vt[c] * W;
            D = phif[c] * W;
            
        if (etas != 0. or zetas != 0.){
            I = pixx[c] * W;
            J = pixy[c] * W;
            L = piyy[c] * W;
            P =  PPi[c] * W;
            }
            
            Wnew[(n*a+b)] += W;
            totA[(n*a+b)] += A;
            totB[(n*a+b)] += B;
            totD[(n*a+b)] += D;
        if (etas != 0. or zetas != 0.){
            totI[(n*a+b)] += I;
            totJ[(n*a+b)] += J;
            totL[(n*a+b)] += L;
            totP[(n*a+b)] += P;
            }
            } else {
              W = 1.;
              A =    d[c] * W;
              B =   vt[c] * W;
              D = phif[c] * W;
            if (etas != 0. or zetas != 0.){
              I = pixx[c] * W;
              J = pixy[c] * W;
              L = piyy[c] * W;
              P =  PPi[c] * W;
			}
              Wnew[(n*a+b)] = W;
              totA[(n*a+b)] = A;
              totB[(n*a+b)] = B;
              totD[(n*a+b)] = D; 
            if (etas != 0. or zetas != 0.){
              totI[(n*a+b)] = I;
              totJ[(n*a+b)] = J;
              totL[(n*a+b)] = L;
              totP[(n*a+b)] = P; 
			  }
        }
 }

    dnew[(n*a+b)] = totA[(n*a+b)] / Wnew[(n*a+b)];
    vtnew[(n*a+b)] = totB[(n*a+b)] / Wnew[(n*a+b)];
    phifnew[(n*a+b)] = totD[(n*a+b)] / Wnew[(n*a+b)];
            if (etas != 0. or zetas != 0.){
    pixxnew[(n*a+b)] = totI[(n*a+b)] / Wnew[(n*a+b)];
    pixynew[(n*a+b)] = totJ[(n*a+b)] / Wnew[(n*a+b)];
    piyynew[(n*a+b)] = totL[(n*a+b)] / Wnew[(n*a+b)];
    Pinew[(n*a+b)] = totP[(n*a+b)] / Wnew[(n*a+b)];
	  }
        Wnew[(n*a+b)]=0.0;      
        totA[(n*a+b)]=0.0;
        totB[(n*a+b)]=0.0;
        totD[(n*a+b)]=0.0;
            if (etas != 0. or zetas != 0.){
        totI[(n*a+b)]=0.0;
        totJ[(n*a+b)]=0.0;
        totL[(n*a+b)]=0.0;
        totP[(n*a+b)]=0.0; 	
      	}
	} 
}

ofstream fon("therm.xml", ios_base::out);
fon << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n"
<< "<HYPERSURFACE filename=\"FURG-QGP\" version=\"1.0\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://therminator2.ifj.edu.pl/hypersurface.xsd\">" << "\n"
<< "<PARAMETERS>"<< "\n"
<<  "  <PARAMETER name=\"Tau_i\" unit=\"fm\">"<< tau0 <<"</PARAMETER>"  << "\n"
<<  "  <PARAMETER name=\"Temperature\" unit=\"MeV\">"<< (Tfo*1000) <<"</PARAMETER>"  << "\n"
<<  "  <PARAMETER name=\"Mu_B\" unit=\"MeV\">0.</PARAMETER>"  << "\n"
<<  "  <PARAMETER name=\"Mu_I\" unit=\"MeV\">0.</PARAMETER>"  << "\n"
<<  "  <PARAMETER name=\"Mu_S\" unit=\"MeV\">0.</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"Mu_C\" unit=\"MeV\">0.</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"device\">LHC</PARAMETER>" << "\n";
if (etas != 0. or zetas != 0.){
fon<<  "  <PARAMETER name=\"cs^2\" unit=\"unit\">"<< cs2 << "</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"eta/s\" unit=\"unit\">"<< etas << "</PARAMETER>" << "\n"
 <<  "  <PARAMETER name=\"zeta/s\" unit=\"unit\">"<< zetas << "</PARAMETER>" << "\n";
};
fon<<  "  <PARAMETER name=\"colliding_system\">AuAu</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"colliding_energy\" unit=\"GeV\">5500.</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"centrality_min\" unit=\"%\">30.</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"centrality_max\" unit=\"%\">40.</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"impact_parameter\" unit=\"fm\">-1.</PARAMETER>" << "\n"
<<  "  <PARAMETER name=\"temperature_at_center\" unit=\"MeV\">500.</PARAMETER>" << "\n"
<<  "</PARAMETERS>" << "\n"
<<  "<DESCRIPTION>" << "\n"
<<  "  <ITEM name=\"initial_nuclear_profile\">" << "\n"
<<  "    <DETAIL name=\"type\">gaussian fit to Glissando</DETAIL>" << "\n"
<<  "    <DETAIL name=\"sigma_x\" unit=\"fm\">2.547</DETAIL>" << "\n"
<<  "    <DETAIL name=\"sigma_y\" unit=\"fm\">1.947</DETAIL>" << "\n"
<<  "  </ITEM>" << "\n"
<<  "  <ITEM name=\"initial_transverse_flow\">not available</ITEM>" << "\n"
<<  "  <ITEM name=\"equatiom_of_state\">" << "\n"
<<  "    <DETAIL name=\"T_gt_TC\">lattice QCD P(T) fit [hep-ph/0607079]</DETAIL>" << "\n"
<<  "    <DETAIL name=\"T_eq_TC\">smooth cross-over</DETAIL>" << "\n"
<<  "    <DETAIL name=\"T_lt_TC\">massive hadron gas, SHARE input [nucl-th/0404083]</DETAIL>" << "\n"
<<  "  </ITEM>" << "\n"
<<  "  <ITEM name=\"eos_to_initial_condition\">entropy density profile</ITEM>" << "\n"
<<  "  <ITEM name=\"hydrodynamic_equation\">" << "\n"
<<  "    <DETAIL name=\"type\">2+1D boost-invariant</DETAIL>" << "\n"
<<  "    <DETAIL name=\"r_min\" unit=\"fm\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"r_max\" unit=\"fm\">15.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"r_pts\">122</DETAIL>" << "\n"
<<  "    <DETAIL name=\"phi_min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"phi_max\" unit=\"fm\">6.28319</DETAIL>" << "\n"
<<  "    <DETAIL name=\"phi_pts\">61</DETAIL>" << "\n"
<<  "    <DETAIL name=\"t_min\" unit=\"fm\">"<< tau0 <<"</DETAIL>" << "\n"
<<  "    <DETAIL name=\"t_max\" unit=\"fm\">20.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"t_pts\">179</DETAIL>"<< "\n"
<<  "  </ITEM>" << "\n"
<<  "</DESCRIPTION>" << "\n"
<<  "<!-- Included VECTOR3D names : Distance, FluidVt, FluidPhi, Viscous -->" << "\n"
<<  "<VECTOR3D name=\"Distance\" unit=\"GeV^-1\">" << "\n"
<<  "  <AXIS name=\"Zeta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.5708</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">"<< n <<"</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Phi\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">6.28319</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">"<< n <<"</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Theta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">1</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<<  "  <DATA points=\""<< pow(n,2) <<"\">" << "\n";
for (a = 0; a < n; a++){
  for (b = 0; b < n; b++){
      fon <<   dnew[(n*a+b)] << endl;
  }
}

fon <<"  </DATA>" << "\n"
<< "</VECTOR3D>" << "\n"
<< "<VECTOR3D name=\"FluidVt\" unit=\"1\">" << "\n"
<< "  <AXIS name=\"Zeta\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">1.5708</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n <<"</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<< "  <AXIS name=\"Phi\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">6.28319</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n <<"</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Theta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">1</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<< "  <DATA points=\""<< pow(n,2) <<"\"> " << "\n";
for (a = 0; a < n; a++){
  for (b = 0; b < n; b++){
      fon <<   vtnew[(n*a+b)] << endl;
  }
}


fon << "  </DATA>" << "\n"
<< "</VECTOR3D>" << "\n"
<< "<VECTOR3D name=\"FluidPhi\" unit=\"1\">" << "\n"
<< "  <AXIS name=\"Zeta\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">1.5708</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n <<"</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<< "  <AXIS name=\"Phi\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">6.28319</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n <<"</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Theta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">1</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<< "  <DATA points=\""<< pow(n,2) <<"\"> " << "\n";
for (a = 0; a < n; a++){
  for (b = 0; b < n; b++){
      fon <<   phifnew[n*a+b] << endl;
  }
}

// Pixx
if (etas != 0. or zetas != 0.){
fon <<"  </DATA>" << "\n"
<< "</VECTOR3D>" << "\n"
<< "<VECTOR3D name=\"Pixx\" unit=\"1\">" << "\n"
<< "  <AXIS name=\"Zeta\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">1.5708</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<< "  <AXIS name=\"Phi\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">6.28319</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Theta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">1</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<< "  <DATA points=\""<< pow(n,2) <<"\"> " << "\n";
for (a = 0; a < n; a++){
  for (b = 0; b < n; b++){
      fon <<   (pixxnew[(n*a+b)]/(efo+pfo)) << endl;
  }
}  

// Pixy
  
fon <<"  </DATA>" << "\n"
<< "</VECTOR3D>" << "\n"
<< "<VECTOR3D name=\"Pixy\" unit=\"1\">" << "\n"
<< "  <AXIS name=\"Zeta\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">1.5708</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<< "  <AXIS name=\"Phi\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">6.28319</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Theta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">1</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<< "  <DATA points=\""<< pow(n,2) <<"\"> " << "\n";
for (a = 0; a < n; a++){
  for (b = 0; b < n; b++){
      fon <<   (pixynew[(n*a+b)]/ (efo+pfo)) << endl;
  }
}  
// Piyy
fon <<"  </DATA>" << "\n"
<< "</VECTOR3D>" << "\n"
<< "<VECTOR3D name=\"Piyy\" unit=\"1\">" << "\n"
<< "  <AXIS name=\"Zeta\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">1.5708</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<< "  <AXIS name=\"Phi\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">6.28319</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Theta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">1</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<< "  <DATA points=\""<< pow(n,2) <<"\"> " << "\n";  
for (a = 0; a < n; a++){
  for (b = 0; b < n; b++){
      fon <<   (piyynew[(n*a+b)]/ (efo+pfo)) << endl;
  }
}  
// PI - bulk
fon <<"  </DATA>" << "\n"
<< "</VECTOR3D>" << "\n"
<< "<VECTOR3D name=\"PI\" unit=\"1\">" << "\n"
<< "  <AXIS name=\"Zeta\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">1.5708</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<< "  <AXIS name=\"Phi\">" << "\n"
<< "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<< "    <DETAIL name=\"max\" unit=\"rad\">6.28319</DETAIL>" << "\n"
<< "    <DETAIL name=\"pts\">"<< n << "</DETAIL>" << "\n"
<< "  </AXIS>" << "\n"
<<  "  <AXIS name=\"Theta\">" << "\n"
<<  "    <DETAIL name=\"min\" unit=\"rad\">0.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"max\" unit=\"rad\">1.</DETAIL>" << "\n"
<<  "    <DETAIL name=\"pts\">1</DETAIL>" << "\n"
<<  "  </AXIS>" << "\n"
<< "  <DATA points=\""<< pow(n,2) <<"\"> " << "\n";
 for (a = 0; a < n; a++){
  for (b = 0; b < n; b++){
      fon <<   Pinew[(n*a+b)] << endl;
  }
 }  
}
fon << "  </DATA>" << "\n"
<< "</VECTOR3D>" << "\n"
<< "</HYPERSURFACE>" << "\n";
// finish file .xml
sec1 = clock();
return 0;
} // end program
