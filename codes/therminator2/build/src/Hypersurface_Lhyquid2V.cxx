/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include <TMath.h>
#include "Configurator.h"
#include "Parser.h"
#include "Hypersurface_Library.h"
#include "Hypersurface_Lhyquid2V.h"
#include "THGlobal.h"

using namespace std;
using namespace TMath;


extern void AddLogEntry(const char* aEntry);
extern Configurator* sMainConfig;
extern TString	sModelINI;
extern TString	sHyperXML;
extern TString	sEventDIR;

Hypersurface_Lhyquid2V::Hypersurface_Lhyquid2V()
{
  mThermo = new Thermodynamics();
  ReadParameters();
}

Hypersurface_Lhyquid2V::~Hypersurface_Lhyquid2V()
{
  delete mThermo;
  delete mDistance;
  delete mDistanceDZeta;
  delete mDistanceDPhi;
  delete mFluidVt;
  delete mFluidPhi;

// ******************Dener**********************

/*
  delete mFluidUx;
  delete mFluidUy;
*/
  delete mpixx;
  delete mpixy;
  delete mpiyy;
  delete mpibulk;
//  delete mpixe;
//  delete mpiye;
//  delete mpiee;
//  delete mpitx;
//  delete mpity;
//  delete mpite;
//  delete mpitt;
//*************************************

}

double Hypersurface_Lhyquid2V::GetDSigmaP(double aMt, double aPt, double aPhiP, double aRapP)
{
  return Dhs * (mTauI + Dhs * Sin(Zeta)) * (
         Dhs     * Cos(Zeta) * ( aMt * Sin(Zeta) * CosH(aRapP - RapS) + aPt * Cos(Zeta) * Cos(PhiS - aPhiP))
       + dDdZeta * Cos(Zeta) * (-aMt * Cos(Zeta) * CosH(aRapP - RapS) + aPt * Sin(Zeta) * Cos(PhiS - aPhiP))
       + dDdPhi  * aPt * Sin(PhiS - aPhiP)
    );
}

double Hypersurface_Lhyquid2V::GetPdotU(double aMt, double aPt, double aPhiP, double aRapP)
{
  return 1.0 / Sqrt(1 - Vt * Vt) * (aMt * CosH(aRapP - RapS) - Vt * aPt * Cos(PhiS - aPhiP + PhiF));
}

/*
double Hypersurface_Lhyquid2V::GetU0()
{
  return Sqrt(1 + Ux * Ux + Uy * Uy);
}

double Hypersurface_Lhyquid2V::GetUx()
{
  return Ux;
}

double Hypersurface_Lhyquid2V::GetUy()
{
  return Uy;
}

double Hypersurface_Lhyquid2V::GetUz()
{
  return 0.0;
}

*/

double Hypersurface_Lhyquid2V::GetStress(double aMt, double aPx, double aPy,
					 double aRapP)
{
 // double RapSS =0.0;
// Dener Lemos
// Corrections in f(x,p)
//  double Peta =  aMt*SinH(aRapP - RapS);
//  double Ptau =  aMt*CosH(aRapP - RapS); // Done!

double Peta =  aMt*SinH(RapS - aRapP);
double Ptau =  aMt*CosH(RapS - aRapP);

double stresspp = spitt*Ptau*Ptau + spiee*Peta*Peta + spixx*aPx*aPx + spiyy*aPy*aPy + 2.*spixy*aPx*aPy - 2*Ptau*(aPx*spitx + aPy*spity);
//double stresspp =  spitt*Ptau*Ptau + spixx*aPx*aPx + spiyy*aPy*aPy + spiee*Peta*Peta +
//       2.*spixy*aPx*aPy + 2.*spixe*aPx*Peta + 2.*spiye*aPy*Peta
//      - 2.*spitx*Ptau*aPx - 2.*spity*Ptau*aPy - 2.*spite*Ptau*Peta;

  return stresspp*.5 / mThermo->GetTemperature()  / mThermo->GetTemperature();
}

double Hypersurface_Lhyquid2V::GetBulk(double Estar, double Pstar)
{
  return (Pstar*Pstar/3./Estar-mcs2freeze*Estar)*spibulk*mfacbulkfreeze;
}

//*********************************************

double Hypersurface_Lhyquid2V::GetHyperCubeSpatialVolume()
{
  return (mDistance->GetXMax() - mDistance->GetXMin()) * (mDistance->GetYMax() - mDistance->GetYMin()) * mRapSRange;
}

void Hypersurface_Lhyquid2V::SetPositionOnHypersurface(TRandom2* aRandom, double* aXt, double* aXx, double* aXy, double* aXz)
{
  Zeta    = mDistance->GetXMin() + (mDistance->GetXMax() - mDistance->GetXMin()) * aRandom->Rndm();
  PhiS    = mDistance->GetYMin() + (mDistance->GetYMax() - mDistance->GetYMin()) * aRandom->Rndm();
  RapS    = mRapSRange * aRandom->Rndm() - 0.5 * mRapSRange;
  Dhs     = mDistance     ->Interpolate(Zeta, PhiS, 0.0);
  dDdPhi  = mDistanceDPhi ->Interpolate(Zeta, PhiS, 0.0);
  dDdZeta = mDistanceDZeta->Interpolate(Zeta, PhiS, 0.0);
//  Vx	  = mFluidVx      ->Interpolate(Zeta, PhiS, 0.0);
  Vt	  = mFluidVt      ->Interpolate(Zeta, PhiS, 0.0);
  PhiF	  = mFluidPhi     ->Interpolate(Zeta, PhiS, 0.0);
  Tau     = mTauI + Dhs * Sin(Zeta);
  Rho     = Dhs * Cos(Zeta);
  (*aXt)  = Tau * CosH(RapS);
  (*aXx)  = Rho * Cos(PhiS);
  (*aXy)  = Rho * Sin(PhiS);
  (*aXz)  = Tau * SinH(RapS);  
  Vx = Vt * Cos(PhiS + PhiF);
  Vy = Vt * Sin(PhiS + PhiF);
  spixx   = mpixx          ->Interpolate(Zeta, PhiS, 0.0);
  spixy   = mpixy          ->Interpolate(Zeta, PhiS, 0.0);
  spiyy   = mpiyy          ->Interpolate(Zeta, PhiS, 0.0);
  spitx   = Vx * spixx + Vy * spixy;
  spity   = Vx * spixy + Vy * spiyy;
  spitt   = Vx * spitx + Vy * spity;
  spiee   = spitt - spixx - spiyy;
  spibulk = mpibulk        ->Interpolate(Zeta, PhiS, 0.0);  



// Dener ***********************************************
 /* Ux	  = mFluidUx       ->Interpolate(Zeta, PhiS, 0.0);
  Uy	  = mFluidUy       ->Interpolate(Zeta, PhiS, 0.0);

  Ueta= SinH(-RapS);
  Utau= sqrt(1.+Ux*Ux+Uy*Uy+Ueta*Ueta);
*/
 
//  spixe   = mpixe          ->Interpolate(Zeta, PhiS, 0.0);
//  spiye   = mpiye          ->Interpolate(Zeta, PhiS, 0.0);
//  spiee   = mpiee          ->Interpolate(Zeta, PhiS, 0.0);
//  spitx   = mpitx          ->Interpolate(Zeta, PhiS, 0.0);
//  spity   = mpity          ->Interpolate(Zeta, PhiS, 0.0);
//  spite   = mpite          ->Interpolate(Zeta, PhiS, 0.0);
//  spitt   = mpitt          ->Interpolate(Zeta, PhiS, 0.0);

}

void Hypersurface_Lhyquid2V::ReadParameters()
{
  Hypersurface_Library*	tLib;
  Configurator*		tModelParam;
  Parser*		tParser;
  
  tModelParam = new Configurator;
  tParser     = new Parser(sModelINI.Data());
  tParser->ReadINI(tModelParam);
  delete tParser;
    
  try {
    mRapPRange   = tModelParam->GetParameter("RapPRange").Atof();						// [1]
    mRapSRange   = tModelParam->GetParameter("RapSRange").Atof();						// [1]
    if(sHyperXML.IsNull()) {
      sHyperXML = sMainConfig->GetParameter("FreezeOutDir");
      sHyperXML.Prepend("./");
      sHyperXML += tModelParam->GetParameter("FreezeFile");
    }
  } catch (TString tError) {
    PRINT_MESSAGE("<Hypersurface_Lhyquid2V::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  
  char tBuff[2*kFileNameMaxChar];
  sprintf(tBuff,"[input]\t%s",sHyperXML.Data()); 
  AddLogEntry(tBuff);
//###############################################
//		XML file with hypersurface	#
//###############################################

  PRINT_DEBUG_1("Hydro based model input. FreezeFile: " << sHyperXML.Data());
  tLib    = new Hypersurface_Library;
  tParser = new Parser(sHyperXML.Data());
  tParser->ReadXML(tLib);
  delete tParser;
  
  try {
    mTauI	     = tLib->GetXMLTag("PARAMETER","name","Tau_i")->GetXMLContent().Atof() / kHbarC;		// [GeV^-1]
    mThermo->SetTemperature(tLib->GetXMLTag("PARAMETER","name","Temperature")->GetXMLContent().Atof() * 0.001);	// [GeV]
    mThermo->SetChemistry(  tLib->GetXMLTag("PARAMETER","name","Mu_B")->GetXMLContent().Atof() * 0.001,
			    tLib->GetXMLTag("PARAMETER","name","Mu_I")->GetXMLContent().Atof() * 0.001,
			    tLib->GetXMLTag("PARAMETER","name","Mu_S")->GetXMLContent().Atof() * 0.001,
			    tLib->GetXMLTag("PARAMETER","name","Mu_C")->GetXMLContent().Atof() * 0.001);	// [GeV]
    mCollidingEnergy = tLib->GetXMLTag("PARAMETER","name","colliding_energy")->GetXMLContent().Atof();		// [GeV]
    mCentralityMin   = tLib->GetXMLTag("PARAMETER","name","centrality_min")->GetXMLContent().Atof();		// [%]
    mCentralityMax   = tLib->GetXMLTag("PARAMETER","name","centrality_max")->GetXMLContent().Atof();		// [%]
    mImpactParameter = tLib->GetXMLTag("PARAMETER","name","impact_parameter")->GetXMLContent().Atof();		// [fm]
    mTempI	     = tLib->GetXMLTag("PARAMETER","name","temperature_at_center")->GetXMLContent().Atof();	// [MeV]
    sprintf(mDeviceName,     "%s",tLib->GetXMLTag("PARAMETER","name","device")->GetXMLContent().Data());
    sprintf(mCollidingSystem,"%s",tLib->GetXMLTag("PARAMETER","name","colliding_system")->GetXMLContent().Data());
//Dener
    mcs2freeze = tLib->GetXMLTag("PARAMETER","name","cs^2")->GetXMLContent().Atof();		// []
    metas = tLib->GetXMLTag("PARAMETER","name","eta/s")->GetXMLContent().Atof();		// []
    mzetas = tLib->GetXMLTag("PARAMETER","name","zeta/s")->GetXMLContent().Atof();		// []
    mfacbulkfreeze = (22.76/ mThermo->GetTemperature() -379.8 +
		      2371.* mThermo->GetTemperature()
		      -4928.* mThermo->GetTemperature()
		      * mThermo->GetTemperature());
      // factor for bulk viscosity corrections for SHARE hadrons 
      // in relaxation time form
    mDistance	     = tLib->GetXMLTag("VECTOR3D", "name","Distance")->GetXMLVector3D();			// [GeV^-1]
    mFluidVt	     = tLib->GetXMLTag("VECTOR3D", "name","FluidVt" )->GetXMLVector3D();			// [c]
//    mFluidVy	     = tLib->GetXMLTag("VECTOR3D", "name","FluidVy" )->GetXMLVector3D();			// [c]
    mFluidPhi	     = tLib->GetXMLTag("VECTOR3D", "name","FluidPhi")->GetXMLVector3D();			// [rad]
//Dener
 /*   mFluidUx	     = tLib->GetXMLTag("VECTOR3D", "name","FluidUx" )->GetXMLVector3D();			// [1]
    mFluidUy	     = tLib->GetXMLTag("VECTOR3D", "name","FluidUy" )->GetXMLVector3D();			// [1]
*/  mpixx	     = tLib->GetXMLTag("VECTOR3D", "name","Pixx" )->GetXMLVector3D();		// [1]
    mpixy	     = tLib->GetXMLTag("VECTOR3D", "name","Pixy" )->GetXMLVector3D();		// [1]
//    mpixe	     = tLib->GetXMLTag("VECTOR3D", "name","Pixe" )->GetXMLVector3D();		// [1]
    mpiyy	     = tLib->GetXMLTag("VECTOR3D", "name","Piyy" )->GetXMLVector3D();		// [1]
//    mpiye	     = tLib->GetXMLTag("VECTOR3D", "name","Piye" )->GetXMLVector3D();		// [1]
//   mpiee	     = tLib->GetXMLTag("VECTOR3D", "name","Piee" )->GetXMLVector3D();		// [1]
//    mpitt	     = tLib->GetXMLTag("VECTOR3D", "name","Pitt" )->GetXMLVector3D();		// [1]
//    mpitx	     = tLib->GetXMLTag("VECTOR3D", "name","Pitx" )->GetXMLVector3D();		// [1]
//    mpity	     = tLib->GetXMLTag("VECTOR3D", "name","Pity" )->GetXMLVector3D();		// [1]
//    mpite	     = tLib->GetXMLTag("VECTOR3D", "name","Pite" )->GetXMLVector3D();		// [1]
    mpibulk	     = tLib->GetXMLTag("VECTOR3D", "name","PI" )->GetXMLVector3D();		// [1]
  } catch (int tError) {
    PRINT_MESSAGE("<Hypersurface_Lhyquid2V::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the XML file.");
    exit(_ERROR_LIBRARY_TAG_NOT_FOUND_);
  }

// mDistance derivatives
  try {
    mDistanceDZeta = tLib->GetXMLTag("VECTOR3D","name","DistanceDZeta")->GetXMLVector3D();			// [GeV^-1/rad]
  } catch (int tError) {
    PRINT_DEBUG_1("<Hypersurface_Lhyquid2V::ReadParameters>\tCalculating derivative Distance->DerivativeX()");
    mDistanceDZeta = mDistance->DerivativeX("DistanceDZeta");   
  }
  try {
    mDistanceDPhi = tLib->GetXMLTag("VECTOR3D","name","DistanceDPhi")->GetXMLVector3D();			// [GeV^-1/rad]
  } catch (int tError) {
    PRINT_DEBUG_1("<Hypersurface_Lhyquid2V::ReadParameters>\tCalculating derivative Distance->DerivativeY()");
    mDistanceDPhi = mDistance->DerivativeY("DistanceDPhi");   
  }
  
// event subdirectory
  try {
    sEventDIR += tModelParam->GetParameter("EventSubDir");
  } catch (TString tError) {
    TString tTemp = tModelParam->GetParameter("FreezeFile");
    tTemp.ReplaceAll("/","-"); // if you comment this the event directory will have the same structure as the XML files inside fomodel directory
    tTemp.ReplaceAll(".xml",TString(sHyperXML)+"/");// Dener Lemos    
    tTemp.ReplaceAll(".xml/","/");// Dener Lemos
    sEventDIR += tTemp;
  }

  delete tLib; 
//  cout<<" del 1"<<endl;
  delete tModelParam;
//  cout<<" del 2"<<endl;  
//  cout<<" \n"<<endl; 
//  cout<<" RUN!"<<endl; 
}
