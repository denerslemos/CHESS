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

#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include "events2chain.C" 
#include "model2legend.C"
#include "hist2xml.C"

#define  _FIGURE_NAME_ "fig_v1chpt2"
#define  _N_HISTOGRAMS_ 3

void figure_v1chptd2(TString aEventDir = "./events/", Int_t aEventFiles = 1)
{
// ##########################################################################
// # READ ROOT FILES
// ##########################################################################
  static ParticleCoor Particle;
  Int_t   Events;
  TChain* Chain = events2chain(aEventDir, aEventFiles, &Particle, &Events);

// ##########################################################################
// HISTOGRAMS
// ##########################################################################
  Int_t   XBins  = 20;
  Float_t XMin   = 0.0;
  Float_t XMax   = 5.0;
  Float_t shift   = -0.0;
  Float_t dX     = (XMax - XMin) / XBins;
  TH1D*   H1D[2*_N_HISTOGRAMS_];
// Create histograms
  H1D[0] = new TH1D("H0", "#pi^{+}", XBins, XMin, XMax);
  H1D[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  H1D[0]->GetYaxis()->SetTitle("v_{1}");
  H1D[0]->Sumw2();
  H1D[1] = (TH1D*) H1D[0]->Clone("H1");  H1D[1]->SetTitle("K^{+}");
  H1D[2] = (TH1D*) H1D[0]->Clone("H2");  H1D[2]->SetTitle("p");
//Create norm histograms
  H1D[3] = (TH1D*) H1D[0]->Clone("N1");  H1D[3]->SetTitle("K^{+}");
  H1D[4] = (TH1D*) H1D[0]->Clone("N1");  H1D[4]->SetTitle("K^{+}");
  H1D[5] = (TH1D*) H1D[0]->Clone("N2");  H1D[5]->SetTitle("p");
// Fill histograms
// Fill histograms
  Float_t Rap, Pt, v2, s2, vv2=0., ss2=0., count=0., angle;
  Int_t   pid;
  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Double_t fang;
    Chain->GetEntry(i);
    if(Particle.e == Particle.pz)
      continue;
    Rap = 0.5 * TMath::Log((Particle.e + Particle.pz) / (Particle.e - Particle.pz));
    if( TMath::Abs(Rap-shift) >= 1.0 )
      continue;
    Pt  = TMath::Sqrt(Particle.px*Particle.px + Particle.py*Particle.py);
    if(Pt == 0.0)
      continue;
    fang=atan2(Particle.py,Particle.px);
    v2=cos(1.*fang);
    s2=sin(1.*fang);
    pid = Particle.pid;

    if((pid == 211) || (pid == -211) || (pid == -321) || (pid == 321) || (pid == -2212) || (pid == 2212)) {
      ss2+=s2;
      vv2+=v2;
      count++;
    } else {
      if((pid==-10321)||(pid==10321)||(pid==-3324)||(pid==3324)||(pid==-323)||(pid==323)||(pid==-3114)||(pid==3114)||(pid==-3224)||(pid==3224)||(pid==-10323)||(pid==10323)||(pid==-3112)||(pid==3112)||(pid==-67719)||(pid==67719)||(pid==-3334)||(pid==3334)||(pid==13214)||(pid==13114)||(pid==-46653)||(pid==46653)||(pid==-3312)||(pid==3312)||(pid==-3222)||(pid==3222)||(pid==-9000211)||(pid==9000211))
      ss2+=s2;
      vv2+=v2;
      count++;

    }    
  }
 
    cout<< " averages" <<ss2/count<<" "<<vv2/count<<endl;
    angle=atan2(ss2,vv2)/1.;
    cout<<angle*180./3.14159265<<endl;

  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Chain->GetEntry(i);
    if(Particle.e == Particle.pz)
      continue;
    Rap = 0.5 * TMath::Log((Particle.e + Particle.pz) / (Particle.e - Particle.pz));
    if( TMath::Abs(Rap-shift) >= 1.0 )
      continue;
    Pt  = TMath::Sqrt(Particle.px*Particle.px + Particle.py*Particle.py);
    if(Pt == 0.0)
      continue;
    v2=atan2(Particle.py,Particle.px);
    v2=v2-angle;
    v2=cos(1.*v2);
    pid = Particle.pid;
     if((pid == 211) || (pid == -211) || (pid == -321) || (pid == 321) || (pid == -2212) || (pid == 2212)) {
      H1D[0]->Fill(Pt, 1.0 * v2);
      H1D[3]->Fill(Pt, 1.0);
    } else {
      if((pid==-10321)||(pid==10321)||(pid==-3324)||(pid==3324)||(pid==-323)||(pid==323)||(pid==-3114)||(pid==3114)||(pid==-3224)||(pid==3224)||(pid==-10323)||(pid==10323)||(pid==-3112)||(pid==3112)||(pid==-67719)||(pid==67719)||(pid==-3334)||(pid==3334)||(pid==13214)||(pid==13114)||(pid==-46653)||(pid==46653)||(pid==-3312)||(pid==3312)||(pid==-3222)||(pid==3222)||(pid==-9000211)||(pid==9000211))
      H1D[0]->Fill(Pt, 1.0 * v2);
  //    H1D[3]->Fill(Pt, 1.0);
    }

  }
// Rescale histograms
  for(Int_t i=0; i<2*_N_HISTOGRAMS_; i++)
    H1D[i]->Scale(1.0 / (Events * 2*1.0 * 2.0*TMath::Pi() * dX ));
  H1D[0]->Divide(H1D[3]);
  //  H1D[1]->Divide(H1D[4]);
  // H1D[2]->Divide(H1D[5]);

// ##########################################################################
// # SAVE HISTOGRAMS TO XML FILE
// ##########################################################################
  hist2xml(aEventDir + _FIGURE_NAME_ + ".xml", H1D, _N_HISTOGRAMS_, Events, Chain->GetEntries());

// ##########################################################################
// # PLOT HISTOGRAMS
// ##########################################################################
  Int_t tMinBin, tMaxBin;
  gStyle->SetOptStat("");
// Canvas
  TCanvas* Canvas	= new TCanvas("canvas", H1D[0]->GetYaxis()->GetTitle(), 800, 600);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetMargin(0.14, 0.02, 0.14, 0.08);




 ofstream wout;
  wout.open(aEventDir + _FIGURE_NAME_ + ".dat");


  for (int i = 1;i< =XBins; i++){
    cout<<H1D[0]->GetBinCenter(i)<<" "<<H1D[0]->GetBinContent(i)<<endl;
    wout<<H1D[0]->GetBinCenter(i)<<" "<<H1D[0]->GetBinContent(i)<<" "<<H1D[0]->GetBinError(i)<<endl;
  }
  wout.close();


// Histograms  
  H1D[0]->GetXaxis()->SetTitleSize(0.06);
  H1D[0]->GetXaxis()->CenterTitle(kTRUE);
  H1D[0]->GetXaxis()->SetLabelSize(0.05); 
  H1D[0]->GetYaxis()->SetTitleSize(0.06);
  H1D[0]->GetYaxis()->CenterTitle(kTRUE); 
  H1D[0]->GetYaxis()->SetLabelSize(0.05);
  H1D[0]->SetMinimum(0.0);
  tMaxBin = H1D[0]->GetMaximumBin();	H1D[0]->SetMaximum((H1D[0]->GetBinContent(tMaxBin) + H1D[0]->GetBinError(tMaxBin)) * 1.05);
  H1D[0]->SetMarkerColor(2);	H1D[0]->SetMarkerStyle(20);
  H1D[1]->SetMarkerColor(4);	H1D[1]->SetMarkerStyle(21);
  H1D[2]->SetMarkerColor(1);	H1D[2]->SetMarkerStyle(22);
// Plot
  H1D[0]->SetTitle("v_{1}(p_{T})");
  H1D[0]->Draw();
  H1D[1]->Draw("SAME");
  H1D[2]->Draw("SAME");
 
// Save to files
  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".eps");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".png");
//    Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".C");
}
