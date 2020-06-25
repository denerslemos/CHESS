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

#include <fstream>
#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include "events2chain.C"
#include "model2legend.C"
#include "hist2xml.C"

#define  _FIGURE_NAME_ "fig_v22eta"
#define  _N_HISTOGRAMS_ 2

void figure_v22etad(TString aEventDir = "./events/", Int_t aEventFiles = 1)
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
  Int_t   XBins  = 10;
  Float_t XMin   = -5.0;
  Float_t XMax   =  5.0;
  Float_t Pmax   =  2.0;
  Float_t dX     = (XMax - XMin) / XBins;
  TH1D*   H1D[2*_N_HISTOGRAMS_];
// Create histograms
  H1D[0] = new TH1D("H0", "charged", XBins, XMin, XMax);
  H1D[0]->GetXaxis()->SetTitle("#eta");
  H1D[0]->GetYaxis()->SetTitle("v_{2}");
  H1D[0]->Sumw2();
  
  H1D[1] = (TH1D*) H1D[0]->Clone("H1");  H1D[1]->SetTitle("charged");

  H1D[2] = new TH1D("H2", "charged", XBins, XMin, XMax);
  H1D[2]->GetXaxis()->SetTitle("#eta");
  H1D[2]->GetYaxis()->SetTitle("s_{2}");
  H1D[2]->Sumw2();
  
  H1D[3] = (TH1D*) H1D[2]->Clone("H2");  H1D[3]->SetTitle("charged sin");

// Fill histograms
  Float_t Eta, Pt, v2, s2, P;
  Int_t   pid;
  
  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Chain->GetEntry(i);
       P   = TMath::Sqrt(Particle.px*Particle.px + Particle.py*Particle.py + Particle.pz*Particle.pz);
    if(P == Particle.pz)
      continue;
    Pt  = TMath::Sqrt(Particle.px*Particle.px + Particle.py*Particle.py);
    if(Pt == 0.0)
      continue;
    v2 = (Particle.px*Particle.px - Particle.py*Particle.py) / (Pt*Pt);
    s2 = (2.*Particle.px*Particle.py) / (Pt*Pt);
    Eta = 0.5 * TMath::Log((P+Particle.pz) / (P-Particle.pz));
    pid = Particle.pid;
    if(((pid == 211) || (pid == -211) || (pid == -321) || (pid == 321) || (pid == -2212) || (pid == 2212))&&(Pt<Pmax)) {
      H1D[0]->Fill(Eta, v2);
      H1D[1]->Fill(Eta, 1.0);
      H1D[2]->Fill(Eta, s2);
      H1D[3]->Fill(Eta, 1.0);
    }
  }
// Rescale histograms
  for(Int_t i=0; i<2*_N_HISTOGRAMS_; i++)
    H1D[i]->Scale(1.0 / (Events * 2.0*TMath::Pi() * dX ));
  H1D[0]->Divide(H1D[1]);
  H1D[2]->Divide(H1D[3]);

  H1D[0]->Multiply(H1D[0]);
  H1D[2]->Multiply(H1D[2]);
  H1D[0]->Add(H1D[2]);


// ##########################################################################
// # SAVE HISTOGRAMS TO XML FILE
// ##########################################################################
  hist2xml(aEventDir + _FIGURE_NAME_ + ".xml", H1D, _N_HISTOGRAMS_, Events, Chain->GetEntries());

  ofstream wout;
  wout.open(aEventDir + _FIGURE_NAME_ + ".dat");


  for (int i = 1;i< =XBins; i++){
    cout<<H1D[0]->GetBinCenter(i)<<" "<<H1D[0]->GetBinContent(i)<<endl;
    wout<<H1D[0]->GetBinCenter(i)<<" "<<H1D[0]->GetBinContent(i)<<" "<<H1D[0]->GetBinError(i)<<endl;
  }
  wout.close();

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
// Histogram legend

// Model legend
//  TLegend* LegendModel	= new TLegend(0.70, 0.40, 0.98, 0.92, "model",	"brNDC");
//  LegendModel->SetFillColor(0);
//  LegendModel->SetFillStyle(4000);
//  LegendModel->SetTextSize(0.035);
//  model2legend(aEventDir, aEventFiles, LegendModel);
// Histograms  
  H1D[0]->GetXaxis()->SetTitleSize(0.06);
  H1D[0]->GetXaxis()->CenterTitle(kTRUE);
  H1D[0]->GetXaxis()->SetLabelSize(0.05); 
  H1D[0]->GetYaxis()->SetTitleSize(0.06);
  H1D[0]->GetYaxis()->CenterTitle(kTRUE); 
  H1D[0]->GetYaxis()->SetLabelSize(0.05);
  H1D[0]->SetMinimum(0.00);
  tMaxBin = H1D[0]->GetMaximumBin();	H1D[0]->SetMaximum((H1D[0]->GetBinContent(tMaxBin) + H1D[0]->GetBinError(tMaxBin)) * 1.05);
  H1D[0]->SetMarkerColor(2);	H1D[0]->SetMarkerStyle(20);
 
// Plot
  H1D[0]->SetTitle("v_{2}(#eta)");
  H1D[0]->Draw();  
  //  LegendPart ->Draw();
//  LegendModel->Draw();
// Save to files

  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".eps");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".png");
    Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".C");
}
