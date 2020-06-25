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

#define  _FIGURE_NAME_ "fig_distpt"
#define  _N_HISTOGRAMS_ 3

void figure_distpt(TString aEventDir = "./events/", Int_t aEventFiles = 1)
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
  Int_t   XBins  = 30;
  Float_t XMin   = 0.0;
  Float_t XMax   = 3.0;
  Float_t dX     = (XMax - XMin) / XBins;
  TH1D*   H1D[_N_HISTOGRAMS_];
// Create histograms
  H1D[0] = new TH1D("H0", "#pi^{+}", XBins, XMin, XMax);
  H1D[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  H1D[0]->GetYaxis()->SetTitle("dN/(2 #pi p_{T} dp_{T} dy)");
  H1D[0]->Sumw2();
  H1D[1] = (TH1D*) H1D[0]->Clone("H1");  H1D[1]->SetTitle("h^{+}");
  H1D[2] = (TH1D*) H1D[0]->Clone("H2");  H1D[2]->SetTitle("h^{#pm}");
// Fill histograms
  Float_t P, Rap, Pt;
  Int_t   pid;
  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Chain->GetEntry(i);
    if(Particle.e == Particle.pz)
      continue;
    Rap = 0.5 * TMath::Log((Particle.e+Particle.pz) / (Particle.e-Particle.pz));
    if( Rap > 1.4  || Rap < 0.2)      continue;
    pid = Particle.pid;
    Pt  = TMath::Sqrt(Particle.px*Particle.px + Particle.py*Particle.py);
    if((pid == 211) || (pid == -211) || (pid == -321) || (pid == 321) || (pid == -2212) || (pid == 2212)) {
      H1D[0]->Fill(Pt, 1.0/Pt);
      if(pid == 211)
	H1D[1]->Fill(Pt, 1.0/Pt);
      if(pid == 321)
	H1D[2]->Fill(Pt, 1.0/Pt);
    } else {
      if((pid==-10321)||(pid==10321)||(pid==-3324)||(pid==3324)||(pid==-323)||(pid==323)||(pid==-3114)||(pid==3114)||(pid==-3224)||(pid==3224)||(pid==-10323)||(pid==10323)||(pid==-3112)||(pid==3112)||(pid==-67719)||(pid==67719)||(pid==-3334)||(pid==3334)||(pid==13214)||(pid==13114)||(pid==-46653)||(pid==46653)||(pid==-3312)||(pid==3312)||(pid==-3222)||(pid==3222)||(pid==-9000211)||(pid==9000211))
        H1D[0]->Fill(Pt, 1.0/Pt);
    }
/*
    if((pid == 211)||(pid == 321)||(pid == 2212)||(pid==10321)||(pid==3324)||(pid==323)||(pid==3114)||(pid==3224)||(pid==10323)||(pid==3112)||(pid==67719)||(pid==3334)||(pid==13214)||(pid==13114)||(pid==46653)||(pid==3312)||(pid==3222)||(pid==9000211)||(pid == -211)||(pid == -321)||(pid == -2212)||(pid==-10321)||(pid==-3324)||(pid==-323)||(pid==-3114)||(pid==-3224)||(pid==-10323)||(pid==-3112)||(pid==-67719)||(pid==-3334)||(pid==-46653)||(pid==-3312)||(pid==-3222)||(pid==-9000211))
      H1D[0]->Fill(Pt, 1.0/Pt);
    if((pid == 211)||(pid == 321)||(pid == 2212)||(pid==10321)||(pid==3324)||(pid==323)||(pid==3114)||(pid==3224)||(pid==10323)||(pid==3112)||(pid==67719)||(pid==3334)||(pid==13214)||(pid==13114)||(pid==46653)||(pid==3312)||(pid==3222)||(pid==9000211))
      H1D[1]->Fill(Pt, 1.0/Pt);
    if((pid == -211)||(pid == -321)||(pid == -2212)||(pid==-10321)||(pid==-3324)||(pid==-323)||(pid==-3114)||(pid==-3224)||(pid==-10323)||(pid==-3112)||(pid==-67719)||(pid==-3334)||(pid==-46653)||(pid==-3312)||(pid==-3222)||(pid==-9000211)) // protons without the ones form (Lm1115zer)
      H1D[2]->Fill(Pt, 1.0/Pt);
*/
  }
// Rescale histograms
// Rescale histograms
  for(Int_t i=0; i<_N_HISTOGRAMS_; i++)
    H1D[i]->Scale(1.0 / (Events * 2.0 * 1.0 * 2.0*TMath::Pi() * dX ));

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


  for (int i = 1;i<=XBins; i++){
//    cout<<H1D[0]->GetBinCenter(i)<<" "<<H1D[0]->GetBinContent(i) << endl;
    wout<<H1D[0]->GetBinCenter(i)<<" "<<H1D[0]->GetBinContent(i)/2.0<<" "<< H1D[0]->GetBinError(i)<<endl;
  }
  wout.close();


  gPad->SetLogy();
// Histogram legend
  TLegend* LegendPart	= new TLegend(0.17, 0.18, 0.34, 0.37, "",	"brNDC");
  LegendPart->SetFillColor(0);
  LegendPart->SetFillStyle(4000);
  LegendPart->SetTextSize(0.05);
  for(Int_t i=0; i<_N_HISTOGRAMS_; i++)
    LegendPart->AddEntry(H1D[i], H1D[i]->GetTitle());
// Model legend
  TLegend* LegendModel	= new TLegend(0.70, 0.40, 0.98, 0.92, "model",	"brNDC");
  LegendModel->SetFillColor(0);
  LegendModel->SetFillStyle(4000);
  LegendModel->SetTextSize(0.035);
  model2legend(aEventDir, aEventFiles, LegendModel);
// Histograms  
  H1D[0]->GetXaxis()->SetTitleSize(0.06);
  H1D[0]->GetXaxis()->CenterTitle(kTRUE);
  H1D[0]->GetXaxis()->SetLabelSize(0.05); 
  H1D[0]->GetYaxis()->SetTitleSize(0.06);
  H1D[0]->GetYaxis()->CenterTitle(kTRUE); 
  H1D[0]->GetYaxis()->SetLabelSize(0.05);
//  tMinBin = H1D[0]->GetMinimumBin();	H1D[0]->SetMinimum((H1D[0]->GetBinContent(tMinBin) - H1D[0]->GetBinError(tMinBin)) * 0.5);
  tMaxBin = H1D[0]->GetMaximumBin();	H1D[0]->SetMaximum((H1D[0]->GetBinContent(tMaxBin) + H1D[0]->GetBinError(tMaxBin)) * 2.0);
  H1D[0]->SetMarkerColor(2);	H1D[0]->SetMarkerStyle(20);
  H1D[1]->SetMarkerColor(4);	H1D[1]->SetMarkerStyle(21);
  H1D[2]->SetMarkerColor(1);	H1D[2]->SetMarkerStyle(22);
// Plot
  H1D[0]->SetTitle("m_{T} distribution");
  H1D[0]->Draw();
  H1D[1]->Draw("SAME");
  H1D[2]->Draw("SAME");
  LegendPart ->Draw();
  LegendModel->Draw();
// Save to files
  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".eps");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".png");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".C");
}
