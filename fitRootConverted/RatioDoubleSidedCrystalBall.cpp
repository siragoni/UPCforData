#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <vector>



//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void Ratio(){

  TDatime d;

  // TFile* file1D_HE      = new TFile(Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* file1D_HE_DSCB = new TFile("pngResults/CosThetaHeFrame_DSCB.root");
  TFile* file1D_HE      = new TFile(Form("pngResults/%d-%2.2d-%2.2d/CosThetaCS/CosThetaCsFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* file1D_HE_DSCB = new TFile("pngResults/CosThetaCsFrame_DSCB.root");

  TH1F* CosThetaHE_CB   = (TH1F*) file1D_HE     ->Get("CosThetaAfterSignalExtractionErrorsH");
  TH1F* CosThetaHE_DSCB = (TH1F*) file1D_HE_DSCB->Get("CosThetaAfterSignalExtractionErrorsH");

  TCanvas* RatioCanvas = new TCanvas("RatioCanvas", "RatioCanvas", 1000, 800);



  TH1F* RatioH = (TH1F*)CosThetaHE_CB->Clone("RatioH");
  RatioH->Sumw2();
  RatioH->Divide(CosThetaHE_DSCB);


  RatioH->SetMarkerStyle(21);
  RatioH->GetYaxis()->SetRangeUser(0.95, 1.05);
  RatioH->SetLineColor(kBlue);
  RatioH->SetLineWidth(2);
  RatioH->Draw("ep");
  // gPad->SaveAs("pngResults/RatioPtDistr.png", "RECREATE");



}
