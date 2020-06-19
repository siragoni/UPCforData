#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TDatime.h"


#include "TH2.h"

//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void PolarisationIntegrationRangeBiasHe1D( ){

  TDatime d;

  TFile* fileDataRawCosTheta[6] = {0,0,0,0,0,0};
  TFile* fileDataRawPhi[6]      = {0,0,0,0,0,0};
  TFile* fileDataRawTildePhi[6] = {0,0,0,0,0,0};

  fileDataRawCosTheta[0] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi[0]      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi[0] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );

  fileDataRawCosTheta[1] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_1.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi[1]      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_1.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi[1] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_1.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );

  fileDataRawCosTheta[2] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi[2]      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_2.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi[2] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );

  fileDataRawCosTheta[3] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi[3]      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_3.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi[3] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );

  fileDataRawCosTheta[4] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi[4]      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_4.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi[4] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );

  fileDataRawCosTheta[5] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi[5]      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_5.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi[5] = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );

  TH1F* CosThetaAfterSignalExtractionErrorsRawH[6] = {0,0,0,0,0,0};
  TH1F* PhiAfterSignalExtractionErrorsRawH[6]      = {0,0,0,0,0,0};
  TH1F* TildePhiAfterSignalExtractionErrorsRawH[6] = {0,0,0,0,0,0};

  for ( Int_t iLoop = 0; iLoop < 6; iLoop++ ) {

    CosThetaAfterSignalExtractionErrorsRawH[iLoop] = (TH1F*)fileDataRawCosTheta[iLoop]->Get("CosThetaAfterSignalExtractionErrorsH");
    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->Sumw2();

    PhiAfterSignalExtractionErrorsRawH[iLoop] = (TH1F*)fileDataRawPhi[iLoop]->Get("PhiAfterSignalExtractionErrorsH");
    PhiAfterSignalExtractionErrorsRawH[iLoop]->Sumw2();

    TildePhiAfterSignalExtractionErrorsRawH[iLoop] = (TH1F*)fileDataRawTildePhi[iLoop]->Get("TildePhiAfterSignalExtractionErrorsH");
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->Sumw2();

  }

  for ( Int_t iLoop = 1; iLoop < 6; iLoop++ ) {

    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->Divide( CosThetaAfterSignalExtractionErrorsRawH[0] );
    PhiAfterSignalExtractionErrorsRawH[iLoop]     ->Divide( PhiAfterSignalExtractionErrorsRawH[0]      );
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->Divide( TildePhiAfterSignalExtractionErrorsRawH[0] );

  }

  TH1F* AllOneCosTheta = new TH1F("AllOneCosTheta", "AllOneCosTheta", 25, -1.,   1.  );
  TH1F* AllOnePhi      = new TH1F("AllOnePhi",      "AllOnePhi",      25, -3.14, 3.14);
  TH1F* AllOneTildePhi = new TH1F("AllOneTildePhi", "AllOneTildePhi", 25,  0.,   6.28);
  for ( Int_t iBin = 1; iBin <= 25; iBin++ ) {
    AllOneCosTheta->SetBinContent(iBin,1.);
    AllOnePhi     ->SetBinContent(iBin,1.);
    AllOneTildePhi->SetBinContent(iBin,1.);
    AllOneCosTheta->SetBinError(iBin,0.0000000001);
    AllOnePhi     ->SetBinError(iBin,0.0000000001);
    AllOneTildePhi->SetBinError(iBin,0.0000000001);
  }


  for ( Int_t iLoop = 1; iLoop < 6; iLoop++ ) {
    for ( Int_t iBin = 1; iBin <= 25; iBin++ ) {
      CosThetaAfterSignalExtractionErrorsRawH[iLoop]->SetBinError(iBin,0);
      PhiAfterSignalExtractionErrorsRawH[iLoop]     ->SetBinError(iBin,0);
      TildePhiAfterSignalExtractionErrorsRawH[iLoop]->SetBinError(iBin,0);
    }
  }



  // for ( Int_t iLoop = 1; iLoop < 6; iLoop++ ) {
  //
  //   CosThetaAfterSignalExtractionErrorsRawH[iLoop]->Add( AllOneCosTheta );
  //   PhiAfterSignalExtractionErrorsRawH[iLoop]     ->Add( AllOnePhi      );
  //   TildePhiAfterSignalExtractionErrorsRawH[iLoop]->Add( AllOneTildePhi );
  //
  //   CosThetaAfterSignalExtractionErrorsRawH[iLoop]->Scale( -1. );
  //   PhiAfterSignalExtractionErrorsRawH[iLoop]     ->Scale( -1. );
  //   TildePhiAfterSignalExtractionErrorsRawH[iLoop]->Scale( -1. );
  //
  // }


  TCanvas* CosThetaCanvas = new TCanvas("CosThetaCanvas","CosThetaCanvas",900,800);


  for ( Int_t iLoop = 1; iLoop < 6; iLoop++ ) {

    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerStyle(20);
    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerColor(iLoop);
    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->SetLineColor(iLoop);
    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerSize(3);
    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->SetLineWidth(2);
    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->GetYaxis()->SetRangeUser(0.5, 1.1);


    PhiAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerStyle(20);
    PhiAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerColor(iLoop);
    PhiAfterSignalExtractionErrorsRawH[iLoop]->SetLineColor(iLoop);
    PhiAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerSize(3);
    PhiAfterSignalExtractionErrorsRawH[iLoop]->SetLineWidth(2);
    PhiAfterSignalExtractionErrorsRawH[iLoop]->GetYaxis()->SetRangeUser(0.5, 1.1);

    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerStyle(20);
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerColor(iLoop);
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->SetLineColor(iLoop);
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->SetMarkerSize(3);
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->SetLineWidth(2);
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->GetYaxis()->SetRangeUser(0.5, 1.1);

  }



  for ( Int_t iLoop = 1; iLoop < 6; iLoop++ ) {
    CosThetaAfterSignalExtractionErrorsRawH[iLoop]->Draw("epsame");
  }

  // gPad->BuildLegend();

  TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
  l->SetMargin(0.1);
  l->SetBorderSize(0);
  l->AddEntry(  CosThetaAfterSignalExtractionErrorsRawH[1], "2.85<M<3.35");
  l->AddEntry(  CosThetaAfterSignalExtractionErrorsRawH[2], "2.80<M<3.35");
  l->AddEntry(  CosThetaAfterSignalExtractionErrorsRawH[3], "2.90<M<3.35");
  l->AddEntry(  CosThetaAfterSignalExtractionErrorsRawH[4], "2.85<M<3.40");
  l->AddEntry(  CosThetaAfterSignalExtractionErrorsRawH[5], "2.85<M<3.30");
  l->Draw("same");



  TCanvas* PhiCanvas = new TCanvas("PhiCanvas","PhiCanvas",900,800);

  for ( Int_t iLoop = 1; iLoop < 6; iLoop++ ) {
    PhiAfterSignalExtractionErrorsRawH[iLoop]->Draw("epsame");
  }

  // gPad->BuildLegend();

  l->Draw("same");


  TCanvas* TildePhiCanvas = new TCanvas("TildePhiCanvas","TildePhiCanvas",900,800);

  for ( Int_t iLoop = 1; iLoop < 6; iLoop++ ) {
    TildePhiAfterSignalExtractionErrorsRawH[iLoop]->Draw("epsame");
  }

  // gPad->BuildLegend();

  l->Draw("same");

}
