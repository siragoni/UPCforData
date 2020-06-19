#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <math.h>
#include "TH2D.h"
#include "TF2.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"

#include "TDatime.h"

#include <vector>
#include <map>



//_____________________________________________________________________________
void AxEcheck(){

  TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);
  /* - We now do the same as before to ascertain if the TList was there and
   * - to try to retrieve the plots. Result:
   *   listings->ls()
   *     OBJ: TList	  MyOutputContainer	          Doubly linked list          : 0
   *     OBJ: TH1F	  fNumberMuonsH	              fNumberMuonsH               : 0 at: 0x5a145f0
   *     OBJ: TH1F	  fCounterH	                  fCounterH                   : 0 at: 0x5a3b570
   *     OBJ: TH1F	  fEtaMuonH	                  fEtaMuonH                   : 0 at: 0x5a3ba80
   *     OBJ: TH1F	  fRAbsMuonH	                fRAbsMuonH                  : 0 at: 0x5a3c0c0
   *     OBJ: TH1F	  fInvariantMassDistributionH	fInvariantMassDistributionH : 0 at: 0x5a3c720
   */
  // TH1F* fReconCosThetaH = (TH1F*)listings->FindObject("fCosThetaCsFrameTwentyfiveBinsH");
  // TH1F* fGenerCosThetaH = (TH1F*)listings->FindObject("fMCCosThetaCsFrameTwentyfiveBinsH");
  TH1F* fReconCosThetaH = (TH1F*)listings->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH = (TH1F*)listings->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  fReconCosThetaH->Sumw2();
  fGenerCosThetaH->Sumw2();
  TH1F* ReconTheta = (TH1F*) fReconCosThetaH->Clone("ReconTheta");

  // TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiCsFrameTwentyfiveBinsH");
  // TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiCsFrameTwentyfiveBinsH");
  TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsH");
  fReconPhiH->Sumw2();
  fGenerPhiH->Sumw2();

  TH1F* fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiCsFrameTwentyfiveBinsH");
  TH1F* fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiCsFrameTwentyfiveBinsH");
  fReconTildePhiH->Sumw2();
  fGenerTildePhiH->Sumw2();

  Double_t Integral_ReconCosTheta = fReconCosThetaH -> Integral();
  Double_t Integral_ReconPhi      = fReconPhiH      -> Integral();
  fReconCosThetaH  -> Scale( 1/Integral_ReconCosTheta  );
  fReconPhiH       -> Scale( 1/Integral_ReconPhi       );


  // TFile* fileDataRawCosTheta = new TFile("pngResults/CosThetaCsFrame.root");
  // TFile* fileDataRawPhi      = new TFile("pngResults/PhiCsFrame.root");
  // TFile* fileDataRawTildePhi = new TFile("pngResults/TildePhiCsFrame.root");
  TDatime d;
  // TFile* fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaCS/CosThetaCSFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiCS/PhiCsFrame.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiCS/TildePhiCsFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );

  TFile* fileDataRawCosTheta = 0x0;
  TFile* fileDataRawPhi      = 0x0;
  TFile* fileDataRawTildePhi = 0x0;
  // fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaCS/CosThetaCSFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiCS/PhiCsFrame.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiCS/TildePhiCsFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );




  TH1F* CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");

  TH1F* PhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  PhiAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawPhiH = (TH1F*) PhiAfterSignalExtractionErrorsRawH->Clone("RawPhiH");

  TH1F* TildePhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
  TildePhiAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawTildePhiH = (TH1F*) TildePhiAfterSignalExtractionErrorsRawH->Clone("RawTildePhiH");


  Double_t Integral_RawCosTheta = RawCosThetaH -> Integral();
  Double_t Integral_RawPhi      = RawPhiH      -> Integral();
  RawCosThetaH  -> Scale( 1/Integral_RawCosTheta  );
  RawPhiH       -> Scale( 1/Integral_RawPhi       );







  new TCanvas;
  RawCosThetaH->Draw();
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  RawCosThetaH->SetLineColor(kRed);
  RawCosThetaH->SetLineWidth(3);
  fReconCosThetaH->SetLineColor(kBlue);
  fReconCosThetaH->SetLineWidth(3);
  RawCosThetaH->GetXaxis()->SetTitleOffset(1.15);
  // RawCosThetaH->GetYaxis()->SetTitleOffset(1.25);
  RawCosThetaH->GetYaxis()->SetTitleOffset(1.);
  RawCosThetaH->GetXaxis()->SetTitleSize(0.045);
  RawCosThetaH->GetYaxis()->SetTitleSize(0.045);
  RawCosThetaH->GetXaxis()->SetLabelSize(0.045);
  RawCosThetaH->GetYaxis()->SetLabelSize(0.045);
  RawCosThetaH->GetXaxis()->SetTitleFont(42);
  RawCosThetaH->GetYaxis()->SetTitleFont(42);
  RawCosThetaH->GetXaxis()->SetLabelFont(42);
  RawCosThetaH->GetYaxis()->SetLabelFont(42);
  RawCosThetaH->GetXaxis()->SetNdivisions(408);
  RawCosThetaH->GetYaxis()->SetRangeUser(0., RawCosThetaH->GetMaximum()*2.0);
  RawCosThetaH->SetTitle( ";#cos(#theta); Counts [a.u.]" );
  RawCosThetaH->Draw("");
  fReconCosThetaH->Draw("same");
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"CosTheta");
  latex->SetTextSize(0.045);
  gPad->BuildLegend();
  gPad->Modified();







  new TCanvas;
  RawPhiH->Draw();
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  RawPhiH->SetLineColor(kRed);
  RawPhiH->SetLineWidth(3);
  fReconPhiH->SetLineColor(kBlue);
  fReconPhiH->SetLineWidth(3);
  RawPhiH->GetXaxis()->SetTitleOffset(1.15);
  // RawPhiH->GetYaxis()->SetTitleOffset(1.25);
  RawPhiH->GetYaxis()->SetTitleOffset(1.);
  RawPhiH->GetXaxis()->SetTitleSize(0.045);
  RawPhiH->GetYaxis()->SetTitleSize(0.045);
  RawPhiH->GetXaxis()->SetLabelSize(0.045);
  RawPhiH->GetYaxis()->SetLabelSize(0.045);
  RawPhiH->GetXaxis()->SetTitleFont(42);
  RawPhiH->GetYaxis()->SetTitleFont(42);
  RawPhiH->GetXaxis()->SetLabelFont(42);
  RawPhiH->GetYaxis()->SetLabelFont(42);
  RawPhiH->GetXaxis()->SetNdivisions(408);
  RawPhiH->GetYaxis()->SetRangeUser(0., RawCosThetaH->GetMaximum()*2.0);
  RawPhiH->SetTitle( ";#phi; Counts [a.u.]" );
  RawPhiH->Draw("");
  fReconPhiH->Draw("same");
  TLatex* latex1 = new TLatex();
  latex1->SetTextSize(0.05);
  latex1->SetTextFont(42);
  latex1->SetTextAlign(11);
  latex1->SetNDC();
  latex1->DrawLatex(0.17,0.94,"Phi");
  latex1->SetTextSize(0.045);
  gPad->BuildLegend();
  gPad->Modified();


  // latex->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  // latex->DrawLatex(0.55,0.78,"Minuit 2D Fit");
  // latex->DrawLatex(0.55,0.70,Form("#lambda_{#theta} = %.3f #pm %.3f",     LambdaTheta,      LambdaThetaErr   ));
  // latex->DrawLatex(0.55,0.62,Form("#lambda_{#phi} = %.3f #pm %.3f",       LambdaPhi,        LambdaPhiErr     ));
  // latex->DrawLatex(0.55,0.54,Form("#lambda_{#theta#phi} = %.3f #pm %.3f", LambdaThetaPhi,   LambdaThetaPhiErr));
  // latex->DrawLatex(0.55,0.44,Form("#tilde{#chi} = %3.3f/%3.3d = %2.2f",      GlobalChi,        ndf, (Double_t)GlobalChi/(Double_t)ndf ));



}
