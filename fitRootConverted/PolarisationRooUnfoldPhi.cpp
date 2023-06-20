#include <TError.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include "TUnfoldDensity.h"
// #define VERBOSE_LCURVE_SCAN
using namespace std;












int testUnfold1()
{

  // do it before running
  // gSystem->Load("../RooUnfold/libRooUnfold");
  RooUnfoldResponse response(10,-1,1);

  // switch on histogram errors
  TH1::SetDefaultSumw2();
  // show fit result
  gStyle->SetOptFit(1111);

  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_19032021.root");
  TFile* fileList = new TFile("AnalysisResultsLHC18l7_coherent_30112021.root");
  // TFile* fileList = new TFile("AnalysisResults.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);
  //============================================
  // generate MC distribution
  //
  TH1F *histMgenMC    = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  // TH1F *histMdetMC    = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  TH1D *histMdetMC    = (TH1D*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  // TH2F *histMdetGenMC = (TH2F*)listings->FindObject("fPhiRecVsGenHelicityH");
  TH2D *histMdetGenMC = (TH2D*)listings->FindObject("fPhiRecVsGenHelicityH");
  // generate data distribution
  //
  // TH1D *histMgenData = new TH1D("MgenData",";mass(gen)",nGen,xminGen,xmaxGen);
  TFile* fileDataRawPhi = new TFile("pngResults/2021-09-21/PhiHEv2/PhiHeFrameV2.root");
  TH1F *histMdetData    = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");



  RooUnfoldResponse *h1ameas_h1atrue = new RooUnfoldResponse(histMdetMC, histMgenMC, histMdetGenMC);
  RooUnfoldBayes unfold1(h1ameas_h1atrue, histMdetMC, 1);
  RooUnfoldBayes unfold2(h1ameas_h1atrue, histMdetMC, 2);
  RooUnfoldBayes unfold3(h1ameas_h1atrue, histMdetMC, 3);
  RooUnfoldBayes unfold4(h1ameas_h1atrue, histMdetMC, 4);

  TH1D* hist_reco1= (TH1D*)unfold1.Hreco();
  TH1D* hist_reco2= (TH1D*)unfold2.Hreco();
  TH1D* hist_reco3= (TH1D*)unfold3.Hreco();
  TH1D* hist_reco4= (TH1D*)unfold4.Hreco();


  hist_reco1->SetLineColor(kRed);
  hist_reco2->SetLineColor(kGreen);
  hist_reco3->SetLineColor(kBlack);
  hist_reco4->SetLineColor(kBlue);

  new TCanvas;
  hist_reco1->Draw();
  hist_reco2->Draw("same");
  hist_reco3->Draw("same");
  hist_reco4->Draw("same");
  // cout<<"1 iteration "<< "2 iterations "<< "3 iterations "<< "4 iterations " <<endl;
  // for(Int_t i=1;i<15;i++){
  //
  // cout<<hist_reco1->GetBinContent(i)<< " "<<hist_reco2->GetBinContent(i)<< " "<<hist_reco3->GetBinContent(i)<< " "<<hist_reco4->GetBinContent(i)<< " "<<endl;
  // }

return 1;
}
