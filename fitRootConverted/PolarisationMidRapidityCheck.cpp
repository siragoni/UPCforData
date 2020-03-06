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
/* - Mid-rapidity check of polarisation feasibility.
 * -
 */
void PolarisationMidRapidityCheck(){

  // TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  // TDirectory* dir = fileList->GetDirectory("MyTask");
  // TList* listings;
  // dir->GetObject("MyOutputContainer", listings);
  // TDatime d;

  TFile* fileDataRawCosTheta = 0x0;
  TFile* fileDataRawPhi      = 0x0;
  TFile* fileDataRawTildePhi = 0x0;
  fileDataRawCosTheta = new TFile("pngResults/2019-12-17/CosThetaHE/CosThetaHeFrame.root");
  fileDataRawPhi      = new TFile("pngResults/2019-12-17/PhiHE/PhiHeFrame.root"          );
  fileDataRawTildePhi = new TFile("pngResults/2019-12-17/TildePhiHE/TildePhiHeFrame.root");

  TH1F* CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");

  TCanvas* CorrCanvasCosTheta = new TCanvas("CorrCanvasCosTheta","CorrCanvasCosTheta",900,800);
  RawCosThetaH->Draw("ep");

  TH1F*  CorrCosThetaH  = (TH1F*) RawCosThetaH->Clone("CorrCosThetaH");
  TFile* correctedPlot  = new TFile("pngResults/2019-12-17/1Dresults/PolarisationCorrectedHe1D.root");
  TH1F*  CorrectedPlotH = (TH1F*)correctedPlot->Get("CorrCosThetaH");

  TFile* midRapidity    = new TFile("polarizationCCUP31-B-NOPF-CENTNOTRD.root");
  TH1F*  midRapidityH   = (TH1F*)midRapidity->Get("fCosThetaHelicityFrameJPsiH");
  // TH1F*  midRapidityH   = (TH1F*)midRapidity->Get("fCosThetaCollinsSoperFrameJPsiH");
  midRapidityH->Rebin(40);
  Double_t pointsMidRapidity[25];
  Double_t pointsForward[25];
  for( Int_t i = 0; i < 25; i++ ){
    pointsMidRapidity[i] = 0;
  }
  for( Int_t i = 0; i < 25; i++ ){
    pointsForward[i] = 0;
  }



  for( Int_t i = 0; i < 25; i++ ){
    pointsMidRapidity[i] = midRapidityH->GetBinContent(i+1);
    cout << "midRapidity  " << i << " = " << pointsMidRapidity[i] << endl;
  }

  for( Int_t i = 0; i < 25; i++ ){
    pointsForward[i] = RawCosThetaH->GetBinContent(i+1);
    cout << "pointsForward  " << i << " = " << pointsForward[i] << endl;
  }


  Double_t weights[25];
  for (size_t i = 0; i < 25; i++) {
    weights[i] = 0;
    if ( pointsMidRapidity[i] != 0 ){
      weights[i] = TMath::Sqrt(pointsForward[i]/pointsMidRapidity[i]);
    } else {
      weights[i] = 0;
    }
    cout << "weights  " << i << " = " << weights[i] << endl;
  }


  Double_t previousErrors[25];
  for (size_t i = 0; i < 25; i++) {
    previousErrors[i] = 0;
    previousErrors[i] = CorrectedPlotH->GetBinError(i+1);
    cout << "previousErrors  " << i << " = " << previousErrors[i] << endl;
  }


  Double_t newErrors[25];
  for (size_t i = 0; i < 25; i++) {
    newErrors[i] = 0;
    newErrors[i] = previousErrors[i] * weights[i];
    cout << "newErrors  " << i << " = " << newErrors[i] << endl;
  }


  for (size_t i = 0; i < 25; i++) {
    newErrors[i] = 0;
    newErrors[i] = previousErrors[i] * weights[i];
    cout << "newErrors  " << i << " = " << newErrors[i] << endl;
  }



  TCanvas* CorrCanvasCosThetaV2 = new TCanvas("CorrCanvasCosThetaV2","CorrCanvasCosThetaV2",900,800);
  for (size_t i = 0; i < 25; i++) {
    CorrectedPlotH->SetBinError(i+1, newErrors[i]);
  }

  CorrectedPlotH->Draw("ep");

  // acceptanceCosTheta->Write();
  // CorrCosThetaH     ->Write();
  // acceptancePhi     ->Write();
  // CorrPhiH          ->Write();
  // acceptanceTildePhi->Write();
  // CorrTildePhiH     ->Write();
  // // AccErrors         ->Write();
  // // EntErrors         ->Write();
  // // ReconTheta        ->Write();
  // f.Close();

}
