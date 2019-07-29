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


#include "TH2.h"

//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitHelicityCorrected1D(){

  TFile* fileList = new TFile("AnalysisResultsMC03072019.root");
  // TFile* fileList = new TFile("MCtrainResults/2019-06-08/kCohJpsiToMu/AnalysisResults.root");
  // TFile* fileList = new TFile("MCtrainResults/2019-06-24/kCohJpsiToMu/AnalysisResults.root");
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
  TH1F* fReconCosThetaH = (TH1F*)listings->FindObject("fCosThetaHelicityFrameMyBinningH");
  TH1F* fGenerCosThetaH = (TH1F*)listings->FindObject("fMCCosThetaHelicityFrameMyBinningH");
  fReconCosThetaH->Sumw2();
  fGenerCosThetaH->Sumw2();
  TH1F* ReconTheta = (TH1F*) fReconCosThetaH->Clone("ReconTheta");

  // TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiHeFrameForSignalExH");
  // TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHeFrameForSignalExH");
  // fReconPhiH->Sumw2();
  // fGenerPhiH->Sumw2();


  // gSystem->cd("pngResults/2019-05-30-18qr15o-NoSPD/SignalExtraction/");

  // TFile* fileDataRawCosTheta = new TFile("pngResults/TH1signalCosThetaEX.root");
  TFile* fileDataRawCosTheta = new TFile("pngResults/TH1functionalCosThetaMyBinningEX.root");
  // TFile* fileDataRawPhi      = new TFile("pngResults/TH1signalPhiEX.root");
  // TFile* fileDataRaw = new TFile("TH2signalEX.root");
  TH1F* CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");

  // TH1F* PhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  // PhiAfterSignalExtractionErrorsRawH->Sumw2();
  // TH1F* RawPhiH = (TH1F*) PhiAfterSignalExtractionErrorsRawH->Clone("RawPhiH");


  TCanvas* AcceptanceCanvasCosTheta = new TCanvas("AcceptanceCanvasCosTheta","AcceptanceCanvasCosTheta",900,800);
  TH1F* acceptanceCosTheta = (TH1F*) fReconCosThetaH->Clone("acceptanceCosTheta");
  acceptanceCosTheta->Divide(fGenerCosThetaH);
  acceptanceCosTheta->Draw("ep");

  TCanvas* CorrCanvasCosTheta = new TCanvas("CorrCanvasCosTheta","CorrCanvasCosTheta",900,800);
  RawCosThetaH->Divide(acceptanceCosTheta);
  RawCosThetaH->Draw("colZ");

  // TCanvas* AcceptanceCanvasPhi = new TCanvas("AcceptanceCanvasPhi","AcceptanceCanvasPhi",900,800);
  // TH1F* acceptancePhi = (TH1F*) fReconPhiH->Clone("acceptancePhi");
  // acceptancePhi->Divide(fGenerPhiH);
  // acceptancePhi->Draw("ep");
  //
  // TCanvas* CorrCanvasPhi = new TCanvas("CorrCanvasPhi","CorrCanvasPhi",900,800);
  // RawPhiH->Divide(acceptancePhi);
  // RawPhiH->Draw("colZ");

  // TH1F* AccErrors = new TH1F("AccErrors", "AccErrors", 40, -1, 1);
  // TH1F* EntErrors = new TH1F("EntErrors", "EntErrors", 40, -1, 1);
  // for( Int_t iLoop = 1; iLoop <= 40; iLoop++ ){
  //   Double_t AcceptanceEntries = ReconTheta        ->GetBinContent(iLoop);
  //   Double_t AcceptanceValue   = acceptanceCosTheta->GetBinContent(iLoop);
  //   if( AcceptanceEntries >= 1  ){
  //     AccErrors->SetBinContent( iLoop, AcceptanceValue/(TMath::Sqrt(AcceptanceEntries)) );
  //     EntErrors->SetBinContent( iLoop, 1/(TMath::Sqrt(AcceptanceEntries)) );
  //   }
  // }


  TFile f("pngResults/TH1corrMyBinning.root", "recreate");
  acceptanceCosTheta->Write();
  RawCosThetaH      ->Write();
  // acceptancePhi     ->Write();
  // RawPhiH           ->Write();
  // AccErrors         ->Write();
  // EntErrors         ->Write();
  // ReconTheta        ->Write();
  f.Close();
}
