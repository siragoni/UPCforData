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
void PolarisationCorrectedDistribHe1D( Int_t selectionFlag = 0 ){

  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_19032021.root");
  // TFile* fileList2 = new TFile("AnalysisResultsLHC18l7_long_19032021.root"); // trial longitudinal
  TFile* fileList = new TFile("AnalysisResultsLHC18l7_trans_30032021.root");
  TFile* fileList2 = new TFile("AnalysisResultsLHC18l7_long_rapidity_30032021.root"); // trial longitudinal
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TDirectory* dir2 = fileList2->GetDirectory("MyTask");
  TList* listings;
  TList* listings2;
  dir->GetObject("MyOutputContainer", listings);
  dir2->GetObject("MyOutputContainer", listings2);
  TFile* SavedFile = new TFile("SavedCosThetaLongitudinalAfterBkg.root");
  // TH1F* fReconCosThetaH_nice = (TH1F*)SavedFile->Get("fCosThetaReconstructed");

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
  TH1F* fReconCosThetaH = 0x0;
  TH1F* fGenerCosThetaH = 0x0;
  if (selectionFlag == 0) {
    fReconCosThetaH = (TH1F*)listings->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
    fGenerCosThetaH = (TH1F*)listings->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  } else {
    fReconCosThetaH = (TH1F*)SavedFile->Get("fCosThetaReconstructed");
    fGenerCosThetaH = (TH1F*)listings2->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  }

  // fReconCosThetaH->Rebin(4);
  // fGenerCosThetaH->Rebin(4);
  fReconCosThetaH->Sumw2();
  fGenerCosThetaH->Sumw2();
  // new TCanvas;
  // fGenerCosThetaH->Draw();
  TH1F* ReconTheta = (TH1F*) fReconCosThetaH->Clone("ReconTheta");

  TH1F* fReconPhiH = 0x0;
  TH1F* fGenerPhiH = 0x0;
  if (selectionFlag == 0) {
    fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
    fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  } else {
    // fReconPhiH = (TH1F*)listings2->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
    // fGenerPhiH = (TH1F*)listings2->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");
    fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
    fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");

  }
  fReconPhiH->Sumw2();
  fGenerPhiH->Sumw2();

  TH1F* fReconTildePhiH = 0x0;
  TH1F* fGenerTildePhiH = 0x0;
  if (selectionFlag == 0) {
    fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
    fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
  } else {
    // fReconTildePhiH = (TH1F*)listings2->FindObject("fTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
    // fGenerTildePhiH = (TH1F*)listings2->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
    fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
    fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");

  }
  fReconTildePhiH->Sumw2();
  fGenerTildePhiH->Sumw2();


  TH1F* CosThetaAfterSignalExtractionErrorsRawH = 0x0;
  if (selectionFlag == 0) {
    CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fReconCosThetaH->Clone("CosThetaAfterSignalExtractionErrorsRawH");
  } else {
    // CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)listings->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
    CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fReconCosThetaH->Clone("CosThetaAfterSignalExtractionErrorsRawH");
  }
  CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");

  TH1F* PhiAfterSignalExtractionErrorsRawH = 0x0;
  if (selectionFlag == 0) {
    PhiAfterSignalExtractionErrorsRawH = (TH1F*) fReconPhiH->Clone("PhiAfterSignalExtractionErrorsRawH");
  } else {
    // PhiAfterSignalExtractionErrorsRawH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
    PhiAfterSignalExtractionErrorsRawH = (TH1F*) fReconPhiH->Clone("PhiAfterSignalExtractionErrorsRawH");
  }
  PhiAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawPhiH = (TH1F*) PhiAfterSignalExtractionErrorsRawH->Clone("RawPhiH");

  TH1F* TildePhiAfterSignalExtractionErrorsRawH = 0x0;
  if (selectionFlag == 0) {
    TildePhiAfterSignalExtractionErrorsRawH = (TH1F*) fReconTildePhiH->Clone("TildePhiAfterSignalExtractionErrorsRawH");
  } else {
    TildePhiAfterSignalExtractionErrorsRawH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
  }
  TildePhiAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawTildePhiH = (TH1F*) TildePhiAfterSignalExtractionErrorsRawH->Clone("RawTildePhiH");


  TCanvas* AcceptanceCanvasCosTheta = new TCanvas("AcceptanceCanvasCosTheta","AcceptanceCanvasCosTheta",900,800);
  TH1F* acceptanceCosTheta = (TH1F*) fReconCosThetaH->Clone("acceptanceCosTheta");
  acceptanceCosTheta->Divide(fGenerCosThetaH);
  acceptanceCosTheta->Draw("ep");

  TCanvas* CorrCanvasCosTheta = new TCanvas("CorrCanvasCosTheta","CorrCanvasCosTheta",900,800);
  RawCosThetaH->Divide(acceptanceCosTheta);
  RawCosThetaH->Draw("ep");

  TCanvas* AcceptanceCanvasPhi = new TCanvas("AcceptanceCanvasPhi","AcceptanceCanvasPhi",900,800);
  TH1F* acceptancePhi = (TH1F*) fReconPhiH->Clone("acceptancePhi");
  acceptancePhi->Divide(fGenerPhiH);
  acceptancePhi->Draw("ep");

  TCanvas* CorrCanvasPhi = new TCanvas("CorrCanvasPhi","CorrCanvasPhi",900,800);
  RawPhiH->Divide(acceptancePhi);
  RawPhiH->Draw("ep");

  TCanvas* AcceptanceCanvasTildePhi = new TCanvas("AcceptanceCanvasTildePhi","AcceptanceCanvasTildePhi",900,800);
  TH1F* acceptanceTildePhi = (TH1F*) fReconTildePhiH->Clone("acceptanceTildePhi");
  acceptanceTildePhi->Divide(fGenerTildePhiH);
  acceptanceTildePhi->Draw("ep");

  TCanvas* CorrCanvasTildePhi = new TCanvas("CorrCanvasTildePhi","CorrCanvasTildePhi",900,800);
  RawTildePhiH->Divide(acceptanceTildePhi);
  RawTildePhiH->Draw("ep");


  TH1F* CorrCosThetaH = (TH1F*) RawCosThetaH->Clone("CorrCosThetaH");
  TH1F* CorrPhiH      = (TH1F*) RawPhiH     ->Clone("CorrPhiH");
  TH1F* CorrTildePhiH = (TH1F*) RawTildePhiH->Clone("CorrTildePhiH");


  if        ( selectionFlag == 0 ) {
    TDatime d;
    TFile f(Form("pngResults/%d-%2.2d-%2.2d/Closure/PolarisationCorrectedHe1Dv2_RecTAxET.root",   d.GetYear(), d.GetMonth(), d.GetDay()), "recreate");
    // TFile f("pngResults/PolarisationCorrectedHe1D.root", "recreate");
    acceptanceCosTheta->Write();
    CorrCosThetaH     ->Write();
    acceptancePhi     ->Write();
    CorrPhiH          ->Write();
    acceptanceTildePhi->Write();
    CorrTildePhiH     ->Write();
    // AccErrors         ->Write();
    // EntErrors         ->Write();
    // ReconTheta        ->Write();
    f.Close();
  } else if ( selectionFlag == 1 ) {
    TDatime d;
    TFile f(Form("pngResults/%d-%2.2d-%2.2d/Closure/PolarisationCorrectedHe1Dv2_RecTAxEL.root",   d.GetYear(), d.GetMonth(), d.GetDay()), "recreate");
    // TFile f("pngResults/PolarisationCorrectedHe1D.root", "recreate");
    acceptanceCosTheta->Write();
    CorrCosThetaH     ->Write();
    acceptancePhi     ->Write();
    CorrPhiH          ->Write();
    acceptanceTildePhi->Write();
    CorrTildePhiH     ->Write();
    // AccErrors         ->Write();
    // EntErrors         ->Write();
    // ReconTheta        ->Write();
    f.Close();
  }
}
