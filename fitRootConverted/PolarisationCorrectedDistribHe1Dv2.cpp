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
void PolarisationCorrectedDistribHe1D( Int_t RangeSelectionMode = 0 ){

  // TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  // TFile* fileList = new TFile("AnalysisResultsLHC18l716b2_nophoton_phi.root");
  TFile* fileList = new TFile("AnalysisResultsLHC18l7_19032021.root");
  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_trans_30032021.root");
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
  TH1F* fReconCosThetaH = (TH1F*)listings->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH = (TH1F*)listings->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  // fReconCosThetaH->Rebin(4);
  // fGenerCosThetaH->Rebin(4);
  fReconCosThetaH->Sumw2();
  fGenerCosThetaH->Sumw2();
  // new TCanvas;
  // fGenerCosThetaH->Draw();
  TH1F* ReconTheta = (TH1F*) fReconCosThetaH->Clone("ReconTheta");

  TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2");
  TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2");
  // TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  // TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  fReconPhiH->Sumw2();
  fGenerPhiH->Sumw2();

  TH1F* fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsHv2");
  TH1F* fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsHv2");
  // TH1F* fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
  // TH1F* fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
  fReconTildePhiH->Sumw2();
  fGenerTildePhiH->Sumw2();

  TDatime d;
  // TFile* fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );

  TFile* fileDataRawCosTheta = 0x0;
  TFile* fileDataRawPhi      = 0x0;
  TFile* fileDataRawTildePhi = 0x0;
  if (        RangeSelectionMode == 0 ) {
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root",       d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHEv2/PhiHeFrameV2.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHEv2/TildePhiHeFrameV2.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 1 ) {
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_1.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_1.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_1.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 2 ) {
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_2.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 3 ) {
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_3.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 4 ) {
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_4.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 5 ) {
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame_5.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else {
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  }

  TH1F* CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");

  TH1F* PhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  PhiAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawPhiH = (TH1F*) PhiAfterSignalExtractionErrorsRawH->Clone("RawPhiH");

  TH1F* TildePhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
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

  if        ( RangeSelectionMode == 0 ) {
    TDatime d;
    TFile f(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2.root",   d.GetYear(), d.GetMonth(), d.GetDay()), "recreate");
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
  } else if ( RangeSelectionMode == 1 ) {
    TFile f("pngResults/PolarisationCorrectedHe1D_1.root", "recreate");
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
  } else if ( RangeSelectionMode == 2 ) {
    TFile f("pngResults/PolarisationCorrectedHe1D_2.root", "recreate");
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
  } else if ( RangeSelectionMode == 3 ) {
    TFile f("pngResults/PolarisationCorrectedHe1D_3.root", "recreate");
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
  } else if ( RangeSelectionMode == 4 ) {
    TFile f("pngResults/PolarisationCorrectedHe1D_4.root", "recreate");
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
  } else if ( RangeSelectionMode == 5 ) {
    TFile f("pngResults/PolarisationCorrectedHe1D_5.root", "recreate");
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
  } else {
    TFile f("pngResults/PolarisationCorrectedHe1D.root", "recreate");
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
