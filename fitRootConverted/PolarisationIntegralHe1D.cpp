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

  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_long_19032021.root"); // trial longitudinal
  TFile* fileList = new TFile("AnalysisResultsLHC18l7_19032021.root"); 
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

  // TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsH");
  // TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsH");
  TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  fReconPhiH->Sumw2();
  fGenerPhiH->Sumw2();

  // TH1F* fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsH");
  // TH1F* fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsH");
  TH1F* fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
  TH1F* fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsHv2_restrict");
  fReconTildePhiH->Sumw2();
  fGenerTildePhiH->Sumw2();

  TDatime d;

  TFile* fileDataRawCosTheta = 0x0;
  TFile* fileDataRawPhi      = 0x0;
  TFile* fileDataRawTildePhi = 0x0;
  fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root",       d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHEv2/PhiHeFrameV2.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
  fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHEv2/TildePhiHeFrameV2.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );

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


  // // TDatime d;
  // TFile f(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2_long.root",   d.GetYear(), d.GetMonth(), d.GetDay()), "recreate");
  // // TFile f("pngResults/PolarisationCorrectedHe1D.root", "recreate");
  // acceptanceCosTheta->Write();
  // CorrCosThetaH     ->Write();
  // acceptancePhi     ->Write();
  // CorrPhiH          ->Write();
  // acceptanceTildePhi->Write();
  // CorrTildePhiH     ->Write();




  Double_t IntegralCosTheta = 0.;
  Double_t IntegralPhi      = 0.;
  Double_t IntegralTildePhi = 0.;
  Double_t ErrorCosTheta    = 0.;
  Double_t ErrorPhi         = 0.;
  Double_t ErrorTildePhi    = 0.;


  for (size_t i = 1; i < 26; i++) {

    if ( CorrCosThetaH->GetBinContent(i) > 0.0001 ) {
      IntegralCosTheta += CorrCosThetaH->GetBinContent(i);
      ErrorCosTheta    += CorrCosThetaH->GetBinError(i);
    }
    if ( CorrPhiH->GetBinContent(i) > 0.0001 ) {
      IntegralPhi += CorrPhiH->GetBinContent(i);
      ErrorPhi    += CorrPhiH->GetBinError(i);
    }
    if ( CorrTildePhiH->GetBinContent(i) > 0.0001 ) {
      IntegralTildePhi += CorrTildePhiH->GetBinContent(i);
      ErrorTildePhi    += CorrTildePhiH->GetBinError(i);
    }


  }



  cout << "IntegralCosTheta = " << IntegralCosTheta << " +/- " << ErrorCosTheta << endl;
  cout << "IntegralPhi      = " << IntegralPhi      << " +/- " << ErrorPhi      << endl;
  cout << "IntegralTildePhi = " << IntegralTildePhi << " +/- " << ErrorTildePhi << endl;
}
