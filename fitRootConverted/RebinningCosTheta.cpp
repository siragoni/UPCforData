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
/* - Histograms to be used for the fit.
 * - What happens is that we will interpolate the many points together...
 * -
 */
TH1F*     fCohJpsiToMuLHC18l7;
TH1F*     fCohJpsiToMuLHC16b2;
TH1F*     fCohJpsiToMuGenLHC18l7;
TH1F*     fCohJpsiToMuGenLHC16b2;
TH1F*     fCosThetaData;
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void RebinningCosTheta( const char* AnalysisName = "AnalysisResultsLHC18qr15o17072019.root" ){

  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);

  TFile* fileMC[2];
  fileMC[0] = new TFile("MCtrainResults/2019-06-24/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("AnalysisResultsLHC16b2MC27062019.root");
  TDirectory* dirMC[2];
  for(Int_t iDirectory = 0; iDirectory < 2; iDirectory++) {
    dirMC[iDirectory] = fileMC[iDirectory]->GetDirectory("MyTask");
  }
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listingsMC[2];
  for(Int_t iDirectory = 0; iDirectory < 2; iDirectory++) {
    dirMC[iDirectory]->GetObject("MyOutputContainer", listingsMC[iDirectory]);
  }
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
  fCohJpsiToMuLHC18l7    = (TH1F*)listingsMC[0]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fCohJpsiToMuLHC16b2    = (TH1F*)listingsMC[1]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fCohJpsiToMuGenLHC18l7 = (TH1F*)listingsMC[0]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fCohJpsiToMuGenLHC16b2 = (TH1F*)listingsMC[1]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fCosThetaData          = (TH1F*)listings     ->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");

  TH1F* Recon = (TH1F*) fCohJpsiToMuLHC18l7->Clone("Recon");
  Recon->Add(fCohJpsiToMuLHC16b2);
  Recon->Sumw2();
  TH1F* Gener = (TH1F*) fCohJpsiToMuGenLHC18l7->Clone("Gener");
  Gener->Add(fCohJpsiToMuGenLHC16b2);
  Gener->Sumw2();

  TH1F* Acceptance = (TH1F*) Recon->Clone("Acceptance");
  Acceptance->Divide(Gener);

  new TCanvas;
  Acceptance->Draw("ep");

  Double_t MyVariableCosThetaBinning1Dv3[] = { -0.65,  -0.5,  -0.425, -0.35,
                                               -0.275, -0.2,  -0.125, -0.075,
                                               -0.025,  0.025, 0.075,  0.125,
                                                0.2,    0.275, 0.35,   0.425,
                                                0.5,    0.65 };

  TH1F* ReconSeventeenBins = new TH1F( "ReconSeventeenBins",
                                       "ReconSeventeenBins",
                                       17, MyVariableCosThetaBinning1Dv3);
  TH1F* GenerSeventeenBins = new TH1F( "GenerSeventeenBins",
                                       "GenerSeventeenBins",
                                       17, MyVariableCosThetaBinning1Dv3);
  Double_t errorRecon[1000];
  Double_t errorGener[1000];
  Double_t errorReconSeventeenBins[17];
  Double_t errorGenerSeventeenBins[17];
  Int_t counter = 0;
  for( Int_t iBins = 1; iBins < 1001; iBins++ ) {
    Double_t ReconValue = Recon->GetBinContent(iBins);
    Double_t GenerValue = Gener->GetBinContent(iBins);
    Double_t x          = Recon->GetBinCenter (iBins);
    if( x > MyVariableCosThetaBinning1Dv3[counter] ) counter+=1;
    ReconSeventeenBins->Fill( x, ReconValue );
    GenerSeventeenBins->Fill( x, GenerValue );
    // errorRecon[iBins - 1] = Recon->GetBinError(iBins);
    // errorGener[iBins - 1] = Gener->GetBinError(iBins);
    // errorReconSeventeenBins[counter - 1] += errorRecon[iBins - 1]*errorRecon[iBins - 1];
    // errorGenerSeventeenBins[counter - 1] += errorGener[iBins - 1]*errorGener[iBins - 1];
    if( counter > 0. ) {
      errorReconSeventeenBins[counter - 1] += Recon->GetBinError(iBins)*Recon->GetBinError(iBins);
      errorGenerSeventeenBins[counter - 1] += Gener->GetBinError(iBins)*Gener->GetBinError(iBins);
    }
  }
  for(Int_t iBins = 0; iBins < 17; iBins++ ) {
    errorReconSeventeenBins[iBins] = TMath::Sqrt(errorReconSeventeenBins[iBins]);
    errorGenerSeventeenBins[iBins] = TMath::Sqrt(errorGenerSeventeenBins[iBins]);
    ReconSeventeenBins->SetBinError( iBins + 1, errorReconSeventeenBins[iBins] );
    GenerSeventeenBins->SetBinError( iBins + 1, errorGenerSeventeenBins[iBins] );
  }
  new TCanvas;
  ReconSeventeenBins->Draw("ep");
  new TCanvas;
  GenerSeventeenBins->Draw("ep");

  ReconSeventeenBins->Sumw2();
  GenerSeventeenBins->Sumw2();
  TH1F* AcceptanceSeventeenBins = (TH1F*) ReconSeventeenBins->Clone("AcceptanceSeventeenBins");
  AcceptanceSeventeenBins->Divide(GenerSeventeenBins);
  new TCanvas;
  AcceptanceSeventeenBins->Draw();

  TFile* fileDataRawCosTheta = new TFile("pngResults/TH1functionalCosThetaMySeventeenBinningEX.root");
  TH1F* CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");

  new TCanvas;
  CosThetaAfterSignalExtractionErrorsRawH->Draw("ep");
  new TCanvas;
  RawCosThetaH->Divide(AcceptanceSeventeenBins);
  RawCosThetaH->Draw("ep");














  Double_t CosThetaBinning[25];
  for( Int_t iBins = 0; iBins < 25; iBins++ ){
    CosThetaBinning[iBins] = -1 + 0.08 * iBins;
  }

  TH1F* ReconTwentyfiveBins = new TH1F( "ReconTwentyfiveBins",
                                        "ReconTwentyfiveBins",
                                        25, -1, 1);
  TH1F* GenerTwentyfiveBins = new TH1F( "GenerTwentyfiveBins",
                                        "GenerTwentyfiveBins",
                                        25, -1, 1);
  Double_t errorRecon2[1000];
  Double_t errorGener2[1000];
  Double_t errorReconTwentyfiveBins[25];
  Double_t errorGenerTwentyfiveBins[25];
  Int_t counter2 = 0;
  for( Int_t iBins = 1; iBins < 1001; iBins++ ) {
    Double_t ReconValue = Recon->GetBinContent(iBins);
    Double_t GenerValue = Gener->GetBinContent(iBins);
    Double_t x          = Recon->GetBinCenter (iBins);
    if( x > CosThetaBinning[counter] ) counter+=1;
    ReconTwentyfiveBins->Fill( x, ReconValue );
    GenerTwentyfiveBins->Fill( x, GenerValue );
    // errorRecon[iBins - 1] = Recon->GetBinError(iBins);
    // errorGener[iBins - 1] = Gener->GetBinError(iBins);
    // errorReconSeventeenBins[counter - 1] += errorRecon[iBins - 1]*errorRecon[iBins - 1];
    // errorGenerSeventeenBins[counter - 1] += errorGener[iBins - 1]*errorGener[iBins - 1];
    if( counter > 0. ) {
      errorReconTwentyfiveBins[counter - 1] += Recon->GetBinError(iBins)*Recon->GetBinError(iBins);
      errorGenerTwentyfiveBins[counter - 1] += Gener->GetBinError(iBins)*Gener->GetBinError(iBins);
    }
  }
  for(Int_t iBins = 0; iBins < 25; iBins++ ) {
    errorReconTwentyfiveBins[iBins] = TMath::Sqrt(errorReconTwentyfiveBins[iBins]);
    errorGenerTwentyfiveBins[iBins] = TMath::Sqrt(errorGenerTwentyfiveBins[iBins]);
    ReconTwentyfiveBins->SetBinError( iBins + 1, errorReconTwentyfiveBins[iBins] );
    GenerTwentyfiveBins->SetBinError( iBins + 1, errorGenerTwentyfiveBins[iBins] );
  }
  new TCanvas;
  ReconTwentyfiveBins->Draw("ep");
  new TCanvas;
  GenerTwentyfiveBins->Draw("ep");

  ReconTwentyfiveBins->Sumw2();
  GenerTwentyfiveBins->Sumw2();
  TH1F* AcceptanceTwentyfiveBins = (TH1F*) ReconTwentyfiveBins->Clone("AcceptanceTwentyfiveBins");
  AcceptanceTwentyfiveBins->Divide(GenerTwentyfiveBins);
  new TCanvas;
  AcceptanceTwentyfiveBins->Draw();

  TFile* fileDataRawCosTheta2 = new TFile("pngResults/TH1functionalCosTheta25binsEX.root");
  TH1F* CosThetaAfterSignalExtractionErrorsRaw2H = (TH1F*)fileDataRawCosTheta2->Get("CosThetaAfterSignalExtractionErrorsH");
  CosThetaAfterSignalExtractionErrorsRaw2H->Sumw2();
  TH1F* RawCosTheta2H = (TH1F*) CosThetaAfterSignalExtractionErrorsRaw2H->Clone("RawCosTheta2H");

  new TCanvas;
  CosThetaAfterSignalExtractionErrorsRaw2H->Draw("ep");
  new TCanvas;
  RawCosTheta2H->Divide(AcceptanceTwentyfiveBins);
  RawCosTheta2H->Draw("ep");

  TFile* outputFile = new TFile("pngResults/TH1corr25bins.root", "recreate");
  RawCosTheta2H->Write();
  outputFile->Close();

  // fCosThetaData->Sumw2();
  // fCosThetaData->Divide(Acceptance);
  // fCosThetaData->Draw("ep");

  // /* - Rebin
  //    -
  //  */
  // // fCohJpsiToMu           -> Rebin(5);
  // // fCohPsi2sToMu          -> Rebin(5);
  // // fCohPsi2sToMuPi        -> Rebin(5);
  // // fIncohJpsiToMu         -> Rebin(5);
  // // fIncohPsi2sToMu        -> Rebin(5);
  // // fIncohPsi2sToMuPi      -> Rebin(5);
  // // fTwoGammaToMuMedium    -> Rebin(5);
  // // fTwoGammaToMuHigh      -> Rebin(5);
  // // fCohJpsiToMuGen        -> Rebin(5);
  // // fCohPsi2sToMuGen       -> Rebin(5);
  // // fCohPsi2sToMuPiGen     -> Rebin(5);
  // // fIncohJpsiToMuGen      -> Rebin(5);
  // // fIncohPsi2sToMuGen     -> Rebin(5);
  // // fIncohPsi2sToMuPiGen   -> Rebin(5);
  // // fTwoGammaToMuMediumGen -> Rebin(5);
  // // fTwoGammaToMuHighGen   -> Rebin(5);
  // // fCohJpsiToMu           -> Rebin(20);
  // // fCohPsi2sToMu          -> Rebin(20);
  // // fCohPsi2sToMuPi        -> Rebin(20);
  // // fIncohJpsiToMu         -> Rebin(20);
  // // fIncohPsi2sToMu        -> Rebin(20);
  // // fIncohPsi2sToMuPi      -> Rebin(20);
  // // fTwoGammaToMuMedium    -> Rebin(20);
  // // fTwoGammaToMuHigh      -> Rebin(20);
  // // fCohJpsiToMuGen        -> Rebin(20);
  // // fCohPsi2sToMuGen       -> Rebin(20);
  // // fCohPsi2sToMuPiGen     -> Rebin(20);
  // // fIncohJpsiToMuGen      -> Rebin(20);
  // // fIncohPsi2sToMuGen     -> Rebin(20);
  // // fIncohPsi2sToMuPiGen   -> Rebin(20);
  // // fTwoGammaToMuMediumGen -> Rebin(20);
  // // fTwoGammaToMuHighGen   -> Rebin(20);
  //
  // /* - Firstly we normalize the histograms.
  //    - Remember to always Sumw2()!!
  //    -
  //  */
  // fCohJpsiToMu           -> Sumw2();
  // fCohPsi2sToMu          -> Sumw2();
  // fCohPsi2sToMuPi        -> Sumw2();
  // fIncohJpsiToMu         -> Sumw2();
  // fIncohPsi2sToMu        -> Sumw2();
  // fIncohPsi2sToMuPi      -> Sumw2();
  // fTwoGammaToMuMedium    -> Sumw2();
  // fTwoGammaToMuHigh      -> Sumw2();
  // fCohJpsiToMuGen        -> Sumw2();
  // fCohPsi2sToMuGen       -> Sumw2();
  // fCohPsi2sToMuPiGen     -> Sumw2();
  // fIncohJpsiToMuGen      -> Sumw2();
  // fIncohPsi2sToMuGen     -> Sumw2();
  // fIncohPsi2sToMuPiGen   -> Sumw2();
  // fTwoGammaToMuMediumGen -> Sumw2();
  // fTwoGammaToMuHighGen   -> Sumw2();
  // // Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral();
  // // Double_t Integral_fCohPsi2sToMu       = fCohPsi2sToMu       -> Integral();
  // // Double_t Integral_fCohPsi2sToMuPi     = fCohPsi2sToMuPi     -> Integral();
  // // Double_t Integral_fIncohJpsiToMu      = fIncohJpsiToMu      -> Integral();
  // // Double_t Integral_fIncohPsi2sToMu     = fIncohPsi2sToMu     -> Integral();
  // // Double_t Integral_fIncohPsi2sToMuPi   = fIncohPsi2sToMuPi   -> Integral();
  // // Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral();
  // // Double_t Integral_fTwoGammaToMuHigh   = fTwoGammaToMuHigh   -> Integral();
  // // fCohJpsiToMu        -> Scale( 1/Integral_fCohJpsiToMu        );
  // // fCohPsi2sToMu       -> Scale( 1/Integral_fCohPsi2sToMu       );
  // // fCohPsi2sToMuPi     -> Scale( 1/Integral_fCohPsi2sToMuPi     );
  // // fIncohJpsiToMu      -> Scale( 1/Integral_fIncohJpsiToMu      );
  // // fIncohPsi2sToMu     -> Scale( 1/Integral_fIncohPsi2sToMu     );
  // // fIncohPsi2sToMuPi   -> Scale( 1/Integral_fIncohPsi2sToMuPi   );
  // // fTwoGammaToMuMedium -> Scale( 1/Integral_fTwoGammaToMuMedium );
  // // fTwoGammaToMuHigh   -> Scale( 1/Integral_fTwoGammaToMuHigh   );
  // fCohJpsiToMu        -> Divide(fCohJpsiToMuGen);
  // fCohPsi2sToMu       -> Divide(fCohPsi2sToMuGen);
  // fCohPsi2sToMuPi     -> Divide(fCohPsi2sToMuPiGen);
  // fIncohJpsiToMu      -> Divide(fIncohJpsiToMuGen);
  // fIncohPsi2sToMu     -> Divide(fIncohPsi2sToMuGen);
  // fIncohPsi2sToMuPi   -> Divide(fIncohPsi2sToMuPiGen);
  // fTwoGammaToMuMedium -> Divide(fTwoGammaToMuMediumGen);
  // fTwoGammaToMuHigh   -> Divide(fTwoGammaToMuHighGen);
  //
  //
  // fCohJpsiToMu        -> SetLineColor(kRed);
  // fCohPsi2sToMu       -> SetLineColor(kMagenta);
  // fCohPsi2sToMuPi     -> SetLineColor(kYellow+1);
  // fIncohJpsiToMu      -> SetLineColor(kCyan);
  // fIncohPsi2sToMu     -> SetLineColor(kYellow);
  // fIncohPsi2sToMuPi   -> SetLineColor(kBlue+2);
  // fTwoGammaToMuMedium -> SetLineColor(kGreen);
  // fTwoGammaToMuHigh   -> SetLineColor(kBlue+3);
  // // fFitInvMassResults  -> SetLineColor(kGreen+1);
  // // fCohJpsiToMu        -> SetLineWidth(3);
  // // fCohPsi2sToMu       -> SetLineWidth(3);
  // // fCohPsi2sToMuPi     -> SetLineWidth(3);
  // // fIncohJpsiToMu      -> SetLineWidth(3);
  // // fIncohPsi2sToMu     -> SetLineWidth(3);
  // // fIncohPsi2sToMuPi   -> SetLineWidth(3);
  // // fTwoGammaToMuMedium -> SetLineWidth(3);
  // // fTwoGammaToMuHigh   -> SetLineWidth(3);
  //
  // fCohJpsiToMu        -> Draw();
  // fCohPsi2sToMu       -> Draw("SAME");
  // // fCohPsi2sToMuPi     -> Draw("SAME");
  // fIncohJpsiToMu      -> Draw("SAME");
  // fIncohPsi2sToMu     -> Draw("SAME");
  // // fIncohPsi2sToMuPi   -> Draw("SAME");
  // fTwoGammaToMuMedium -> Draw("SAME");
  // fTwoGammaToMuHigh   -> Draw("SAME");
  //
  //
  // TLegend* l = new TLegend(0.2,0.55,0.5,0.85);
  // l->SetMargin(0.1);
  // l->SetBorderSize(0);
  // l->AddEntry(  fCohJpsiToMu,        "fCohJpsiToMu"       );
  // l->AddEntry(  fCohPsi2sToMu,       "fCohPsi2sToMu"      );
  // // l->AddEntry(  fCohPsi2sToMuPi,     "fCohPsi2sToMuPi"    );
  // l->AddEntry(  fIncohJpsiToMu,      "fIncohJpsiToMu"     );
  // l->AddEntry(  fIncohPsi2sToMu,     "fIncohPsi2sToMu"    );
  // // l->AddEntry(  fIncohPsi2sToMuPi,   "fIncohPsi2sToMuPi"  );
  // l->AddEntry(  fTwoGammaToMuHigh,   "fTwoGammaToMuHigh"  );
  // l->AddEntry(  fTwoGammaToMuMedium, "fTwoGammaToMuMedium");
  // l->Draw("same");
  //
  // // gPad->SaveAs("pngResults/Coh0N0Nleg.png",      "RECREATE");


}
