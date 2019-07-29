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
TH1F*       fCohJpsiToMu;
TH1F*       fCohPsi2sToMu;
TH1F*       fCohPsi2sToMuPi;
TH1F*       fIncohJpsiToMu;
TH1F*       fIncohPsi2sToMu;
TH1F*       fIncohPsi2sToMuPi;
TH1F*       fTwoGammaToMuMedium;
TH1F*       fTwoGammaToMuHigh;
TH1F*       fHighPtTail;
Double_t    ptrSelectionFlag = 0;
Double_t    JPsiPeakValue    = 0;
Double_t    JPsiPeakValueErr = 0;
Double_t    BkgValue         = 0;
Double_t    BkgValueError    = 0;
TFile*      fileList;
TDirectory* dir;
TFile*      fileMC[8];
TDirectory* dirMC[8];
TList*      listings;
TList*      listingsMC[8];


//_____________________________________________________________________________
/* - Fit function for the Pt plots.
 * - I am using simple ROOT to make gaussian fits to the plot.
 */
Double_t fInvariantMass(Double_t* x,Double_t* par)
{
  /* - Par 0, 1, 2:   coherent.
     - Par 3, 4, 5:   incoherent.
     - Par 6      :   gamma+gamma.
     -
   */
  Double_t val = 0;
  if ( ptrSelectionFlag == 2 ) {
    val += par[0]* ( fIncohJpsiToMu     ->Interpolate(x[0]) );
    val += par[1]* ( fIncohPsi2sToMu    ->Interpolate(x[0]) );
  } else {
    val += par[0]* ( fCohJpsiToMu       ->Interpolate(x[0]) );
    val += par[1]* ( fCohPsi2sToMu      ->Interpolate(x[0]) );
  }
  // val   += par[2]* ( fTwoGammaToMuMedium->Interpolate(x[0]) );
  val   += par[2]* ( fTwoGammaToMuHigh->Interpolate(x[0]) );

  return val;
}
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitJPsiTemplate(const char* AnalysisName, const int selectionFlag){
  /* - There are three cases for the selectionFlag:
     - 1) = 0 ; this implies the traditional pt-integrated plot;
     - 2) = 1 ; this is instead the coherent component;
     - 3) = 2 ; this is the incoherent component;
     - 4) = 3 ; ******************* ;
     -
   */
  ptrSelectionFlag = 0;
  TF1* fFitInvMass = new TF1( "fFitInvMass", fInvariantMass , 2., 6.5, 3 );
  fFitInvMass->SetParameter(0, 1);
  fFitInvMass->SetParameter(1, 1);
  fFitInvMass->SetParameter(2, 1);
  fFitInvMass->SetNpx( fCohJpsiToMu->GetNbinsX() );
  // fFitInvMassResults = (TH1F*) fFitInvMass->GetHistogram()->Clone("fFitInvMassResults");
  // for (Int_t ibin=1; ibin<=fFitInvMassResults->GetNbinsX(); ibin++) {
  //   fFitInvMassResults->SetBinError(ibin,0);
  // }

  TH1F *fInvariantMassDistributionH = 0x0;
  fInvariantMassDistributionH = (TH1F*)listings->FindObject( Form("fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH_%d", selectionFlag) );
  fInvariantMassDistributionH->Rebin(5);
  fInvariantMassDistributionH->Draw("PE");

  fInvariantMassDistributionH->SetLineColor(kBlue);
  fInvariantMassDistributionH->SetLineStyle(kSolid);
  fInvariantMassDistributionH->SetLineWidth(3);
  fInvariantMassDistributionH->SetMarkerStyle(kFullCircle);
  fInvariantMassDistributionH->SetMarkerColor(kBlue);
  fInvariantMassDistributionH->SetMarkerSize(1);
  fInvariantMassDistributionH->GetXaxis()->SetTitle("M_{#mu#mu} [GeV/#it{c}^{2}");
  fInvariantMassDistributionH->GetYaxis()->SetTitle( Form( "Counts / (%.3f GeV/#it{c})",
                                                          fInvariantMassDistributionH->GetXaxis()->GetBinWidth(1)
                                                        )
                                                    );
  fInvariantMassDistributionH->SetTitle("");
  fInvariantMassDistributionH->Fit( fFitInvMass,"","", 2.2, 6. );
  TCanvas* PtDistrCanvas = new TCanvas( "InvariantMassDimuonFit", "InvariantMassDimuonFit", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gPad->SetTitle(  Form(  ";M_{#mu#mu} (GeV/c^{2});Counts / (%.0f MeV/c^{2})",
                           fInvariantMassDistributionH->GetXaxis()->GetBinWidth(1)*1000.  )  );
  /* - Beautifying is starting now.
     -
   */
  fInvariantMassDistributionH->GetXaxis()->SetTitleOffset(1.25);
  // fInvariantMassDistributionH->GetYaxis()->SetTitleOffset(1.25);
  fInvariantMassDistributionH->GetYaxis()->SetTitleOffset(1.45);
  fInvariantMassDistributionH->GetXaxis()->SetTitleSize(0.045);
  fInvariantMassDistributionH->GetYaxis()->SetTitleSize(0.045);
  fInvariantMassDistributionH->GetXaxis()->SetLabelSize(0.045);
  fInvariantMassDistributionH->GetYaxis()->SetLabelSize(0.045);
  fInvariantMassDistributionH->GetXaxis()->SetTitleFont(42);
  fInvariantMassDistributionH->GetYaxis()->SetTitleFont(42);
  fInvariantMassDistributionH->GetXaxis()->SetLabelFont(42);
  fInvariantMassDistributionH->GetYaxis()->SetLabelFont(42);
  fInvariantMassDistributionH->GetXaxis()->SetNdivisions(408);
  fInvariantMassDistributionH->GetYaxis()->SetRangeUser(0.0000000000001, fInvariantMassDistributionH->GetMaximum()*1.3);
  fInvariantMassDistributionH->GetXaxis()->SetRangeUser(2, 6);
  // gPad ->SetLogy();
  fInvariantMassDistributionH->Draw("PEsame");
  fCohJpsiToMu        -> SetLineColor(kRed);
  fCohPsi2sToMu       -> SetLineColor(kMagenta);
  fCohPsi2sToMuPi     -> SetLineColor(kYellow+1);
  fIncohJpsiToMu      -> SetLineColor(kCyan);
  fIncohPsi2sToMu     -> SetLineColor(kYellow);
  fIncohPsi2sToMuPi   -> SetLineColor(kBlue+2);
  fTwoGammaToMuMedium -> SetLineColor(kGreen);
  fTwoGammaToMuHigh   -> SetLineColor(kBlue+3);
  // fFitInvMassResults  -> SetLineColor(kGreen+1);
  fCohJpsiToMu        -> SetLineWidth(3);
  fCohPsi2sToMu       -> SetLineWidth(3);
  fCohPsi2sToMuPi     -> SetLineWidth(3);
  fIncohJpsiToMu      -> SetLineWidth(3);
  fIncohPsi2sToMu     -> SetLineWidth(3);
  fIncohPsi2sToMuPi   -> SetLineWidth(3);
  fTwoGammaToMuMedium -> SetLineWidth(3);
  fTwoGammaToMuHigh   -> SetLineWidth(3);
  // fFitInvMassResults  -> SetLineWidth(3);
  TH1F* fCohJpsiToMuC        = (TH1F*) fCohJpsiToMu        -> Clone("fCohJpsiToMuC");
  TH1F* fCohPsi2sToMuC       = (TH1F*) fCohPsi2sToMu       -> Clone("fCohPsi2sToMuC");
  TH1F* fCohPsi2sToMuPiC     = (TH1F*) fCohPsi2sToMuPi     -> Clone("fCohPsi2sToMuPiC");
  TH1F* fIncohJpsiToMuC      = (TH1F*) fIncohJpsiToMu      -> Clone("fIncohJpsiToMuC");
  TH1F* fIncohPsi2sToMuC     = (TH1F*) fIncohPsi2sToMu     -> Clone("fIncohPsi2sToMuC");
  TH1F* fIncohPsi2sToMuPiC   = (TH1F*) fIncohPsi2sToMuPi   -> Clone("fIncohPsi2sToMuPiC");
  TH1F* fTwoGammaToMuMediumC = (TH1F*) fTwoGammaToMuMedium -> Clone("fTwoGammaToMuMediumC");
  TH1F* fTwoGammaToMuHighC   = (TH1F*) fTwoGammaToMuHigh   -> Clone("fTwoGammaToMuHighC");
  // TH1F* fFitInvMassResultsC  = (TH1F*) fFitInvMassResults  -> Clone("fFitInvMassResultsC");
  fCohJpsiToMuC        -> Scale( fFitInvMass->GetParameter(0) );
  fCohPsi2sToMuC       -> Scale( fFitInvMass->GetParameter(1) );
  // fCohPsi2sToMuPiC     -> Scale( fFitInvMass->GetParameter(2) );
  fIncohJpsiToMuC      -> Scale( fFitInvMass->GetParameter(0) );
  fIncohPsi2sToMuC     -> Scale( fFitInvMass->GetParameter(1) );
  // fIncohPsi2sToMuPiC   -> Scale( fFitInvMass->GetParameter(5) );
  fTwoGammaToMuMediumC -> Scale( fFitInvMass->GetParameter(2) );
  fTwoGammaToMuHighC   -> Scale( fFitInvMass->GetParameter(2) );
  // fFitInvMassResultsC  -> Scale( fFitInvMass->GetParameter(7) );
  if ( ptrSelectionFlag == 2 ) {
    fIncohJpsiToMuC      -> Draw("lsame");
    fIncohPsi2sToMuC     -> Draw("lsame");
  } else {
    fCohJpsiToMuC        -> Draw("lsame");
    fCohPsi2sToMuC       -> Draw("lsame");
  }
  // fTwoGammaToMuMediumC -> Draw("same");
  fTwoGammaToMuHighC   -> Draw("Csame");
  // fFitInvMassResultsC  -> Draw("Esame");


  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, LHC18q+LHC18r");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T} < 0.3 GeV/#it{c}");
  latex->DrawLatex(0.55,0.78,Form("SignalPhi_%d", selectionFlag));
  latex->DrawLatex(0.55,0.72,Form("%.1f < y < %.1f",-4.0,-2.5));

  /* - This is the part where we obtain the actual number of J/Psi, PsiPrime
     - and the background. This is still Kay's original code. I will modify it.
     - Hopefully if everything goes alright, I should have been able to complete
     - it by the time you are reading this.
     -
   */
  Double_t numberOfTotalJPsi  = 0;
  Double_t numberOfTotalPsi2s = 0;
  Double_t numberOfTotalBkg   = 0;
  if ( ptrSelectionFlag == 2 ) {
    numberOfTotalJPsi  = fIncohJpsiToMuC  -> Integral();
    numberOfTotalPsi2s = fIncohPsi2sToMuC -> Integral();
  } else {
    numberOfTotalJPsi  = fCohJpsiToMuC    -> Integral();
    numberOfTotalPsi2s = fCohPsi2sToMuC   -> Integral();
  }
  numberOfTotalBkg = fTwoGammaToMuMediumC -> Integral();
  latex->DrawLatex(0.55,0.66,Form("N_{J/#psi} = %.0f #pm %.0f",        numberOfTotalJPsi,  fFitInvMass->GetParError(0) ) );
  latex->DrawLatex(0.55,0.60,Form("N_{#psi(2S)} = %.0f #pm %.0f",      numberOfTotalPsi2s, fFitInvMass->GetParError(1) ) );
  latex->DrawLatex(0.55,0.54,Form("N_{#gamma#gamma} = %.0f #pm %.0f",  numberOfTotalBkg,   fFitInvMass->GetParError(2) ) );
  // latex->DrawLatex(0.55,0.54,Form("#sigma_{J/#psi} = %.0f #pm %.0f MeV/c^{2}", sigma.getVal()*1000,    sigma.getError()*1000));
  // latex->DrawLatex(0.55,0.48,Form("#sigma_{#psi(2S)} = %.0f MeV/c^{2} fixed",  sigma2.getVal()*1000));

  /* - This part concerns the background of the two signals.
     - Here, we extrapolate the background and compute the significance maybe?
     -
   */
  Double_t JPsiPeakBkg        = 0;
  Double_t Psi2JPsiPeakBkg    = 0;
  Double_t JPsiPeakSignal     = 0;
  Double_t Psi2JPsiPeakSignal = 0;
  if ( ptrSelectionFlag == 2 ) {
    JPsiPeakSignal     = fIncohJpsiToMuC  -> Integral( fIncohJpsiToMuC ->GetXaxis()->FindBin(2.75), fIncohJpsiToMuC ->GetXaxis()->FindBin(3.45) );
    Psi2JPsiPeakSignal = fIncohPsi2sToMuC -> Integral( fIncohPsi2sToMuC->GetXaxis()->FindBin(3.45), fIncohPsi2sToMuC->GetXaxis()->FindBin(3.90) );
  } else {
    JPsiPeakSignal     = fCohJpsiToMuC    -> Integral( fCohJpsiToMuC ->GetXaxis()->FindBin(2.75), fCohJpsiToMuC ->GetXaxis()->FindBin(3.45) );
    Psi2JPsiPeakSignal = fCohPsi2sToMuC   -> Integral( fCohPsi2sToMuC->GetXaxis()->FindBin(3.45), fCohPsi2sToMuC->GetXaxis()->FindBin(3.90) );
  }
  JPsiPeakBkg     = fTwoGammaToMuHighC -> Integral( fTwoGammaToMuHighC->GetXaxis()->FindBin(2.75), fTwoGammaToMuHighC->GetXaxis()->FindBin(3.45) );
  Psi2JPsiPeakBkg = fTwoGammaToMuHighC -> Integral( fTwoGammaToMuHighC->GetXaxis()->FindBin(3.45), fTwoGammaToMuHighC->GetXaxis()->FindBin(3.90) );
  latex->DrawLatex(0.55,0.42,Form("N_{BG J/#psi} = %.0f #pm %.0f",   JPsiPeakBkg,     JPsiPeakBkg     * fFitInvMass->GetParError(2) / numberOfTotalJPsi ));
  latex->DrawLatex(0.55,0.36,Form("N_{BG #psi(2s)} = %.0f #pm %.0f", Psi2JPsiPeakBkg, Psi2JPsiPeakBkg * fFitInvMass->GetParError(2) / numberOfTotalPsi2s));
  latex->DrawLatex(0.55,0.18,Form("      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                     fFitInvMass->GetChisquare(),
                                     fFitInvMass->GetNDF(),
                                     fFitInvMass->GetChisquare()/fFitInvMass->GetNDF()
                                     )
                                    );


  //Signal to background ratios
  Double_t errSB1 = sqrt(   fFitInvMass->GetParError(0) * fFitInvMass->GetParError(0)/( numberOfTotalJPsi  *numberOfTotalJPsi )  +
                          ( fFitInvMass->GetParError(2) * fFitInvMass->GetParError(2)/( JPsiPeakBkg        * JPsiPeakBkg)) );
  Double_t errSB2 = sqrt(   fFitInvMass->GetParError(1) * fFitInvMass->GetParError(1)/( numberOfTotalPsi2s * numberOfTotalPsi2s )  +
                            fFitInvMass->GetParError(2) * fFitInvMass->GetParError(2)/( Psi2JPsiPeakBkg    * Psi2JPsiPeakBkg ) );
  latex->DrawLatex(0.55,0.30,Form("S/B for J/#psi = %.2f #pm %.2f",   JPsiPeakSignal/(Double_t)JPsiPeakBkg,         errSB1*JPsiPeakSignal/(Double_t)JPsiPeakBkg));
  latex->DrawLatex(0.55,0.24,Form("S/B for #psi(2s) = %.2f #pm %.2f", Psi2JPsiPeakSignal/(Double_t)Psi2JPsiPeakBkg, errSB2*Psi2JPsiPeakSignal/(Double_t)Psi2JPsiPeakBkg));
  gPad ->SaveAs(Form("pngResults/fitSignalExtractionPhi_%d.png", selectionFlag),     "RECREATE");



  JPsiPeakValue    = fFitInvMass->GetParameter(0);
  JPsiPeakValueErr = fFitInvMass->GetParError(0);
  BkgValue         = fFitInvMass->GetParameter(2);
  BkgValueError    = fFitInvMass->GetParError(2);

  // err = errorTotalBackGroundJPsiPeak;

}
//_____________________________________________________________________________
/* - Here I create the new TH2 for the after the signal extraction.
 * - Basically I run the fit function many times and then I memorise
 * - the values each time. After that I fill with a setbincontent
 * - and a setbinerror the
 * -
 */
void CreateTH2(const char* AnalysisName){

  fileList = new TFile(AnalysisName);
  dir = fileList->GetDirectory("MyTask");
  fileMC[0] = new TFile("MCtrainResults/2019-05-17/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2019-05-17/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2019-05-17/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2019-05-17/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2019-05-17/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2019-05-17/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2019-05-17/kTwoGammaToMuHigh/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2019-05-17/kTwoGammaToMuMedium/AnalysisResults.root");
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory] = fileMC[iDirectory]->GetDirectory("MyTask");
  }
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  dir->GetObject("MyOutputContainer", listings);
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
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
  fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fInvariantMassDistributionH");
  fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fInvariantMassDistributionH");
  fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fInvariantMassDistributionH");
  fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fInvariantMassDistributionH");
  fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fInvariantMassDistributionH");
  fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fInvariantMassDistributionH");
  fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fInvariantMassDistributionH");
  fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fInvariantMassDistributionH");

  fCohJpsiToMu        -> Rebin(5);
  fCohPsi2sToMu       -> Rebin(5);
  fCohPsi2sToMuPi     -> Rebin(5);
  fIncohJpsiToMu      -> Rebin(5);
  fIncohPsi2sToMu     -> Rebin(5);
  fIncohPsi2sToMuPi   -> Rebin(5);
  fTwoGammaToMuMedium -> Rebin(5);
  fTwoGammaToMuHigh   -> Rebin(5);
  /* - Firstly we normalize the histograms.
     - Remember to always Sumw2()!!
     -
   */
  fCohJpsiToMu        -> Sumw2();
  fCohPsi2sToMu       -> Sumw2();
  fCohPsi2sToMuPi     -> Sumw2();
  fIncohJpsiToMu      -> Sumw2();
  fIncohPsi2sToMu     -> Sumw2();
  fIncohPsi2sToMuPi   -> Sumw2();
  fTwoGammaToMuMedium -> Sumw2();
  fTwoGammaToMuHigh   -> Sumw2();
  Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral();
  Double_t Integral_fCohPsi2sToMu       = fCohPsi2sToMu       -> Integral();
  Double_t Integral_fCohPsi2sToMuPi     = fCohPsi2sToMuPi     -> Integral();
  Double_t Integral_fIncohJpsiToMu      = fIncohJpsiToMu      -> Integral();
  Double_t Integral_fIncohPsi2sToMu     = fIncohPsi2sToMu     -> Integral();
  Double_t Integral_fIncohPsi2sToMuPi   = fIncohPsi2sToMuPi   -> Integral();
  Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral();
  Double_t Integral_fTwoGammaToMuHigh   = fTwoGammaToMuHigh   -> Integral();
  fCohJpsiToMu        -> Scale( 1/Integral_fCohJpsiToMu        );
  fCohPsi2sToMu       -> Scale( 1/Integral_fCohPsi2sToMu       );
  fCohPsi2sToMuPi     -> Scale( 1/Integral_fCohPsi2sToMuPi     );
  fIncohJpsiToMu      -> Scale( 1/Integral_fIncohJpsiToMu      );
  fIncohPsi2sToMu     -> Scale( 1/Integral_fIncohPsi2sToMu     );
  fIncohPsi2sToMuPi   -> Scale( 1/Integral_fIncohPsi2sToMuPi   );
  fTwoGammaToMuMedium -> Scale( 1/Integral_fTwoGammaToMuMedium );
  fTwoGammaToMuHigh   -> Scale( 1/Integral_fTwoGammaToMuHigh   );

  // Double_t CosThetaBinCenters[] = { -0.9, -0.7, -0.5, -0.3, -0.1,
  //                                    0.1,  0.3,  0.5,  0.7,  0.9 };
  // Double_t PhiBinCenters[]      = { -2.826, -2.198, -1.57, -0.942, -0.314,
  //                                    0.314,  0.942,  1.57,  2.198,  2.826};
  Double_t CosThetaBinCenters[40];
  Double_t PhiBinCenters     [50];
  for( Int_t iFillingCenters = 0; iFillingCenters < 50; iFillingCenters++ ){
    if( iFillingCenters < 40 ) {
      CosThetaBinCenters[iFillingCenters] = -0.975 + (Double_t)iFillingCenters * 0.025;
    }
    PhiBinCenters[iFillingCenters] = -3.14 + (Double_t)iFillingCenters * 0.1256;
  }
  for( Int_t looping = 0; looping < 40; looping++ ) {cout << "CosThetaBinCenters  " << CosThetaBinCenters[looping] << endl; }
  for( Int_t looping = 0; looping < 50; looping++ ) {cout << "PhiBinCenters       " << PhiBinCenters[looping]      << endl; }

  TH1F* PhiAfterSignalExtractionH =
            new TH1F( "PhiAfterSignalExtractionH",
                      "PhiAfterSignalExtractionH",
                      50, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      );
  TH1F* PhiGammaGammaH =
            new TH1F( "PhiGammaGammaH",
                      "PhiGammaGammaH",
                      50, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      );

  TH1F* PhiAfterSignalExtractionErrorsH =
            new TH1F( "PhiAfterSignalExtractionErrorsH",
                      "PhiAfterSignalExtractionErrorsH",
                      50, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      );
  TH1F* PhiGammaGammaErrorsH =
            new TH1F( "PhiGammaGammaErrorsH",
                      "PhiGammaGammaErrorsH",
                      50, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      );

  // for (size_t iCosThetaBins = 2; iCosThetaBins < 40; iCosThetaBins++) {
    for (size_t iPhiBins = 0; iPhiBins < 50; iPhiBins++) {
      JPsiPeakValue    = 0;
      JPsiPeakValueErr = 0;
      BkgValue         = 0;
      BkgValueError    = 0;
      fitJPsiTemplate(AnalysisName, iPhiBins);

      PhiAfterSignalExtractionH->Fill( -3.14 + 0.0628 + (Double_t)iPhiBins * 0.1256,
                                            // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                            // -2.826 + (Double_t)iPhiBins      * 0.628,
                                            // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                            JPsiPeakValue
                                            );
      PhiGammaGammaH->Fill( -3.14 + 0.0628  + (Double_t)iPhiBins * 0.1256,
                                 // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                 // -2.826 + (Double_t)iPhiBins      * 0.628,
                                 // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                 BkgValue
                                 );

      PhiAfterSignalExtractionErrorsH->Fill(  -3.14 + 0.0628  + (Double_t)iPhiBins * 0.1256,
                                                   // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                                   // -2.826 + (Double_t)iPhiBins      * 0.628,
                                                   // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                                   JPsiPeakValue
                                                   );
      PhiGammaGammaErrorsH->Fill(   -3.14 + 0.0628  + (Double_t)iPhiBins * 0.1256,
                                        // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                        // -2.826 + (Double_t)iPhiBins      * 0.628,
                                        // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                        BkgValue
                                        );
      PhiAfterSignalExtractionErrorsH->SetBinError(  // iCosThetaBins + 1 ,
                                                          iPhiBins      + 1 ,
                                                          JPsiPeakValueErr
                                                          );
      PhiGammaGammaErrorsH->SetBinError(  // iCosThetaBins + 1 ,
                                               iPhiBins      + 1 ,
                                               BkgValueError
                                               );

    }
  // }

  TFile f("pngResults/TH1signalPhiEX.root", "new");
  PhiAfterSignalExtractionH      ->Write();
  PhiGammaGammaH                 ->Write();
  PhiAfterSignalExtractionErrorsH->Write();
  PhiGammaGammaErrorsH           ->Write();
  f.Close();
}
