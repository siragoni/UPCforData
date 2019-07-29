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
TH1F*     fCohJpsiToMu;
TH1F*     fCohPsi2sToMu;
TH1F*     fCohPsi2sToMuPi;
TH1F*     fIncohJpsiToMu;
TH1F*     fIncohPsi2sToMu;
TH1F*     fIncohPsi2sToMuPi;
TH1F*     fTwoGammaToMuMedium;
TH1F*     fTwoGammaToMuHigh;
TH1F*     fCohJpsiToMuGen;
TH1F*     fCohPsi2sToMuGen;
TH1F*     fCohPsi2sToMuPiGen;
TH1F*     fIncohJpsiToMuGen;
TH1F*     fIncohPsi2sToMuGen;
TH1F*     fIncohPsi2sToMuPiGen;
TH1F*     fTwoGammaToMuMediumGen;
TH1F*     fTwoGammaToMuHighGen;
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void ShowAcceptancesTogether(){

  TFile* fileMC[8];
  fileMC[0] = new TFile("MCtrainResults/2019-06-24/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2019-06-24/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2019-06-24/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2019-06-24/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2019-06-24/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2019-06-24/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2019-06-24/kTwoGammaToMuHigh/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2019-06-24/kTwoGammaToMuMedium/AnalysisResults.root");
  TDirectory* dirMC[8];
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory] = fileMC[iDirectory]->GetDirectory("MyTask");
  }
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listingsMC[8];
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
  fCohJpsiToMu           = (TH1F*)listingsMC[0]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fCohPsi2sToMu          = (TH1F*)listingsMC[1]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fCohPsi2sToMuPi        = (TH1F*)listingsMC[2]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fIncohJpsiToMu         = (TH1F*)listingsMC[3]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fIncohPsi2sToMu        = (TH1F*)listingsMC[4]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fIncohPsi2sToMuPi      = (TH1F*)listingsMC[5]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fTwoGammaToMuMedium    = (TH1F*)listingsMC[6]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fTwoGammaToMuHigh      = (TH1F*)listingsMC[7]->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fCohJpsiToMuGen        = (TH1F*)listingsMC[0]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fCohPsi2sToMuGen       = (TH1F*)listingsMC[1]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fCohPsi2sToMuPiGen     = (TH1F*)listingsMC[2]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fIncohJpsiToMuGen      = (TH1F*)listingsMC[3]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fIncohPsi2sToMuGen     = (TH1F*)listingsMC[4]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fIncohPsi2sToMuPiGen   = (TH1F*)listingsMC[5]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fTwoGammaToMuMediumGen = (TH1F*)listingsMC[6]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fTwoGammaToMuHighGen   = (TH1F*)listingsMC[7]->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  /* - Rebin
     -
   */
  // fCohJpsiToMu           -> Rebin(5);
  // fCohPsi2sToMu          -> Rebin(5);
  // fCohPsi2sToMuPi        -> Rebin(5);
  // fIncohJpsiToMu         -> Rebin(5);
  // fIncohPsi2sToMu        -> Rebin(5);
  // fIncohPsi2sToMuPi      -> Rebin(5);
  // fTwoGammaToMuMedium    -> Rebin(5);
  // fTwoGammaToMuHigh      -> Rebin(5);
  // fCohJpsiToMuGen        -> Rebin(5);
  // fCohPsi2sToMuGen       -> Rebin(5);
  // fCohPsi2sToMuPiGen     -> Rebin(5);
  // fIncohJpsiToMuGen      -> Rebin(5);
  // fIncohPsi2sToMuGen     -> Rebin(5);
  // fIncohPsi2sToMuPiGen   -> Rebin(5);
  // fTwoGammaToMuMediumGen -> Rebin(5);
  // fTwoGammaToMuHighGen   -> Rebin(5);
  // fCohJpsiToMu           -> Rebin(20);
  // fCohPsi2sToMu          -> Rebin(20);
  // fCohPsi2sToMuPi        -> Rebin(20);
  // fIncohJpsiToMu         -> Rebin(20);
  // fIncohPsi2sToMu        -> Rebin(20);
  // fIncohPsi2sToMuPi      -> Rebin(20);
  // fTwoGammaToMuMedium    -> Rebin(20);
  // fTwoGammaToMuHigh      -> Rebin(20);
  // fCohJpsiToMuGen        -> Rebin(20);
  // fCohPsi2sToMuGen       -> Rebin(20);
  // fCohPsi2sToMuPiGen     -> Rebin(20);
  // fIncohJpsiToMuGen      -> Rebin(20);
  // fIncohPsi2sToMuGen     -> Rebin(20);
  // fIncohPsi2sToMuPiGen   -> Rebin(20);
  // fTwoGammaToMuMediumGen -> Rebin(20);
  // fTwoGammaToMuHighGen   -> Rebin(20);

  /* - Firstly we normalize the histograms.
     - Remember to always Sumw2()!!
     -
   */
  fCohJpsiToMu           -> Sumw2();
  fCohPsi2sToMu          -> Sumw2();
  fCohPsi2sToMuPi        -> Sumw2();
  fIncohJpsiToMu         -> Sumw2();
  fIncohPsi2sToMu        -> Sumw2();
  fIncohPsi2sToMuPi      -> Sumw2();
  fTwoGammaToMuMedium    -> Sumw2();
  fTwoGammaToMuHigh      -> Sumw2();
  fCohJpsiToMuGen        -> Sumw2();
  fCohPsi2sToMuGen       -> Sumw2();
  fCohPsi2sToMuPiGen     -> Sumw2();
  fIncohJpsiToMuGen      -> Sumw2();
  fIncohPsi2sToMuGen     -> Sumw2();
  fIncohPsi2sToMuPiGen   -> Sumw2();
  fTwoGammaToMuMediumGen -> Sumw2();
  fTwoGammaToMuHighGen   -> Sumw2();
  // Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral();
  // Double_t Integral_fCohPsi2sToMu       = fCohPsi2sToMu       -> Integral();
  // Double_t Integral_fCohPsi2sToMuPi     = fCohPsi2sToMuPi     -> Integral();
  // Double_t Integral_fIncohJpsiToMu      = fIncohJpsiToMu      -> Integral();
  // Double_t Integral_fIncohPsi2sToMu     = fIncohPsi2sToMu     -> Integral();
  // Double_t Integral_fIncohPsi2sToMuPi   = fIncohPsi2sToMuPi   -> Integral();
  // Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral();
  // Double_t Integral_fTwoGammaToMuHigh   = fTwoGammaToMuHigh   -> Integral();
  // fCohJpsiToMu        -> Scale( 1/Integral_fCohJpsiToMu        );
  // fCohPsi2sToMu       -> Scale( 1/Integral_fCohPsi2sToMu       );
  // fCohPsi2sToMuPi     -> Scale( 1/Integral_fCohPsi2sToMuPi     );
  // fIncohJpsiToMu      -> Scale( 1/Integral_fIncohJpsiToMu      );
  // fIncohPsi2sToMu     -> Scale( 1/Integral_fIncohPsi2sToMu     );
  // fIncohPsi2sToMuPi   -> Scale( 1/Integral_fIncohPsi2sToMuPi   );
  // fTwoGammaToMuMedium -> Scale( 1/Integral_fTwoGammaToMuMedium );
  // fTwoGammaToMuHigh   -> Scale( 1/Integral_fTwoGammaToMuHigh   );
  fCohJpsiToMu        -> Divide(fCohJpsiToMuGen);
  fCohPsi2sToMu       -> Divide(fCohPsi2sToMuGen);
  fCohPsi2sToMuPi     -> Divide(fCohPsi2sToMuPiGen);
  fIncohJpsiToMu      -> Divide(fIncohJpsiToMuGen);
  fIncohPsi2sToMu     -> Divide(fIncohPsi2sToMuGen);
  fIncohPsi2sToMuPi   -> Divide(fIncohPsi2sToMuPiGen);
  fTwoGammaToMuMedium -> Divide(fTwoGammaToMuMediumGen);
  fTwoGammaToMuHigh   -> Divide(fTwoGammaToMuHighGen);


  fCohJpsiToMu        -> SetLineColor(kRed);
  fCohPsi2sToMu       -> SetLineColor(kMagenta);
  fCohPsi2sToMuPi     -> SetLineColor(kYellow+1);
  fIncohJpsiToMu      -> SetLineColor(kCyan);
  fIncohPsi2sToMu     -> SetLineColor(kYellow);
  fIncohPsi2sToMuPi   -> SetLineColor(kBlue+2);
  fTwoGammaToMuMedium -> SetLineColor(kGreen);
  fTwoGammaToMuHigh   -> SetLineColor(kBlue+3);
  // fFitInvMassResults  -> SetLineColor(kGreen+1);
  // fCohJpsiToMu        -> SetLineWidth(3);
  // fCohPsi2sToMu       -> SetLineWidth(3);
  // fCohPsi2sToMuPi     -> SetLineWidth(3);
  // fIncohJpsiToMu      -> SetLineWidth(3);
  // fIncohPsi2sToMu     -> SetLineWidth(3);
  // fIncohPsi2sToMuPi   -> SetLineWidth(3);
  // fTwoGammaToMuMedium -> SetLineWidth(3);
  // fTwoGammaToMuHigh   -> SetLineWidth(3);

  fCohJpsiToMu        -> Draw();
  fCohPsi2sToMu       -> Draw("SAME");
  // fCohPsi2sToMuPi     -> Draw("SAME");
  fIncohJpsiToMu      -> Draw("SAME");
  fIncohPsi2sToMu     -> Draw("SAME");
  // fIncohPsi2sToMuPi   -> Draw("SAME");
  fTwoGammaToMuMedium -> Draw("SAME");
  fTwoGammaToMuHigh   -> Draw("SAME");


  TLegend* l = new TLegend(0.2,0.55,0.5,0.85);
  l->SetMargin(0.1);
  l->SetBorderSize(0);
  l->AddEntry(  fCohJpsiToMu,        "fCohJpsiToMu"       );
  l->AddEntry(  fCohPsi2sToMu,       "fCohPsi2sToMu"      );
  // l->AddEntry(  fCohPsi2sToMuPi,     "fCohPsi2sToMuPi"    );
  l->AddEntry(  fIncohJpsiToMu,      "fIncohJpsiToMu"     );
  l->AddEntry(  fIncohPsi2sToMu,     "fIncohPsi2sToMu"    );
  // l->AddEntry(  fIncohPsi2sToMuPi,   "fIncohPsi2sToMuPi"  );
  l->AddEntry(  fTwoGammaToMuHigh,   "fTwoGammaToMuHigh"  );
  l->AddEntry(  fTwoGammaToMuMedium, "fTwoGammaToMuMedium");
  l->Draw("same");

  // gPad->SaveAs("pngResults/Coh0N0Nleg.png",      "RECREATE");


}
