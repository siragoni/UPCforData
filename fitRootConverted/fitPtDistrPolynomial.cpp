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
TH1F*     fHighPtTail;
Double_t  ptrSelectionFlag = 0;
Double_t* CBparameters[9];
Int_t     counter = 0;



/* - Functions used for the fit!!
 * -
 */
TF1* fCohJpsiToMuPolynomial        = new TF1( "fCohJpsiToMuPolynomial",        "pol8", 0, 1.5 );
TF1* fCohPsi2sToMuPolynomial       = new TF1( "fCohPsi2sToMuPolynomial",       "pol8", 0, 1.5 );
TF1* fCohPsi2sToMuPiPolynomial     = new TF1( "fCohPsi2sToMuPiPolynomial",     "pol8", 0, 1.5 );
TF1* fIncohJpsiToMuPolynomial      = new TF1( "fIncohJpsiToMuPolynomial",      "pol8", 0, 1.5 );
TF1* fIncohPsi2sToMuPolynomial     = new TF1( "fIncohPsi2sToMuPolynomial",     "pol8", 0, 1.5 );
TF1* fIncohPsi2sToMuPiPolynomial   = new TF1( "fIncohPsi2sToMuPiPolynomial",   "pol8", 0, 1.5 );
TF1* fTwoGammaToMuMediumPolynomial = new TF1( "fTwoGammaToMuMediumPolynomial", "pol8", 0, 1.5 );
TF1* fTwoGammaToMuHighPolynomial   = new TF1( "fTwoGammaToMuHighPolynomial",   "pol8", 0, 1.5 );
TF1* fHighPtTailSameAsMacro        = new TF1( "fHighPtTailSameAsMacro",        "[0]*x*(1+[1]/[2]*x*x)^(-[2])",0.00000001,6);

//_____________________________________________________________________________
/* - My fit function, taken from the template macro.
 * -
 */
double fFitPtDistr(double *x, double *par){
  Double_t val = 0;
  val += par[0]* ( fCohJpsiToMuPolynomial->Eval(x[0]) );   //needed
  // // val += par[1]* ( fCohPsi2sToMuPolynomial   ->Eval(x[0]) );
  // val += par[0]* ( fCohPsi2sToMuPiPolynomial    ->Eval(x[0]) ) * par[3];   //needed
  // val += par[1]* ( fIncohJpsiToMuPolynomial     ->Eval(x[0]) );   //needed
  // // val += par[4]* ( fIncohPsi2sToMuPolynomial ->Eval(x[0]) );
  // val += par[1]* ( fIncohPsi2sToMuPiPolynomial  ->Eval(x[0]) ) * par[3];   //needed
  // val += par[2]* ( fTwoGammaToMuMediumPolynomial->Eval(x[0]) );   //needed
  // Double_t parHighPt[3];
  // for( Int_t i=0 ; i<3 ; i++ ) parHighPt[i]=par[i+4];
  // val += fHighPtTailSameAsMacro->EvalPar( x[0], parHighPt );   //needed

  return val;
}
//_____________________________________________________________________________
/* - Fit function for the single templates of the MC.
 * -
 */
void fPolynomialFitToMC(TH1F* histoToBeFit, Double_t* &bookKeeping)
{
  TF1* Polynomial8 = 0x0;
  if ( counter == 0 || counter == 7 ) {
    Polynomial8 = new TF1("Polynomial8","pol8",0,0.4);
    // Polynomial8 = new TF1("Polynomial8","pol8",0,0.5);
  }
  else if ( counter == 3 || counter == 5 ) {
    Polynomial8 = new TF1("Polynomial8","pol8",0,0.9);
  }
  else { Polynomial8 = new TF1("Polynomial8","pol8",0,0.6); }
  // TF1* Polynomial8 = new TF1("Polynomial8","pol8",1.8,7);
  // Polynomial8 ->SetParameter(0,1);
  // Polynomial8 ->SetParameter(3,1.08);
  // Polynomial8 ->SetParameter(4,3689197);
  // // Polynomial8 ->SetParameter(4,6);
  // Polynomial8 ->SetParLimits(4,1,99999999);
  // // Polynomial8 ->SetParLimits(4,1,12);
  // // Polynomial8 ->SetParLimits(4,12,99999999);
  // Polynomial8 ->SetParameter(2,0.090);
  // Polynomial8 ->SetParameter(1,3.1);
  for(Int_t i = 0; i < 9; i++){
    Polynomial8->SetParLimits(i,0,2000);
  }
  // if ( counter == 0 || counter == 7 || counter == 6 ) {
  //   Polynomial8->FixParameter(4,1000);
  //   Polynomial8->FixParameter(5,1000);
  //   Polynomial8->FixParameter(6,1000);
  //   Polynomial8->FixParameter(7,1000);
  //   Polynomial8->FixParameter(8,1000);
  // }

  Polynomial8 ->SetNpx(1000);
  TCanvas*    CanvasRec = new TCanvas( Form("MC_%d", counter), Form("MC_%d", counter), 900, 800 );
  Polynomial8 ->Draw();
  histoToBeFit->Fit(Polynomial8, "R");
  // Polynomial8       ->SetParameter(0,1/CBfit->Integral(2,15));
  bookKeeping = new Double_t[9];
  for(Int_t i = 0; i < 9; i++){
    bookKeeping[i] = Polynomial8->GetParameter(i);
  }
}
//_____________________________________________________________________________
/* - Fit function for the MC plots. => AUTO mode
 * -
 */
void fitJPsiTemplateMC(){
  TFile* fileMC[8];
  fileMC[0] = new TFile("MCtrainResults/2019-05-17/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2019-05-17/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2019-05-17/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2019-05-17/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2019-05-17/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2019-05-17/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2019-05-17/kTwoGammaToMuHigh/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2019-05-17/kTwoGammaToMuMedium/AnalysisResults.root");
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
   *     OBJ: TH1F	  fDimuonPtDistributionDataH	fDimuonPtDistributionDataH : 0 at: 0x5a3c720
   */
  fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fTemplatePtDistributionH");
  fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fTemplatePtDistributionH");
  fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fTemplatePtDistributionH");
  fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fTemplatePtDistributionH");
  fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fTemplatePtDistributionH");
  fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fTemplatePtDistributionH");
  fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fTemplatePtDistributionH");
  fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fTemplatePtDistributionH");
  /* - Rebin
     -
   */
  // fCohJpsiToMu        -> Rebin(5);
  // fCohPsi2sToMu       -> Rebin(5);
  // fCohPsi2sToMuPi     -> Rebin(5);
  // fIncohJpsiToMu      -> Rebin(5);
  // fIncohPsi2sToMu     -> Rebin(5);
  // fIncohPsi2sToMuPi   -> Rebin(5);
  // fTwoGammaToMuMedium -> Rebin(5);
  // fTwoGammaToMuHigh   -> Rebin(5);
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


  fTwoGammaToMuMedium->Rebin(5); // mmmh

  fPolynomialFitToMC(fCohJpsiToMu,        CBparameters[0]);
  counter++;
  fPolynomialFitToMC(fCohPsi2sToMu,       CBparameters[1]);
  counter++;
  fPolynomialFitToMC(fCohPsi2sToMuPi,     CBparameters[2]);
  counter++;
  fPolynomialFitToMC(fIncohJpsiToMu,      CBparameters[3]);
  counter++;
  fPolynomialFitToMC(fIncohPsi2sToMu,     CBparameters[4]);
  counter++;
  fPolynomialFitToMC(fIncohPsi2sToMuPi,   CBparameters[5]);
  counter++;
  fPolynomialFitToMC(fTwoGammaToMuMedium, CBparameters[6]);
  counter++;
  fPolynomialFitToMC(fTwoGammaToMuHigh,   CBparameters[7]);

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
  ptrSelectionFlag = selectionFlag;
  TFile* fileList = new TFile(AnalysisName);
  // TFile* fileList = new TFile("AnalysisResultsLHC18qr15o21052019noSPD.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);

  TH1F *fDimuonPtDistributionDataH = 0x0;
  if      ( selectionFlag == 0 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
  else if ( selectionFlag == 1 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAzeroH");
  else if ( selectionFlag == 2 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCzeroZNAanyH");
  else if ( selectionFlag == 3 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAzeroH");
  else if ( selectionFlag == 4 ) fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionZNCanyZNAanyH");
  else                           fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
  fDimuonPtDistributionDataH->Rebin(5);
  fDimuonPtDistributionDataH->Draw("PE");

  fDimuonPtDistributionDataH->SetLineColor(kBlue);
  fDimuonPtDistributionDataH->SetLineStyle(kSolid);
  fDimuonPtDistributionDataH->SetLineWidth(3);
  fDimuonPtDistributionDataH->SetMarkerStyle(kFullCircle);
  fDimuonPtDistributionDataH->SetMarkerColor(kBlue);
  fDimuonPtDistributionDataH->SetMarkerSize(1);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  fDimuonPtDistributionDataH->GetYaxis()->SetTitle( Form( "Counts / (%.3f GeV/#it{c})",
                                                          fDimuonPtDistributionDataH->GetXaxis()->GetBinWidth(1)
                                                        )
                                                    );
  fDimuonPtDistributionDataH->SetTitle("");


  fitJPsiTemplateMC();
  new TCanvas;
  fDimuonPtDistributionDataH->Draw();

  /* - High Pt-tail, with HERA's data.
     -
   */
  TF1* fModelForHighPtTail = new TF1("fModelForHighPtTail","[0]*x*(1+[1]/[2]*x*x)^(-[2])",1.,3);
  fModelForHighPtTail->SetParameter(0,1000);
  fModelForHighPtTail->SetParLimits(0,0,999999999999);
  // fModelForHighPtTail->SetParameter(1,debug==4 ? 1.25 : 1.);
  // fModelForHighPtTail->SetParameter(2,debug==4 ? 6.1 : 1.);
  fModelForHighPtTail->SetParameter(1, 1.6/*1.79*/);
  fModelForHighPtTail->SetParameter(2, 3.58);
  // fModelForHighPtTail->SetNpx( fDimuonPtDistributionDataH->GetNbinsX() );
  fModelForHighPtTail->SetNpx( 1000 );
  fModelForHighPtTail->Draw("same");
  // fHighPtTail = (TH1F*) fModelForHighPtTail->GetHistogram()->Clone("fHighPtTail");
  // for (Int_t ibin=1; ibin<=fHighPtTail->GetNbinsX(); ibin++) {
  //   fHighPtTail->SetBinError(ibin,0);
  // }
  fDimuonPtDistributionDataH->Fit( fModelForHighPtTail,"LR","", 1.5, 2. );
  // fDimuonPtDistributionDataH->Fit( fModelForHighPtTail );

  for(Int_t i = 0; i < 9; i++) fCohJpsiToMuPolynomial       ->SetParameter(i, CBparameters[0][i]);
  for(Int_t i = 0; i < 9; i++) fCohPsi2sToMuPolynomial      ->SetParameter(i, CBparameters[1][i]);
  for(Int_t i = 0; i < 9; i++) fCohPsi2sToMuPiPolynomial    ->SetParameter(i, CBparameters[2][i]);
  for(Int_t i = 0; i < 9; i++) fIncohJpsiToMuPolynomial     ->SetParameter(i, CBparameters[3][i]);
  for(Int_t i = 0; i < 9; i++) fIncohPsi2sToMuPolynomial    ->SetParameter(i, CBparameters[4][i]);
  for(Int_t i = 0; i < 9; i++) fIncohPsi2sToMuPiPolynomial  ->SetParameter(i, CBparameters[5][i]);
  for(Int_t i = 0; i < 9; i++) fTwoGammaToMuMediumPolynomial->SetParameter(i, CBparameters[6][i]);
  for(Int_t i = 0; i < 9; i++) fTwoGammaToMuHighPolynomial  ->SetParameter(i, CBparameters[7][i]);


  for(Int_t i = 0; i < 9; i++) cout << fCohJpsiToMuPolynomial->GetParameter(i) << endl << flush;
  for(Int_t i = 0; i < 9; i++) cout << fCohPsi2sToMuPolynomial->GetParameter(i) << endl << flush;
  for(Int_t i = 0; i < 9; i++) cout << fCohPsi2sToMuPiPolynomial->GetParameter(i) << endl << flush;
  for(Int_t i = 0; i < 9; i++) cout << fIncohJpsiToMuPolynomial->GetParameter(i) << endl << flush;
  for(Int_t i = 0; i < 9; i++) cout << fIncohPsi2sToMuPolynomial->GetParameter(i) << endl << flush;
  for(Int_t i = 0; i < 9; i++) cout << fIncohPsi2sToMuPiPolynomial->GetParameter(i) << endl << flush;
  for(Int_t i = 0; i < 9; i++) cout << fTwoGammaToMuMediumPolynomial->GetParameter(i) << endl << flush;
  for(Int_t i = 0; i < 9; i++) cout << fTwoGammaToMuHighPolynomial->GetParameter(i) << endl << flush;

  new TCanvas;
  fCohJpsiToMuPolynomial       ->Draw();
  new TCanvas;
  fCohPsi2sToMuPolynomial      ->Draw();
  new TCanvas;
  fCohPsi2sToMuPiPolynomial    ->Draw();
  new TCanvas;
  fIncohJpsiToMuPolynomial     ->Draw();
  new TCanvas;
  fIncohPsi2sToMuPolynomial    ->Draw();
  new TCanvas;
  fIncohPsi2sToMuPiPolynomial  ->Draw();
  new TCanvas;
  fTwoGammaToMuMediumPolynomial->Draw();
  new TCanvas;
  fTwoGammaToMuHighPolynomial  ->Draw();
  new TCanvas;

  TF1 *fFitPitDistrFunc = new TF1("fFitPitDistrFunc",fFitPtDistr,0.0001,3,7);
  fFitPitDistrFunc->SetNpx(1000000);
  fFitPitDistrFunc->SetParLimits(0, 0.0000001, 9999999999);
  fFitPitDistrFunc->SetParameter(0, 11000);
  fFitPitDistrFunc->SetParLimits(1, 0.0000001, 9999999999);
  fFitPitDistrFunc->SetParameter(2, 7000);
  fFitPitDistrFunc->SetParLimits(2, 0.0000001, 9999999999);
  fFitPitDistrFunc->SetParLimits(3, 0.0000001, 9999999999);
  // fFitPitDistrFunc->SetParLimits(4, 0.0000001, 9999999999);
  fFitPitDistrFunc->SetParameter(3, 0.05);
  // fFitPitDistrFunc->FixParameter(0, 0);
  // fFitPitDistrFunc->FixParameter(1, 0);
  // fFitPitDistrFunc->FixParameter(2, 0);
  // fFitPitDistrFunc->FixParameter(3, 0);
  fFitPitDistrFunc->FixParameter(4, fModelForHighPtTail->GetParameter(0) );
  fFitPitDistrFunc->FixParameter(5, fModelForHighPtTail->GetParameter(1) );
  fFitPitDistrFunc->FixParameter(6, fModelForHighPtTail->GetParameter(2) );
  fFitPitDistrFunc->Print();
  for(Int_t i = 0; i < 7; i++) cout << fFitPitDistrFunc->GetParameter(i) << endl << flush;


  fDimuonPtDistributionDataH->Fit( fFitPitDistrFunc,"LR","", 0.00001, 3 );
  TCanvas* PtDistrCanvas = new TCanvas( "PtDimuonFit", "PtDimuonFit", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gPad->SetTitle(  Form(  ";p_{T} (GeV/c);Counts / (%.0f MeV/c)",
                           fDimuonPtDistributionDataH->GetXaxis()->GetBinWidth(1)*1000.  )  );
  /* - Beautifying is starting now.
     -
   */
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleOffset(1.25);
  // fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.25);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.45);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleSize(0.045);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleSize(0.045);
  fDimuonPtDistributionDataH->GetXaxis()->SetLabelSize(0.045);
  fDimuonPtDistributionDataH->GetYaxis()->SetLabelSize(0.045);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleFont(42);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleFont(42);
  fDimuonPtDistributionDataH->GetXaxis()->SetLabelFont(42);
  fDimuonPtDistributionDataH->GetYaxis()->SetLabelFont(42);
  fDimuonPtDistributionDataH->GetXaxis()->SetNdivisions(408);
  fDimuonPtDistributionDataH->GetYaxis()->SetRangeUser(0.0000000000001, fDimuonPtDistributionDataH->GetMaximum()*1.3);
  fDimuonPtDistributionDataH->GetXaxis()->SetRangeUser(0, 4);
  // gPad ->SetLogy();
  fDimuonPtDistributionDataH->Draw("PEsame");
  // JPsiPeakFit    ->SetLineColor(kRed);
  // PsiPrimePeakFit->SetLineColor(kMagenta);
  // GammaGammaFit  ->SetLineColor(kGreen);
  // GammaGammaFit  ->SetLineStyle(kDashed);
  // JPsiPeakFit    -> SetLineWidth(3);
  // PsiPrimePeakFit-> SetLineWidth(3);
  // GammaGammaFit  -> SetLineWidth(3);
  // JPsiPeakFit    ->SetNpx(fDimuonPtDistributionDataH->GetNbinsX());
  // PsiPrimePeakFit->SetNpx(fDimuonPtDistributionDataH->GetNbinsX());
  // GammaGammaFit  ->SetNpx(fDimuonPtDistributionDataH->GetNbinsX());
  // JPsiPeakFit->FixParameter( 0, fFitPitDistrFunc->GetParameter(0)*fFitPitDistrFunc->GetParameter(15) );
  // JPsiPeakFit->FixParameter( 1, fFitPitDistrFunc->GetParameter(1) );
  // JPsiPeakFit->FixParameter( 2, fFitPitDistrFunc->GetParameter(2) );
  // JPsiPeakFit->FixParameter( 3, fFitPitDistrFunc->GetParameter(3) );
  // JPsiPeakFit->FixParameter( 4, fFitPitDistrFunc->GetParameter(4) );
  // PsiPrimePeakFit->FixParameter( 0, fFitPitDistrFunc->GetParameter(5)*fFitPitDistrFunc->GetParameter(16) );
  // PsiPrimePeakFit->FixParameter( 1, fFitPitDistrFunc->GetParameter(1+5) );
  // PsiPrimePeakFit->FixParameter( 2, fFitPitDistrFunc->GetParameter(2+5) );
  // PsiPrimePeakFit->FixParameter( 3, fFitPitDistrFunc->GetParameter(3+5) );
  // PsiPrimePeakFit->FixParameter( 4, fFitPitDistrFunc->GetParameter(4+5) );
  // GammaGammaFit->FixParameter( 0, fFitPitDistrFunc->GetParameter(10)*fFitPitDistrFunc->GetParameter(17) );
  // GammaGammaFit->FixParameter( 1, fFitPitDistrFunc->GetParameter(1+10) );
  // GammaGammaFit->FixParameter( 2, fFitPitDistrFunc->GetParameter(2+10) );
  // GammaGammaFit->FixParameter( 3, fFitPitDistrFunc->GetParameter(3+10) );
  // GammaGammaFit->FixParameter( 4, fFitPitDistrFunc->GetParameter(4+10) );
  // JPsiPeakFit    ->Draw("SAME");
  // PsiPrimePeakFit->Draw("SAME");
  // GammaGammaFit  ->Draw("SAME");
  // // GammaGammaFit  ->Draw("SAME");
  // // JPsiPeakFit    ->SetNpx(fDimuonPtDistributionDataH->GetNbinsX()/5);
  // // PsiPrimePeakFit->SetNpx(fDimuonPtDistributionDataH->GetNbinsX()/5);
  // // GammaGammaFit  ->SetNpx(fDimuonPtDistributionDataH->GetNbinsX()/5);



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
  if      ( selectionFlag == 0 ) latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  else if ( selectionFlag == 1 ) latex->DrawLatex(0.55,0.78,"#it{p}_{T} < 0.25 GeV/#it{c}");
  else if ( selectionFlag == 2 ) latex->DrawLatex(0.55,0.78,"#it{p}_{T} > 0.25 GeV/#it{c}");
  else                           latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  latex->DrawLatex(0.55,0.72,Form("%.1f < y < %.1f",-4.0,-2.5));

  /* - This is the part where we obtain the actual number of J/Psi, PsiPrime
     - and the background. This is still Kay's original code. I will modify it.
     - Hopefully if everything goes alright, I should have been able to complete
     - it by the time you are reading this.
     -
   */
  // TH1F* fCohJpsiToMuFromModelH = (TH1F*) JPsiPeakFit->GetHistogram()->Clone("fCohJpsiToMuFromModelH");
  // for (Int_t ibin=1; ibin<=fCohJpsiToMuFromModelH->GetNbinsX(); ibin++) {
  //   fCohJpsiToMuFromModelH->SetBinError(ibin,0);
  // }
  // fCohJpsiToMuFromModelH->Rebin(5);
  // fCohJpsiToMuFromModelH->Scale(0.2);
  // TH1F* fCohPsi2sToMuFromModelH = (TH1F*) PsiPrimePeakFit->GetHistogram()->Clone("fCohPsi2sToMuFromModelH");
  // for (Int_t ibin=1; ibin<=fCohPsi2sToMuFromModelH->GetNbinsX(); ibin++) {
  //   fCohPsi2sToMuFromModelH->SetBinError(ibin,0);
  // }
  // fCohPsi2sToMuFromModelH->Scale(0.2);
  // fCohPsi2sToMuFromModelH->Rebin(5);
  // TH1F* fTwoGammaFromModelH = (TH1F*) GammaGammaFit->GetHistogram()->Clone("fTwoGammaFromModelH");
  // for (Int_t ibin=1; ibin<=fTwoGammaFromModelH->GetNbinsX(); ibin++) {
  //   fTwoGammaFromModelH->SetBinError(ibin,0);
  // }
  // fTwoGammaFromModelH->Scale(0.2);
  // fTwoGammaFromModelH->Rebin(5);
  //
  //
  // Double_t numberOfTotalJPsi     = 0;
  // Double_t numberOfTotalPsi2s    = 0;
  // Double_t numberOfTotalBkg      = 0;
  // Double_t numberOfTotalJPsiErr  = 0;
  // Double_t numberOfTotalPsi2sErr = 0;
  // Double_t numberOfTotalBkgErr   = 0;
  // // if ( ptrSelectionFlag == 2 ) {
  // //   numberOfTotalJPsi  = fIncohJpsiToMuC  -> Integral();
  // //   numberOfTotalPsi2s = fIncohPsi2sToMuC -> Integral();
  // // } else {
  //   // numberOfTotalJPsi  = fCohJpsiToMuFromModelH -> Integral();
  //   // numberOfTotalPsi2s = fCohPsi2sToMuFromModelH-> Integral();
  //   // numberOfTotalJPsi  = JPsiPeakFit    -> Integral(2.2,6,1.E-15);
  //   // numberOfTotalPsi2s = PsiPrimePeakFit-> Integral(2.2,6,1.E-15);
  //   numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.05;
  //   numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.05;
  //   numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitPitDistrFunc->GetParError(15)/fFitPitDistrFunc->GetParameter(15);
  //   numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitPitDistrFunc->GetParError(16)/fFitPitDistrFunc->GetParameter(16);
  //
  // // }
  // // numberOfTotalBkg = fTwoGammaFromModelH      -> Integral();
  // // numberOfTotalBkg = GammaGammaFit-> Integral(2.2,6,1.E-15);
  // numberOfTotalBkg    = (GammaGammaFit-> Integral(2.2,6))/0.05;
  // numberOfTotalBkgErr = numberOfTotalBkg*fFitPitDistrFunc->GetParError(17)/fFitPitDistrFunc->GetParameter(17);
  // latex->DrawLatex(0.55,0.66,Form("N_{J/#psi} = %.0f #pm %.0f",        numberOfTotalJPsi,  numberOfTotalJPsiErr ));//fFitPitDistrFunc->GetParameter(0) *fFitPitDistrFunc->GetParError(15)/0.05 ) );
  // latex->DrawLatex(0.55,0.60,Form("N_{#psi(2S)} = %.0f #pm %.0f",      numberOfTotalPsi2s, numberOfTotalPsi2sErr));//fFitPitDistrFunc->GetParameter(5) *fFitPitDistrFunc->GetParError(16)/0.05 ) );
  // latex->DrawLatex(0.55,0.54,Form("N_{#gamma#gamma} = %.0f #pm %.0f",  numberOfTotalBkg,   numberOfTotalBkgErr  ));//fFitPitDistrFunc->GetParameter(10)*fFitPitDistrFunc->GetParError(17)/0.05 ) );
  // // latex->DrawLatex(0.55,0.54,Form("#sigma_{J/#psi} = %.0f #pm %.0f MeV/c^{2}", sigma.getVal()*1000,    sigma.getError()*1000));
  // // latex->DrawLatex(0.55,0.48,Form("#sigma_{#psi(2S)} = %.0f MeV/c^{2} fixed",  sigma2.getVal()*1000));
  //
  // /* - This part concerns the background of the two signals.
  //    - Here, we extrapolate the background and compute the significance maybe?
  //    -
  //  */
  // Double_t JPsiPeakBkg        = 0;
  // Double_t Psi2JPsiPeakBkg    = 0;
  // Double_t JPsiPeakSignal     = 0;
  // Double_t Psi2JPsiPeakSignal = 0;
  // // if ( ptrSelectionFlag == 2 ) {
  // //   JPsiPeakSignal     = fIncohJpsiToMuC  -> Integral( fIncohJpsiToMuC ->GetXaxis()->FindBin(2.75), fIncohJpsiToMuC ->GetXaxis()->FindBin(3.45) );
  // //   Psi2JPsiPeakSignal = fIncohPsi2sToMuC -> Integral( fIncohPsi2sToMuC->GetXaxis()->FindBin(3.45), fIncohPsi2sToMuC->GetXaxis()->FindBin(3.90) );
  // // } else {
  //   // JPsiPeakSignal     = fCohJpsiToMuFromModelH -> Integral(fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(2.75), fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(3.45));
  //   // Psi2JPsiPeakSignal = fCohPsi2sToMuFromModelH-> Integral(fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.45), fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.90));
  //   // JPsiPeakSignal     = fCohJpsiToMuFromModelH -> Integral(fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(2.75), fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(3.45));
  //   // Psi2JPsiPeakSignal = fCohPsi2sToMuFromModelH-> Integral(fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.45), fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.90));
  // // }
  // JPsiPeakBkg     = GammaGammaFit->Integral(2.75,3.45);
  // Psi2JPsiPeakBkg = GammaGammaFit->Integral(3.45,3.90);
  // // latex->DrawLatex(0.55,0.42,Form("N_{BG J/#psi} = %.0f #pm %.0f",   JPsiPeakBkg,     JPsiPeakBkg     * fFitPitDistrFunc->GetParError(17) / numberOfTotalJPsi ));
  // // latex->DrawLatex(0.55,0.36,Form("N_{BG #psi(2s)} = %.0f #pm %.0f", Psi2JPsiPeakBkg, Psi2JPsiPeakBkg * fFitPitDistrFunc->GetParError(17) / numberOfTotalPsi2s));
  // latex->DrawLatex(0.55,0.18,Form("      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
  //                                    fFitPitDistrFunc->GetChisquare(),
  //                                    fFitPitDistrFunc->GetNDF(),
  //                                    fFitPitDistrFunc->GetChisquare()/fFitPitDistrFunc->GetNDF()
  //                                    )
  //                                   );
  //
  //
  // // //Signal to background ratios
  // // Double_t errSB1 = sqrt(   fFitPitDistrFunc->GetParError(15) * fFitPitDistrFunc->GetParError(15)/( numberOfTotalJPsi  *numberOfTotalJPsi )  +
  // //                         ( fFitPitDistrFunc->GetParError(17) * fFitPitDistrFunc->GetParError(17)/( JPsiPeakBkg        * JPsiPeakBkg)) );
  // // Double_t errSB2 = sqrt(   fFitPitDistrFunc->GetParError(16) * fFitPitDistrFunc->GetParError(16)/( numberOfTotalPsi2s * numberOfTotalPsi2s )  +
  // //                           fFitPitDistrFunc->GetParError(17) * fFitPitDistrFunc->GetParError(17)/( Psi2JPsiPeakBkg    * Psi2JPsiPeakBkg ) );
  // // latex->DrawLatex(0.55,0.30,Form("S/B for J/#psi = %.2f #pm %.2f",   JPsiPeakSignal/(Double_t)JPsiPeakBkg,         errSB1*JPsiPeakSignal/(Double_t)JPsiPeakBkg));
  // // latex->DrawLatex(0.55,0.24,Form("S/B for #psi(2s) = %.2f #pm %.2f", Psi2JPsiPeakSignal/(Double_t)Psi2JPsiPeakBkg, errSB2*Psi2JPsiPeakSignal/(Double_t)Psi2JPsiPeakBkg));
  // // if      ( selectionFlag == 0 ) gPad->SaveAs("pngResults/PtIntJPsiTF1.png",      "RECREATE");
  // // else if ( selectionFlag == 1 ) gPad->SaveAs("pngResults/CoherentJPsiTF1.png",   "RECREATE");
  // // else if ( selectionFlag == 2 ) gPad->SaveAs("pngResults/IncoherentJPsiTF1.png", "RECREATE");
  // // else                           gPad->SaveAs("pngResults/ops.png",               "RECREATE");
  // //
  // // TLegend* l = new TLegend(0.2,0.55,0.5,0.85);
  // // l->SetMargin(0.1);
  // // l->SetBorderSize(0);
  // // l->AddEntry(  fDimuonPtDistributionDataH, "ALICE data 2018");
  // // l->AddEntry(  fFitPitDistrFunc, Form( "Fit: #chi^{2}/NDF = %.2f / %.2d = %.2f  ",
  // //                                   fFitPitDistrFunc->GetChisquare(),
  // //                                   fFitPitDistrFunc->GetNDF(),
  // //                                   fFitPitDistrFunc->GetChisquare()/fFitPitDistrFunc->GetNDF()
  // //                                   )
  // //                                  );
  // // if ( selectionFlag == 2 ) {
  // //   l->AddEntry(  JPsiPeakFit,        "Incoherent J/#psi");
  // //   l->AddEntry(  PsiPrimePeakFit,    "Incoherent #psi(2S)");
  // // } else {
  // //   l->AddEntry(  JPsiPeakFit,        "Coherent   J/#psi");
  // //   l->AddEntry(  PsiPrimePeakFit,    "Coherent   #psi(2S)");
  // // }
  // // // l->AddEntry(  fTwoGammaToMuMediumC, "Continuum  #gamma#gamma #rightarrow #mu#mu");
  // // l->AddEntry(  GammaGammaFit,        "Continuum  #gamma#gamma #rightarrow #mu#mu");
  // // l->Draw();
  //
  // if      ( selectionFlag == 0 ) gPad->SaveAs("pngResults/PtIntJPsiLegTF1.png",      "RECREATE");
  // else if ( selectionFlag == 1 ) gPad->SaveAs("pngResults/CoherentJPsiLegTF1.png",   "RECREATE");
  // else if ( selectionFlag == 2 ) gPad->SaveAs("pngResults/IncoherentJPsiLegTF1.png", "RECREATE");
  // else                           gPad->SaveAs("pngResults/ops.png",                  "RECREATE");
  //
  //
  // // cout << "nBackGround:           " << nBackGround.getVal() << endl;
  // // cout << "fValueBkgJPsiPeak:     " << fValueBkgJPsiPeak << endl;
  // // cout << "fValueBkgPsiPrimePeak: " << fValueBkgPsiPrimePeak << endl;
  //
  // // exp = TotalBackGroundJPsiPeak;
  // // err = errorTotalBackGroundJPsiPeak;
  //
  // new TCanvas;
  // fCohJpsiToMuFromModelH->Draw();
  //
  //
  // TCanvas* ChiSquareCanvas = new TCanvas( "ChiSquareCanvas", "ChiSquareCanvas", 900, 800 );
  // TH1F* ChiSquareH  = new TH1F( "ChiSquareH" ,"ChiSquareH" ,76, 2.2, 6);
  // TH1F* ChiSquareH2 = new TH1F( "ChiSquareH2","ChiSquareH2",76, 2.2, 6);
  // // TH1F* fTwoGammaFromModelH = (TH1F*) GammaGammaFit->GetHistogram()->Clone("fTwoGammaFromModelH");
  // // for (Int_t ibin=1; ibin<=fTwoGammaFromModelH->GetNbinsX(); ibin++) {
  // //   fTwoGammaFromModelH->SetBinError(ibin,0);
  // // }
  // Int_t whichBin  = fDimuonPtDistributionDataH->GetXaxis()->FindBin(2.225);
  // Int_t whichBin2 = ChiSquareH                 ->GetXaxis()->FindBin(2.225);
  // for( Int_t iLoop = 0; iLoop < 78; iLoop++ ) {
  //   // Double_t binContent = fDimuonPtDistributionDataH->GetBinContent(fDimuonPtDistributionDataH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   Double_t binContent = fDimuonPtDistributionDataH->GetBinContent(whichBin + iLoop);
  //   Double_t binError   = fDimuonPtDistributionDataH->GetBinError(whichBin + iLoop);
  //   Double_t binCenter  = fDimuonPtDistributionDataH->GetXaxis()->GetBinCenter(whichBin + iLoop);
  //   // Double_t fitValue   = fFitPitDistrFunc->Eval(2.225+iLoop*0.5);
  //   Double_t fitValue   = fFitPitDistrFunc->Eval(binCenter);
  //   cout << "binContent = " << binContent << endl;
  //   cout << "fitValue   = " << fitValue   << endl;
  //   // Double_t pull = fDimuonPtDistributionDataH->GetBinContent(fDimuonPtDistributionDataH->GetXaxis()->FindBin(2.225+iLoop*0.5)) - fFitPitDistrFunc->Eval(2.225+iLoop*0.5);
  //   Double_t pull = binContent - fitValue;
  //   cout << pull << endl;
  //   // Double_t help = fDimuonPtDistributionDataH->GetBinContent(fDimuonPtDistributionDataH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   // Double_t help = fDimuonPtDistributionDataH->GetBinContent(fDimuonPtDistributionDataH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   Double_t pullHelp;
  //   if( pull != 0 && binContent != 0 ) pullHelp = pull /binContent;
  //   else pullHelp = 0;
  //   cout << "pullHelp = " << pullHelp << endl;
  //   Double_t pullHelp2;
  //   if( pull != 0 && binError != 0 ) pullHelp2 = pull /binError;
  //   else pullHelp2 = 0;
  //   cout << "pullHelp2 = " << pullHelp2 << endl;
  //   cout << "binCenter = " << ChiSquareH                ->GetBinCenter( whichBin2 + iLoop ) << endl;
  //   cout << "binCenter = " << fDimuonPtDistributionDataH->GetBinCenter( whichBin  + iLoop ) << endl;
  //
  //   ChiSquareH->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fDimuonPtDistributionDataH->GetBinContent(fDimuonPtDistributionDataH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitPitDistrFunc->Eval(2.2+iLoop*0.5)
  //                              pullHelp
  //                              // pull
  //                              );
  //   ChiSquareH2->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fDimuonPtDistributionDataH->GetBinContent(fDimuonPtDistributionDataH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitPitDistrFunc->Eval(2.2+iLoop*0.5)
  //                              pullHelp2
  //                              // pull
  //                              );
  //
  // }
  // ChiSquareH->Draw();
  // TCanvas* ChiSquareCanvas2 = new TCanvas( "ChiSquareCanvas2", "ChiSquareCanvas2", 900, 800 );
  // ChiSquareH2->Draw();

}
