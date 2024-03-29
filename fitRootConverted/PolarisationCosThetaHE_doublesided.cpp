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
// Double_t* CBparameters[8];
Double_t  CBparameters[8][7];
Double_t  helpppppp = 0;

//____________________________
/* - Too many fits...
 * - The pointers have to be
 * - declared global!
 * -
 */
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
TH1F *fInvariantMassDistributionH = 0x0;


/* - Functions used for the fit!!
 * -
 */
// TF1* JPsiPeakFit     = new TF1( "JPsiPeakFit",    "crystalball",2.2,6);
// TF1* JPsiPeakFit     = new TF1( "JPsiPeakFit",    DoubleSidedCrystalballFunction,2.2,6, 7);
TF1* PsiPrimePeakFit = new TF1( "PsiPrimePeakFit","crystalball",2.2,6);
TF1* GammaGammaFit   = new TF1( "GammaGammaFit",
                                "[0]*TMath::Exp(-[1]*x)*( (x > 4) ? 1 : 1 + [2]*(x-4)*(x-4) + [3]*(x-4)*(x-4)*(x-4) + [4]*(x-4)*(x-4)*(x-4)*(x-4) )",
                                2.2,6
                                );


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
/* - Evgeny's own fit function.
 * -
 */
// double fsum2(double *x, double *par){
//   Double_t parPsiPrime[5];
//   Double_t parBkg[5];
//   for( Int_t i=0 ; i<5 ; i++ ) parPsiPrime[i]=par[i+7];
//   for( Int_t i=0 ; i<5 ; i++ ) parBkg[i]     =par[i+12];
//   return par[17]*JPsiPeakFit->EvalPar(x,par)+par[18]*PsiPrimePeakFit->EvalPar(x,parPsiPrime)+par[19]*GammaGammaFit->EvalPar(x,parBkg);
// }
//
// double fsum(double *x, double *par){
//   Double_t parPsiPrime[5];
//   Double_t parBkg[5];
//   for( Int_t i=0 ; i<5 ; i++ ) parPsiPrime[i]=par[i+5];
//   for( Int_t i=0 ; i<5 ; i++ ) parBkg[i]     =par[i+10];
//   return par[15]*JPsiPeakFit->EvalPar(x,par)+par[16]*PsiPrimePeakFit->EvalPar(x,parPsiPrime)+par[17]*GammaGammaFit->EvalPar(x,parBkg);
// }

Double_t DoubleSidedCrystalballFunction(double *x, double *par)
{
  Double_t alpha_l = par[0];
  Double_t alpha_h = par[1];
  Double_t n_l     = par[2];
  Double_t n_h     = par[3];
  Double_t mean	   = par[4];
  Double_t sigma   = par[5];
  Double_t N	     = par[6];
  Double_t t       = (x[0]-mean)/sigma;
  Double_t result;
  Double_t fact1TLessMinosAlphaL = alpha_l/n_l;
  Double_t fact2TLessMinosAlphaL = (n_l/alpha_l) - alpha_l -t;
  Double_t fact1THihgerAlphaH    = alpha_h/n_h;
  Double_t fact2THigherAlphaH    = (n_h/alpha_h) - alpha_h +t;

  if (-alpha_l <= t && alpha_h >= t){
      result = exp(-0.5*t*t);
  } else if (t < -alpha_l) {
      result = exp(-0.5*alpha_l*alpha_l)*pow(fact1TLessMinosAlphaL*fact2TLessMinosAlphaL, -n_l);
  } else if (t > alpha_h) {
      result = exp(-0.5*alpha_h*alpha_h)*pow(fact1THihgerAlphaH*fact2THigherAlphaH, -n_h);
  }

  return N*result;
}




TF1* JPsiPeakFit     = new TF1( "JPsiPeakFit",    DoubleSidedCrystalballFunction,2.2,6, 7);

double fsum2(double *x, double *par){
  Double_t parPsiPrime[5];
  Double_t parBkg[5];
  for( Int_t i=0 ; i<5 ; i++ ) parPsiPrime[i]=par[i+7];
  for( Int_t i=0 ; i<5 ; i++ ) parBkg[i]     =par[i+12];
  return par[17]*JPsiPeakFit->EvalPar(x,par)+par[18]*PsiPrimePeakFit->EvalPar(x,parPsiPrime)+par[19]*GammaGammaFit->EvalPar(x,parBkg);
}

double fsum(double *x, double *par){
  Double_t parPsiPrime[5];
  Double_t parBkg[5];
  for( Int_t i=0 ; i<5 ; i++ ) parPsiPrime[i]=par[i+5];
  for( Int_t i=0 ; i<5 ; i++ ) parBkg[i]     =par[i+10];
  return par[15]*JPsiPeakFit->EvalPar(x,par)+par[16]*PsiPrimePeakFit->EvalPar(x,parPsiPrime)+par[17]*GammaGammaFit->EvalPar(x,parBkg);
}

//_____________________________________________________________________________
/* - Fit function for the templates of the JPsi.
 * -
 */
// void fCrystalBallPDFFitJPsi(TH1F* histoToBeFit, Double_t &bookKeeping[5])
// {
//   TF1* CBfit     = new TF1("CBfit","[4]*ROOT::Math::crystalball_pdf(x, [0], [1], [2], [3])",2,15);
//   CBfit       ->SetParameter(4,1);
//   CBfit       ->SetParameter(0,1.08);
//   CBfit       ->SetParameter(1,3689197);
//   CBfit       ->SetParameter(2,0.090);
//   CBfit       ->SetParameter(3,3.1);
//   CBfit       ->SetNpx(1000);
//   CBfit       ->Draw();
//   histoToBeFit->Fit(CBfit, "R");
//   CBfit       ->SetParameter(0,1/CBfit->Integral(2,15));
//   // bookKeeping = new Double_t[5];
//   for(Int_t i = 0; i < 5; i++){
//     bookKeeping[i] = CBfit->GetParameter(i);
//   }
// }
//_____________________________________________________________________________
/* - Fit function for the templates of the PsiPrime.
 * -
 */
void fCrystalBallFitJPsi(TH1F* histoToBeFit)//, Double_t &bookKeeping[5])
{
  TF1* CBfit     = new TF1("CBfit",DoubleSidedCrystalballFunction,2,15,7);
  // TF1* CBfit     = new TF1("CBfit","crystalball",2,15);
  // TF1* CBfit     = new TF1("CBfit","crystalball",1.8,7);
  // CBfit       ->SetParameter(0,1);
  // // CBfit       ->SetParameter(3,1.08);
  // CBfit       ->FixParameter(3,1.08);
  // CBfit       ->SetParameter(4,115);
  // // CBfit       ->SetParameter(4,6);
  // CBfit       ->SetParLimits(4,110,120);
  // // CBfit       ->SetParLimits(4,1,12);
  // // CBfit       ->SetParLimits(4,12,99999999);
  // CBfit       ->SetParameter(2,0.090);
  // CBfit       ->SetParameter(1,3.115);
  // CBfit       ->SetParLimits(1,3.113,3.17);
  // CBfit       ->SetNpx(1000);
  CBfit       ->FixParameter(6,5.3e-02); // Normalisation
  // CBfit       ->SetParameter(3,1.08);
  CBfit       ->SetParameter(1,1.08); //alpha_high
  CBfit       ->SetParameter(0,1.08); //alpha_low
  CBfit       ->SetParameter(2,10);
  CBfit       ->SetParLimits(2,2,10000);
  CBfit       ->SetParameter(3,120);
  CBfit       ->SetParLimits(3,100,1000);
  CBfit       ->SetParameter(4,3.115);
  // CBfit       ->SetParLimits(4,3.113,3.17);
  CBfit       ->SetParLimits(4,3.,3.17);
  CBfit       ->SetParameter(5,0.090);
  CBfit       ->SetNpx(1000);

  TCanvas*    JPsiCanvas = new TCanvas( "JPsiCanvas", "JPsiCanvas", 900, 800 );
  CBfit       ->Draw();
  histoToBeFit->Fit(CBfit, "R");
  // CBfit       ->SetParameter(0,1/CBfit->Integral(2,15));
  // bookKeeping = new Double_t[5];
  for(Int_t i = 0; i < 7; i++){
    CBparameters[0][i] = CBfit->GetParameter(i);
  }
}

void trialfunc(Int_t select = 15){
  fileMC[0] = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  dirMC[0]  = fileMC[0]->GetDirectory("MyTask");
  dirMC[0]->GetObject("MyOutputContainer", listingsMC[0]);
  fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject(Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH_%d", select));
  fCohJpsiToMu        -> Sumw2();
  Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral();
  fCohJpsiToMu        -> Scale( 1/Integral_fCohJpsiToMu        );
  fCrystalBallFitJPsi(fCohJpsiToMu);

}


//_____________________________________________________________________________
/* - Fit function for the templates of the PsiPrime.
 * -
 */
void fCrystalBallFitPsiPrime(TH1F* histoToBeFit)//, Double_t &bookKeeping[5])
{
  TF1* CBfit     = new TF1("CBfit","crystalball",2,15);
  CBfit       ->SetParameter(0,1);
  CBfit       ->SetParameter(3,1.08);
  // CBfit       ->SetParameter(4,3689197);
  // CBfit       ->SetParLimits(4,1,99999999);
  CBfit       ->SetParameter(4,20);
  CBfit       ->SetParLimits(4,1.5,100);
  CBfit       ->SetParameter(2,0.070);
  CBfit       ->SetParameter(1,3.67);
  CBfit       ->SetNpx(1000);
  TCanvas*    PsiPrimeCanvas = new TCanvas( "PsiPrimeCanvas", "PsiPrimeCanvas", 900, 800 );
  CBfit       ->Draw();
  cout << "histo integral = " << histoToBeFit->Integral() << endl << flush;
  histoToBeFit->Fit(CBfit, "R");
  cout << "CBfit integral = " << CBfit->Integral(2.1,9) << endl << flush;
  // CBfit       ->SetParameter(0,1/(CBfit->Integral(2.1,9)));
  cout << "CBfit integral / histo width = " << CBfit->Integral(2.1,9)/histoToBeFit->GetXaxis()->GetBinWidth(1) << endl << flush;
  // CBfit       ->SetParameter(0,1/CBfit->Integral(3,14));
  // bookKeeping = new Double_t[5];
  for(Int_t i = 0; i < 5; i++){
    CBparameters[1][i] = CBfit->GetParameter(i);
  }
}
//_____________________________________________________________________________
/* - Fit function for the templates of the PsiPrime.
 * -
 */
void fBkgPolFit(TH1F* histoToBeFit)//, Double_t &bookKeeping[5])
{
  TF1* PolBkg     = new TF1("PolBkg","[0]*TMath::Exp(-[1]*x)*( (x > 4) ? 1 : 1 + [2]*(x-4)*(x-4) + [3]*(x-4)*(x-4)*(x-4) + [4]*(x-4)*(x-4)*(x-4)*(x-4) )",2,8);
  // TF1* PolBkg     = new TF1("PolBkg","[0]*TMath::Exp(-[1]*x)*( (x > 4) ? 1 : 1 + [2]*(x-4)*(x-4) + [3]*(x-4)*(x-4)*(x-4) + [4]*(x-4)*(x-4)*(x-4)*(x-4) )",1.8,12);
  // PolBkg       ->SetParameter(0,0.025);   // best fit
  // PolBkg       ->SetParameter(3,0);       // best fit
  // PolBkg       ->SetParameter(4,0);       // best fit
  // PolBkg       ->SetParameter(2,-0.25);   // best fit
  PolBkg       ->SetParameter(0,0.025);
  PolBkg       ->SetParameter(3,0.0001);
  PolBkg       ->SetParameter(4,0.0001);
  PolBkg       ->SetParameter(2,0.0001);
  PolBkg       ->SetParLimits(0,0.0001,1);
  PolBkg       ->SetParLimits(3,0.00000001,1);
  PolBkg       ->SetParLimits(4,0.00000001,1);
  PolBkg       ->SetParLimits(2,0.00000001,1);

  // PolBkg       ->SetParameter(0,0.07);
  // PolBkg       ->SetParameter(3,0);
  // PolBkg       ->SetParameter(4,0.25);
  // PolBkg       ->SetParameter(2,0.65);
  PolBkg       ->SetParameter(1,0.9);
  PolBkg       ->SetParLimits(1,0.8,1);
  PolBkg       ->SetNpx(1000);
  TCanvas*      BkgCanvas = new TCanvas( "BkgCanvas", "BkgCanvas", 900, 800 );
  PolBkg       ->Draw();
  histoToBeFit ->Fit(PolBkg, "R");
  // PolBkg       ->SetParameter(0,1/PolBkg->Integral(2,8));
  // bookKeeping = new Double_t[5];
  for(Int_t i = 0; i < 5; i++){
    CBparameters[4][i] = PolBkg->GetParameter(i);
  }
}
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitJPsiTemplateMC(const int selectionFlag = 0){

  // if ( selectionFlag != 0 ) {
  //   fCohJpsiToMu = (TH1F*)listingsMC[0]->FindObject( Form( "fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH_%d", selectionFlag ) );
  // } else {
  //   fCohJpsiToMu = (TH1F*)listingsMC[0]->FindObject("fInvariantMassDistributionH");
  // }
  fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject(Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH_%d", selectionFlag));
  // fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject(Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH_%d", selectionFlag));
  fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fInvariantMassDistributionH");
  fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fInvariantMassDistributionH");
  fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fInvariantMassDistributionH");
  fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fInvariantMassDistributionH");
  fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fInvariantMassDistributionH");
  // fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject(Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH_%d", selectionFlag));
  fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fInvariantMassDistributionH");
  fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fInvariantMassDistributionH");
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

  fCrystalBallFitJPsi    (fCohJpsiToMu);//,      CBparameters[0]);
  // if( selectionFlag < 10 ){
  fCrystalBallFitPsiPrime(fCohPsi2sToMu);//,     CBparameters[1]);
  // fCrystalBallFitJPsi    (fIncohJpsiToMu,    CBparameters[2]);
  // fCrystalBallFitPsiPrime(fIncohPsi2sToMu,   CBparameters[3]);
  fBkgPolFit             (fTwoGammaToMuHigh);//, CBparameters[4]);
  // }
}
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitJPsiTemplate(const int selectionFlag, Int_t SignalRangeModeFlag = 0){

  // TH1F *fInvariantMassDistributionH = 0x0;
  fInvariantMassDistributionH = 0x0;
  fInvariantMassDistributionH = (TH1F*)listings->FindObject( Form("fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH_%d", selectionFlag) );
  fInvariantMassDistributionH->Rebin(5);
  // fInvariantMassDistributionH->Rebin(4);
  // if( selectionFlag < 16 || selectionFlag > 23 ) {
  //   fInvariantMassDistributionH->Rebin(2);
  // }
  // fInvariantMassDistributionH->Rebin(8);
  fInvariantMassDistributionH->Draw("PE");

  fInvariantMassDistributionH->SetLineColor(kBlue);
  fInvariantMassDistributionH->SetLineStyle(kSolid);
  fInvariantMassDistributionH->SetLineWidth(3);
  fInvariantMassDistributionH->SetMarkerStyle(kFullCircle);
  fInvariantMassDistributionH->SetMarkerColor(kBlue);
  fInvariantMassDistributionH->SetMarkerSize(1);
  fInvariantMassDistributionH->GetXaxis()->SetTitle("M_{#mu#mu} [GeV/#it{c}^{2}]");
  fInvariantMassDistributionH->GetYaxis()->SetTitle( Form( "Counts / (%.3f GeV/#it{c})",
                                                          fInvariantMassDistributionH->GetXaxis()->GetBinWidth(1)
                                                        )
                                                    );
  fInvariantMassDistributionH->SetTitle("");


  // fitJPsiTemplateMC();
  new TCanvas;



  TF1 *fFitInvMass = new TF1("fFitInvMass",fsum2,1.8,8,20);
  fFitInvMass->SetNpx(1000000);
  fFitInvMass->FixParameter(0, CBparameters[0][0]);    // best
  fFitInvMass->FixParameter(1, CBparameters[0][1]);    // best
  fFitInvMass->FixParameter(2, CBparameters[0][2]);    // best
  fFitInvMass->FixParameter(3, CBparameters[0][3]);    // best
  fFitInvMass->SetParameter(4, CBparameters[0][4]);
  fFitInvMass->SetParLimits(4, CBparameters[0][4]*0.9, CBparameters[0][4]*1.1);
  fFitInvMass->SetParameter(5, CBparameters[0][5]);
  fFitInvMass->SetParLimits(5, CBparameters[0][5]*0.8, CBparameters[0][5]*1.2);
  fFitInvMass->FixParameter(6, CBparameters[0][6]);    // best
  fFitInvMass->FixParameter(0+7, CBparameters[1][0]);
  fFitInvMass->FixParameter(3+7, CBparameters[1][3]);
  fFitInvMass->FixParameter(4+7, CBparameters[1][4]);
  fFitInvMass->FixParameter(1+7, CBparameters[1][1]);
  fFitInvMass->FixParameter(2+7, fFitInvMass->GetParameter(2)*CBparameters[1][2]/CBparameters[0][2]);
  fFitInvMass->FixParameter(0+12, CBparameters[4][0]);    // best
  fFitInvMass->FixParameter(2+12, CBparameters[4][2]);    // best
  fFitInvMass->FixParameter(3+12, CBparameters[4][3]);    // best
  fFitInvMass->FixParameter(4+12, CBparameters[4][4]);    // best
  // fFitInvMass->FixParameter(2+10, 0.509393 );    // Fix to previous fit...
  // fFitInvMass->FixParameter(3+10, 0.833773 );    // Fix to previous fit...
  // fFitInvMass->FixParameter(4+10, 0.237383 );    // Fix to previous fit...
  // fFitInvMass->SetParameter(0+10, CBparameters[4][0]);    // mmmh
  // fFitInvMass->SetParameter(2+10, CBparameters[4][2]);    // mmmh
  // fFitInvMass->SetParameter(3+10, CBparameters[4][3]);    // mmmh
  // fFitInvMass->SetParameter(4+10, CBparameters[4][4]);    // mmmh
  fFitInvMass->SetParameter(1+12, CBparameters[4][1]);
  // fFitInvMass->FixParameter(1+10, CBparameters[4][1]);
  fFitInvMass->SetParLimits(1+12, CBparameters[4][1]*0.9, CBparameters[4][1]*1.1);
  // fFitInvMass->SetParameter(15, 1);
  // fFitInvMass->SetParameter(16, 1);
  // fFitInvMass->SetParameter(17, 1);
  // fFitInvMass->SetParLimits(15, 0.001, 9999999999);
  // fFitInvMass->SetParLimits(16, 0.001, 9999999999);
  // fFitInvMass->SetParLimits(17, 0.001, 9999999999);
  // if(selectionFlag == 19) {
  //   fFitInvMass->SetParameter(15, helpppppp * 2);
  // } else {
  //   fFitInvMass->SetParameter(15, 1);
  // }
  fFitInvMass->SetParameter(17, 1);
  fFitInvMass->SetParameter(18, 1);
  fFitInvMass->SetParameter(19, 1);
  // if(selectionFlag == 19) {
  //   fFitInvMass->SetParLimits(15, helpppppp * 1.2, 9999999999);
  // } else {
  //   fFitInvMass->SetParLimits(15, 0.001, 9999999999);
  // }
  fFitInvMass->SetParLimits(17, 0.001, 9999999999);
  fFitInvMass->SetParLimits(18, 0.001, 9999999999);
  fFitInvMass->SetParLimits(19, 0.001, 9999999999);
  fFitInvMass->Print();
  for(Int_t i = 0; i < 20; i++) cout << fFitInvMass->GetParameter(i) << endl << flush;


  fInvariantMassDistributionH->Fit( fFitInvMass,"LR","", 2.2, 6. );
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
  JPsiPeakFit    ->SetLineColor(kRed);
  PsiPrimePeakFit->SetLineColor(kMagenta);
  GammaGammaFit  ->SetLineColor(kGreen);
  GammaGammaFit  ->SetLineStyle(kDashed);
  JPsiPeakFit    -> SetLineWidth(3);
  PsiPrimePeakFit-> SetLineWidth(3);
  GammaGammaFit  -> SetLineWidth(3);
  JPsiPeakFit    ->SetNpx(fInvariantMassDistributionH->GetNbinsX());
  PsiPrimePeakFit->SetNpx(fInvariantMassDistributionH->GetNbinsX());
  GammaGammaFit  ->SetNpx(fInvariantMassDistributionH->GetNbinsX());
  JPsiPeakFit->FixParameter( 0, fFitInvMass->GetParameter(0) );
  JPsiPeakFit->FixParameter( 1, fFitInvMass->GetParameter(1) );
  JPsiPeakFit->FixParameter( 2, fFitInvMass->GetParameter(2) );
  JPsiPeakFit->FixParameter( 3, fFitInvMass->GetParameter(3) );
  JPsiPeakFit->FixParameter( 4, fFitInvMass->GetParameter(4) );
  JPsiPeakFit->FixParameter( 5, fFitInvMass->GetParameter(5) );
  JPsiPeakFit->FixParameter( 6, fFitInvMass->GetParameter(6)*fFitInvMass->GetParameter(17) );
  PsiPrimePeakFit->FixParameter( 0, fFitInvMass->GetParameter(7)*fFitInvMass->GetParameter(18) );
  PsiPrimePeakFit->FixParameter( 1, fFitInvMass->GetParameter(1+7) );
  PsiPrimePeakFit->FixParameter( 2, fFitInvMass->GetParameter(2+7) );
  PsiPrimePeakFit->FixParameter( 3, fFitInvMass->GetParameter(3+7) );
  PsiPrimePeakFit->FixParameter( 4, fFitInvMass->GetParameter(4+7) );
  GammaGammaFit->FixParameter( 0, fFitInvMass->GetParameter(12)*fFitInvMass->GetParameter(19) );
  GammaGammaFit->FixParameter( 1, fFitInvMass->GetParameter(1+12) );
  GammaGammaFit->FixParameter( 2, fFitInvMass->GetParameter(2+12) );
  GammaGammaFit->FixParameter( 3, fFitInvMass->GetParameter(3+12) );
  GammaGammaFit->FixParameter( 4, fFitInvMass->GetParameter(4+12) );
  JPsiPeakFit    ->Draw("SAME");
  PsiPrimePeakFit->Draw("SAME");
  GammaGammaFit  ->Draw("SAME");
  if(selectionFlag == 18) helpppppp = fFitInvMass->GetParameter(17);
  // GammaGammaFit  ->Draw("SAME");
  // JPsiPeakFit    ->SetNpx(fInvariantMassDistributionH->GetNbinsX()/5);
  // PsiPrimePeakFit->SetNpx(fInvariantMassDistributionH->GetNbinsX()/5);
  // GammaGammaFit  ->SetNpx(fInvariantMassDistributionH->GetNbinsX()/5);



  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  latex->DrawLatex(0.55,0.78,Form("#color[2]{%.3f < cos(#theta) < %.3f}", -1 + (Double_t)selectionFlag * 0.08, -1 + ((Double_t)selectionFlag + 1.00) * 0.08));
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T} < 0.25 GeV/#it{c}");
  latex->DrawLatex(0.55,0.72,Form("%.1f < y < %.1f",-4.0,-2.5));

  /* - This is the part where we obtain the actual number of J/Psi, PsiPrime
     - and the background. This is still Kay's original code. I will modify it.
     - Hopefully if everything goes alright, I should have been able to complete
     - it by the time you are reading this.
     -
   */


  Double_t numberOfTotalJPsi     = 0;
  Double_t numberOfTotalPsi2s    = 0;
  Double_t numberOfTotalBkg      = 0;
  Double_t numberOfTotalJPsiErr  = 0;
  Double_t numberOfTotalPsi2sErr = 0;
  Double_t numberOfTotalBkgErr   = 0;
  // numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.05;  // USUAL FIT
  // numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.05;  // USUAL FIT
  // numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(15)/fFitInvMass->GetParameter(15);
  // numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(16)/fFitInvMass->GetParameter(16);

  if (        SignalRangeModeFlag == 0 ) {
    numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.05;  // USUAL FIT
    numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(18)/fFitInvMass->GetParameter(18);
  } else if ( SignalRangeModeFlag == 1 ) {
    numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.85,3.35))/0.05;  // USUAL FIT
    numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.85,3.35))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(18)/fFitInvMass->GetParameter(18);
  } else if ( SignalRangeModeFlag == 2 ) {
    numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.8,3.35))/0.05;  // USUAL FIT
    numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.8,3.35))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(18)/fFitInvMass->GetParameter(18);
  } else if ( SignalRangeModeFlag == 3 ) {
    numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.9,3.35))/0.05;  // USUAL FIT
    numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.9,3.35))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(18)/fFitInvMass->GetParameter(18);
  } else if ( SignalRangeModeFlag == 4 ) {
    numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.85,3.4))/0.05;  // USUAL FIT
    numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.85,3.4))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(18)/fFitInvMass->GetParameter(18);
  } else if ( SignalRangeModeFlag == 5 ) {
    numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.85,3.3))/0.05;  // USUAL FIT
    numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.85,3.3))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(18)/fFitInvMass->GetParameter(18);
  } else {
    numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.05;  // USUAL FIT
    numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(18)/fFitInvMass->GetParameter(18);
  }


  numberOfTotalBkg    = (GammaGammaFit-> Integral(2.2,6))/0.05;
  numberOfTotalBkgErr = numberOfTotalBkg*fFitInvMass->GetParError(19)/fFitInvMass->GetParameter(19);
  latex->DrawLatex(0.55,0.66,Form("N_{J/#psi} = %.0f #pm %.0f",        numberOfTotalJPsi,  numberOfTotalJPsiErr ));//fFitInvMass->GetParameter(0) *fFitInvMass->GetParError(15)/0.05 ) );
  latex->DrawLatex(0.55,0.60,Form("N_{#psi(2S)} = %.0f #pm %.0f",      numberOfTotalPsi2s, numberOfTotalPsi2sErr));//fFitInvMass->GetParameter(5) *fFitInvMass->GetParError(16)/0.05 ) );
  latex->DrawLatex(0.55,0.54,Form("N_{#gamma#gamma} = %.0f #pm %.0f",  numberOfTotalBkg,   numberOfTotalBkgErr  ));//fFitInvMass->GetParameter(10)*fFitInvMass->GetParError(17)/0.05 ) );

  /* - This part concerns the background of the two signals.
     - Here, we extrapolate the background and compute the significance maybe?
     -
   */
  Double_t JPsiPeakBkg        = 0;
  Double_t Psi2JPsiPeakBkg    = 0;
  Double_t JPsiPeakSignal     = 0;
  Double_t Psi2JPsiPeakSignal = 0;
  // JPsiPeakBkg     = GammaGammaFit->Integral(2.75,3.45);
  JPsiPeakBkg     = GammaGammaFit->Integral(2.4,2.8)/0.05;
  Psi2JPsiPeakBkg = GammaGammaFit->Integral(3.45,3.90)/0.05;

  JPsiPeakValue    = numberOfTotalJPsi;
  JPsiPeakValueErr = numberOfTotalJPsiErr;
  BkgValue         = JPsiPeakBkg;
  BkgValueError    = JPsiPeakBkg * fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);


  // latex->DrawLatex(0.55,0.42,Form("N_{BG J/#psi} = %.0f #pm %.0f",   JPsiPeakBkg,     JPsiPeakBkg     * fFitInvMass->GetParError(17) / numberOfTotalJPsi ));
  // latex->DrawLatex(0.55,0.36,Form("N_{BG #psi(2s)} = %.0f #pm %.0f", Psi2JPsiPeakBkg, Psi2JPsiPeakBkg * fFitInvMass->GetParError(17) / numberOfTotalPsi2s));
  latex->DrawLatex(0.55,0.18,Form("      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                     fFitInvMass->GetChisquare(),
                                     fFitInvMass->GetNDF(),
                                     fFitInvMass->GetChisquare()/fFitInvMass->GetNDF()
                                     )
                                    );



  if      ( SignalRangeModeFlag == 0 ) { gPad->SaveAs(Form("pngResults/CosThetaHe_DSCB_%d.png",        selectionFlag), "recreate"); }
  else if ( SignalRangeModeFlag == 1 ) { gPad->SaveAs(Form("pngResults/CosThetaHe_DSCB_%d_Range1.png", selectionFlag), "recreate"); }
  else if ( SignalRangeModeFlag == 2 ) { gPad->SaveAs(Form("pngResults/CosThetaHe_DSCB_%d_Range2.png", selectionFlag), "recreate"); }
  else if ( SignalRangeModeFlag == 3 ) { gPad->SaveAs(Form("pngResults/CosThetaHe_DSCB_%d_Range3.png", selectionFlag), "recreate"); }
  else if ( SignalRangeModeFlag == 4 ) { gPad->SaveAs(Form("pngResults/CosThetaHe_DSCB_%d_Range4.png", selectionFlag), "recreate"); }
  else if ( SignalRangeModeFlag == 5 ) { gPad->SaveAs(Form("pngResults/CosThetaHe_DSCB_%d_Range5.png", selectionFlag), "recreate"); }
  else                                 { gPad->SaveAs(Form("pngResults/CosThetaHe_DSCB_%d.png",        selectionFlag), "recreate"); }



}
//_____________________________________________________________________________
/* - Here I create the new TH1 for the after the signal extraction.
 * - Basically I run the fit function many times and then I memorise
 * - the values each time. After that I fill with a setbincontent
 * - and a setbinerror the
 * -
 * - SignalRangeMode, selects the signal extraction range:
 * - 0: standard
 * - 1: ehm
 */
void CreateCosThetaTh1(const char* AnalysisName, Int_t SignalRangeMode = 0){
  fileMC[0] = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2019-09-17/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2019-09-17/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2019-09-17/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2019-09-17/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2019-09-17/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2019-09-17/kTwoGammaToMuHigh/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2019-09-17/kTwoGammaToMuMedium/AnalysisResults.root");
  // TDirectory* dirMC[8];
  cout << "CHECKPOINT 4 " << endl << flush;

  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory] = fileMC[iDirectory]->GetDirectory("MyTask");
  }
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  // TList* listingsMC[8];
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory]->GetObject("MyOutputContainer", listingsMC[iDirectory]);
  }


  // fitJPsiTemplateMC();
  fileList = new TFile(AnalysisName);
  dir      = fileList->GetDirectory("MyTask");
  dir->GetObject("MyOutputContainer", listings);


  TH1F* CosThetaAfterSignalExtractionH =
            new TH1F( "CosThetaAfterSignalExtractionH",
                      "CosThetaAfterSignalExtractionH",
                      25, -1, 1
                      );
  TH1F* CosThetaGammaGammaH =
            new TH1F( "CosThetaGammaGammaH",
                      "CosThetaGammaGammaH",
                      25, -1, 1
                      );

  TH1F* CosThetaAfterSignalExtractionErrorsH =
            new TH1F( "CosThetaAfterSignalExtractionErrorsH",
                      "CosThetaAfterSignalExtractionErrorsH",
                      25, -1, 1
                      );
  TH1F* CosThetaGammaGammaErrorsH =
            new TH1F( "CosThetaGammaGammaErrorsH",
                      "CosThetaGammaGammaErrorsH",
                      25, -1, 1
                      );

  for (size_t iCosThetaBins = 5; iCosThetaBins < 20; iCosThetaBins++) {
    // for (size_t iPhiBins = 0; iPhiBins < 10; iPhiBins++) {
      JPsiPeakValue    = 0;
      JPsiPeakValueErr = 0;
      BkgValue         = 0;
      BkgValueError    = 0;
      // if( iCosThetaBins > 17 && iCosThetaBins < 23 ) {
      cout << "CHECKPOINT 1 " << endl << flush;
        fitJPsiTemplateMC(iCosThetaBins);
        cout << "CHECKPOINT 2 " << endl << flush;

      // } else {
      //   fitJPsiTemplateMC();
      // }
      fitJPsiTemplate(iCosThetaBins, SignalRangeMode);
      cout << "CHECKPOINT 3 " << endl << flush;


      CosThetaAfterSignalExtractionH->Fill(  -0.96 + (Double_t)iCosThetaBins * 0.08,
                                            // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                            JPsiPeakValue
                                            );

      CosThetaAfterSignalExtractionErrorsH->Fill(  -0.96 + (Double_t)iCosThetaBins * 0.08,
                                                   JPsiPeakValue
                                                   );
      CosThetaGammaGammaErrorsH->Fill(  -0.96 + (Double_t)iCosThetaBins * 0.08,
                                        BkgValue
                                        );
      CosThetaAfterSignalExtractionErrorsH->SetBinError(  iCosThetaBins + 1 ,
                                                          // iPhiBins      + 1 ,
                                                          JPsiPeakValueErr
                                                          );
      CosThetaGammaGammaErrorsH->SetBinError(  iCosThetaBins + 1 ,
                                               // iPhiBins      + 1 ,
                                               BkgValueError
                                               );

    // }
  }

  // TFile f("pngResults/CosThetaHeFrame.root", "recreate");
  if      ( SignalRangeMode == 0 ) {
    TFile f("pngResults/CosThetaHeFrame_DSCB.root",   "recreate");
    CosThetaAfterSignalExtractionH      ->Write();
    CosThetaGammaGammaH                 ->Write();
    CosThetaAfterSignalExtractionErrorsH->Write();
    CosThetaGammaGammaErrorsH           ->Write();
    f.Close();
  }
  else if ( SignalRangeMode == 1 ) {
    TFile f("pngResults/CosThetaHeFrame_DSCB_1.root", "recreate");
    CosThetaAfterSignalExtractionH      ->Write();
    CosThetaGammaGammaH                 ->Write();
    CosThetaAfterSignalExtractionErrorsH->Write();
    CosThetaGammaGammaErrorsH           ->Write();
    f.Close();
  }
  else if ( SignalRangeMode == 2 ) {
    TFile f("pngResults/CosThetaHeFrame_DSCB_2.root", "recreate");
    CosThetaAfterSignalExtractionH      ->Write();
    CosThetaGammaGammaH                 ->Write();
    CosThetaAfterSignalExtractionErrorsH->Write();
    CosThetaGammaGammaErrorsH           ->Write();
    f.Close();
  }
  else if ( SignalRangeMode == 3 ) {
    TFile f("pngResults/CosThetaHeFrame_DSCB_3.root", "recreate");
    CosThetaAfterSignalExtractionH      ->Write();
    CosThetaGammaGammaH                 ->Write();
    CosThetaAfterSignalExtractionErrorsH->Write();
    CosThetaGammaGammaErrorsH           ->Write();
    f.Close();
  }
  else if ( SignalRangeMode == 4 ) {
    TFile f("pngResults/CosThetaHeFrame_DSCB_4.root", "recreate");
    CosThetaAfterSignalExtractionH      ->Write();
    CosThetaGammaGammaH                 ->Write();
    CosThetaAfterSignalExtractionErrorsH->Write();
    CosThetaGammaGammaErrorsH           ->Write();
    f.Close();
  }
  else if ( SignalRangeMode == 5 ) {
    TFile f("pngResults/CosThetaHeFrame_DSCB_5.root", "recreate");
    CosThetaAfterSignalExtractionH      ->Write();
    CosThetaGammaGammaH                 ->Write();
    CosThetaAfterSignalExtractionErrorsH->Write();
    CosThetaGammaGammaErrorsH           ->Write();
    f.Close();
  }
  else                             {
    TFile f("pngResults/CosThetaHeFrame_DSCB.root",   "recreate");
    CosThetaAfterSignalExtractionH      ->Write();
    CosThetaGammaGammaH                 ->Write();
    CosThetaAfterSignalExtractionErrorsH->Write();
    CosThetaGammaGammaErrorsH           ->Write();
    f.Close();
  }
  // CosThetaAfterSignalExtractionH      ->Write();
  // CosThetaGammaGammaH                 ->Write();
  // CosThetaAfterSignalExtractionErrorsH->Write();
  // CosThetaGammaGammaErrorsH           ->Write();
  // f.Close();
}
