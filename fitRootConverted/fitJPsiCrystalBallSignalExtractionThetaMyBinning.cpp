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
Double_t  CBparameters[8][5];
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
TF1* JPsiPeakFit     = new TF1( "JPsiPeakFit",    "crystalball",2.2,6);
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
  TF1* CBfit     = new TF1("CBfit","crystalball",2,15);
  // TF1* CBfit     = new TF1("CBfit","crystalball",1.8,7);
  CBfit       ->SetParameter(0,1);
  CBfit       ->SetParameter(3,1.08);
  CBfit       ->SetParameter(4,10);
  // CBfit       ->SetParameter(4,6);
  CBfit       ->SetParLimits(4,2,20);
  // CBfit       ->SetParLimits(4,1,12);
  // CBfit       ->SetParLimits(4,12,99999999);
  CBfit       ->SetParameter(2,0.090);
  CBfit       ->SetParameter(1,3.1);
  CBfit       ->SetNpx(1000);
  TCanvas*    JPsiCanvas = new TCanvas( "JPsiCanvas", "JPsiCanvas", 900, 800 );
  CBfit       ->Draw();
  histoToBeFit->Fit(CBfit, "R");
  // CBfit       ->SetParameter(0,1/CBfit->Integral(2,15));
  // bookKeeping = new Double_t[5];
  for(Int_t i = 0; i < 5; i++){
    CBparameters[0][i] = CBfit->GetParameter(i);
  }
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

  if ( selectionFlag != 0 ) {
    fCohJpsiToMu = (TH1F*)listingsMC[0]->FindObject( Form( "fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH_%d", selectionFlag ) );
  } else {
    fCohJpsiToMu = (TH1F*)listingsMC[0]->FindObject("fInvariantMassDistributionH");
  }
  fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fInvariantMassDistributionH");
  fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fInvariantMassDistributionH");
  fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fInvariantMassDistributionH");
  fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fInvariantMassDistributionH");
  fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fInvariantMassDistributionH");
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
  if( selectionFlag < 10 ){
  fCrystalBallFitPsiPrime(fCohPsi2sToMu);//,     CBparameters[1]);
  // fCrystalBallFitJPsi    (fIncohJpsiToMu,    CBparameters[2]);
  // fCrystalBallFitPsiPrime(fIncohPsi2sToMu,   CBparameters[3]);
  fBkgPolFit             (fTwoGammaToMuHigh);//, CBparameters[4]);
  }
}
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitJPsiTemplate(const int selectionFlag){

  // TH1F *fInvariantMassDistributionH = 0x0;
  fInvariantMassDistributionH = 0x0;
  fInvariantMassDistributionH = (TH1F*)listings->FindObject( Form("fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH_%d", selectionFlag) );
  // fInvariantMassDistributionH->Rebin(2);
  fInvariantMassDistributionH->Rebin(5);
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
  // TF1* JPsiPeakFit     = new TF1( "JPsiPeakFit",    "crystalball",2.2,6);
  // TF1* PsiPrimePeakFit = new TF1( "PsiPrimePeakFit","crystalball",2.2,6);
  // TF1* GammaGammaFit   = new TF1( "GammaGammaFit",
  //                                 "[0]*TMath::Exp(-[1]*x)*( (x > 4) ? 1 : 1 + [2]*(x-4)*(x-4) + [3]*(x-4)*(x-4)*(x-4) + [4]*(x-4)*(x-4)*(x-4)*(x-4) )",
  //                                 2.2,6
  //                                 );

  // for(Int_t i=0; i<5; i++){
  //   for(Int_t j=0; j<5; j++) cout << CBparameters[i][j] << " ";
  //   cout << endl << flush;
  // }

  // JPsiPeakFit->FixParameter(0, CBparameters[0][0]);
  // JPsiPeakFit->FixParameter(3, CBparameters[0][3]);
  // JPsiPeakFit->FixParameter(4, CBparameters[0][4]);
  // JPsiPeakFit->SetParameter(1, CBparameters[0][1]);
  // JPsiPeakFit->SetParLimits(1, CBparameters[0][1]*0.8, CBparameters[0][1]*1.2);
  // JPsiPeakFit->SetParameter(2, CBparameters[0][2]);
  // JPsiPeakFit->SetParLimits(2, CBparameters[0][2]*0.8, CBparameters[0][2]*1.2);
  // PsiPrimePeakFit->FixParameter(0, CBparameters[1][0]);
  // PsiPrimePeakFit->FixParameter(3, CBparameters[1][3]);
  // PsiPrimePeakFit->FixParameter(4, CBparameters[1][4]);
  // PsiPrimePeakFit->FixParameter(1, CBparameters[1][1]);
  // PsiPrimePeakFit->FixParameter(2, JPsiPeakFit->GetParameter(2)*CBparameters[1][2]/CBparameters[0][2]);
  // GammaGammaFit->FixParameter(0, CBparameters[4][0]);
  // GammaGammaFit->FixParameter(2, CBparameters[4][2]);
  // GammaGammaFit->FixParameter(3, CBparameters[4][3]);
  // GammaGammaFit->FixParameter(4, CBparameters[4][4]);
  // GammaGammaFit->SetParameter(1, CBparameters[4][1]);
  // GammaGammaFit->SetParLimits(1, CBparameters[4][1]*0.8, CBparameters[4][1]*1.2);
  // TF1 *fFitInvMass = new TF1("fFitInvMass","[0]*JPsiPeakFit+[1]*PsiPrimePeakFit+[2]*GammaGammaFit",2.2001,5.9999);



  TF1 *fFitInvMass = new TF1("fFitInvMass",fsum,1.8,8,18);
  fFitInvMass->SetNpx(1000000);
  fFitInvMass->FixParameter(0, CBparameters[0][0]);    // best
  // fFitInvMass->FixParameter(3, CBparameters[0][3]);    // best
  // fFitInvMass->FixParameter(4, CBparameters[0][4]);    // best
  fFitInvMass->SetParameter(3, CBparameters[0][3]);
  fFitInvMass->SetParLimits(3, CBparameters[0][3]*0.8, CBparameters[0][3]*1.2);
  fFitInvMass->SetParameter(4, CBparameters[0][4]);
  fFitInvMass->SetParLimits(4, CBparameters[0][4]*0.8, CBparameters[0][4]*1.2);
  // fFitInvMass->FixParameter(1, CBparameters[0][1]);
  // fFitInvMass->FixParameter(2, CBparameters[0][2]);
  fFitInvMass->SetParameter(1, CBparameters[0][1]);
  fFitInvMass->SetParLimits(1, CBparameters[0][1]*0.8, CBparameters[0][1]*1.2);
  fFitInvMass->SetParameter(2, CBparameters[0][2]);
  fFitInvMass->SetParLimits(2, CBparameters[0][2]*0.8, CBparameters[0][2]*1.2);
  fFitInvMass->FixParameter(0+5, CBparameters[1][0]);
  fFitInvMass->FixParameter(3+5, CBparameters[1][3]);
  fFitInvMass->FixParameter(4+5, CBparameters[1][4]);
  fFitInvMass->FixParameter(1+5, CBparameters[1][1]);
  fFitInvMass->FixParameter(2+5, fFitInvMass->GetParameter(2)*CBparameters[1][2]/CBparameters[0][2]);
  // fFitInvMass->FixParameter(0+10, CBparameters[4][0]);    // best
  // fFitInvMass->FixParameter(2+10, CBparameters[4][2]);    // best
  // fFitInvMass->FixParameter(3+10, CBparameters[4][3]);    // best
  // fFitInvMass->FixParameter(4+10, CBparameters[4][4]);    // best
  fFitInvMass->SetParameter(0+10, CBparameters[4][0]);    // mmmh
  fFitInvMass->SetParameter(2+10, CBparameters[4][2]);    // mmmh
  fFitInvMass->SetParameter(3+10, CBparameters[4][3]);    // mmmh
  fFitInvMass->SetParameter(4+10, CBparameters[4][4]);    // mmmh
  fFitInvMass->SetParameter(1+10, CBparameters[4][1]);
  // fFitInvMass->FixParameter(1+10, CBparameters[4][1]);
  fFitInvMass->SetParLimits(1+10, CBparameters[4][1]*0.9, CBparameters[4][1]*1.1);
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
  fFitInvMass->SetParameter(15, 1);
  fFitInvMass->SetParameter(16, 1);
  fFitInvMass->SetParameter(17, 1);
  // if(selectionFlag == 19) {
  //   fFitInvMass->SetParLimits(15, helpppppp * 1.2, 9999999999);
  // } else {
  //   fFitInvMass->SetParLimits(15, 0.001, 9999999999);
  // }
  fFitInvMass->SetParLimits(15, 0.001, 9999999999);
  fFitInvMass->SetParLimits(16, 0.001, 9999999999);
  fFitInvMass->SetParLimits(17, 0.001, 9999999999);
  fFitInvMass->Print();
  for(Int_t i = 0; i < 18; i++) cout << fFitInvMass->GetParameter(i) << endl << flush;


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
  JPsiPeakFit->FixParameter( 0, fFitInvMass->GetParameter(0)*fFitInvMass->GetParameter(15) );
  JPsiPeakFit->FixParameter( 1, fFitInvMass->GetParameter(1) );
  JPsiPeakFit->FixParameter( 2, fFitInvMass->GetParameter(2) );
  JPsiPeakFit->FixParameter( 3, fFitInvMass->GetParameter(3) );
  JPsiPeakFit->FixParameter( 4, fFitInvMass->GetParameter(4) );
  PsiPrimePeakFit->FixParameter( 0, fFitInvMass->GetParameter(5)*fFitInvMass->GetParameter(16) );
  PsiPrimePeakFit->FixParameter( 1, fFitInvMass->GetParameter(1+5) );
  PsiPrimePeakFit->FixParameter( 2, fFitInvMass->GetParameter(2+5) );
  PsiPrimePeakFit->FixParameter( 3, fFitInvMass->GetParameter(3+5) );
  PsiPrimePeakFit->FixParameter( 4, fFitInvMass->GetParameter(4+5) );
  GammaGammaFit->FixParameter( 0, fFitInvMass->GetParameter(10)*fFitInvMass->GetParameter(17) );
  GammaGammaFit->FixParameter( 1, fFitInvMass->GetParameter(1+10) );
  GammaGammaFit->FixParameter( 2, fFitInvMass->GetParameter(2+10) );
  GammaGammaFit->FixParameter( 3, fFitInvMass->GetParameter(3+10) );
  GammaGammaFit->FixParameter( 4, fFitInvMass->GetParameter(4+10) );
  JPsiPeakFit    ->Draw("SAME");
  PsiPrimePeakFit->Draw("SAME");
  GammaGammaFit  ->Draw("SAME");
  if(selectionFlag == 18) helpppppp = fFitInvMass->GetParameter(15);
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
  TH1F* fCohJpsiToMuFromModelH = (TH1F*) JPsiPeakFit->GetHistogram()->Clone("fCohJpsiToMuFromModelH");
  for (Int_t ibin=1; ibin<=fCohJpsiToMuFromModelH->GetNbinsX(); ibin++) {
    fCohJpsiToMuFromModelH->SetBinError(ibin,0);
  }
  // fCohJpsiToMuFromModelH->Rebin(5);
  // fCohJpsiToMuFromModelH->Scale(0.2);
  // fCohJpsiToMuFromModelH->Rebin(4);
  // fCohJpsiToMuFromModelH->Scale(0.25);
  if( selectionFlag < 16 || selectionFlag > 23 ) {
    fCohJpsiToMuFromModelH->Rebin(2);
    fCohJpsiToMuFromModelH->Scale(0.5);
  }
  // fCohJpsiToMuFromModelH->Rebin(8);
  // fCohJpsiToMuFromModelH->Scale(0.125);
  TH1F* fCohPsi2sToMuFromModelH = (TH1F*) PsiPrimePeakFit->GetHistogram()->Clone("fCohPsi2sToMuFromModelH");
  for (Int_t ibin=1; ibin<=fCohPsi2sToMuFromModelH->GetNbinsX(); ibin++) {
    fCohPsi2sToMuFromModelH->SetBinError(ibin,0);
  }
  // fCohPsi2sToMuFromModelH->Scale(0.2);
  // fCohPsi2sToMuFromModelH->Rebin(5);
  // fCohPsi2sToMuFromModelH->Scale(0.25);
  // fCohPsi2sToMuFromModelH->Rebin(4);
  if( selectionFlag < 16 || selectionFlag > 23 ) {
    fCohPsi2sToMuFromModelH->Scale(0.5);
    fCohPsi2sToMuFromModelH->Rebin(2);
  }
  // fCohPsi2sToMuFromModelH->Scale(0.125);
  // fCohPsi2sToMuFromModelH->Rebin(8);
  TH1F* fTwoGammaFromModelH = (TH1F*) GammaGammaFit->GetHistogram()->Clone("fTwoGammaFromModelH");
  for (Int_t ibin=1; ibin<=fTwoGammaFromModelH->GetNbinsX(); ibin++) {
    fTwoGammaFromModelH->SetBinError(ibin,0);
  }
  // fTwoGammaFromModelH->Scale(0.2);
  // fTwoGammaFromModelH->Rebin(5);
  // fTwoGammaFromModelH->Scale(0.25);
  // fTwoGammaFromModelH->Rebin(4);
  if( selectionFlag < 16 || selectionFlag > 23 ) {
    fTwoGammaFromModelH->Scale(0.5);
    fTwoGammaFromModelH->Rebin(2);
  }
  // fTwoGammaFromModelH->Scale(0.125);
  // fTwoGammaFromModelH->Rebin(8);


  Double_t numberOfTotalJPsi     = 0;
  Double_t numberOfTotalPsi2s    = 0;
  Double_t numberOfTotalBkg      = 0;
  Double_t numberOfTotalJPsiErr  = 0;
  Double_t numberOfTotalPsi2sErr = 0;
  Double_t numberOfTotalBkgErr   = 0;
  // if ( ptrSelectionFlag == 2 ) {
  //   numberOfTotalJPsi  = fIncohJpsiToMuC  -> Integral();
  //   numberOfTotalPsi2s = fIncohPsi2sToMuC -> Integral();
  // } else {
    // numberOfTotalJPsi  = fCohJpsiToMuFromModelH -> Integral();
    // numberOfTotalPsi2s = fCohPsi2sToMuFromModelH-> Integral();
    // numberOfTotalJPsi  = JPsiPeakFit    -> Integral(2.2,6,1.E-15);
    // numberOfTotalPsi2s = PsiPrimePeakFit-> Integral(2.2,6,1.E-15);
    // numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.04;
    // numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.04;
    // if( selectionFlag < 16 || selectionFlag > 23 ) {
      numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.05;
      numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.05;
    // } else {
    //   numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.01;
    //   numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.01;
    // }
    // numberOfTotalJPsi     = (JPsiPeakFit    -> Integral(2.2,6))/0.05;  // USUAL FIT
    // numberOfTotalPsi2s    = (PsiPrimePeakFit-> Integral(2.2,6))/0.05;  // USUAL FIT
    numberOfTotalJPsiErr  = numberOfTotalJPsi *fFitInvMass->GetParError(15)/fFitInvMass->GetParameter(15);
    numberOfTotalPsi2sErr = numberOfTotalPsi2s*fFitInvMass->GetParError(16)/fFitInvMass->GetParameter(16);

  // }
  // numberOfTotalBkg = fTwoGammaFromModelH      -> Integral();
  // numberOfTotalBkg = GammaGammaFit-> Integral(2.2,6,1.E-15);
  numberOfTotalBkg    = (GammaGammaFit-> Integral(2.2,6))/0.05;
  numberOfTotalBkgErr = numberOfTotalBkg*fFitInvMass->GetParError(17)/fFitInvMass->GetParameter(17);
  latex->DrawLatex(0.55,0.66,Form("N_{J/#psi} = %.0f #pm %.0f",        numberOfTotalJPsi,  numberOfTotalJPsiErr ));//fFitInvMass->GetParameter(0) *fFitInvMass->GetParError(15)/0.05 ) );
  latex->DrawLatex(0.55,0.60,Form("N_{#psi(2S)} = %.0f #pm %.0f",      numberOfTotalPsi2s, numberOfTotalPsi2sErr));//fFitInvMass->GetParameter(5) *fFitInvMass->GetParError(16)/0.05 ) );
  latex->DrawLatex(0.55,0.54,Form("N_{#gamma#gamma} = %.0f #pm %.0f",  numberOfTotalBkg,   numberOfTotalBkgErr  ));//fFitInvMass->GetParameter(10)*fFitInvMass->GetParError(17)/0.05 ) );
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
  // if ( ptrSelectionFlag == 2 ) {
  //   JPsiPeakSignal     = fIncohJpsiToMuC  -> Integral( fIncohJpsiToMuC ->GetXaxis()->FindBin(2.75), fIncohJpsiToMuC ->GetXaxis()->FindBin(3.45) );
  //   Psi2JPsiPeakSignal = fIncohPsi2sToMuC -> Integral( fIncohPsi2sToMuC->GetXaxis()->FindBin(3.45), fIncohPsi2sToMuC->GetXaxis()->FindBin(3.90) );
  // } else {
    // JPsiPeakSignal     = fCohJpsiToMuFromModelH -> Integral(fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(2.75), fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(3.45));
    // Psi2JPsiPeakSignal = fCohPsi2sToMuFromModelH-> Integral(fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.45), fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.90));
    // JPsiPeakSignal     = fCohJpsiToMuFromModelH -> Integral(fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(2.75), fCohJpsiToMuFromModelH ->GetXaxis()->FindBin(3.45));
    // Psi2JPsiPeakSignal = fCohPsi2sToMuFromModelH-> Integral(fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.45), fCohPsi2sToMuFromModelH->GetXaxis()->FindBin(3.90));
  // }
  JPsiPeakBkg     = GammaGammaFit->Integral(2.75,3.45);
  Psi2JPsiPeakBkg = GammaGammaFit->Integral(3.45,3.90);

  JPsiPeakValue    = numberOfTotalJPsi;
  JPsiPeakValueErr = numberOfTotalJPsiErr;
  BkgValue         = JPsiPeakBkg;
  BkgValueError    = JPsiPeakBkg * fFitInvMass->GetParError(15)/fFitInvMass->GetParameter(15);


  // latex->DrawLatex(0.55,0.42,Form("N_{BG J/#psi} = %.0f #pm %.0f",   JPsiPeakBkg,     JPsiPeakBkg     * fFitInvMass->GetParError(17) / numberOfTotalJPsi ));
  // latex->DrawLatex(0.55,0.36,Form("N_{BG #psi(2s)} = %.0f #pm %.0f", Psi2JPsiPeakBkg, Psi2JPsiPeakBkg * fFitInvMass->GetParError(17) / numberOfTotalPsi2s));
  latex->DrawLatex(0.55,0.18,Form("      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                     fFitInvMass->GetChisquare(),
                                     fFitInvMass->GetNDF(),
                                     fFitInvMass->GetChisquare()/fFitInvMass->GetNDF()
                                     )
                                    );


  // //Signal to background ratios
  // Double_t errSB1 = sqrt(   fFitInvMass->GetParError(15) * fFitInvMass->GetParError(15)/( numberOfTotalJPsi  *numberOfTotalJPsi )  +
  //                         ( fFitInvMass->GetParError(17) * fFitInvMass->GetParError(17)/( JPsiPeakBkg        * JPsiPeakBkg)) );
  // Double_t errSB2 = sqrt(   fFitInvMass->GetParError(16) * fFitInvMass->GetParError(16)/( numberOfTotalPsi2s * numberOfTotalPsi2s )  +
  //                           fFitInvMass->GetParError(17) * fFitInvMass->GetParError(17)/( Psi2JPsiPeakBkg    * Psi2JPsiPeakBkg ) );
  // latex->DrawLatex(0.55,0.30,Form("S/B for J/#psi = %.2f #pm %.2f",   JPsiPeakSignal/(Double_t)JPsiPeakBkg,         errSB1*JPsiPeakSignal/(Double_t)JPsiPeakBkg));
  // latex->DrawLatex(0.55,0.24,Form("S/B for #psi(2s) = %.2f #pm %.2f", Psi2JPsiPeakSignal/(Double_t)Psi2JPsiPeakBkg, errSB2*Psi2JPsiPeakSignal/(Double_t)Psi2JPsiPeakBkg));
  // if      ( selectionFlag == 0 ) gPad->SaveAs("pngResults/PtIntJPsiTF1.png",      "RECREATE");
  // else if ( selectionFlag == 1 ) gPad->SaveAs("pngResults/CoherentJPsiTF1.png",   "RECREATE");
  // else if ( selectionFlag == 2 ) gPad->SaveAs("pngResults/IncoherentJPsiTF1.png", "RECREATE");
  // else                           gPad->SaveAs("pngResults/ops.png",               "RECREATE");
  //
  // TLegend* l = new TLegend(0.2,0.55,0.5,0.85);
  // l->SetMargin(0.1);
  // l->SetBorderSize(0);
  // l->AddEntry(  fInvariantMassDistributionH, "ALICE data 2018");
  // l->AddEntry(  fFitInvMass, Form( "Fit: #chi^{2}/NDF = %.2f / %.2d = %.2f  ",
  //                                   fFitInvMass->GetChisquare(),
  //                                   fFitInvMass->GetNDF(),
  //                                   fFitInvMass->GetChisquare()/fFitInvMass->GetNDF()
  //                                   )
  //                                  );
  // if ( selectionFlag == 2 ) {
  //   l->AddEntry(  JPsiPeakFit,        "Incoherent J/#psi");
  //   l->AddEntry(  PsiPrimePeakFit,    "Incoherent #psi(2S)");
  // } else {
  //   l->AddEntry(  JPsiPeakFit,        "Coherent   J/#psi");
  //   l->AddEntry(  PsiPrimePeakFit,    "Coherent   #psi(2S)");
  // }
  // // l->AddEntry(  fTwoGammaToMuMediumC, "Continuum  #gamma#gamma #rightarrow #mu#mu");
  // l->AddEntry(  GammaGammaFit,        "Continuum  #gamma#gamma #rightarrow #mu#mu");
  // l->Draw();

  gPad->SaveAs(Form("pngResults/FunctionalSignalCosThetaMyBinning_%d.png", selectionFlag), "recreate");


  // cout << "nBackGround:           " << nBackGround.getVal() << endl;
  // cout << "fValueBkgJPsiPeak:     " << fValueBkgJPsiPeak << endl;
  // cout << "fValueBkgPsiPrimePeak: " << fValueBkgPsiPrimePeak << endl;

  // exp = TotalBackGroundJPsiPeak;
  // err = errorTotalBackGroundJPsiPeak;

  new TCanvas;
  fCohJpsiToMuFromModelH->Draw();

  // Double_t ChiSquareSingleContribution[76];
  // TCanvas* ChiSquareCanvas = new TCanvas( "ChiSquareCanvas", "ChiSquareCanvas", 900, 800 );
  // TH1F* ChiSquareH  = new TH1F( "ChiSquareH" ,"ChiSquareH" ,76, 2.2, 6);
  // TH1F* ChiSquareH2 = new TH1F( "ChiSquareH2","ChiSquareH2",76, 2.2, 6);
  // TH1F* ChiSquareH3 = new TH1F( "ChiSquareH3","ChiSquareH3",76, 2.2, 6);
  // // TH1F* fTwoGammaFromModelH = (TH1F*) GammaGammaFit->GetHistogram()->Clone("fTwoGammaFromModelH");
  // // for (Int_t ibin=1; ibin<=fTwoGammaFromModelH->GetNbinsX(); ibin++) {
  // //   fTwoGammaFromModelH->SetBinError(ibin,0);
  // // }
  // Int_t whichBin  = fInvariantMassDistributionH->GetXaxis()->FindBin(2.225);
  // Int_t whichBin2 = ChiSquareH                 ->GetXaxis()->FindBin(2.225);
  // for( Int_t iLoop = 0; iLoop < 76; iLoop++ ) {
  //   // Double_t binContent = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   Double_t binContent = fInvariantMassDistributionH->GetBinContent(whichBin + iLoop);
  //   Double_t binError   = fInvariantMassDistributionH->GetBinError(whichBin + iLoop);
  //   Double_t binCenter  = fInvariantMassDistributionH->GetXaxis()->GetBinCenter(whichBin + iLoop);
  //   // Double_t fitValue   = fFitInvMass->Eval(2.225+iLoop*0.5);
  //   Double_t fitValue   = fFitInvMass->Eval(binCenter);
  //   cout << "binContent = " << binContent << endl;
  //   cout << "fitValue   = " << fitValue   << endl;
  //   // Double_t pull = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5)) - fFitInvMass->Eval(2.225+iLoop*0.5);
  //   Double_t pull = binContent - fitValue;
  //   cout << pull << endl;
  //   // Double_t help = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   // Double_t help = fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.225+iLoop*0.5));
  //   Double_t pullHelp;
  //   if( pull != 0 && binContent != 0 ) pullHelp = pull /binContent;
  //   else pullHelp = 0;
  //   cout << "pullHelp = " << pullHelp << endl;
  //   Double_t pullHelp2;
  //   if( pull != 0 && binError != 0 ) pullHelp2 = pull /binError;
  //   else pullHelp2 = 0;
  //   cout << "pullHelp2 = " << pullHelp2 << endl;
  //   cout << "binCenter = " << ChiSquareH                 ->GetBinCenter( whichBin2 + iLoop ) << endl;
  //   cout << "binCenter = " << fInvariantMassDistributionH->GetBinCenter( whichBin  + iLoop ) << endl;
  //
  //   ChiSquareH->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitInvMass->Eval(2.2+iLoop*0.5)
  //                              pullHelp
  //                              // pull
  //                              );
  //   ChiSquareH2->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitInvMass->Eval(2.2+iLoop*0.5)
  //                              pullHelp2
  //                              // pull
  //                              );
  //   ChiSquareH3->SetBinContent(
  //                              // ChiSquareH->GetXaxis()->FindBin(2.2+iLoop*0.5),
  //                              whichBin2 + iLoop,
  //                              // fInvariantMassDistributionH->GetBinContent(fInvariantMassDistributionH->GetXaxis()->FindBin(2.2+iLoop*0.5)) - fFitInvMass->Eval(2.2+iLoop*0.5)
  //                              pull*pull/binContent
  //                              // pull
  //                              );
  //   ChiSquareSingleContribution[iLoop] = pull*pull/binContent;
  //
  //
  // }
  // Double_t AllChiSquare = 0;
  // for( Int_t i = 0; i < 76; i++){
  //   cout << "ChiSquareSingleContribution = " << ChiSquareSingleContribution[i] << endl;
  //   AllChiSquare += ChiSquareSingleContribution[i];
  // }
  // cout << "AllChiSquare = " << AllChiSquare << endl;
  // // ChiSquareH->Draw();
  // // TCanvas* ChiSquareCanvas2 = new TCanvas( "ChiSquareCanvas2", "ChiSquareCanvas2", 900, 800 );
  // // ChiSquareH2->Draw();
  // TCanvas* ChiSquareCanvas3 = new TCanvas( "ChiSquareCanvas3", "ChiSquareCanvas3", 900, 800 );
  // ChiSquareH3->Draw();
  // ChiSquareCanvas3->SaveAs(Form("pngResults/FunctionalSignalCosThetaChi_%d.png", selectionFlag), "recreate");

}
//_____________________________________________________________________________
/* - Here I create the new TH1 for the after the signal extraction.
 * - Basically I run the fit function many times and then I memorise
 * - the values each time. After that I fill with a setbincontent
 * - and a setbinerror the
 * -
 */
void CreateCosThetaTh1(const char* AnalysisName){
  // fileMC[0] = new TFile("MCtrainResults/2019-05-17/kCohJpsiToMu/AnalysisResults.root");
  // fileMC[0] = new TFile("AnalysisResultsMC20062019.root");
  fileMC[0] = new TFile("MCtrainResults/2019-06-24/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2019-06-24/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2019-06-24/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2019-06-24/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2019-06-24/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2019-06-24/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2019-06-24/kTwoGammaToMuHigh/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2019-06-24/kTwoGammaToMuMedium/AnalysisResults.root");
  // fileMC[0] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kCohJpsiToMu/AnalysisResults.root");
  // fileMC[1] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kCohPsi2sToMu/AnalysisResults.root");
  // fileMC[2] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kCohPsi2sToMuPi/AnalysisResults.root");
  // fileMC[3] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kIncohJpsiToMu/AnalysisResults.root");
  // fileMC[4] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kIncohPsi2sToMu/AnalysisResults.root");
  // fileMC[5] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kIncohPsi2sToMuPi/AnalysisResults.root");
  // fileMC[6] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kTwoGammaToMuHigh/AnalysisResults.root");
  // fileMC[7] = new TFile("MCtrainResults/2019-06-16-LHC18qr/kTwoGammaToMuMedium/AnalysisResults.root");
  // TDirectory* dirMC[8];
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

  fitJPsiTemplateMC();
  fileList = new TFile(AnalysisName);
  dir      = fileList->GetDirectory("MyTask");
  dir->GetObject("MyOutputContainer", listings);


  // Double_t CosThetaBinCenters[] = { -0.9, -0.7, -0.5, -0.3, -0.1,
  //                                    0.1,  0.3,  0.5,  0.7,  0.9 };
  // Double_t PhiBinCenters[]      = { -2.826, -2.198, -1.57, -0.942, -0.314,
  //                                    0.314,  0.942,  1.57,  2.198,  2.826};
  Double_t CosThetaBinCenters[26] = { -0.575,  -0.45,   -0.375,  -0.325,  -0.275,
                                      -0.225,  -0.175,  -0.1375, -0.1125, -0.0875,
                                      -0.0625, -0.0375, -0.0125,  0.0125,  0.0375,
                                       0.0625,  0.0875,  0.1125,  0.1375,  0.175,
                                       0.225,   0.275,   0.325,   0.375,   0.45,
                                       0.575 };
  Double_t PhiBinCenters[30] = {        -3.14*14.5/15,  -3.14*13.5/15,  -3.14*12.5/15,
                                        -3.14*11.5/15,  -3.14*10.5/15,  -3.14*9.5/15,   -3.14*8.5/15,
                                        -3.14*7.5/15,   -3.14*6.5/15,   -3.14*5.5/15,   -3.14*4.5/15,
                                        -3.14*3.5/15,   -3.14*2.5/15,   -3.14*1.5/15,   -3.14*0.5/15,
                                        +3.14*0.5/15,   +3.14*1.5/15,   +3.14*2.5/15,   +3.14*3.5/15,
                                        +3.14*4.5/15,   +3.14*5.5/15,   +3.14*6.5/15,   +3.14*7.5/15,
                                        +3.14*8.5/15,   +3.14*9.5/15,   +3.14*10.5/15,  +3.14*11.5/15,
                                        +3.14*12.5/15,  +3.14*13.5/15,  +3.14*14.5/15 };
  Double_t MyVariableCosThetaBinning1D[] = { -0.65,  -0.5,  -0.4,  -0.35,  -0.3,
                                             -0.25,  -0.2,  -0.15, -0.125, -0.1,
                                             -0.075, -0.05, -0.025, 0,      0.025,
                                              0.05,   0.075, 0.1,   0.125,  0.15,
                                              0.2,    0.25,  0.3,   0.35,   0.4,
                                              0.5,    0.65 };
  Double_t MyVariablePhiBinning1D[] = { -3.14*1,      -3.14*14/15,  -3.14*13/15,  -3.14*12/15,
                                        -3.14*11/15,  -3.14*10/15,  -3.14*9/15,   -3.14*8/15,
                                        -3.14*7/15,   -3.14*6/15,   -3.14*5/15,   -3.14*4/15,
                                        -3.14*3/15,   -3.14*2/15,   -3.14*1/15,    0,
                                        +3.14*1/15,   +3.14*2/15,   +3.14*3/15,   +3.14*4/15,
                                        +3.14*5/15,   +3.14*6/15,   +3.14*7/15,   +3.14*8/15,
                                        +3.14*9/15,   +3.14*10/15,  +3.14*11/15,  +3.14*12/15,
                                        +3.14*13/15,  +3.14*14/15,  +3.14*1 };
  Double_t Widths[] = { 0.15,  0.1,   0.05,  0.05,  0.05,  0.05,  0.05,  0.025, 0.025,
                        0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025,
                        0.025, 0.05,  0.05,  0.05,  0.05,  0.05,  0.1,   0.15 };
  // for( Int_t iFillingCenters = 0; iFillingCenters < 50; iFillingCenters++ ){
  //   if( iFillingCenters < 40 ) {
  //     CosThetaBinCenters[iFillingCenters] = -0.975 + (Double_t)iFillingCenters * 0.025;
  //   }
  //   PhiBinCenters[iFillingCenters] = -3.14 + (Double_t)iFillingCenters * 0.1256;
  // }
  for( Int_t looping = 0; looping < 26; looping++ ) {cout << "CosThetaBinCenters  " << CosThetaBinCenters[looping] << endl; }
  for( Int_t looping = 0; looping < 30; looping++ ) {cout << "PhiBinCenters       " << PhiBinCenters[looping]      << endl; }

  TH1F* CosThetaAfterSignalExtractionH =
            new TH1F( "CosThetaAfterSignalExtractionH",
                      "CosThetaAfterSignalExtractionH",
                      // 40, -1, 1//, 10, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      26, MyVariableCosThetaBinning1D
                      );
  TH1F* CosThetaGammaGammaH =
            new TH1F( "CosThetaGammaGammaH",
                      "CosThetaGammaGammaH",
                      // 40, -1, 1//, 10, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      26, MyVariableCosThetaBinning1D
                      );

  TH1F* CosThetaAfterSignalExtractionErrorsH =
            new TH1F( "CosThetaAfterSignalExtractionErrorsH",
                      "CosThetaAfterSignalExtractionErrorsH",
                      // 40, -1, 1//, 10, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      26, MyVariableCosThetaBinning1D
                      );
  TH1F* CosThetaGammaGammaErrorsH =
            new TH1F( "CosThetaGammaGammaErrorsH",
                      "CosThetaGammaGammaErrorsH",
                      // 40, -1, 1//, 10, -3.14, 3.14
                      // 10, -1, 1, 10, -4, 4
                      26, MyVariableCosThetaBinning1D
                      );

  for (size_t iCosThetaBins = 0; iCosThetaBins < 26; iCosThetaBins++) {
    // for (size_t iPhiBins = 0; iPhiBins < 10; iPhiBins++) {
      JPsiPeakValue    = 0;
      JPsiPeakValueErr = 0;
      BkgValue         = 0;
      BkgValueError    = 0;
      // if( iCosThetaBins > 17 && iCosThetaBins < 23 ) {
      //   fitJPsiTemplateMC(iCosThetaBins);
      // } else {
      //   fitJPsiTemplateMC();
      // }
      fitJPsiTemplate(iCosThetaBins);

      CosThetaAfterSignalExtractionH->Fill( // -0.975 + (Double_t)iCosThetaBins * 0.05,
                                            // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                            // -2.826 + (Double_t)iPhiBins      * 0.628,
                                            // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                            CosThetaBinCenters[iCosThetaBins],
                                            // JPsiPeakValue/( MyVariableCosThetaBinning1D[iCosThetaBins+1] - MyVariableCosThetaBinning1D[iCosThetaBins] )
                                            JPsiPeakValue/Widths[iCosThetaBins]
                                            );
      CosThetaGammaGammaH->Fill( // -0.975 + (Double_t)iCosThetaBins * 0.05,
                                 // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                 // -2.826 + (Double_t)iPhiBins      * 0.628,
                                 // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                 CosThetaBinCenters[iCosThetaBins],
                                 BkgValue/Widths[iCosThetaBins]
                                 );

      CosThetaAfterSignalExtractionErrorsH->Fill(  // -0.975 + (Double_t)iCosThetaBins * 0.05,
                                                   // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                                   // -2.826 + (Double_t)iPhiBins      * 0.628,
                                                   // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                                   CosThetaBinCenters[iCosThetaBins],
                                                   // JPsiPeakValue/( MyVariableCosThetaBinning1D[iCosThetaBins+1] - MyVariableCosThetaBinning1D[iCosThetaBins] )
                                                   JPsiPeakValue/Widths[iCosThetaBins]
                                                   );
      CosThetaGammaGammaErrorsH->Fill(  // -0.975 + (Double_t)iCosThetaBins * 0.05,
                                        // -3.6 + (Double_t)iPhiBins      * 8.0 / 10.0,
                                        // -2.826 + (Double_t)iPhiBins      * 0.628,
                                        // -3.14 + 0.314 * ( 2.0 * (Double_t)iPhiBins + 1.0) ,
                                        CosThetaBinCenters[iCosThetaBins],
                                        // BkgValue/( MyVariableCosThetaBinning1D[iCosThetaBins+1] - MyVariableCosThetaBinning1D[iCosThetaBins] )
                                        BkgValue/Widths[iCosThetaBins]
                                        );
      CosThetaAfterSignalExtractionErrorsH->SetBinError(  iCosThetaBins + 1 ,
                                                          // iPhiBins      + 1 ,
                                                          // JPsiPeakValueErr/( MyVariableCosThetaBinning1D[iCosThetaBins+1] - MyVariableCosThetaBinning1D[iCosThetaBins] )
                                                          JPsiPeakValueErr/Widths[iCosThetaBins]
                                                          );
      CosThetaGammaGammaErrorsH->SetBinError(  iCosThetaBins + 1 ,
                                               // iPhiBins      + 1 ,
                                               // BkgValueError/( MyVariableCosThetaBinning1D[iCosThetaBins+1] - MyVariableCosThetaBinning1D[iCosThetaBins] )
                                               BkgValueError/Widths[iCosThetaBins]
                                               );

    // }
  }

  TFile f("pngResults/TH1functionalCosThetaMyBinningEX.root", "recreate");
  CosThetaAfterSignalExtractionH      ->Write();
  CosThetaGammaGammaH                 ->Write();
  CosThetaAfterSignalExtractionErrorsH->Write();
  CosThetaGammaGammaErrorsH           ->Write();
  f.Close();
}
