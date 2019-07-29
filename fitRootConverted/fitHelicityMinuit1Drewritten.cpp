#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <math.h>
#include "TH2D.h"
#include "TF2.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"

#include <vector>
#include <map>

Int_t switchFlag = 0;


//_____________________________________________________________________________
/* - Coding in the fit functions.
   - The fit is modelled as the sum of 3 different functions:
   - 1) CosTheta only
   - 2) Phi      only
   - 3) Mix of the two
   -
 */
Double_t CosTheta(Double_t *x, Double_t *par) {
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t returnValue     = 1 + par[0] * CosSquaredTheta;
  returnValue              = par[1] * returnValue / ( 3 + par[0] );
  return   returnValue;
}
//______________________________________________
Double_t Phi(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi = TMath::Cos( 2 * x[1] );
  Double_t returnValue = par[2] * ( 1 + 2 * par[3] * CosOfTwoPhi / ( 3 + par[0] ) );
  return   returnValue;
}
//______________________________________________
Double_t PhiV2(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi = TMath::Cos( 2 * x[0] );
  Double_t returnValue = par[2] * ( 1 + 2 * par[3] * CosOfTwoPhi / ( 3 + par[0] ) );
  return   returnValue;
}
//______________________________________________
Double_t DummyPhi(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi = TMath::Cos( 2 * x[1] );
  Double_t returnValue = par[0] * ( 1 + 2 * par[1] * CosOfTwoPhi / ( 3 + 1.13220 ) );
  return   returnValue;
}
//______________________________________________
Double_t MixWithPositiveCosTheta(Double_t *x, Double_t *par) {
  Double_t CosTwoTildePhi = TMath::Cos( 2 * x[1] - 0.50 * TMath::Pi() );
  Double_t returnValue    = par[4] * ( 1 + TMath::Sqrt(2) * par[5] * CosTwoTildePhi / ( 3 + par[0] ) );
  return   returnValue;
}
//______________________________________________
Double_t MixWithNegativeCosTheta(Double_t *x, Double_t *par) {
  Double_t CosTwoTildePhi = TMath::Cos( 2 * x[1] - 1.50 * TMath::Pi() );
  Double_t returnValue    = par[4] * ( 1 + TMath::Sqrt(2) * par[5] * CosTwoTildePhi / ( 3 + par[0] ) );
  return   returnValue;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t SimultaneousFitComplete(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  if ( x[0] < 0 ) {
    sumOfTheSubFits = CosTheta( x, par ) + Phi( x, par ) + MixWithNegativeCosTheta( x, par );
  } else {
    sumOfTheSubFits = CosTheta( x, par ) + Phi( x, par ) + MixWithPositiveCosTheta( x, par );
  }
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t SimultaneousFit(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = CosTheta( x, par ) + Phi( x, par );
  // if( switchFlag == 0 ) {
  //   sumOfTheSubFits = CosTheta( x, par );
  //   return sumOfTheSubFits;
  // } else {
  //   sumOfTheSubFits = Phi( x, par );
  //   return sumOfTheSubFits;
  // }
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t SimultaneousFitLastHope(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = ( x[0] < 2*TMath::Pi() ) ? CosTheta( x, par ) : PhiV2( x, par );
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Dummy fit CosTheta only.
 * -
 */
Double_t DummyFitCosTheta(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = CosTheta( x, par );
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t DummyFitPhi(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = DummyPhi( x, par );
  // if( switchFlag == 0 ) {
  //   sumOfTheSubFits = CosTheta( x, par );
  //   return sumOfTheSubFits;
  // } else {
  //   sumOfTheSubFits = Phi( x, par );
  //   return sumOfTheSubFits;
  // }
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Data needed for the global fit.
 * - They have to be global to be visible.
 */
// std::vector< std::pair< Double_t, Double_t > > coords;
std::vector< Double_t > coords;
std::vector< Double_t > values;
std::vector< Double_t > errors;

// void FcnForMinimisation(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
void FcnForMinimisation(Int_t &npar, Double_t *gin, Double_t &f, Double_t *p, Int_t iflag)
{
  // cout << "HI" << flush << endl;
  Int_t n = coords.size();
  // cout << "n is " << n << flush << endl;

  // Int_t n = 40;
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  for ( Int_t i = 0; i < n; ++i ) {
    // // if ( i < 40 ) {
    if ( i < 12 ) {
      x[0] = coords[i];
      x[1] = 0;
      // // cout << "x[0] = " << x[0] << endl;
      // // switchFlag = 0;
    } else {
      x[0] = 0;
      x[1] = coords[i];
      // switchFlag = 1;
    }
    // x[0] = 0;
    // x[1] = coords[i];

    // cout << "HI2" << flush << endl;
    if ( values[i] != 0 ) {
      // tmp = ( values[i] - SimultaneousFit( x, p ) ) / errors[i];
      tmp = ( values[i] - DummyFitCosTheta( x, p ) ) / errors[i];
      // tmp = ( values[i] - DummyFitPhi( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  f = chi2;
  cout << "ChiSquared = " << chi2 << endl;
}
//_____________________________________________________________________________
void FcnForMinimisationV2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *p, Int_t iflag)
{
  Int_t n = coords.size();
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  for ( Int_t i = 0; i < n; ++i ) {
    if ( i < 12 ) {
      x[0] = coords[i];
      x[1] = 0;
    } else {
      x[0] = coords[i] + 4*TMath::Pi();
      x[1] = 0;
    }
    if ( values[i] != 0 ) {
      tmp = ( values[i] - SimultaneousFitLastHope( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  f = chi2;
  cout << "ChiSquared = " << chi2 << endl;
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void fitHelicitySimultaneously1d(){

  TFile* file1D = new TFile("pngResults/TH1corr.root");
  // TFile* file2D = new TFile("pngResults/TH2corr.root");

  // TFile* fileNew = new TFile("pngResults/TH1corrMyBinning.root");
  TFile* fileNew = new TFile("pngResults/TH1corr25bins.root");
  //
  TH1F* CorrectedCosTheta = (TH1F*) fileNew->Get("RawCosTheta2H");
  // TH1F* CorrectedCosTheta = (TH1F*) file1D->Get("RawCosThetaH");
  TH1F* CorrectedPhi      = (TH1F*) file1D->Get("RawPhiH");


  Double_t CosThetaLowLimit   = -1;
  Double_t CosThetaUpperLimit = +1;
  Double_t PhiLowLimit        = -3.14;
  Double_t PhiUpperLimit      = +3.14;


  // TF2 * helicitySimultaneously1d = new TF2( "helicitySimultaneously1d",
  //                                           my2Dfunc,xlow2,xup2,ylow2,yup2, 10);


  Int_t nBinsCosTheta = CorrectedCosTheta->GetNbinsX();
  Int_t nBinsPhi      = CorrectedPhi     ->GetNbinsX();


  /// reset data structure
  // coords = std::vector<std::pair<double,double> >();
  coords = std::vector<Double_t>();
  values = std::vector<Double_t>();
  errors = std::vector<Double_t>();
  /// fill data structure
  for (Int_t ix = 8; ix <= nBinsCosTheta-6; ++ix) {
    // coords.push_back( std::make_pair(xaxis1->GetBinCenter(ix), yaxis1->GetBinCenter(iy) ) );
    coords.push_back( CorrectedCosTheta->GetXaxis()->GetBinCenter(ix) );
    values.push_back( CorrectedCosTheta->GetBinContent(ix)            );
    errors.push_back( CorrectedCosTheta->GetBinError(ix)              );
  }
  for (Int_t iy = 1; iy <= nBinsPhi; ++iy) {
    coords.push_back( CorrectedPhi     ->GetXaxis()->GetBinCenter(iy) );
    values.push_back( CorrectedPhi     ->GetBinContent(iy)            );
    errors.push_back( CorrectedPhi     ->GetBinError(iy)              );
  }

  for(int i=0; i < 50; i++){
    cout << i << "  " << coords[i] << "  " << values[i] << endl;
  }


  TMinuit *gMinuit = new TMinuit(4);
  // gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->SetFCN(FcnForMinimisationV2);
  gMinuit->DefineParameter(0, "LambdaTheta", 1., 0.1, -2, 2);
  gMinuit->DefineParameter(1, "NormalTheta", 2.60e+04, 100,  2.58e+04, 2.8e+04);
  // gMinuit->DefineParameter(1, "NormalTheta", 2.52497e+04, 100,  2.50e+04, 2.55e+04);
  gMinuit->DefineParameter(2, "NormalisPhi", 4137, 100,  4000, 4400);
  gMinuit->DefineParameter(3, "LambdaPhi",   0, 0.1, -2, 2);
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MINOS");
  Double_t LambdaTheta,LambdaPhi,NormalTheta,NormalisPhi;
  Double_t LambdaThetaErr,LambdaPhiErr,NormalThetaErr,NormalisPhiErr;
  gMinuit->GetParameter(0,LambdaTheta,LambdaThetaErr);
  gMinuit->GetParameter(1,NormalTheta,NormalThetaErr);
  gMinuit->GetParameter(2,NormalisPhi,NormalisPhiErr);
  gMinuit->GetParameter(3,LambdaPhi,  LambdaPhiErr);
  printf("LambdaTheta: %+.7f +- %.7f\n",LambdaTheta,LambdaThetaErr);
  printf("NormalTheta: %+.7f +- %.7f\n",NormalTheta,NormalThetaErr);
  printf("NormalisPhi: %+.7f +- %.7f\n",NormalisPhi,NormalisPhiErr);
  printf("LambdaPhi  : %+.7f +- %.7f\n",LambdaPhi,  LambdaPhiErr);
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);




  TF1* Model = new TF1("Model", "[1]*(1+[0]*x*x)/(3+[0])", -1 ,1 );
  new TCanvas;
  CorrectedCosTheta->Draw();
  Model->SetParameter( 0, LambdaTheta );
  Model->SetParameter( 1, NormalTheta );
  Model->SetNpx(500);
  Model->Draw("same");

  TF1* Model2 = new TF1("Model2", "[1]*(1+2*[2]*cos(2*x)/(3+[0]))", -3.1 ,3.1 );
  new TCanvas;
  CorrectedPhi->Draw();
  Model2->SetParameter( 0, LambdaTheta );
  Model2->SetParameter( 2, LambdaPhi );
  Model2->SetParameter( 1, NormalisPhi );
  Model2->SetNpx(500);
  Model2->Draw("same");


}
//_____________________________________________________________________________
/* - Fit function for the helicity case.
 * - It is basically a parabolic fit...
 * - Only CosTheta Fit.
 */
void fitOnlyCosTheta(){

  TFile* fileNew           = new TFile("pngResults/TH1corr25bins.root");
  TH1F*  CorrectedCosTheta = (TH1F*) fileNew->Get("RawCosTheta2H");

  Double_t CosThetaLowLimit   = -1;
  Double_t CosThetaUpperLimit = +1;

  Int_t nBinsCosTheta = CorrectedCosTheta->GetNbinsX();

  /// reset data structure
  // coords = std::vector<std::pair<double,double> >();
  coords = std::vector<Double_t>();
  values = std::vector<Double_t>();
  errors = std::vector<Double_t>();
  /// fill data structure
  for (Int_t ix = 7; ix <= nBinsCosTheta-6; ++ix) {
    // coords.push_back( std::make_pair(xaxis1->GetBinCenter(ix), yaxis1->GetBinCenter(iy) ) );
    coords.push_back( CorrectedCosTheta->GetXaxis()->GetBinCenter(ix) );
    values.push_back( CorrectedCosTheta->GetBinContent(ix)            );
    errors.push_back( CorrectedCosTheta->GetBinError(ix)              );
  }

  for(int i=0; i < 18; i++){
    cout << i << "  " << coords[i] << "  " << values[i] << endl;
  }

  TMinuit *gMinuit = new TMinuit(2);
  gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->DefineParameter(0, "LambdaTheta", 1., 0.1, -2, 2);
  gMinuit->DefineParameter(1, "NormalTheta", 2.60e+04, 100,  2.58e+04, 2.62e+04);
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  // gMinuit->Command("MIGRAD");
  gMinuit->Command("MINOS");
  Double_t LambdaTheta,LambdaPhi,NormalTheta,NormalisPhi;
  Double_t LambdaThetaErr,LambdaPhiErr,NormalThetaErr,NormalisPhiErr;
  gMinuit->GetParameter(0,LambdaTheta,LambdaThetaErr);
  gMinuit->GetParameter(1,NormalTheta,NormalThetaErr);
  printf("LambdaTheta: %+.7f +- %.7f\n",LambdaTheta,LambdaThetaErr);
  printf("NormalTheta: %+.7f +- %.7f\n",NormalTheta,NormalThetaErr);
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);

  TF1* Model = new TF1("Model", "[1]*(1+[0]*x*x)/(3+[0])", -1 ,1 );
  new TCanvas;
  CorrectedCosTheta->Draw();
  Model->SetParameter( 0, LambdaTheta );
  Model->SetParameter( 1, NormalTheta );
  Model->SetNpx(500);
  Model->Draw("same");

}
//_____________________________________________________________________________
/* - Fit function for the helicity case.
 * - Flat Phi ONLY...
 * -
 */
void fitOnlyPhi(){

  TFile* file1D = new TFile("pngResults/TH1corr.root");
  TH1F* CorrectedPhi      = (TH1F*) file1D->Get("RawPhiH");

  Double_t PhiLowLimit        = -3.14;
  Double_t PhiUpperLimit      = +3.14;

  Int_t nBinsPhi      = CorrectedPhi     ->GetNbinsX();


  /// reset data structure
  // coords = std::vector<std::pair<double,double> >();
  coords = std::vector<Double_t>();
  values = std::vector<Double_t>();
  errors = std::vector<Double_t>();
  /// fill data structure
  for (Int_t iy = 1; iy <= nBinsPhi; ++iy) {
    coords.push_back( CorrectedPhi     ->GetXaxis()->GetBinCenter(iy) );
    values.push_back( CorrectedPhi     ->GetBinContent(iy)            );
    errors.push_back( CorrectedPhi     ->GetBinError(iy)              );
  }

  for(int i=0; i < 50; i++){
    cout << i << "  " << coords[i] << "  " << values[i] << endl;
  }


  TMinuit *gMinuit = new TMinuit(2);
  gMinuit->SetFCN(FcnForMinimisation);
  // gMinuit->DefineParameter(0, "LambdaTheta", 1., 0.1, -2, 2);
  // gMinuit->DefineParameter(1, "NormalTheta", 2.60e+04, 100,  2.58e+04, 2.62e+04);
  // gMinuit->DefineParameter(1, "NormalTheta", 2.52497e+04, 100,  2.50e+04, 2.55e+04);
  gMinuit->DefineParameter(0, "NormalisPhi", 4000, 100,  3500, 5000);
  gMinuit->DefineParameter(1, "LambdaPhi",   0, 0.1, -2, 2);
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  // gMinuit->Command("MIGRAD");
  // gMinuit->Command("MINOS");
  Double_t LambdaTheta,LambdaPhi,NormalTheta,NormalisPhi;
  Double_t LambdaThetaErr,LambdaPhiErr,NormalThetaErr,NormalisPhiErr;
  // gMinuit->GetParameter(0,LambdaTheta,LambdaThetaErr);
  // gMinuit->GetParameter(1,NormalTheta,NormalThetaErr);
  gMinuit->GetParameter(0,NormalisPhi,NormalisPhiErr);
  gMinuit->GetParameter(1,LambdaPhi,  LambdaPhiErr);
  printf("LambdaTheta: %+.7f +- %.7f\n",LambdaTheta,LambdaThetaErr);
  printf("NormalTheta: %+.7f +- %.7f\n",NormalTheta,NormalThetaErr);
  printf("NormalisPhi: %+.7f +- %.7f\n",NormalisPhi,NormalisPhiErr);
  printf("LambdaPhi  : %+.7f +- %.7f\n",LambdaPhi,  LambdaPhiErr);
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);



  TF1* Model2 = new TF1("Model2", "[1]*(1+2*[2]*cos(2*x)/(3+[0]))", -3.1 ,3.1 );
  new TCanvas;
  CorrectedPhi->Draw();
  Model2->SetParameter( 0, 1.13220 );
  Model2->SetParameter( 1, NormalisPhi );
  Model2->SetParameter( 2, LambdaPhi );
  Model2->SetNpx(500);
  Model2->Draw("same");


}
