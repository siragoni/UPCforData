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



#include "TROOT.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include "TH1F.h"




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
/* - Data needed for the global fit.
 * - They have to be global to be visible.
 */
// std::vector< std::pair< Double_t, Double_t > > coords;
std::vector< Double_t > coords;
std::vector< Double_t > values;
std::vector< Double_t > errors;

void FcnForMinimisation(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  // cout << "HI" << flush << endl;
  Int_t n = coords.size();
  // cout << "n is " << n << flush << endl;

  // Int_t n = 40;
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  for ( Int_t i = 0; i < n; ++i ) {
    // if ( i < 40 ) {
    if ( i < 13 ) {
      x[0] = coords[i];
      x[1] = 0;
      // switchFlag = 0;
    } else {
      x[0] = 0;
      x[1] = coords[i];
      // switchFlag = 1;
    }
    // cout << "HI2" << flush << endl;
    if ( values[i] != 0 ) {
      tmp = ( values[i] - SimultaneousFit( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  fval = chi2;
  cout << "chisquare = " << chi2 << endl;
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
  for (Int_t ix = 7; ix <= nBinsCosTheta-6; ++ix) {
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



  const Int_t npars = 4;
  TMinuit *gMinuit = new TMinuit(npars); //initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(FcnForMinimisation);
  Double_t arglist[300];
  //  Double_t * arglist = new Double_t [size];
  Int_t ierflg = 0;
  arglist[0] = 0.1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  gMinuit->mnparm(0, "LambdaTheta",    1., 0.1,     -10,     10, ierflg ); // LambdaTheta
  // gMinuit->mnparm(1, "NormalTheta", 4000, 100,  0, 100000, ierflg ); // normalisation THETA
  gMinuit->mnparm(1, "NormalTheta", 25500, 100,  1000, 27000, ierflg ); // normalisation THETA
  gMinuit->mnparm(2, "NormalisPhi",  4000, 100,   3500,  5000, ierflg ); // normalisation PHI
  gMinuit->mnparm(3, "LambdaPhi"  ,     0, 0.1,     -2,     2, ierflg ); // LambdaTheta



  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1;//1
  gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
  gMinuit->mnexcm("MIGRAD",  arglist ,2,ierflg);


  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);

  //_____________________________________________________________________
     Int_t iuext;
     TString chnam;   // The name of the parameter
     Double_t val;    // The current (external) value of the parameter
     Double_t errl;   // The current estimate of the parameter uncertainty
     Double_t xlolim; // The lower bound (or zero if no limits)
     Double_t xuplim; // The upper bound (or zero if no limits)
     Int_t iuint;     // The internal parameter number


     Double_t currentPar[npars] = {0};
     for (Int_t i=0; i< npars;i++) {
        gMinuit->mnpout(i, chnam, currentPar[i], errl, xlolim, xuplim, iuint);
     }



  // TVirtualFitter::SetDefaultFitter("Minuit");
  // TVirtualFitter * minuit = TVirtualFitter::Fitter(0,4);
  // // minuit->SetParameter(0,  1   ); // LambdaTheta
  // // minuit->SetParLimits(0, -2, 2);
  // // minuit->SetParameter(1, 4000 ); // normalisation THETA
  // // // minuit->SetParLimits(0, -2, 2);
  // // minuit->SetParameter(2, 4000 ); // normalisation PHI
  // // minuit->SetParLimits(3, -1, 1);
  // // minuit->SetParameter(3,  0   ); // LambdaPhi
  // minuit->SetParameter(0, "LambdaTheta",    1., 0.1, -2, 2 ); // LambdaTheta
  // // minuit->SetParameter(1, "NormalTheta", 4000, 100,  0, 100000 ); // normalisation THETA
  // minuit->SetParameter(1, "NormalTheta", 25500, 100,  25000, 27000 ); // normalisation THETA
  // minuit->SetParameter(2, "NormalisPhi", 4000, 100,  3500, 5000 ); // normalisation PHI
  // minuit->SetParameter(3, "LambdaPhi"  ,    0, 0.1, -2, 2 ); // LambdaTheta
  //
  // // cout << "AHM" << flush << endl;
  //
  //
  // minuit->SetFCN(FcnForMinimisation);
  //
  // double arglist[100];
  // arglist[0] = 0;
  // // set print level
  // minuit->ExecuteCommand("SET PRINT",arglist,2);
  // // cout << "UHM" << flush << endl;
  //
  //
  // // minimize
  // arglist[0] = 5000; // number of function calls
  // arglist[1] = 0.1; // tolerance
  // minuit->ExecuteCommand("MIGRAD",arglist,2);
  //
  // // cout << "EHM" << flush << endl;
  //
  // // get result
  // double minParams[10];
  // double parErrors[10];
  // for (int i = 0; i < 4; ++i) {
  //   minParams[i] = minuit->GetParameter(i);
  //   parErrors[i] = minuit->GetParError(i);
  // }
  // double chi2, edm, errdef;
  // int nvpar, nparx;
  // minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
  //
  // // func->SetParameters(minParams);
  // // func->SetParErrors(parErrors);
  // // func->SetChisquare(chi2);
  // int ndf = coords.size()-nvpar;
  // // func->SetNDF(ndf);
  //
  // std::cout << "Chi2 Fit = " << chi2 << " ndf = " << ndf << "  " << endl; //func->GetNDF() << std::endl;

  Double_t par[4], error[4];
  gMinuit->GetParameter(0,par[0],error[0]);
  gMinuit->GetParameter(1,par[1],error[1]);
  gMinuit->GetParameter(2,par[2],error[2]);
  gMinuit->GetParameter(3,par[3],error[3]);

  TF1* Model = new TF1("Model", "[1]*(1+[0]*x*x)/(3+[0])", -1 ,1 );
  new TCanvas;
  CorrectedCosTheta->Draw();
  Model->SetParameter( 0, par[0] );
  Model->SetParameter( 1, par[1] );
  Model->SetNpx(500);
  Model->Draw("same");

  TF1* Model2 = new TF1("Model2", "[1]*(1+2*[2]*cos(2*x)/(3+[0]))", -3.1 ,3.1 );
  new TCanvas;
  CorrectedPhi->Draw();
  Model2->SetParameter( 0, par[0] );
  Model2->SetParameter( 1, par[2] );
  Model2->SetParameter( 2, par[3] );
  Model2->SetNpx(500);
  Model2->Draw("same");






}
