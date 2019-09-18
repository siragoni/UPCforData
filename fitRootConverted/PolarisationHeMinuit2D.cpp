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




//_____________________________________________________________________________
/* - Codign in the fit functions.
   - The fit is modelled as the sum of 3 different functions:
   - 1) CosTheta only
   - 2) Phi      only
   - 3) Mix of the two
   -
 */
Double_t CosTheta(Double_t *x, Double_t *par) {
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t returnValue     = 1 + par[0] * CosSquaredTheta;
  return   returnValue;
}
//______________________________________________
Double_t Phi(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi     = TMath::Cos( 2 * x[1] );
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t SinSquaredTheta = 1 - CosSquaredTheta;

  Double_t returnValue = par[0] * SinSquaredTheta * CosOfTwoPhi;
  return   returnValue;
}
//______________________________________________
Double_t Mix(Double_t *x, Double_t *par) {
  Double_t CosPhi               = TMath::Cos( x[1] );
  Double_t CosSquaredTheta      = x[0]*x[0];
  Double_t SinSquaredTheta      = 1 - CosSquaredTheta;
  Double_t SinSquaredOfTwoTheta = 4 * CosSquaredTheta * SinSquaredTheta;
  Double_t SinOfTwoTheta        = TMath::Sqrt(SinSquaredOfTwoTheta);
  Double_t returnValue          = par[0] * SinOfTwoTheta * CosPhi;
  return   returnValue;
}
//______________________________________________
Double_t helicity2D(Double_t *x, Double_t *par) {
   Double_t *lambdaTheta    = &par[0];
   Double_t *lambdaPhi      = &par[1];
   Double_t *lambdaThetaPhi = &par[2];
   Double_t sumOfTheSubFits = CosTheta( x, lambdaTheta) + Phi( x, lambdaPhi ) + Mix( x, lambdaThetaPhi );
   Double_t FinalResult     = par[3] * sumOfTheSubFits / ( 3 + par[0] );
   return   FinalResult;
}
//_____________________________________________________________________________
/* - Data needed for the global fit.
 * - They have to be global to be visible.
 */
// std::vector< std::pair< Double_t, Double_t > > coords;
std::vector< Double_t > coordsX;
std::vector< Double_t > coordsY;
std::vector< Double_t > values;
std::vector< Double_t > errors;

void FcnForMinimisation(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  // cout << "HI" << flush << endl;
  Int_t n = coordsX.size();
  // cout << "HI2" << flush << endl;

  // Int_t n = 40;
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  for ( Int_t i = 0; i < n; ++i ) {
    x[0] = coordsX[i];
    x[1] = coordsY[i];
    if ( values[i] != 0 ) {
      tmp = ( values[i] - helicity2D( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  fval = chi2;
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void fitHelicitySimultaneously2D(){

  // TFile* file2D = new TFile("pngResults/TH2corrPerfect.root");
  TFile* file2D = new TFile("pngResults/TH2corrPerfect2.root");
  TH2F* Distr2D = (TH2F*) file2D->Get("RawH");


  // Double_t CosThetaLowLimit   = -1;
  // Double_t CosThetaUpperLimit = +1;
  // Double_t PhiLowLimit        = -3.14;
  // Double_t PhiUpperLimit      = +3.14;

  Int_t nBinsCosTheta = Distr2D->GetNbinsX();
  Int_t nBinsPhi      = Distr2D->GetNbinsY();


  /// reset data structure
  coordsX = std::vector<Double_t>();
  coordsY = std::vector<Double_t>();
  values  = std::vector<Double_t>();
  errors  = std::vector<Double_t>();
  /// fill data structure
  for (Int_t ix = 1; ix <= nBinsCosTheta; ++ix) {
    for (Int_t iy = 1; iy <= nBinsPhi; ++iy) {
      if( ((ix==1)&&(iy==18)) || ((ix==7)&&(iy==9)) || ((ix==7)&&(iy==10)) ) continue;
      coordsX.push_back( ((TAxis*) Distr2D->GetXaxis())->GetBinCenter( ix ) );
      coordsY.push_back( Distr2D->GetYaxis()->GetBinCenter( iy ) );
      values.push_back(  Distr2D->GetBinContent( ix, iy )        );
      errors.push_back(  Distr2D->GetBinError(   ix, iy )        );
    }
  }



  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * minuit = TVirtualFitter::Fitter(0,4);
  minuit->SetParameter(0, "LambdaTheta",         1, 0.1, -2, 2 ); // LambdaTheta
  minuit->SetParameter(1, "LambdaPhi"  ,         0, 0.1, -2, 2 ); // LambdaPhi
  minuit->SetParameter(2, "LambdaThetaPhi"  ,    0, 0.1, -2, 2 ); // LambdaThetaPhi
  // minuit->SetParameter(3, "Norm",            52000, 100,  0, 100000 ); // normalisation
  minuit->SetParameter(3, "Norm",            48200, 100,  0, 100000 ); // normalisation

  // cout << "AHM" << flush << endl;


  minuit->SetFCN(FcnForMinimisation);

  double arglist[100];
  arglist[0] = 0;
  // set print level
  minuit->ExecuteCommand("SET PRINT",arglist,2);
  // cout << "UHM" << flush << endl;


  // minimize
  arglist[0] = 5000; // number of function calls
  arglist[1] = 1; // tolerance
  minuit->ExecuteCommand("MIGRAD",arglist, 2);
  minuit->ExecuteCommand("HESSE", arglist, 2);
	minuit->ExecuteCommand("MINOS", arglist, 2);

  // cout << "EHM" << flush << endl;

  // get result
  double minParams[10];
  double parErrors[10];
  for (int i = 0; i < 4; ++i) {
    minParams[i] = minuit->GetParameter(i);
    parErrors[i] = minuit->GetParError(i);
  }
  double chi2, edm, errdef;
  int nvpar, nparx;
  minuit->GetStats(chi2,edm,errdef,nvpar,nparx);

  // func->SetParameters(minParams);
  // func->SetParErrors(parErrors);
  // func->SetChisquare(chi2);
  int ndf = coordsX.size()-nvpar;
  // func->SetNDF(ndf);

  std::cout << "Chi2 Fit = " << chi2 << " ndf = " << ndf << "  " << endl; //func->GetNDF() << std::endl;

}
