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
#include "TF2.h"
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

Double_t GlobalChi = 0;


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
  fval      = chi2;
  GlobalChi = chi2;
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void PolarisationCsMinuit2D_0N0N(){

  TDatime d;
  TFile* file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DCS_0N0N/PolarisationCorrectedCs2D_0N0N.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TH2F* Distr2D = (TH2F*) file2D->Get("RawH");
  new TCanvas;
  Distr2D->Draw("ep");
  new TCanvas;


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
      // if( ((ix==0)&&(iy==18)) || ((ix==6)&&(iy==8))/*|| ((ix==7)&&(iy==9)) || ((ix==7)&&(iy==10))*/ ) continue;
      coordsX.push_back( ((TAxis*) Distr2D->GetXaxis())->GetBinCenter( ix ) );
      coordsY.push_back( Distr2D->GetYaxis()->GetBinCenter( iy ) );
      values.push_back(  Distr2D->GetBinContent( ix, iy )        );
      errors.push_back(  Distr2D->GetBinError(   ix, iy )        );
    }
  }





  /* - NB:
   * - COMPUTING THE
   * - CHI SQUARE
   */
  int ndf = coordsX.size()-4;

  TMinuit *gMinuit = new TMinuit(4);
  gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.1,    -2, 2        );
  gMinuit->DefineParameter(1, "LambdaPhi",           0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(2, "LambdaThetaPhi",      0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(3, "Normalisation",   50000, 100, 45000, 55000    );
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MINOS");
  Double_t LambdaTheta,    LambdaPhi,        LambdaThetaPhi,    Normalisation;
  Double_t LambdaThetaErr, LambdaPhiErr,     LambdaThetaPhiErr, NormalisationErr;
  gMinuit->GetParameter(0, LambdaTheta,      LambdaThetaErr     );
  gMinuit->GetParameter(1, LambdaPhi,        LambdaPhiErr       );
  gMinuit->GetParameter(2, LambdaThetaPhi,   LambdaThetaPhiErr  );
  gMinuit->GetParameter(3, Normalisation,    NormalisationErr   );
  printf("LambdaTheta     : %+.7f +- %.7f\n",LambdaTheta,     LambdaThetaErr     );
  printf("LambdaPhi       : %+.7f +- %.7f\n",LambdaPhi,       LambdaPhiErr       );
  printf("LambdaThetaPhi  : %+.7f +- %.7f\n",LambdaThetaPhi,  LambdaThetaPhiErr  );
  printf("Normalisation   : %+.7f +- %.7f\n",Normalisation,   NormalisationErr   );
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);

  gStyle->SetOptStat(0);

  cout << "OK1" << endl << flush;

  TF2* Model = new TF2("Model", helicity2D, -0.6 ,0.6, -3.14, 3.14, 4 );
  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  Distr2D->GetXaxis()->SetTitleOffset(1.15);
  // Distr2D->GetYaxis()->SetTitleOffset(1.25);
  Distr2D->GetYaxis()->SetTitleOffset(1.);
  Distr2D->GetXaxis()->SetTitleSize(0.045);
  Distr2D->GetYaxis()->SetTitleSize(0.045);
  Distr2D->GetXaxis()->SetLabelSize(0.045);
  Distr2D->GetYaxis()->SetLabelSize(0.045);
  Distr2D->GetXaxis()->SetTitleFont(42);
  Distr2D->GetYaxis()->SetTitleFont(42);
  Distr2D->GetXaxis()->SetLabelFont(42);
  Distr2D->GetYaxis()->SetLabelFont(42);
  Distr2D->GetXaxis()->SetNdivisions(408);
  cout << "OK2" << endl << flush;
  // Distr2D->GetYaxis()->SetRangeUser(0., Distr2D->GetMaximum()*2);
  // Distr2D->GetXaxis()->SetRangeUser(2, 6);
  Distr2D->SetTitle( ";cos(#theta); #phi" );
  Distr2D->Draw("colZ");
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  latex->DrawLatex(0.55,0.84,"UPC, Run 2, 0N0N, CS");
  latex->DrawLatex(0.55,0.78,"Minuit 2D Fit");
  latex->DrawLatex(0.55,0.70,Form("#lambda_{#theta} = %.3f #pm %.3f",     LambdaTheta,      LambdaThetaErr   ));
  latex->DrawLatex(0.55,0.62,Form("#lambda_{#phi} = %.3f #pm %.3f",       LambdaPhi,        LambdaPhiErr     ));
  latex->DrawLatex(0.55,0.54,Form("#lambda_{#theta#phi} = %.3f #pm %.3f", LambdaThetaPhi,   LambdaThetaPhiErr));
  latex->DrawLatex(0.55,0.44,Form("#tilde{#chi} = %3.3f/%3.3d = %2.2f",      GlobalChi,        ndf, (Double_t)GlobalChi/(Double_t)ndf ));
  Model->SetParameter( 0, LambdaTheta    );
  Model->SetParameter( 1, LambdaPhi      );
  Model->SetParameter( 2, LambdaThetaPhi );
  Model->SetParameter( 3, Normalisation  );
  cout << "OK3" << endl << flush;
  // Model->SetNpx(500);
  Model->Draw("same");
  gPad->SaveAs("pngResults/Cs2DMinuit_0N0N.png", "recreate");
  cout << "OK4" << endl << flush;

  new TCanvas;
  Model->Draw("surf1");


}
