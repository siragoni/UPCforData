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
void PolarisationHeMinuit2D( Int_t SignalRangeSelectionMode = 0, Int_t FitRangeMode = 0 ){

  TDatime d;
  TFile* file2D = 0x0;
  // TFile* file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  if        ( SignalRangeSelectionMode == 0 ) {
    file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 1 ) {
    file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D_1.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 2 ) {
    file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 3 ) {
    file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 4 ) {
    file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 5 ) {
    file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else {
    file2D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE/PolarisationCorrectedHe2D.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  }

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
      // if( ((ix==0)&&(iy==18)) || ((ix==6)&&(iy==8))/*|| ((ix==7)&&(iy==9)) || ((ix==7)&&(iy==10))*/ ) continue;

      if        ( FitRangeMode == 1 ) {
        if (ix == 1)   continue;
      } else if ( FitRangeMode == 2 ) {
        if (ix == 7)   continue;
      } else if ( FitRangeMode == 3 ) {
        if ( (ix == 1) || (ix == 7) )  continue;
      } else {
      }
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













  Double_t CosThetaCenters[] = { -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5 };
  Double_t PhiCenters[] = { -3.14*19.5*0.05, -3.14*18.5*0.05, -3.14*17.5*0.05, -3.14*15*0.05,
                            -3.14*11*0.05,   -3.14*7.5*0.05,  -3.14*5*0.05,    -3.14*3*0.05,
                            -3.14*1.5*0.05,  -3.14*0.5*0.05,  +3.14*0.5*0.05,  +3.14*1.5*0.05,
                            +3.14*3*0.05,    +3.14*5*0.05,    +3.14*7.5*0.05,  +3.14*11*0.05,
                            +3.14*15*0.05,   +3.14*17.5*0.05, +3.14*18.5*0.05, +3.14*19.5*0.05 };


  Double_t MyVariableCosThetaBinning[] = { -0.65, -0.35, -0.15, -0.05,
                                            0.05,  0.15,  0.35,  0.65 };
  Double_t MyVariablePhiBinning[] = { -3.14*1,       -3.14*19*0.05, -3.14*18*0.05, -3.14*17*0.05,
                                      -3.14*13*0.05, -3.14*9*0.05,  -3.14*6*0.05,  -3.14*4*0.05,
                                      -3.14*2*0.05,  -3.14*1*0.05,   0,            +3.14*1*0.05,
                                      +3.14*2*0.05,  +3.14*4*0.05,  +3.14*6*0.05,  +3.14*9*0.05,
                                      +3.14*13*0.05, +3.14*17*0.05, +3.14*18*0.05, +3.14*19*0.05,
                                      +3.14*1 };

  TH2F* helicity2DafterSignalExtractionErrors =
            new TH2F( "helicity2DafterSignalExtractionErrors",
                      "helicity2DafterSignalExtractionErrors",
                      20, MyVariablePhiBinning, 7, MyVariableCosThetaBinning
                      // 10, -1, 1, 10, -4, 4
                      );
  for (size_t iCosThetaBins = 0; iCosThetaBins < 7; iCosThetaBins++) {
    for (size_t iPhiBins = 0; iPhiBins < 20; iPhiBins++) {

      Int_t binx = Distr2D->GetXaxis()->FindBin(CosThetaCenters[iCosThetaBins]);
      Int_t biny = Distr2D->GetYaxis()->FindBin(PhiCenters[iPhiBins]);

      helicity2DafterSignalExtractionErrors->Fill( PhiCenters[iPhiBins],
                                                   CosThetaCenters[iCosThetaBins],
                                                   Distr2D->GetBinContent(binx, biny)
                                                   );
      helicity2DafterSignalExtractionErrors->SetBinError( iPhiBins      + 1 ,
                                                          iCosThetaBins + 1 ,
                                                          Distr2D->GetBinError(binx, biny)
                                                          );

    }
  }










  TF2* Model = new TF2("Model", helicity2D, -0.6 ,0.6, -3.14, 3.14, 4 );
  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  // Distr2D->GetXaxis()->SetTitleOffset(1.15);
  // // Distr2D->GetYaxis()->SetTitleOffset(1.25);
  // Distr2D->GetYaxis()->SetTitleOffset(1.);
  // Distr2D->GetXaxis()->SetTitleSize(0.045);
  // Distr2D->GetYaxis()->SetTitleSize(0.045);
  // Distr2D->GetXaxis()->SetLabelSize(0.045);
  // Distr2D->GetYaxis()->SetLabelSize(0.045);
  // Distr2D->GetXaxis()->SetTitleFont(42);
  // Distr2D->GetYaxis()->SetTitleFont(42);
  // Distr2D->GetXaxis()->SetLabelFont(42);
  // Distr2D->GetYaxis()->SetLabelFont(42);
  // Distr2D->GetXaxis()->SetNdivisions(408);

  helicity2DafterSignalExtractionErrors->GetXaxis()->SetTitleOffset(1.15);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.25);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.);
  helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(0.9);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetLabelSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetLabelSize(0.045);
  helicity2DafterSignalExtractionErrors->GetXaxis()->SetTitleSize(0.055);
  helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleSize(0.055);
  helicity2DafterSignalExtractionErrors->GetXaxis()->SetLabelSize(0.05);
  helicity2DafterSignalExtractionErrors->GetYaxis()->SetLabelSize(0.05);
  helicity2DafterSignalExtractionErrors->GetXaxis()->SetTitleFont(42);
  helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleFont(42);
  helicity2DafterSignalExtractionErrors->GetXaxis()->SetLabelFont(42);
  helicity2DafterSignalExtractionErrors->GetYaxis()->SetLabelFont(42);
  helicity2DafterSignalExtractionErrors->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
  // // Distr2D->GetYaxis()->SetRangeUser(0., Distr2D->GetMaximum()*2);
  // // Distr2D->GetXaxis()->SetRangeUser(2, 6);
  // Distr2D->SetTitle( ";cos(#theta); #phi" );
  // Distr2D->Draw("colZ");
  helicity2DafterSignalExtractionErrors->SetTitle( "; #it{#varphi} ;cos(#it{#theta})" );
  // Int_t Number = 5;
  // Double_t Red[5] = { 0.00, 0.09, 0.18, 0.09, 0.00 };
  // Double_t Green[5] = { 0.01, 0.02, 0.39, 0.68, 0.97 };
  // Double_t Blue[5] = { 0.17, 0.39, 0.62, 0.79, 0.97 };
  // Double_t Stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  // Double_t Red[13]   = { 0.00, 0.030, 0.060, 0.09, 0.120, 0.150, 0.18, 0.15, 0.12, 0.09, 0.06, 0.03, 0.00 };
  // Double_t Green[13] = { 0.01, 0.013, 0.017, 0.02, 0.026, 0.032, 0.39, 0.49, 0.58, 0.68, 0.78, 0.88, 0.97 };
  // Double_t Blue[13]  = { 0.17, 0.24,  0.31,  0.39, 0.46,  0.53,  0.62, 0.68, 0.74, 0.79, 0.86, 0.92, 0.97 };
  // Double_t Stops[13] = { 10000, 12000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 25000 };
  //
  //
  // Int_t nb=1000;
  // TColor::CreateGradientColorTable(13,Stops,Red,Green,Blue,nb);
  // helicity2DafterSignalExtractionErrors->SetContour(nb);
  // helicity2DafterSignalExtractionErrors->GetZaxis()->SetMoreLogLabels();
  // helicity2DafterSignalExtractionErrors->Draw("colZ");
  // helicity2DafterSignalExtractionErrors->Draw("BOX");









    const Int_t NRGBs = 6;
    const Int_t NCont = 999;

    Double_t stops[NRGBs] = { 0.00, 0.1, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.0, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };


    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetOptStat(0);

    //here the actually interesting code starts
    // const Double_t min = 8000;
    const Double_t min = 6000;
    const Double_t max = 26000;

    const Int_t nLevels = 999;
    Double_t levels[nLevels];


    for(int i = 1; i < nLevels; i++) {
      levels[i] = min + (max - min) / (nLevels - 1) * (i);
    }
    levels[0] = 0.01;
  //  levels[0] = -1; //Interesting, but this also works as I want!

    // TCanvas * c = new TCanvas();
    // TH2D *h  = new TH2D("h", "", 10, 0, 10, 10, 0, 10);
    helicity2DafterSignalExtractionErrors->SetContour((sizeof(levels)/sizeof(Double_t)), levels);

    // h->SetBinContent(5, 7, 1.20);
    // h->SetBinContent(5, 6, 1.05);
    // h->SetBinContent(5, 5, 1.00);
    // h->SetBinContent(5, 4, 0.95);
    // h->SetBinContent(5, 3, 0.80);

    helicity2DafterSignalExtractionErrors->DrawClone("col");// draw "axes", "contents", "statistics box"

    helicity2DafterSignalExtractionErrors->GetZaxis()->SetRangeUser(min, max); // ... set the range ...

    helicity2DafterSignalExtractionErrors->Draw("colz same"); // draw the "color palette"




  TLatex* latex = new TLatex();
  latex->SetTextSize(0.055);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  // latex->DrawLatex(0.17,0.94,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.7,0.94,"Helicity");
  // latex->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");

  gPad->SaveAs("pngResults/Inverted2Dmaps.pdf",  "RECREATE");




}
