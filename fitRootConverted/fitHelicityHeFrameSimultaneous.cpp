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
/* - Coding in the fit functions.
   - The fit is modelled as the sum of 3 different functions:
   - 1) CosTheta only
   - 2) Phi      only
   - 3) Mix of the two
   -
 */
// Double_t CosTheta(Double_t *x, Double_t *par) {
//   Double_t CosSquaredTheta = x[0] * x[0];
//   Double_t returnValue     = 1 + par[0] * CosSquaredTheta;
//   returnValue              = par[1] * returnValue / ( 3 + par[0] );
//   return   returnValue;
// }
// //______________________________________________
// Double_t Phi(Double_t *x, Double_t *par) {
//   Double_t CosOfTwoPhi = TMath::Cos( 2 * x[0] );
//   Double_t returnValue = par[2] * ( 1 + 2 * par[1] * CosOfTwoPhi / ( 3 + par[0] ) );
//   return   returnValue;
// }
// //______________________________________________
// Double_t MixWithPositiveCosTheta(Double_t *x, Double_t *par) {
//   Double_t CosTwoTildePhi = TMath::Cos( 2 * x[0] - 0.50 * TMath::Pi() );
//   Double_t returnValue = par[2] * ( 1 + TMath::Sqrt(2) * par[1] * CosTwoTildePhi / ( 3 + par[0] ) );
//   return   returnValue;
// }
// //______________________________________________
// Double_t MixWithNegativeCosTheta(Double_t *x, Double_t *par) {
//   Double_t CosTwoTildePhi = TMath::Cos( 2 * x[0] - 1.50 * TMath::Pi() );
//   Double_t returnValue = par[2] * ( 1 + TMath::Sqrt(2) * par[1] * CosTwoTildePhi / ( 3 + par[0] ) );
//   return   returnValue;
// }
//_____________________________________________________________________________
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
double SimultaneousFitComplete(double *x, double *par) {
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
double SimultaneousFit(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = CosTheta( x, par ) + Phi( x, par );
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
  // cout << "HI2" << flush << endl;

  // Int_t n = 40;
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  for ( Int_t i = 0; i < n; ++i ) {
    // if ( i < 40 ) {
    if ( i < 22 ) {
      x[0] = coords[i];
      x[1] = 0;
    } else {
      x[0] = 0;
      x[1] = coords[i];
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
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void fitHelicitySimultaneously1d(){

  TFile* file1D = new TFile("pngResults/TH1corr.root");
  // TFile* file2D = new TFile("pngResults/TH2corr.root");

  TFile* fileNew = new TFile("pngResults/TH1corrMyBinning.root");
  //
  TH1F* CorrectedCosTheta = (TH1F*) fileNew->Get("RawCosThetaH");
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
  for (Int_t ix = 1; ix <= nBinsCosTheta-4; ++ix) {
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

  // for(int i=0; i < 74; i++){
  //   cout << i << "  " << coords[i] << "  " << values[i] << endl;
  // }



  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * minuit = TVirtualFitter::Fitter(0,4);
  // minuit->SetParameter(0,  1   ); // LambdaTheta
  // minuit->SetParLimits(0, -2, 2);
  // minuit->SetParameter(1, 4000 ); // normalisation THETA
  // // minuit->SetParLimits(0, -2, 2);
  // minuit->SetParameter(2, 4000 ); // normalisation PHI
  // minuit->SetParLimits(3, -1, 1);
  // minuit->SetParameter(3,  0   ); // LambdaPhi
  minuit->SetParameter(0, "LambdaTheta",    1, 0.1, -2, 2 ); // LambdaTheta
  // minuit->SetParameter(1, "NormalTheta", 4000, 100,  0, 100000 ); // normalisation THETA
  minuit->SetParameter(1, "NormalTheta", 280000, 1000,  250000, 3000000 ); // normalisation THETA
  minuit->SetParameter(2, "NormalisPhi", 4000, 100,  3500, 5000 ); // normalisation PHI
  minuit->SetParameter(3, "LambdaPhi"  ,    0, 0.1, -2, 2 ); // LambdaTheta

  // cout << "AHM" << flush << endl;


  minuit->SetFCN(FcnForMinimisation);

  double arglist[100];
  arglist[0] = 0;
  // set print level
  minuit->ExecuteCommand("SET PRINT",arglist,2);
  // cout << "UHM" << flush << endl;


  // minimize
  arglist[0] = 5000; // number of function calls
  arglist[1] = 0.1; // tolerance
  minuit->ExecuteCommand("MIGRAD",arglist,2);

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
  int ndf = coords.size()-nvpar;
  // func->SetNDF(ndf);

  std::cout << "Chi2 Fit = " << chi2 << " ndf = " << ndf << "  " << endl; //func->GetNDF() << std::endl;



  TF1* Model = new TF1("Model", "[1]*(1+[0]*x*x)/(3+[0])", -1 ,1 );
  new TCanvas;
  CorrectedCosTheta->Draw();
  Model->SetParameter( 0, minuit->GetParameter(0) );
  Model->SetParameter( 1, minuit->GetParameter(1) );
  Model->SetNpx(500);
  Model->Draw("same");

  TF1* Model2 = new TF1("Model2", "[1]*(1+2*[2]*cos(2*x)/(3+[0]))", -3.1 ,3.1 );
  new TCanvas;
  CorrectedPhi->Draw();
  Model2->SetParameter( 0, minuit->GetParameter(0) );
  Model2->SetParameter( 1, minuit->GetParameter(2) );
  Model2->SetParameter( 2, minuit->GetParameter(3) );
  Model2->SetNpx(500);
  Model2->Draw("same");

  // /* - There are three cases for the selectionFlag:
  //    - 1) = 0 ; this implies the traditional pt-integrated plot;
  //    - 2) = 1 ; this is instead the coherent component;
  //    - 3) = 2 ; this is the incoherent component;
  //    - 4) = 3 ; ******************* ;
  //    -
  //  */
  // TFile* dataList = new TFile(AnalysisName);
  // TDirectory* dirData = dataList->GetDirectory("MyTask");
  // TFile* mcList = new TFile(MonteCarloName);
  // TDirectory* dirMC = mcList->GetDirectory("MyTask");
  // /* - At this level you could check if everything was okay.
  //  * - We do a dir->ls() to find out! We get:
  //  *   dir->ls();
  //  *   TDirectoryFile*		MyTask	MyTask
  //  *   KEY: TList	MyOutputContainer;1	Doubly linked list
  //  */
  // TList* listingsData;
  // dirData->GetObject("MyOutputContainer", listingsData);
  // TList* listingsMC;
  // dirMC  ->GetObject("MyOutputContainer", listingsMC);
  // /* - We now do the same as before to ascertain if the TList was there and
  //  * - to try to retrieve the plots. Result:
  //  *   listings->ls()
  //  *     OBJ: TList	  MyOutputContainer	          Doubly linked list          : 0
  //  *     OBJ: TH1F	  fNumberMuonsH	              fNumberMuonsH               : 0 at: 0x5a145f0
  //  *     OBJ: TH1F	  fCounterH	                  fCounterH                   : 0 at: 0x5a3b570
  //  *     OBJ: TH1F	  fEtaMuonH	                  fEtaMuonH                   : 0 at: 0x5a3ba80
  //  *     OBJ: TH1F	  fRAbsMuonH	                fRAbsMuonH                  : 0 at: 0x5a3c0c0
  //  *     OBJ: TH1F	  fInvariantMassDistributionH	fInvariantMassDistributionH : 0 at: 0x5a3c720
  //  */
  // new TCanvas;
  // TH2F *fAngularDistribOfPositiveMuonRestFrameJPsiH = (TH2F*)listingsMC->FindObject("fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH");
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin2D(4,8);
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin2D(4,4);
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin2D(2,2);
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(16);
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("colZ");
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("cont1 Z");
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("SURF1 Z");
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->SetFillColor(kBlue);
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("LEGO1 0");
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->SaveAs("rawDistr.png", "RECREATE");
  //
  //
  //
  //
  // new TCanvas;
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  // // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(5);
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->SetMarkerStyle(21);
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->SetLineColor(kBlue+1);
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->SetLineWidth(2);
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("ep");
  // TH2F *fCorrectedShape = fAngularDistribOfPositiveMuonRestFrameJPsiH;
  //
  // /* - Do the parabolic fit with Root only.
  //    - Retrieve the parameters later for Roofit plotting.
  //    -
  //    - dN/dPhi = 1/2pi [ 1 - 4*rho_{1,-1}*Cos2Phi ]
  // */
  // new TCanvas;
  // // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-1, 1, -4, 4, 3);
  // // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-0.9, 0.9, -3, 3, 3);
  // /* - The global normalisation SHOULD be a parameter too!!!
  //    -
  //  */
  // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-1, 1, -2.9, 2.9, 4);
  // // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-0.9, 0.9, -2.7, 2.7, 4);
  // ParabolicFit->SetNpx(1000);
  // ParabolicFit->SetParameter(0, 1);
  // ParabolicFit->SetParameter(1, 0);
  // ParabolicFit->SetParameter(2, 0);
  // // ParabolicFit->FixParameter(0, 1);
  // // ParabolicFit->FixParameter(1, 0);
  // // ParabolicFit->FixParameter(2, 0);
  // // ParabolicFit->SetParameter(3, 12000);
  // ParabolicFit->SetParameter(3, 760);
  // ParabolicFit->SetParLimits(0, -2, +2);
  // ParabolicFit->SetParLimits(1, -2, +2);
  // ParabolicFit->SetParLimits(2, -2, +2);
  // ParabolicFit->SetParLimits(3, 700, 800);
  // // ParabolicFit->SetParLimits(0,  0.75, 1.25);
  // // ParabolicFit->SetParLimits(1, -0.25, 0.25);
  // // ParabolicFit->SetParLimits(2, -0.25, 0.25);
  // fCorrectedShape->Fit( "ParabolicFit" , "R" );
  // // fCorrectedShape->Fit( ParabolicFit,"","", -0.5, 0.5, -2, 2 );
  //
  // fCorrectedShape->SetLineColor(kBlue);
  // fCorrectedShape->SetLineStyle(kSolid);
  // fCorrectedShape->SetLineWidth(3);
  // fCorrectedShape->SetMarkerStyle(kFullCircle);
  // fCorrectedShape->SetMarkerSize(1);
  // fCorrectedShape->GetXaxis()->SetTitle("cos(#theta)");
  // fCorrectedShape->GetYaxis()->SetTitle("#phi");
  // // fCorrectedShape->GetYaxis()->SetTitle( Form( "Counts / (%.3f a.u.)",
  // //                                              fCorrectedShape->GetXaxis()->GetBinWidth(1)
  // //                                              )
  // //                                             );
  // fCorrectedShape->SetTitle("");
  // new TCanvas;
  // TCanvas* ZNAEnergy = new TCanvas( "CosTheta", "CosTheta", 900, 800 );
  // gPad->SetMargin(0.13,0.01,0.12,0.01);
  // /* - Beautifying is starting now.
  //    -
  //  */
  // fCorrectedShape->GetXaxis()->SetTitleOffset(1.25);
  // fCorrectedShape->GetYaxis()->SetTitleOffset(1.25);
  // // fCorrectedShape->GetYaxis()->SetTitleOffset(1.45);
  // fCorrectedShape->GetXaxis()->SetTitleSize(0.045);
  // fCorrectedShape->GetYaxis()->SetTitleSize(0.045);
  // fCorrectedShape->GetXaxis()->SetLabelSize(0.045);
  // fCorrectedShape->GetYaxis()->SetLabelSize(0.045);
  // fCorrectedShape->GetXaxis()->SetTitleFont(42);
  // fCorrectedShape->GetYaxis()->SetTitleFont(42);
  // fCorrectedShape->GetXaxis()->SetLabelFont(42);
  // fCorrectedShape->GetYaxis()->SetLabelFont(42);
  // fCorrectedShape->GetXaxis()->SetNdivisions(408);
  // fCorrectedShape->GetYaxis()->SetRangeUser(-3.2, 3.2);
  // fCorrectedShape->GetXaxis()->SetRangeUser(-1, 1);
  // // fCorrectedShape->GetXaxis()->SetRangeUser(-0.6, 0.6);
  // // fCorrectedShape->GetYaxis()->SetRangeUser(5, fCorrectedShape->GetMaximum()*10.);
  // // gPad ->SetLogy();
  // gStyle->SetOptFit(0);
  // gStyle->SetOptStat(0);
  // fCorrectedShape->Draw("colZ  SAME");
  // // fCorrectedShape->Draw("surf3  ");
  // ParabolicFit   ->Draw("cont1 SAME");
  // TLatex* latex = new TLatex();
  // latex->SetTextSize(0.05);
  // latex->SetTextFont(42);
  // latex->SetTextAlign(11);
  // latex->SetNDC();
  // latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  // latex->SetTextSize(0.045);
  // // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  // latex->DrawLatex(0.55,0.84,"UPC, LHC18qr");
  // // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  // latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));
  // latex->DrawLatex(0.55,0.72,Form("#lambda_{#theta} = %.4f #pm %.4f",
  //                                 ParabolicFit->GetParameter(0),
  //                                 ParabolicFit->GetParError(0)
  //                                 )
  //                               );
  // latex->DrawLatex(0.55,0.66,Form("#lambda_{#phi} = %.4f #pm %.4f",
  //                                 ParabolicFit->GetParameter(1),
  //                                 ParabolicFit->GetParError(1)
  //                                 )
  //                               );
  // latex->DrawLatex(0.55,0.60,Form("#lambda_{#theta#phi} = %.4f #pm %.4f",
  //                                 ParabolicFit->GetParameter(2),
  //                                 ParabolicFit->GetParError(2)
  //                                 )
  //                               );
  // latex->DrawLatex(0.55,0.18,Form( "      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
  //                                  ParabolicFit->GetChisquare(),
  //                                  ParabolicFit->GetNDF(),
  //                                  ParabolicFit->GetChisquare()/ParabolicFit->GetNDF()
  //                                 )
  //                                );
  //
  //
  // gPad->SaveAs("pngResults/fit2DonlyGeneratedLevel.png",         "RECREATE");

}
