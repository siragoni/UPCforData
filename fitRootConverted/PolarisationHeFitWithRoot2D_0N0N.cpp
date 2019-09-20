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
#include <vector>


#include "TDatime.h"

#include "TH2.h"


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
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void PolarisationHeFitWithRoot2D_0N0N(){

  TDatime d;
  TFile* data = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DHE_0N0N/PolarisationCorrectedHe2D_0N0N.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* data = new TFile("pngResults/2019-09-18/2DHE/PolarisationCorrectedHe2D.root");
  TH2F *fCorrectedShape = (TH2F*) data->Get("RawH");
  TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-0.6, 0.6, -3.1, 3.1, 4);
  // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-0.9, 0.9, -2.7, 2.7, 4);
  ParabolicFit->SetNpx(60);
  ParabolicFit->SetParameter(0, 1);
  ParabolicFit->SetParameter(1, 0);
  ParabolicFit->SetParameter(2, 0);
  // ParabolicFit->FixParameter(0, 1);
  // ParabolicFit->FixParameter(1, 0);
  // ParabolicFit->FixParameter(2, 0);
  ParabolicFit->SetParameter(3, 52000);
  // ParabolicFit->SetParameter(3, 4200);
  ParabolicFit->SetParLimits(0, -2, +2);
  ParabolicFit->SetParLimits(1, -2, +2);
  ParabolicFit->SetParLimits(2, -2, +2);
  // ParabolicFit->SetParLimits(3, 4000, 4500);
  // ParabolicFit->SetParLimits(3, 700, 800);
  // ParabolicFit->SetParLimits(0,  0.75, 1.25);
  // ParabolicFit->SetParLimits(1, -0.25, 0.25);
  // ParabolicFit->SetParLimits(2, -0.25, 0.25);
  // fCorrectedShape->Fit( "ParabolicFit" , "RL" );
  fCorrectedShape->Fit( "ParabolicFit" , "R" );
  // fCorrectedShape->Fit( ParabolicFit,"","", -0.5, 0.5, -2, 2 );

  fCorrectedShape->SetLineColor(kBlue);
  fCorrectedShape->SetLineStyle(kSolid);
  fCorrectedShape->SetLineWidth(3);
  fCorrectedShape->SetMarkerStyle(kFullCircle);
  fCorrectedShape->SetMarkerSize(1);
  fCorrectedShape->GetXaxis()->SetTitle("cos(#theta)");
  fCorrectedShape->GetYaxis()->SetTitle("#phi");
  // fCorrectedShape->GetYaxis()->SetTitle( Form( "Counts / (%.3f a.u.)",
  //                                              fCorrectedShape->GetXaxis()->GetBinWidth(1)
  //                                              )
  //                                             );
  fCorrectedShape->SetTitle("");
  new TCanvas;
  TCanvas* ZNAEnergy = new TCanvas( "CosTheta", "CosTheta", 900, 800 );
  gPad->SetMargin(0.13,0.125,0.12,0.13);
  /* - Beautifying is starting now.
     -
   */
  fCorrectedShape->GetXaxis()->SetTitleOffset(1.25);
  fCorrectedShape->GetYaxis()->SetTitleOffset(0.5);
  // fCorrectedShape->GetYaxis()->SetTitleOffset(1.45);
  fCorrectedShape->GetXaxis()->SetTitleSize(0.045);
  fCorrectedShape->GetYaxis()->SetTitleSize(0.045);
  fCorrectedShape->GetXaxis()->SetLabelSize(0.045);
  fCorrectedShape->GetYaxis()->SetLabelSize(0.045);
  fCorrectedShape->GetXaxis()->SetTitleFont(42);
  fCorrectedShape->GetYaxis()->SetTitleFont(42);
  fCorrectedShape->GetXaxis()->SetLabelFont(42);
  fCorrectedShape->GetYaxis()->SetLabelFont(42);
  fCorrectedShape->GetXaxis()->SetNdivisions(408);
  fCorrectedShape->GetYaxis()->SetRangeUser(-3.14, 3.14);
  fCorrectedShape->GetXaxis()->SetRangeUser(-1, 1);
  // fCorrectedShape->GetXaxis()->SetRangeUser(-0.6, 0.6);
  // fCorrectedShape->GetYaxis()->SetRangeUser(5, fCorrectedShape->GetMaximum()*10.);
  // gPad ->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fCorrectedShape->Draw("colZ  SAME");
  // fCorrectedShape->Draw("surf3  ");
  ParabolicFit   ->Draw("cont1 SAME");
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, Run 2, 0N0N, HE");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));
  latex->DrawLatex(0.55,0.72,Form("#lambda_{#theta} = %.4f #pm %.4f",
                                  ParabolicFit->GetParameter(0),
                                  ParabolicFit->GetParError(0)
                                  )
                                );
  latex->DrawLatex(0.55,0.66,Form("#lambda_{#phi} = %.4f #pm %.4f",
                                  ParabolicFit->GetParameter(1),
                                  ParabolicFit->GetParError(1)
                                  )
                                );
  latex->DrawLatex(0.55,0.60,Form("#lambda_{#theta#phi} = %.4f #pm %.4f",
                                  ParabolicFit->GetParameter(2),
                                  ParabolicFit->GetParError(2)
                                  )
                                );
  latex->DrawLatex(0.5,0.18,Form( "#tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                   ParabolicFit->GetChisquare(),
                                   ParabolicFit->GetNDF(),
                                   ParabolicFit->GetChisquare()/ParabolicFit->GetNDF()
                                  )
                                 );

  gPad->SaveAs("pngResults/PolFitWithRootHe2D_0N0N.png",         "RECREATE");

}
