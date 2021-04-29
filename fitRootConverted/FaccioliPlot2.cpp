#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"


//_____________________________________________________________________________
/* -
 * - Original macro:
 * - https://root.cern/doc/v610/graphShade_8C.html
 */
void FaccioliPlot()
{
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  // gPad->SetLeftMargin(0.2);
  // // gPad->SetLeftMargin(0.12);
  // // gPad->SetRightMargin(0.1);
  // gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.1);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();


  TF1 *LambdaThetaUpper = new TF1("LambdaThetaUpper","1.",-1.,1.);
  LambdaThetaUpper->SetTitle("");
  LambdaThetaUpper->GetXaxis()->SetTitleOffset(1.15);
  // LambdaThetaUpper->GetYaxis()->SetTitleOffset(1.25);
  LambdaThetaUpper->GetYaxis()->SetTitleOffset(1.25);
  LambdaThetaUpper->GetXaxis()->SetTitleSize(0.045);
  LambdaThetaUpper->GetYaxis()->SetTitleSize(0.045);
  LambdaThetaUpper->GetXaxis()->SetLabelSize(0.045);
  LambdaThetaUpper->GetYaxis()->SetLabelSize(0.045);
  LambdaThetaUpper->GetXaxis()->SetTitleFont(42);
  LambdaThetaUpper->GetYaxis()->SetTitleFont(42);
  LambdaThetaUpper->GetXaxis()->SetLabelFont(42);
  LambdaThetaUpper->GetYaxis()->SetLabelFont(42);

  LambdaThetaUpper->GetYaxis()->SetTitle("#lambda_{#theta}");
  LambdaThetaUpper->GetXaxis()->SetTitle("#lambda_{#varphi}");
  LambdaThetaUpper->GetYaxis()->SetRangeUser(-1.,1.);
  LambdaThetaUpper->GetXaxis()->SetRangeUser(-1.,1.);
  // LambdaThetaUpper->GetYaxis()->SetRangeUser(-1.4,1.4);
  // LambdaThetaUpper->GetXaxis()->SetRangeUser(-1.4,1.4);
  LambdaThetaUpper->SetLineWidth(5);
  LambdaThetaUpper->SetLineColor(2);
  // LambdaThetaUpper->SetFillColor(kRed-3);
  // LambdaThetaUpper->SetFillStyle(1001);
  LambdaThetaUpper->Draw();

  TF1 *LambdaThetaLowerPlus = new TF1("LambdaThetaLowerPlus","-1.+2*x",0.,1.);
  LambdaThetaLowerPlus->SetLineWidth(5);
  LambdaThetaLowerPlus->SetLineColor(2);
  // LambdaThetaLowerPlus->SetFillColor(kRed-3);
  // LambdaThetaLowerPlus->SetFillStyle(1001);
  LambdaThetaLowerPlus->Draw("same");

  TF1 *LambdaThetaLowerMinus = new TF1("LambdaThetaLowerMinus","-1.-2*x",-1.,0.);
  LambdaThetaLowerMinus->SetLineWidth(5);
  LambdaThetaLowerMinus->SetLineColor(2);
  LambdaThetaLowerMinus->Draw("same");

  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");

  gPad->SaveAs("pngResults/Faccioli1.pdf", "recreate");



  new TCanvas;
  // TF2 *LambdaThetaPhiUpper = new TF2("LambdaThetaPhiUpper", "0.5*TMath::Sqrt(1-2*x-2*x*y-y*y)", -1., 1., -1., 1.);
  // LambdaThetaPhiUpper->Draw("colZ");
  TF1 *LambdaThetaPhiUpper = new TF1("LambdaThetaPhiUpper","0.5*TMath::Sqrt(1-x*x)",-1.,1.);
  LambdaThetaPhiUpper->SetTitle("");

  gPad->SetMargin(0.13,0.10,0.12,0.10);
  // gPad->SetLeftMargin(0.2);
  // // gPad->SetLeftMargin(0.12);
  // // gPad->SetRightMargin(0.1);
  // gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.1);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();


  LambdaThetaPhiUpper->GetXaxis()->SetTitleOffset(1.15);
  // LambdaThetaPhiUpper->GetYaxis()->SetTitleOffset(1.25);
  LambdaThetaPhiUpper->GetYaxis()->SetTitleOffset(1.25);
  LambdaThetaPhiUpper->GetXaxis()->SetTitleSize(0.045);
  LambdaThetaPhiUpper->GetYaxis()->SetTitleSize(0.045);
  LambdaThetaPhiUpper->GetXaxis()->SetLabelSize(0.045);
  LambdaThetaPhiUpper->GetYaxis()->SetLabelSize(0.045);
  LambdaThetaPhiUpper->GetXaxis()->SetTitleFont(42);
  LambdaThetaPhiUpper->GetYaxis()->SetTitleFont(42);
  LambdaThetaPhiUpper->GetXaxis()->SetLabelFont(42);
  LambdaThetaPhiUpper->GetYaxis()->SetLabelFont(42);
  LambdaThetaPhiUpper->SetLineWidth(5);
  LambdaThetaPhiUpper->SetLineColor(2);


  LambdaThetaPhiUpper->Draw();
  LambdaThetaPhiUpper->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LambdaThetaPhiUpper->GetXaxis()->SetTitle("#lambda_{#theta}");
  LambdaThetaPhiUpper->GetYaxis()->SetRangeUser(-1.,1.);
  // LambdaThetaPhiUpper->GetXaxis()->SetRangeUser(-1.,1.);
  // LambdaThetaPhiUpper->GetYaxis()->SetRangeUser(-1.4,1.4);
  // LambdaThetaPhiUpper->GetXaxis()->SetRangeUser(-1.4,1.4);
  // LambdaThetaPhiUpper->GetYaxis()->SetRangeUser(-2.,2.);
  LambdaThetaPhiUpper->GetXaxis()->SetRangeUser(-2.,2.);
  TF1 *LambdaThetaPhiLower = new TF1("LambdaThetaPhiLower","-0.5*TMath::Sqrt(1-x*x)",-1.,1.);
  LambdaThetaPhiLower->SetLineWidth(5);
  LambdaThetaPhiLower->SetLineColor(2);
  // LambdaThetaPhiLower->GetYaxis()->SetRangeUser(-2.,2.);
  LambdaThetaPhiLower->GetXaxis()->SetRangeUser(-2.,2.);
  LambdaThetaPhiLower->Draw("same");
  TLatex* latex6 = new TLatex();
  latex6->SetTextSize(0.045);
  latex6->SetTextFont(42);
  latex6->SetTextAlign(11);
  latex6->SetNDC();
  latex6->DrawLatex(0.31,0.94,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");


  gPad->SaveAs("pngResults/Faccioli2.pdf", "recreate");

  // TF2 *LambdaThetaPhiLower = new TF2("LambdaThetaPhiLower", "-0.5*TMath::Sqrt(1-2*x-2*x*y-y*y)", -1., 1., -1., 1.);
  // LambdaThetaPhiLower->Draw("same");







  // // LambdaPhi vs LambdaTheta
  // Int_t n = 50;
  // Double_t x0[n], y0[n], x1[n], y1[n];
  // for (size_t i = 0; i < 20; i++) {
  //   x0[i] = -1. + 2.*(Double_t)i/20.;
  //   if ( x0[i] > 0.) {
  //     y0[i] = LambdaThetaLowerPlus->Eval(x0[i]);
  //   } else {
  //     y0[i] = LambdaThetaLowerMinus->Eval(x0[i]);
  //   }
  // }
  // for (size_t i = 0; i < 10; i++) {
  //   x0[i+20] = 1.;
  //   y0[i+20] = 1. - 2.*(Double_t)i/10.;
  // }
  // for (size_t i = 0; i < 20; i++) {
  //   x0[i+30] = -1. + 2.*(Double_t)i/20.;
  //   if ( x0[i+30] > 0.) {
  //     y0[i+30] = LambdaThetaLowerPlus->Eval(x0[i+30]);
  //   } else {
  //     y0[i+30] = LambdaThetaLowerMinus->Eval(x0[i+30]);
  //   }
  // }
  //
  //
  //
  //
  // // LambdaThetaPhi vs LambdaTheta
  // for (size_t i = 0; i < 25; i++) {
  //   x1[i] = -1. + 2.*(Double_t)i/25.;
  //   x1[n-i] = -1. + 2.*(Double_t)i/25.;
  //   y1[i] = LambdaThetaPhiUpper->Eval(x1[i]);
  //   y1[n-1] = LambdaThetaPhiLower->Eval(x1[i]);
  // }
  // // for (size_t i = 0; i < 25; i++) {
  // //   x0[i+25] = 1. - 2.*(Double_t)i/25.;
  // //   y0[i+25] = LambdaThetaPhiLower->Eval(x1[i+25]);
  // // }
  //
  //
  // new TCanvas;
  // TGraph *LPvsLTgr  = new TGraph(n,x0,y0);
  // TGraph *LTPvsLTgr = new TGraph(n,x1,y1);
  // LPvsLTgr->SetFillStyle(3013);
  // LPvsLTgr->SetFillColor(16);
  // // LPvsLTgr->Draw("f");
  // LPvsLTgr->Draw();
  //
  // new TCanvas;
  // LTPvsLTgr->SetFillStyle(3013);
  // LTPvsLTgr->SetFillColor(16);
  // // LTPvsLTgr->Draw("f");
  // LTPvsLTgr->Draw();





}
