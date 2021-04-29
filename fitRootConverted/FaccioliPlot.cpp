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


#include "TH2.h"


//_____________________________________________________________________________
/* -
 * - Original macro:
 * - https://root.cern/doc/v610/graphShade_8C.html
 */
void FaccioliPlot()
{

  TCanvas *c1 = new TCanvas("c1","LambdaPhi vs LambdaTheta",200,10,700,500);
  // TCanvas *c1 = new TCanvas("c1","LambdaPhi vs LambdaTheta",100,100,400,500);
  // TCanvas *c_super = new TCanvas("c_super","comparing modes",100,100,1200,700);
  // c1->SetGrid();
  // c1->DrawFrame(-1.,-1.,1.,1.);
  c1->DrawFrame(-10.,-10.,10.,10.);
  // c1->cd();
  // gPad->SetLeftMargin(0.2);
  // // gPad->SetLeftMargin(0.12);
  // // gPad->SetRightMargin(0.1);
  // gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.1);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();

  // const Int_t n = 20;
  // Double_t x[n], y[n],ymin[n], ymax[n];
  // Int_t i;
  // for (i=0;i<n;i++) {
  //   x[i] = 0.1+i*0.1;
  //   ymax[i] = 10*sin(x[i]+0.2);
  //   ymin[i] = 8*sin(x[i]+0.1);
  //   y[i] = 9*sin(x[i]+0.15);
  // }




  const Int_t n = 40;
  Double_t x[2*n+10], y[2*n+10];
  for (Int_t i=0;i<n;i++) {
    y[i]    = -1.+(Double_t)i*(2.)/(Double_t)n;
    if( y[i] > 0. ) {
      x[i]  = -1. + 2.*y[i];
    } else {
      x[i]  = -1. - 2.*y[i];
    }
  }

  for (Int_t i=n;i<2*n;i++) {
    y[i]    = -1.+(Double_t)(i-n)*(2.)/(Double_t)n;
    if( y[i] > 0. ) {
      x[i]  = -1. + 2.*y[i];
    } else {
      x[i]  = -1. - 2.*y[i];
    }
    // y[i]    = -1.+(Double_t)(i-n)*(2.)/(Double_t)n;
    // x[i]    = 1.;
    // y[i]    = 0;
    // x[i]    = 0;
  }


  for (Int_t i=2*n;i<2*n+5;i++) {
    y[i]    = -1.+(Double_t)(i-2*n)*(2.)/5.;
    x[i]    = 1.;
  }


  // TGraph *grmin = new TGraph(n,x,ymin);
  // TGraph *grmax = new TGraph(n,x,ymax);
  // TGraph *gr    = new TGraph(n,x,y);
  // TH1* help = new TH1F("help", "help", 100, -1., 1.);
  // gStyle->SetOptStat(0);
  // help->SetTitle("");
  // help->GetXaxis()->SetRangeUser(-1., 1.0);
  // help->GetXaxis()->SetRangeUser(-1., 1.0);
  // help->GetXaxis()->SetTitle("#lambda_{#theta}");
  // help->GetYaxis()->SetTitle("#lambda_{#varphi}");
  // help->Draw();

  TGraph *grshade = new TGraph(2*n+5);
  for (Int_t i=0;i<2*n+5;i++) {
  // for (Int_t i=0;i<2*n;i++) {
    grshade->SetPoint(i,x[i],y[i]);
  }
  grshade->SetFillStyle(3013);
  grshade->SetFillColor(16);
  // grshade->Draw("f");



  grshade->SetTitle("");
  grshade->GetXaxis()->SetTitle("#lambda_{#theta}");
  grshade->GetYaxis()->SetTitle("#lambda_{#varphi}");
  grshade->GetXaxis()->CenterTitle(true);
  grshade->GetYaxis()->CenterTitle(true);
  // grshade->Draw("f");

  grshade->GetYaxis()->SetTitleOffset(2.2);
  grshade->GetXaxis()->SetTitleOffset(1.1);
  // grshade->GetYaxis()->SetTitleOffset(0.8);
  // grshade->GetXaxis()->SetTitleOffset(0.8);
  // grshade->GetYaxis()->SetTitleSize(0.07);
  // grshade->GetXaxis()->SetTitleSize(0.07);
  grshade->Draw("fsame");

  TAxis *axis1 = grshade->GetXaxis();
  axis1->SetTitle("#lambda_{#theta}");
  axis1->SetLimits(-1.,1.);                 // along X
  grshade->GetHistogram()->SetMaximum( 1.);   // along
  grshade->GetHistogram()->SetMinimum(-1.);  //   Y
  TAxis *axis2 = grshade->GetYaxis();
  axis2->SetTitle("#lambda_{#varphi}");




  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");

  // grmin->Draw("l");
  // grmax->Draw("l");
  // gr->SetLineWidth(4);
  // gr->SetMarkerColor(4);
  // gr->SetMarkerStyle(21);
  // gr->Draw("CP");
}
