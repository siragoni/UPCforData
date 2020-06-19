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

#include "TDatime.h"

#include "TH2.h"
#include "TGraphErrors.h"
#include "TGraph.h"

//_____________________________________________________________________________
/* - Compute the Signal extraction systematics.
 * - Compute the Fit Range Variation systematics.
 * -
 */
void PolarisationTriggerSysErrorCS(){

  /* - Open all files.
   * -
   * - 6 Signal Extraction bins.
   * - 3 Range Variation bins.
   * -
   */
  TFile*    FitResultFile[7];
  TH1F*     SavedHisto[7];
  // Double_t  LambdaTheta[6][3];
  // Double_t  LambdaPhi[6][3];
  // Double_t  LambdaThetaPhi[6][3];
  // Double_t  LambdaThetaErr[6][3];
  // Double_t  LambdaPhiErr[6][3];
  // Double_t  LambdaThetaPhiErr[6][3];
  Double_t  LambdaTheta2[7];
  Double_t  LambdaPhi2[7];
  Double_t  LambdaThetaPhi2[7];
  Double_t  LambdaThetaErr2[7];
  Double_t  LambdaPhiErr2[7];
  Double_t  LambdaThetaPhiErr2[7];
  Double_t  Zeroes[7];
  Double_t  Xentries[7];
  Double_t  LambdaThetaPlusError[7];
  Double_t  LambdaPhiPlusError[7];
  Double_t  LambdaThetaPhiPlusError[7];
  Double_t  LambdaThetaMinError[7];
  Double_t  LambdaPhiMinError[7];
  Double_t  LambdaThetaPhiMinError[7];
  for ( Int_t SigExBin = 0; SigExBin < 7; SigExBin++ ) {
      FitResultFile[SigExBin]           = 0x0;
      SavedHisto[SigExBin]              = 0x0;
      // LambdaTheta[SigExBin][FitRangeMode]              = 0x0;
      // LambdaPhi[SigExBin][FitRangeMode]                = 0x0;
      // LambdaThetaPhi[SigExBin][FitRangeMode]           = 0x0;
      // LambdaThetaErr[SigExBin][FitRangeMode]           = 0x0;
      // LambdaPhiErr[SigExBin][FitRangeMode]             = 0x0;
      // LambdaThetaPhiErr[SigExBin][FitRangeMode]        = 0x0;
      LambdaTheta2[SigExBin]            = 0x0;
      LambdaPhi2[SigExBin]              = 0x0;
      LambdaThetaPhi2[SigExBin]         = 0x0;
      LambdaThetaErr2[SigExBin]         = 0x0;
      LambdaPhiErr2[SigExBin]           = 0x0;
      LambdaThetaPhiErr2[SigExBin]      = 0x0;
      Zeroes[SigExBin]                  = 0x0;
      Xentries[SigExBin]                = 0x0;
      LambdaThetaPlusError[SigExBin]    = 0x0;
      LambdaPhiPlusError[SigExBin]      = 0x0;
      LambdaThetaPhiPlusError[SigExBin] = 0x0;
      LambdaThetaMinError[SigExBin]     = 0x0;
      LambdaPhiMinError[SigExBin]       = 0x0;
      LambdaThetaPhiMinError[SigExBin]  = 0x0;
  }
  for ( Int_t SigExBin = 0; SigExBin < 7; SigExBin++ ) {
      FitResultFile[SigExBin] = new TFile(
            Form( "pngResults/PolTrigger%d/2DCS/Parameters_SigEx_%d_FitRange_0_CS.root",
                  SigExBin+1, SigExBin+1
                  )
            );
      SavedHisto[SigExBin]               = (TH1F*) FitResultFile[SigExBin]->Get("SavingParamH");
      // LambdaTheta[SigExBin][FitRangeMode]              = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(1);
      // LambdaPhi[SigExBin][FitRangeMode]                = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(2);
      // LambdaThetaPhi[SigExBin][FitRangeMode]           = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(3);
      // LambdaThetaErr[SigExBin][FitRangeMode]           = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(6);
      // LambdaPhiErr[SigExBin][FitRangeMode]             = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(7);
      // LambdaThetaPhiErr[SigExBin][FitRangeMode]        = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(8);
      LambdaTheta2[SigExBin]            = SavedHisto[SigExBin]->GetBinContent(1);
      LambdaPhi2[SigExBin]              = SavedHisto[SigExBin]->GetBinContent(2);
      LambdaThetaPhi2[SigExBin]         = SavedHisto[SigExBin]->GetBinContent(3);
      LambdaThetaErr2[SigExBin]         = SavedHisto[SigExBin]->GetBinContent(6);
      LambdaPhiErr2[SigExBin]           = SavedHisto[SigExBin]->GetBinContent(7);
      LambdaThetaPhiErr2[SigExBin]      = SavedHisto[SigExBin]->GetBinContent(8);
      LambdaThetaPlusError[SigExBin]    = SavedHisto[SigExBin]->GetBinContent(1) + SavedHisto[SigExBin]->GetBinContent(6);
      LambdaPhiPlusError[SigExBin]      = SavedHisto[SigExBin]->GetBinContent(2) + SavedHisto[SigExBin]->GetBinContent(7);
      LambdaThetaPhiPlusError[SigExBin] = SavedHisto[SigExBin]->GetBinContent(3) + SavedHisto[SigExBin]->GetBinContent(8);
      LambdaThetaMinError[SigExBin]     = SavedHisto[SigExBin]->GetBinContent(1) - SavedHisto[SigExBin]->GetBinContent(6);
      LambdaPhiMinError[SigExBin]       = SavedHisto[SigExBin]->GetBinContent(2) - SavedHisto[SigExBin]->GetBinContent(7);
      LambdaThetaPhiMinError[SigExBin]  = SavedHisto[SigExBin]->GetBinContent(3) - SavedHisto[SigExBin]->GetBinContent(8);
      Zeroes[SigExBin]                  = 0;
      Xentries[SigExBin]                = SigExBin+0.5;
      FitResultFile[SigExBin]->Close();
  }

  cout << " max element LambdaThetaPlusError    is: " << *max_element( LambdaThetaPlusError,    LambdaThetaPlusError    + 7 ) << endl;
  cout << " max element LambdaPhiPlusError      is: " << *max_element( LambdaPhiPlusError,      LambdaPhiPlusError      + 7 ) << endl;
  cout << " max element LambdaThetaPhiPlusError is: " << *max_element( LambdaThetaPhiPlusError, LambdaThetaPhiPlusError + 7 ) << endl;
  cout << " min element LambdaThetaMinError     is: " << *min_element( LambdaThetaMinError,     LambdaThetaMinError     + 7 ) << endl;
  cout << " min element LambdaPhiMinError       is: " << *min_element( LambdaPhiMinError,       LambdaPhiMinError       + 7 ) << endl;
  cout << " min element LambdaThetaPhiMinError  is: " << *min_element( LambdaThetaPhiMinError,  LambdaThetaPhiMinError  + 7 ) << endl;

  // Double_t MaxLambdaThetaPlusError     = *max_element( LambdaThetaPlusError,    LambdaThetaPlusError    + 7 );
  // Double_t MaxLambdaPhiPlusError       = *max_element( LambdaPhiPlusError,      LambdaPhiPlusError      + 7 );
  // Double_t MaxLambdaThetaPhiPlusError  = *max_element( LambdaThetaPhiPlusError, LambdaThetaPhiPlusError + 7 );
  // Double_t MinLambdaThetaMinError      = *max_element( LambdaThetaMinError,     LambdaThetaMinError     + 7 );
  // Double_t MinLambdaPhiMinError        = *max_element( LambdaPhiMinError,       LambdaPhiMinError       + 7 );
  // Double_t MinLambdaThetaPhiMinError   = *max_element( LambdaThetaPhiMinError,  LambdaThetaPhiMinError  + 7 );

  Double_t MaxLambdaTheta      = *max_element( LambdaTheta2,     LambdaTheta2     + 7 );
  Double_t MaxLambdaPhi        = *max_element( LambdaPhi2,       LambdaPhi2       + 7 );
  Double_t MaxLambdaThetaPhi   = *max_element( LambdaThetaPhi2,  LambdaThetaPhi2  + 7 );
  Double_t MinLambdaTheta      = *min_element( LambdaTheta2,     LambdaTheta2     + 7 );
  Double_t MinLambdaPhi        = *min_element( LambdaPhi2,       LambdaPhi2       + 7 );
  Double_t MinLambdaThetaPhi   = *min_element( LambdaThetaPhi2,  LambdaThetaPhi2  + 7 );


  // Double_t LambdaThetaSysErr           = TMath::Abs( MaxLambdaThetaPlusError    - MinLambdaThetaMinError    )*0.500;
  // Double_t LambdaPhiSysErr             = TMath::Abs( MaxLambdaPhiPlusError      - MinLambdaPhiMinError      )*0.500;
  // Double_t LambdaThetaPhiSysErr        = TMath::Abs( MaxLambdaThetaPhiPlusError - MinLambdaThetaPhiMinError )*0.500;

  Double_t LambdaThetaSysErr           = TMath::Abs( MaxLambdaTheta    - MinLambdaTheta    )*0.500;
  Double_t LambdaPhiSysErr             = TMath::Abs( MaxLambdaPhi      - MinLambdaPhi      )*0.500;
  Double_t LambdaThetaPhiSysErr        = TMath::Abs( MaxLambdaThetaPhi - MinLambdaThetaPhi )*0.500;


  // Double_t PercentLambdaThetaSysErr    = LambdaTheta2[0]    != 0 ? LambdaThetaSysErr    / TMath::Abs(LambdaTheta2[0])    : 0 ;
  // Double_t PercentLambdaPhiSysErr      = LambdaPhi2[0]      != 0 ? LambdaPhiSysErr      / TMath::Abs(LambdaPhi2[0])      : 0 ;
  // Double_t PercentLambdaThetaPhiSysErr = LambdaThetaPhi2[0] != 0 ? LambdaThetaPhiSysErr / TMath::Abs(LambdaThetaPhi2[0]) : 0 ;

  Double_t PercentLambdaThetaSysErr    = LambdaThetaSysErr    / 1.187 ;
  Double_t PercentLambdaPhiSysErr      = LambdaPhiSysErr      / 0.011 ;
  Double_t PercentLambdaThetaPhiSysErr = LambdaThetaPhiSysErr / 0.096 ;


  // new TCanvas;
  TCanvas* c1 = new TCanvas("c1","c1",1200,1100);
  TGraphErrors* graph1 = new TGraphErrors( 7, Xentries, LambdaTheta2, Zeroes, LambdaThetaErr2 );

  //____________________________
  /* -
   * - COSMETICS
   */
  TString labels[7] = { "0.85", "0.90", "0.95", "1.00", "1.05", "1.10", "1.15" };
  TH1F *h = new TH1F("h",Form("LambdaThetaSys = %f / %f = %f", LambdaThetaSysErr, 1.187, PercentLambdaThetaSysErr),7,Xentries[0]-0.5,Xentries[6]+0.5);
  for (Int_t i=1;i<=7;i++) h->GetXaxis()->SetBinLabel(i,labels[i-1].Data());
  h->SetMaximum(1.8);
  h->SetMinimum(0.6);
  gStyle->SetOptStat(0);

  // Define xAxis and yAxis. We're going to re-use these variables later.
  TAxis* xAxis = 0;
  TAxis* yAxis = 0;

  // Set the axis labels. Note the use of TLatex on the y-axis title.
  xAxis = h->GetXaxis();
  xAxis->SetTitle("p_{T} threshold [GeV/#it{c}]");
  // xAxis->CenterTitle( kTRUE );
  xAxis->SetTitleOffset( 1.2 );

  yAxis = h->GetYaxis();
  yAxis->SetTitle("#lambda_{#theta}");
  yAxis->CenterTitle( kTRUE );
  yAxis->SetTitleOffset( 1.2 );

  h->Draw();
  // TGraph *gr = new TGraph(n,x,y);
  // gr->Draw("CP");


  graph1->SetName("LambdaTheta");
  graph1->SetTitle(Form("LambdaThetaSys = %f / %f = %f", LambdaThetaSysErr, 1.187, PercentLambdaThetaSysErr));
  graph1->SetFillColor(1);
  graph1->SetMarkerColor(4);
  graph1->SetMarkerStyle(21);
  graph1->SetMarkerSize(1.3);


  // Draw the graph on the canvas.
  // graph1->Draw("AP");
  graph1->Draw("Psame");
  gPad->SaveAs("pngResults/CosThetaTriggerCS_2D.png", "recreate");
  // canvas1->Update();


  // new TCanvas;
  TCanvas* c2 = new TCanvas("c2","c2",1200,1100);
  TGraphErrors* graph2 = new TGraphErrors( 7, Xentries, LambdaPhi2, Zeroes, LambdaPhiErr2 );

  TH1F *h2 = new TH1F("h2",Form("LambdaPhiSys = %f / %f = %f", LambdaPhiSysErr, 0.011, PercentLambdaPhiSysErr),7,Xentries[0]-0.5,Xentries[6]+0.5);
  for (Int_t i=1;i<=7;i++) h2->GetXaxis()->SetBinLabel(i,labels[i-1].Data());
  h2->SetMaximum(0.06);
  h2->SetMinimum(-0.04);
  gStyle->SetOptStat(0);

  // Define xAxis and yAxis. We're going to re-use these variables later.
  TAxis* xAxis2 = 0;
  TAxis* yAxis2 = 0;

  // Set the axis labels. Note the use of TLatex on the y-axis title.
  xAxis2 = h2->GetXaxis();
  xAxis2->SetTitle("p_{T} threshold [GeV/#it{c}]");
  // xAxis2->CenterTitle( kTRUE );
  xAxis2->SetTitleOffset( 1.2 );

  yAxis2 = h2->GetYaxis();
  yAxis2->SetTitle("#lambda_{#phi}");
  yAxis2->CenterTitle( kTRUE );
  yAxis2->SetTitleOffset( 1.2 );

  h2->Draw();


  graph2->SetName("LambdaPhi");
  graph2->SetTitle(Form("LambdaPhiSys = %f / %f = %f", LambdaPhiSysErr, 0.011, PercentLambdaPhiSysErr));
  graph2->SetFillColor(1);
  graph2->SetMarkerColor(4);
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerSize(1.3);

  // // Define xAxis and yAxis. We're going to re-use these variables later.
  // TAxis* xAxis2 = 0;
  // TAxis* yAxis2 = 0;
  //
  // // Set the axis labels. Note the use of TLatex on the y-axis title.
  // xAxis2 = graph2->GetXaxis();
  // xAxis2->SetTitle("Different modes");
  // xAxis2->CenterTitle( kTRUE );
  // xAxis2->SetTitleOffset( 1.2 );
  //
  // yAxis2 = graph2->GetYaxis();
  // yAxis2->SetTitle("#lambda{#phi}");
  // yAxis2->CenterTitle( kTRUE );
  // yAxis2->SetTitleOffset( 1.2 );

  // Draw the graph on the canvas.
  // graph2->Draw("AP");
  graph2->Draw("Psame");
  gPad->SaveAs("pngResults/PhiTriggerSysCS_2D.png", "recreate");
  // canvas1->Update();


  // new TCanvas;
  TCanvas* c3 = new TCanvas("c3","c3",1200,1100);
  TGraphErrors* graph3 = new TGraphErrors( 7, Xentries, LambdaThetaPhi2, Zeroes, LambdaThetaPhiErr2 );

  TH1F *h3 = new TH1F("h3",Form("LambdaThetaPhiSys = %f / %f = %f", LambdaThetaPhiSysErr, 0.096, PercentLambdaThetaPhiSysErr),7,Xentries[0]-0.5,Xentries[6]+0.5);
  for (Int_t i=1;i<=7;i++) h3->GetXaxis()->SetBinLabel(i,labels[i-1].Data());
  h3->SetMaximum(0.20);
  h3->SetMinimum(0.02);
  h3->SetMarkerColor(0);
  h3->SetLineColor(0);
  gStyle->SetOptStat(0);

  // Define xAxis and yAxis. We're going to re-use these variables later.
  TAxis* xAxis3 = 0;
  TAxis* yAxis3 = 0;

  // Set the axis labels. Note the use of TLatex on the y-axis title.
  xAxis3 = h3->GetXaxis();
  // xAxis3->SetTitleSize(1.1);
  xAxis3->SetTitle("p_{T} threshold [GeV/#it{c}]");
  // xAxis3->SetTitleSize(0.15);

  // xAxis3->CenterTitle( kTRUE );
  xAxis3->SetTitleOffset( 1.2 );

  yAxis3 = h3->GetYaxis();
  yAxis3->SetTitle("#lambda_{#theta#phi}");
  yAxis3->CenterTitle( kTRUE );
  yAxis3->SetTitleOffset( 1.2 );

  h3->Draw();



  graph3->SetName("LambdaThetaPhi");
  graph3->SetTitle(Form("LambdaThetaPhiSys = %f / %f = %f", LambdaThetaPhiSysErr, 0.096, PercentLambdaThetaPhiSysErr));
  graph3->SetFillColor(1);
  graph3->SetMarkerColor(4);
  graph3->SetMarkerStyle(21);
  graph3->SetMarkerSize(1.3);

  // Define xAxis and yAxis. We're going to re-use these variables later.
  // TAxis* xAxis3 = 0;
  // TAxis* yAxis3 = 0;
  //
  // // Set the axis labels. Note the use of TLatex on the y-axis title.
  // xAxis3 = graph3->GetXaxis();
  // xAxis3->SetTitle("Different modes");
  // xAxis3->CenterTitle( kTRUE );
  // xAxis3->SetTitleOffset( 1.2 );
  //
  // yAxis3 = graph3->GetYaxis();
  // yAxis3->SetTitle("#lambda{#theta#phi}");
  // yAxis3->CenterTitle( kTRUE );
  // yAxis3->SetTitleOffset( 1.2 );

  // Draw the graph on the canvas.
  graph3->Draw("Psame");
  gPad->SaveAs("pngResults/TildeTriggerSysCS_2D.png", "recreate");
  // canvas1->Update();


}
