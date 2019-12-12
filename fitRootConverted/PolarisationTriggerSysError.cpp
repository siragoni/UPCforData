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
void PolarisationTriggerSysError(){

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
            Form( "pngResults/PolTrigger%d/1Dresults/Parameters_SigEx_0_FitRange_0_HE.root",
                  SigExBin+1
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

  Double_t MaxLambdaThetaPlusError     = *max_element( LambdaThetaPlusError,    LambdaThetaPlusError    + 7 );
  Double_t MaxLambdaPhiPlusError       = *max_element( LambdaPhiPlusError,      LambdaPhiPlusError      + 7 );
  Double_t MaxLambdaThetaPhiPlusError  = *max_element( LambdaThetaPhiPlusError, LambdaThetaPhiPlusError + 7 );
  Double_t MinLambdaThetaMinError      = *max_element( LambdaThetaMinError,     LambdaThetaMinError     + 7 );
  Double_t MinLambdaPhiMinError        = *max_element( LambdaPhiMinError,       LambdaPhiMinError       + 7 );
  Double_t MinLambdaThetaPhiMinError   = *max_element( LambdaThetaPhiMinError,  LambdaThetaPhiMinError  + 7 );

  Double_t LambdaThetaSysErr           = TMath::Abs( MaxLambdaThetaPlusError    - MinLambdaThetaMinError    )*0.500;
  Double_t LambdaPhiSysErr             = TMath::Abs( MaxLambdaPhiPlusError      - MinLambdaPhiMinError      )*0.500;
  Double_t LambdaThetaPhiSysErr        = TMath::Abs( MaxLambdaThetaPhiPlusError - MinLambdaThetaPhiMinError )*0.500;


  Double_t PercentLambdaThetaSysErr    = LambdaTheta2[0]    != 0 ? LambdaThetaSysErr    / TMath::Abs(LambdaTheta2[0])    : 0 ;
  Double_t PercentLambdaPhiSysErr      = LambdaPhi2[0]      != 0 ? LambdaPhiSysErr      / TMath::Abs(LambdaPhi2[0])      : 0 ;
  Double_t PercentLambdaThetaPhiSysErr = LambdaThetaPhi2[0] != 0 ? LambdaThetaPhiSysErr / TMath::Abs(LambdaThetaPhi2[0]) : 0 ;

  new TCanvas;
  TGraphErrors* graph1 = new TGraphErrors( 7, Xentries, LambdaTheta2, Zeroes, LambdaThetaErr2 );

  graph1->SetName("LambdaTheta");
  graph1->SetTitle(Form("LambdaThetaSys = %f / %f = %f", LambdaThetaSysErr, LambdaTheta2[0], PercentLambdaThetaSysErr));
  graph1->SetFillColor(1);
  graph1->SetMarkerColor(4);
  graph1->SetMarkerStyle(21);
  graph1->SetMarkerSize(1.3);

  // Define xAxis and yAxis. We're going to re-use these variables later.
  TAxis* xAxis = 0;
  TAxis* yAxis = 0;

  // Set the axis labels. Note the use of TLatex on the y-axis title.
  xAxis = graph1->GetXaxis();
  xAxis->SetTitle("Different modes");
  xAxis->CenterTitle( kTRUE );
  xAxis->SetTitleOffset( 1.2 );

  yAxis = graph1->GetYaxis();
  yAxis->SetTitle("#lambda{#theta}");
  yAxis->CenterTitle( kTRUE );
  yAxis->SetTitleOffset( 1.2 );

  // Draw the graph on the canvas.
  graph1->Draw("AP");
  // canvas1->Update();


  new TCanvas;
  TGraphErrors* graph2 = new TGraphErrors( 7, Xentries, LambdaPhi2, Zeroes, LambdaPhiErr2 );

  graph2->SetName("LambdaPhi");
  graph2->SetTitle(Form("LambdaPhiSys = %f / %f = %f", LambdaPhiSysErr, LambdaPhi2[0], PercentLambdaPhiSysErr));
  graph2->SetFillColor(1);
  graph2->SetMarkerColor(4);
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerSize(1.3);

  // Define xAxis and yAxis. We're going to re-use these variables later.
  TAxis* xAxis2 = 0;
  TAxis* yAxis2 = 0;

  // Set the axis labels. Note the use of TLatex on the y-axis title.
  xAxis2 = graph2->GetXaxis();
  xAxis2->SetTitle("Different modes");
  xAxis2->CenterTitle( kTRUE );
  xAxis2->SetTitleOffset( 1.2 );

  yAxis2 = graph2->GetYaxis();
  yAxis2->SetTitle("#lambda{#phi}");
  yAxis2->CenterTitle( kTRUE );
  yAxis2->SetTitleOffset( 1.2 );

  // Draw the graph on the canvas.
  graph2->Draw("AP");
  // canvas1->Update();


  new TCanvas;
  TGraphErrors* graph3 = new TGraphErrors( 7, Xentries, LambdaThetaPhi2, Zeroes, LambdaThetaPhiErr2 );

  graph3->SetName("LambdaThetaPhi");
  graph3->SetTitle(Form("LambdaThetaPhiSys = %f / %f = %f", LambdaThetaPhiSysErr, LambdaThetaPhi2[0], PercentLambdaThetaPhiSysErr));
  graph3->SetFillColor(1);
  graph3->SetMarkerColor(4);
  graph3->SetMarkerStyle(21);
  graph3->SetMarkerSize(1.3);

  // Define xAxis and yAxis. We're going to re-use these variables later.
  TAxis* xAxis3 = 0;
  TAxis* yAxis3 = 0;

  // Set the axis labels. Note the use of TLatex on the y-axis title.
  xAxis3 = graph3->GetXaxis();
  xAxis3->SetTitle("Different modes");
  xAxis3->CenterTitle( kTRUE );
  xAxis3->SetTitleOffset( 1.2 );

  yAxis3 = graph3->GetYaxis();
  yAxis3->SetTitle("#lambda{#theta#phi}");
  yAxis3->CenterTitle( kTRUE );
  yAxis3->SetTitleOffset( 1.2 );

  // Draw the graph on the canvas.
  graph3->Draw("AP");
  // canvas1->Update();


}
