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
void PolarisationSigExtrSystematicsCS(){

  /* - Open all files.
   * -
   * - 6 Signal Extraction bins.
   * - 3 Range Variation bins.
   * -
   */
  TDatime d;
  TFile*    FitResultFile[6][3];
  TH1F*     SavedHisto[6][3];
  Double_t  LambdaTheta[6][3];
  Double_t  LambdaPhi[6][3];
  Double_t  LambdaThetaPhi[6][3];
  Double_t  LambdaThetaErr[6][3];
  Double_t  LambdaPhiErr[6][3];
  Double_t  LambdaThetaPhiErr[6][3];
  Double_t  LambdaTheta2[6*3];
  Double_t  LambdaPhi2[6*3];
  Double_t  LambdaThetaPhi2[6*3];
  Double_t  LambdaThetaErr2[6*3];
  Double_t  LambdaPhiErr2[6*3];
  Double_t  LambdaThetaPhiErr2[6*3];
  Double_t  Zeroes[6*3];
  Double_t  Xentries[6*3];
  Double_t  LambdaThetaPlusError[6*3];
  Double_t  LambdaPhiPlusError[6*3];
  Double_t  LambdaThetaPhiPlusError[6*3];
  Double_t  LambdaThetaMinError[6*3];
  Double_t  LambdaPhiMinError[6*3];
  Double_t  LambdaThetaPhiMinError[6*3];
  for ( Int_t SigExBin = 0; SigExBin < 6; SigExBin++ ) {
    for ( Int_t FitRangeMode = 0; FitRangeMode < 3; FitRangeMode++ ) {
      FitResultFile[SigExBin][FitRangeMode]            = 0x0;
      SavedHisto[SigExBin][FitRangeMode]               = 0x0;
      LambdaTheta[SigExBin][FitRangeMode]              = 0x0;
      LambdaPhi[SigExBin][FitRangeMode]                = 0x0;
      LambdaThetaPhi[SigExBin][FitRangeMode]           = 0x0;
      LambdaThetaErr[SigExBin][FitRangeMode]           = 0x0;
      LambdaPhiErr[SigExBin][FitRangeMode]             = 0x0;
      LambdaThetaPhiErr[SigExBin][FitRangeMode]        = 0x0;
      LambdaTheta2[SigExBin*3+FitRangeMode]            = 0x0;
      LambdaPhi2[SigExBin*3+FitRangeMode]              = 0x0;
      LambdaThetaPhi2[SigExBin*3+FitRangeMode]         = 0x0;
      LambdaThetaErr2[SigExBin*3+FitRangeMode]         = 0x0;
      LambdaPhiErr2[SigExBin*3+FitRangeMode]           = 0x0;
      LambdaThetaPhiErr2[SigExBin*3+FitRangeMode]      = 0x0;
      Zeroes[SigExBin*3+FitRangeMode]                  = 0x0;
      Xentries[SigExBin*3+FitRangeMode]                = 0x0;
      LambdaThetaPlusError[SigExBin*3+FitRangeMode]    = 0x0;
      LambdaPhiPlusError[SigExBin*3+FitRangeMode]      = 0x0;
      LambdaThetaPhiPlusError[SigExBin*3+FitRangeMode] = 0x0;
      LambdaThetaMinError[SigExBin*3+FitRangeMode]     = 0x0;
      LambdaPhiMinError[SigExBin*3+FitRangeMode]       = 0x0;
      LambdaThetaPhiMinError[SigExBin*3+FitRangeMode]  = 0x0;

    }
  }
  for ( Int_t SigExBin = 0; SigExBin < 6; SigExBin++ ) {
    for ( Int_t FitRangeMode = 0; FitRangeMode < 3; FitRangeMode++ ) {
      FitResultFile[SigExBin][FitRangeMode] = new TFile(
            Form( "pngResults/%d-%2.2d-%2.2d/1Dresults/Parameters_SigEx_%d_FitRange_%d_CS.root",
                  d.GetYear(), d.GetMonth(), d.GetDay(),
                  SigExBin, FitRangeMode
                  )
            );
      SavedHisto[SigExBin][FitRangeMode]               = (TH1F*) FitResultFile[SigExBin][FitRangeMode]->Get("SavingParamH");
      LambdaTheta[SigExBin][FitRangeMode]              = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(1);
      LambdaPhi[SigExBin][FitRangeMode]                = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(2);
      LambdaThetaPhi[SigExBin][FitRangeMode]           = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(3);
      LambdaThetaErr[SigExBin][FitRangeMode]           = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(6);
      LambdaPhiErr[SigExBin][FitRangeMode]             = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(7);
      LambdaThetaPhiErr[SigExBin][FitRangeMode]        = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(8);
      LambdaTheta2[SigExBin*3+FitRangeMode]            = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(1);
      LambdaPhi2[SigExBin*3+FitRangeMode]              = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(2);
      LambdaThetaPhi2[SigExBin*3+FitRangeMode]         = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(3);
      LambdaThetaErr2[SigExBin*3+FitRangeMode]         = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(6);
      LambdaPhiErr2[SigExBin*3+FitRangeMode]           = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(7);
      LambdaThetaPhiErr2[SigExBin*3+FitRangeMode]      = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(8);
      LambdaThetaPlusError[SigExBin*3+FitRangeMode]    = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(1) + SavedHisto[SigExBin][FitRangeMode]->GetBinContent(6);
      LambdaPhiPlusError[SigExBin*3+FitRangeMode]      = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(2) + SavedHisto[SigExBin][FitRangeMode]->GetBinContent(7);
      LambdaThetaPhiPlusError[SigExBin*3+FitRangeMode] = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(3) + SavedHisto[SigExBin][FitRangeMode]->GetBinContent(8);
      LambdaThetaMinError[SigExBin*3+FitRangeMode]     = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(1) - SavedHisto[SigExBin][FitRangeMode]->GetBinContent(6);
      LambdaPhiMinError[SigExBin*3+FitRangeMode]       = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(2) - SavedHisto[SigExBin][FitRangeMode]->GetBinContent(7);
      LambdaThetaPhiMinError[SigExBin*3+FitRangeMode]  = SavedHisto[SigExBin][FitRangeMode]->GetBinContent(3) - SavedHisto[SigExBin][FitRangeMode]->GetBinContent(8);
      Zeroes[SigExBin*3+FitRangeMode]             = 0;
      Xentries[SigExBin*3+FitRangeMode]           = SigExBin*3+FitRangeMode+0.5;
      FitResultFile[SigExBin][FitRangeMode]->Close();
    }
  }

  cout << " max element LambdaThetaPlusError    is: " << *max_element( LambdaThetaPlusError,    LambdaThetaPlusError    + 18 ) << endl;
  cout << " max element LambdaPhiPlusError      is: " << *max_element( LambdaPhiPlusError,      LambdaPhiPlusError      + 18 ) << endl;
  cout << " max element LambdaThetaPhiPlusError is: " << *max_element( LambdaThetaPhiPlusError, LambdaThetaPhiPlusError + 18 ) << endl;
  cout << " min element LambdaThetaMinError     is: " << *min_element( LambdaThetaMinError,     LambdaThetaMinError     + 18 ) << endl;
  cout << " min element LambdaPhiMinError       is: " << *min_element( LambdaPhiMinError,       LambdaPhiMinError       + 18 ) << endl;
  cout << " min element LambdaThetaPhiMinError  is: " << *min_element( LambdaThetaPhiMinError,  LambdaThetaPhiMinError  + 18 ) << endl;

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


  Double_t PercentLambdaThetaSysErr    = LambdaTheta[0][0]    != 0 ? LambdaThetaSysErr    / TMath::Abs(LambdaTheta[0][0])    : 0 ;
  Double_t PercentLambdaPhiSysErr      = LambdaPhi[0][0]      != 0 ? LambdaPhiSysErr      / TMath::Abs(LambdaPhi[0][0])      : 0 ;
  Double_t PercentLambdaThetaPhiSysErr = LambdaThetaPhi[0][0] != 0 ? LambdaThetaPhiSysErr / TMath::Abs(LambdaThetaPhi[0][0]) : 0 ;

  new TCanvas;
  TGraphErrors* graph1 = new TGraphErrors( 18, Xentries, LambdaTheta2, Zeroes, LambdaThetaErr2 );

  graph1->SetName("LambdaTheta");
  graph1->SetTitle(Form("LambdaThetaSys = %f / %f = %f", LambdaThetaSysErr, LambdaTheta[0][0], PercentLambdaThetaSysErr));
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
  gPad->SaveAs("pngResults/CosThetaSysExtractCS.png", "recreate");
  // canvas1->Update();


  new TCanvas;
  TGraphErrors* graph2 = new TGraphErrors( 18, Xentries, LambdaPhi2, Zeroes, LambdaPhiErr2 );

  graph2->SetName("LambdaPhi");
  graph2->SetTitle(Form("LambdaPhiSys = %f / %f = %f", LambdaPhiSysErr, LambdaPhi[0][0], PercentLambdaPhiSysErr));
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
  gPad->SaveAs("pngResults/PhiSysExtractCS.png", "recreate");
  // canvas1->Update();


  new TCanvas;
  TGraphErrors* graph3 = new TGraphErrors( 18, Xentries, LambdaThetaPhi2, Zeroes, LambdaThetaPhiErr2 );

  graph3->SetName("LambdaThetaPhi");
  graph3->SetTitle(Form("LambdaThetaPhiSys = %f / %f = %f", LambdaThetaPhiSysErr, LambdaThetaPhi[0][0], PercentLambdaThetaPhiSysErr));
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
  gPad->SaveAs("pngResults/TildeSysExtractCS.png", "recreate");
  // canvas1->Update();


}
