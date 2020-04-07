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


#include "TH2.h"


//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void fitBinMigration( const char* MonteCarloName = "MCtrainResults/2019-06-24/kCohJpsiToMu/AnalysisResults.root" ){
  TFile* mcList = new TFile(MonteCarloName);
  // TFile* mcList = new TFile("AnalysisResultsLHC1815o15072019.root");
  TDirectory* dirMC = mcList->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listingsMC;
  dirMC  ->GetObject("MyOutputContainer", listingsMC);
  /* - We now do the same as before to ascertain if the TList was there and
   * - to try to retrieve the plots. Result:
   *   listings->ls()
   *     OBJ: TList	  MyOutputContainer	          Doubly linked list          : 0
   *     OBJ: TH1F	  fNumberMuonsH	              fNumberMuonsH               : 0 at: 0x5a145f0
   *     OBJ: TH1F	  fCounterH	                  fCounterH                   : 0 at: 0x5a3b570
   *     OBJ: TH1F	  fEtaMuonH	                  fEtaMuonH                   : 0 at: 0x5a3ba80
   *     OBJ: TH1F	  fRAbsMuonH	                fRAbsMuonH                  : 0 at: 0x5a3c0c0
   *     OBJ: TH1F	  fInvariantMassDistributionH	fInvariantMassDistributionH : 0 at: 0x5a3c720
   */
  TH2F *fBinMigration = (TH2F*)listingsMC->FindObject("fBinMigrationHelicityH");
  fBinMigration->Sumw2();
  TProfile *fBinMigrationProfile = fBinMigration->ProfileY();
  // fBinMigrationProfile->Rebin(4);
  // fCorrectedShape->SetMarkerStyle(21);
  // fCorrectedShape->SetLineColor(kBlue+1);
  // fCorrectedShape->SetLineWidth(2);
  fBinMigrationProfile->Draw();

  /* - Do the linear fit with Root only.
     - Retrieve the parameters later for Roofit plotting.
     -
  */
  TF1* LinearFit = new TF1("LinearFit","[0]*x+[1]",-1., 1.);
  LinearFit->SetNpx(1000);
  // LinearFit->SetParameter(1, 1);
  // LinearFit->SetParLimits(1, -3, +3);
  fBinMigrationProfile->Fit( LinearFit,"","", -0.5, 0.5 );
  fBinMigrationProfile->SetLineColor(kBlue);
  fBinMigrationProfile->SetLineStyle(kSolid);
  fBinMigrationProfile->SetLineWidth(3);
  fBinMigrationProfile->SetMarkerStyle(kFullCircle);
  fBinMigrationProfile->SetMarkerSize(1);
  fBinMigrationProfile->GetXaxis()->SetTitle("cos(#theta_{generated})");
  fBinMigrationProfile->GetYaxis()->SetTitle("cos(#theta_{recon})");
  fBinMigrationProfile->SetTitle("");
  TCanvas* ZNAEnergy = new TCanvas( "BinMigration", "BinMigration", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  /* - Beautifying is starting now.
     -
   */
  fBinMigrationProfile->GetXaxis()->SetTitleOffset(1.25);
  // fBinMigrationProfile->GetYaxis()->SetTitleOffset(1.25);
  fBinMigrationProfile->GetYaxis()->SetTitleOffset(1.45);
  fBinMigrationProfile->GetXaxis()->SetTitleSize(0.045);
  fBinMigrationProfile->GetYaxis()->SetTitleSize(0.045);
  fBinMigrationProfile->GetXaxis()->SetLabelSize(0.045);
  fBinMigrationProfile->GetYaxis()->SetLabelSize(0.045);
  fBinMigrationProfile->GetXaxis()->SetTitleFont(42);
  fBinMigrationProfile->GetYaxis()->SetTitleFont(42);
  fBinMigrationProfile->GetXaxis()->SetLabelFont(42);
  fBinMigrationProfile->GetYaxis()->SetLabelFont(42);
  fBinMigrationProfile->GetXaxis()->SetNdivisions(408);
  fBinMigrationProfile->GetYaxis()->SetRangeUser(-1.0, fBinMigrationProfile->GetMaximum()*10.);
  // gPad ->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fBinMigrationProfile->Draw("SAME");
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, LHC18q+LHC18r data");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));
  latex->DrawLatex(0.55,0.72,"a*x+b");
  latex->DrawLatex(0.55,0.66,Form("a = %.4f #pm %.4f",
                                  LinearFit->GetParameter(0),
                                  LinearFit->GetParError(0)
                                  )
                                );
  latex->DrawLatex(0.55,0.60,Form("b = %.4f #pm %.4f",
                                  LinearFit->GetParameter(1),
                                  LinearFit->GetParError(1)
                                  )
                                );
  latex->DrawLatex(0.55,0.18,Form( "      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                   LinearFit->GetChisquare(),
                                   LinearFit->GetNDF(),
                                   LinearFit->GetChisquare()/LinearFit->GetNDF()
                                  )
                                 );


  gPad->SaveAs("pngResults/BinMigration.png", "RECREATE");

  new TCanvas;
  // gPad->SetMargin(0.13,0.1,0.12,0.1);
  fBinMigration->Rebin2D(10, 10);
  // fBinMigration->GetXaxis()->SetTitle("cos(#theta_{generated})");
  // fBinMigration->GetYaxis()->SetTitle("cos(#theta_{reconstructed})");
  // fBinMigration->GetXaxis()->SetTitleOffset(1.25);
  // // fBinMigration->GetYaxis()->SetTitleOffset(1.25);
  // fBinMigration->GetYaxis()->SetTitleOffset(1.45);
  // fBinMigration->GetXaxis()->SetTitleSize(0.045);
  // fBinMigration->GetYaxis()->SetTitleSize(0.045);
  // fBinMigration->GetXaxis()->SetLabelSize(0.045);
  // fBinMigration->GetYaxis()->SetLabelSize(0.045);
  // fBinMigration->GetXaxis()->SetTitleFont(42);
  // fBinMigration->GetYaxis()->SetTitleFont(42);
  // fBinMigration->GetXaxis()->SetLabelFont(42);
  // fBinMigration->GetYaxis()->SetLabelFont(42);
  // // fBinMigration->GetXaxis()->SetNdivisions(408);
  // fBinMigration->GetYaxis()->SetRangeUser(-1, 1);
  // fBinMigration->GetXaxis()->SetRangeUser(-1, 1);
  fBinMigration->Draw("colZ");
  gPad->SetGrid();
  gPad->SaveAs("pngResults/BinMigration2D.png", "RECREATE");



  TH1F* PurityH    = new TH1F("PurityH",    "PurityH",    100, -1, 1);
  TH1F* StabilityH = new TH1F("StabilityH", "StabilityH", 100, -1, 1);
  /* - NB:
   * - x-axis = generated;
   * - y-axis = reconstructed;
   * -
   * - Hence:
   * - ProjX => Recon distribution;
   * - ProjY => Gener distribution.
   */
  Double_t GenAndReconInBin = 0;
  Double_t GenInBin         = 0;
  Double_t ReconInBin       = 0;
  for( Int_t iLoop = 0; iLoop < 100; iLoop++ ) {
    GenAndReconInBin = 0;
    GenInBin         = 0;
    ReconInBin       = 0;
    for( Int_t iLoop2  = 0; iLoop2 < 100; iLoop2++ ) {
      GenInBin   += fBinMigration->GetBinContent( iLoop+1,  iLoop2+1 );
      ReconInBin += fBinMigration->GetBinContent( iLoop2+1, iLoop+1  );

    }
    GenAndReconInBin = fBinMigration->GetBinContent( iLoop+1, iLoop+1 );
    if ( GenAndReconInBin != 0 && GenInBin != 0  && ReconInBin != 0 ) {
      PurityH   ->SetBinContent( iLoop+1, GenAndReconInBin/ReconInBin );
      StabilityH->SetBinContent( iLoop+1, GenAndReconInBin/GenInBin   );
      PurityH   ->SetBinError  ( iLoop+1, TMath::Sqrt( 1/GenAndReconInBin + 1/ReconInBin ) );
      StabilityH->SetBinError  ( iLoop+1, TMath::Sqrt( 1/GenAndReconInBin + 1/GenInBin   ) );
    } else {
      PurityH   ->SetBinContent( iLoop+1, 0 );
      StabilityH->SetBinContent( iLoop+1, 0 );
      PurityH   ->SetBinError  ( iLoop+1, 0 );
      StabilityH->SetBinError  ( iLoop+1, 0 );
    }
  }
  new TCanvas;
  gPad->SetMargin(0.13,0.1,0.12,0.1);
  // gPad->SetTitle(  ";#cos(#theta_{generated}); #cos(#theta_{gen && rec in the same bin})/ #cos(#theta_{generated})");
  PurityH->GetXaxis()->SetTitle("cos(#theta_{generated})");
  PurityH->GetYaxis()->SetTitle("cos(#theta_{gen && rec in the same bin})/ cos(#theta_{generated})");
  PurityH->GetXaxis()->SetTitleOffset(1.25);
  // PurityH->GetYaxis()->SetTitleOffset(1.25);
  PurityH->GetYaxis()->SetTitleOffset(1.45);
  PurityH->GetXaxis()->SetTitleSize(0.045);
  PurityH->GetYaxis()->SetTitleSize(0.045);
  PurityH->GetXaxis()->SetLabelSize(0.045);
  PurityH->GetYaxis()->SetLabelSize(0.045);
  PurityH->GetXaxis()->SetTitleFont(42);
  PurityH->GetYaxis()->SetTitleFont(42);
  PurityH->GetXaxis()->SetLabelFont(42);
  PurityH->GetYaxis()->SetLabelFont(42);
  // PurityH->GetXaxis()->SetNdivisions(408);
  PurityH->GetYaxis()->SetRangeUser(-0.15, 0.55);
  PurityH->GetXaxis()->SetRangeUser(-0.66, 0.66);
  PurityH->Draw();
  new TCanvas;
  StabilityH->Draw();



  fBinMigration->Rebin2D(2, 2);
  TH1F* PurityHv2    = new TH1F("PurityHv2",    "PurityHv2",    50, -1, 1);
  TH1F* StabilityHv2 = new TH1F("StabilityHv2", "StabilityHv2", 50, -1, 1);
  /* - NB:
   * - x-axis = generated;
   * - y-axis = reconstructed;
   * -
   * - Hence:
   * - ProjX => Recon distribution;
   * - ProjY => Gener distribution.
   */
  Double_t GenAndReconInBinv2 = 0;
  Double_t GenInBinv2         = 0;
  Double_t ReconInBinv2       = 0;
  for( Int_t iLoop = 0; iLoop < 50; iLoop++ ) {
    GenAndReconInBinv2 = 0;
    GenInBinv2         = 0;
    ReconInBinv2       = 0;
    for( Int_t iLoop2  = 0; iLoop2 < 50; iLoop2++ ) {
      GenInBinv2   += fBinMigration->GetBinContent( iLoop+1,  iLoop2+1 );
      ReconInBinv2 += fBinMigration->GetBinContent( iLoop2+1, iLoop+1  );

    }
    GenAndReconInBinv2  = fBinMigration->GetBinContent( iLoop+1, iLoop+1 );
    if ( GenAndReconInBinv2 != 0 && GenInBinv2 != 0  && ReconInBinv2 != 0 ) {
      PurityHv2   ->SetBinContent( iLoop+1, GenAndReconInBinv2/ReconInBinv2 );
      StabilityHv2->SetBinContent( iLoop+1, GenAndReconInBinv2/GenInBinv2   );
      PurityHv2   ->SetBinError  ( iLoop+1, TMath::Sqrt( 1/GenAndReconInBinv2 + 1/ReconInBinv2 ) );
      StabilityHv2->SetBinError  ( iLoop+1, TMath::Sqrt( 1/GenAndReconInBinv2 + 1/GenInBinv2   ) );
    } else {
      PurityHv2   ->SetBinContent( iLoop+1, 0 );
      StabilityHv2->SetBinContent( iLoop+1, 0 );
      PurityHv2   ->SetBinError  ( iLoop+1, 0 );
      StabilityHv2->SetBinError  ( iLoop+1, 0 );
    }
  }
  new TCanvas;
  PurityHv2->Draw();
  new TCanvas;
  StabilityHv2->Draw();



  fBinMigration->Rebin2D(2, 2);
  TH1F* PurityHv3    = new TH1F("PurityHv3",    "PurityHv3",    25, -1, 1);
  TH1F* StabilityHv3 = new TH1F("StabilityHv3", "StabilityHv3", 25, -1, 1);
  TH1F* GenHv3       = new TH1F("GenHv3",       "GenHv3",       25, -1, 1);
  TH1F* GenAndRecHv3 = new TH1F("GenAndRecHv3", "GenAndRecHv3", 25, -1, 1);
  /* - NB:
   * - x-axis = generated;
   * - y-axis = reconstructed;
   * -
   * - Hence:
   * - ProjX => Recon distribution;
   * - ProjY => Gener distribution.
   */
  Double_t GenAndReconInBinv3 = 0;
  Double_t GenInBinv3         = 0;
  Double_t ReconInBinv3       = 0;
  for( Int_t iLoop = 0; iLoop < 25; iLoop++ ) {
    GenAndReconInBinv3 = 0;
    GenInBinv3         = 0;
    ReconInBinv3       = 0;
    for( Int_t iLoop2 = 0; iLoop2 < 25; iLoop2++ ) {
      GenInBinv3   += fBinMigration->GetBinContent( iLoop+1,  iLoop2+1 );
      ReconInBinv3 += fBinMigration->GetBinContent( iLoop2+1, iLoop+1  );

    }
    GenAndReconInBinv3  = fBinMigration->GetBinContent( iLoop+1, iLoop+1 );
    if ( GenAndReconInBinv3 != 0 && GenInBinv3 != 0  && ReconInBinv3 != 0 ) {
      PurityHv3   ->SetBinContent( iLoop+1, GenAndReconInBinv3/ReconInBinv3 );
      StabilityHv3->SetBinContent( iLoop+1, GenAndReconInBinv3/GenInBinv3   );
      PurityHv3   ->SetBinError  ( iLoop+1, TMath::Sqrt( 1/GenAndReconInBinv3 + 1/ReconInBinv3 ) );
      StabilityHv3->SetBinError  ( iLoop+1, TMath::Sqrt( 1/GenAndReconInBinv3 + 1/GenInBinv3   ) );
      GenAndRecHv3->SetBinContent( iLoop+1, GenAndReconInBinv3 );
      GenHv3      ->SetBinContent( iLoop+1, GenInBinv3 );
      GenAndRecHv3->SetBinError  ( iLoop+1, TMath::Sqrt( GenAndReconInBinv3 ) );
      GenHv3      ->SetBinError  ( iLoop+1, TMath::Sqrt( GenInBinv3 ) );
    } else {
      PurityHv3   ->SetBinContent( iLoop+1, 0 );
      StabilityHv3->SetBinContent( iLoop+1, 0 );
      PurityHv3   ->SetBinError  ( iLoop+1, 0 );
      StabilityHv3->SetBinError  ( iLoop+1, 0 );
      GenAndRecHv3->SetBinContent( iLoop+1, GenAndReconInBinv3 );
      GenHv3      ->SetBinContent( iLoop+1, GenInBinv3 );
      GenAndRecHv3->SetBinError  ( iLoop+1, TMath::Sqrt( GenAndReconInBinv3 ) );
      GenHv3      ->SetBinError  ( iLoop+1, TMath::Sqrt( GenInBinv3 ) );
    }
  }
  new TCanvas;
  // gPad->SetTitle(  ";#cos(#theta_{generated}); #cos(#theta_{gen && rec in the same bin})/ #cos(#theta_{generated})");
  gPad->SetMargin(0.13,0.1,0.12,0.1);
  PurityHv3->GetXaxis()->SetTitle("cos(#theta_{generated})");
  PurityHv3->GetYaxis()->SetTitle("cos(#theta_{gen && rec in the same bin})/ cos(#theta_{generated})");
  PurityHv3->GetXaxis()->SetTitleOffset(1.25);
  // PurityHv3->GetYaxis()->SetTitleOffset(1.25);
  PurityHv3->GetYaxis()->SetTitleOffset(1.45);
  PurityHv3->GetXaxis()->SetTitleSize(0.045);
  PurityHv3->GetYaxis()->SetTitleSize(0.045);
  PurityHv3->GetXaxis()->SetLabelSize(0.045);
  PurityHv3->GetYaxis()->SetLabelSize(0.045);
  PurityHv3->GetXaxis()->SetTitleFont(42);
  PurityHv3->GetYaxis()->SetTitleFont(42);
  PurityHv3->GetXaxis()->SetLabelFont(42);
  PurityHv3->GetYaxis()->SetLabelFont(42);
  // PurityHv3->GetXaxis()->SetNdivisions(408);
  PurityHv3->GetYaxis()->SetRangeUser( 0.60, 0.85);
  PurityHv3->GetXaxis()->SetRangeUser(-0.7, 0.7);
  PurityHv3->Draw();
  PurityHv3->Draw();
  new TCanvas;
  StabilityHv3->Draw();
  // new TCanvas;
  // GenAndRecHv3->Draw();
  // new TCanvas;
  // GenHv3->Draw();


  new TCanvas;
  gStyle->SetOptStat("nemr");
  gPad->SetMargin(0.13,0.1,0.12,0.1);
  // fBinMigration->Rebin2D(10, 10);
  fBinMigration->GetXaxis()->SetTitle("cos(#theta_{generated})");
  fBinMigration->GetYaxis()->SetTitle("cos(#theta_{reconstructed})");
  fBinMigration->GetXaxis()->SetTitleOffset(1.25);
  // fBinMigration->GetYaxis()->SetTitleOffset(1.25);
  fBinMigration->GetYaxis()->SetTitleOffset(1.45);
  fBinMigration->GetXaxis()->SetTitleSize(0.045);
  fBinMigration->GetYaxis()->SetTitleSize(0.045);
  fBinMigration->GetXaxis()->SetLabelSize(0.045);
  fBinMigration->GetYaxis()->SetLabelSize(0.045);
  fBinMigration->GetXaxis()->SetTitleFont(42);
  fBinMigration->GetYaxis()->SetTitleFont(42);
  fBinMigration->GetXaxis()->SetLabelFont(42);
  fBinMigration->GetYaxis()->SetLabelFont(42);
  // fBinMigration->GetXaxis()->SetNdivisions(408);
  fBinMigration->GetYaxis()->SetRangeUser(-1, 1);
  fBinMigration->GetXaxis()->SetRangeUser(-1, 1);
  // fBinMigration->Draw("text colZ");
  fBinMigration->Draw("text colZ");
}
