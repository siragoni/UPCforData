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





double gauss2D(double *x, double *par) {
  double z1 = double((x[0]-par[1])/par[2]);
  double z2 = double((x[1]-par[3])/par[4]);
  return par[0]*exp(-0.5*(z1*z1+z2*z2));
}
double my2Dfunc(double *x, double *par) {
  return gauss2D(x,&par[0]) + gauss2D(x,&par[5]);
}




//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void fitBinMigration(){
  TFile* mcList = new TFile("AnalysisResultsLHC18l7_long_faultyevents_31032021.root");
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
  TFile* mcList2 = new TFile("AnalysisResultsLHC18l7_long_rapidity_30032021.root");
  TDirectory* dirMC2 = mcList2->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listingsMC2;
  dirMC2  ->GetObject("MyOutputContainer", listingsMC2);
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







  TH2F *fBinMigrationPhi = (TH2F*)listingsMC->FindObject("fBinMigrationForPhiHelicityH");
  fBinMigrationPhi->Sumw2();
  TCanvas* ZNAEnergy2 = new TCanvas( "BinMigration2", "BinMigration2", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad ->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);



  new TCanvas;
  fBinMigrationPhi->Rebin2D(10, 10);
  fBinMigrationPhi->Draw("colZ");
  gPad->SetGrid();
  gPad->SaveAs("pngResults/BinMigration2D.png", "RECREATE");





  fBinMigrationPhi->Rebin2D(2, 2);
  fBinMigrationPhi->Rebin2D(2, 2);



  // fBinMigrationPhi->Rebin2D(2, 2);


  new TCanvas;
  gStyle->SetOptStat("nemr");
  gPad->SetMargin(0.13,0.1,0.12,0.1);
  // fBinMigrationPhi->Rebin2D(10, 10);
  fBinMigrationPhi->GetXaxis()->SetTitle("cos(#theta_{generated})");
  fBinMigrationPhi->GetYaxis()->SetTitle("cos(#theta_{reconstructed})");
  fBinMigrationPhi->GetXaxis()->SetTitleOffset(1.25);
  // fBinMigrationPhi->GetYaxis()->SetTitleOffset(1.25);
  fBinMigrationPhi->GetYaxis()->SetTitleOffset(1.45);
  fBinMigrationPhi->GetXaxis()->SetTitleSize(0.045);
  fBinMigrationPhi->GetYaxis()->SetTitleSize(0.045);
  fBinMigrationPhi->GetXaxis()->SetLabelSize(0.045);
  fBinMigrationPhi->GetYaxis()->SetLabelSize(0.045);
  fBinMigrationPhi->GetXaxis()->SetTitleFont(42);
  fBinMigrationPhi->GetYaxis()->SetTitleFont(42);
  fBinMigrationPhi->GetXaxis()->SetLabelFont(42);
  fBinMigrationPhi->GetYaxis()->SetLabelFont(42);
  // fBinMigrationPhi->GetXaxis()->SetNdivisions(408);
  // fBinMigrationPhi->GetYaxis()->SetRangeUser(-1, 1);
  // fBinMigrationPhi->GetXaxis()->SetRangeUser(-1, 1);
  // fBinMigrationPhi->Draw("text colZ");
  fBinMigrationPhi->Draw("text colZ");








  TH2F *fBinMigrationPhi2 = (TH2F*)listingsMC2->FindObject("fBinMigrationForPhiHelicityH");
  fBinMigrationPhi2->Sumw2();
  fBinMigrationPhi2->Rebin2D(10, 10);
  fBinMigrationPhi2->Rebin2D(2, 2);
  fBinMigrationPhi2->Rebin2D(2, 2);


  new TCanvas;
  fBinMigrationPhi2->Draw("colZ text");
  new TCanvas;
  fBinMigrationPhi2->Draw("surf");


  TH2F *fBinMigrationPhi2_clone = (TH2F*)fBinMigrationPhi2->Clone("fBinMigrationPhi2_clone");

  new TCanvas;
  fBinMigrationPhi2_clone->Add(fBinMigrationPhi, -1);
  fBinMigrationPhi2_clone->Draw("surf");
  new TCanvas;
  fBinMigrationPhi2_clone->Draw("colZ text");



  new TCanvas;
  TH1F* fPhiReconstructed = (TH1F*)  fBinMigrationPhi2_clone->ProjectionY("fPhiReconstructed");
  fPhiReconstructed->Draw();


  TFile* file2 = new TFile("SavedPhiLongitudinalAfterBkg.root", "recreate");
  fPhiReconstructed->Write();
  file2->Close();


}
