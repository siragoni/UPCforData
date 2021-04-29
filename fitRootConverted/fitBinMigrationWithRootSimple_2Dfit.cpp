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
  TH2F *fBinMigration = (TH2F*)listingsMC->FindObject("fBinMigrationHelicityH");
  fBinMigration->Sumw2();
  TCanvas* ZNAEnergy = new TCanvas( "BinMigration", "BinMigration", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  // gPad ->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);



  new TCanvas;
  fBinMigration->Rebin2D(10, 10);
  fBinMigration->Draw("colZ");
  gPad->SetGrid();
  gPad->SaveAs("pngResults/BinMigration2D.png", "RECREATE");





  fBinMigration->Rebin2D(2, 2);
  fBinMigration->Rebin2D(2, 2);



  // fBinMigration->Rebin2D(2, 2);


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






  // TF2 *fitfunction = new TF2("fitfunction",gauss2D,-1,1,-1,1, 5);
  // // fitfunction->SetNpx(50);
  // // fitfunction->SetNpy(50);
  // fitfunction->SetNpx(25);
  // fitfunction->SetNpy(25);
  // fitfunction->FixParameter(0, 3100. );
  // fitfunction->FixParameter(1, 0. );
  // // fitfunction->SetParameter(2, 0.3);
  // fitfunction->FixParameter(2, 5.10690e-01);
  // fitfunction->FixParameter(3, 0. );
  // fitfunction->FixParameter(4, 2.10074e-01);
  // // fitfunction->SetParLimits(2, 0.01, 1);
  // // fitfunction->SetParLimits(4, 0.01, 1);
  //
  // fBinMigration->Fit(fitfunction);
  // new TCanvas;
  // gPad->SetMargin(0.13,0.1,0.12,0.1);
  // fBinMigration->Draw("surf");
  // new TCanvas;
  // gPad->SetMargin(0.13,0.1,0.12,0.1);
  // fitfunction->Draw("surf");

















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
  TH2F *fBinMigration2 = (TH2F*)listingsMC2->FindObject("fBinMigrationHelicityH");
  fBinMigration2->Sumw2();
  fBinMigration2->Rebin2D(10, 10);
  fBinMigration2->Rebin2D(2, 2);
  fBinMigration2->Rebin2D(2, 2);




  // TH2F* CorrectionH = (TH2F*) fitfunction->GetHistogram();
  // new TCanvas;
  // CorrectionH->Draw("surf");

  new TCanvas;
  fBinMigration2->Draw("colZ text");
  new TCanvas;
  fBinMigration2->Draw("surf");
  // fBinMigration2->Add(CorrectionH, -1);
  // new TCanvas;
  // fBinMigration2->Draw("colZ text");
  // new TCanvas;
  // fBinMigration2->Draw("surf");




  // TH2F* binmap = new TH2F("binmap", "binmap", 25, -1, 1, 25, -1, 1);
  // Int_t ibin = 1;
  // Int_t nXbins = binmap->GetNbinsX();
  // Int_t nYbins = binmap->GetNbinsY();
  // for( Int_t ix = 1; ix < nXbins+1; ix++){
  //   for( Int_t iy = 1; iy < nYbins+1; iy++){
  //     binmap->SetBinContent(ibin, ibin);
  //     ibin++;
  //   }
  // }
  // new TCanvas;
  // binmap->Draw("text");



  TH2F* CorrectionMatrix = (TH2F*) fBinMigration->Clone("CorrectionMatrix");
  new TCanvas;
  CorrectionMatrix->Draw("surf");


  // Int_t binx = h2.GetXaxis()->FindBin(x);
  // Int_t biny = h2.GetYaxis()->FindBin(y);
  // Int_t bin = h2->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number

  Int_t binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.625);
  Int_t biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.625);
  Int_t binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.625);
  Int_t biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.625);
  Int_t bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  Int_t bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.625);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.575);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.625);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.575);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.575);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.625);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.575);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.625);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.575);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.575);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.575);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.575);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.575);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.48);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.575);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.48);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.48);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.575);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.48);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.575);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.48);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.48);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.48);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.48);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.48);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.4);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.48);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.4);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.4);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.48);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.4);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.48);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.4);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.4);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.4);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.4);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.4);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.33);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.4);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.33);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.33);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.4);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.33);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.4);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.33);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.33);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.33);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.33);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.33);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.24);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.33);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.24);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.24);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.33);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.24);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.33);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.24);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.24);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.24);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.24);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.24);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.16);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.24);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.16);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.16);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.24);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.24);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.16);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.16);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.16);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.16);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.08);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.08);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.08);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.16);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.08);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.16);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.08);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.08);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.08);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.08);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.08);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.0);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.0);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.08);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.08);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.08);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(-0.);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(-0.);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(-0.);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.08);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.08);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.08);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.08);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.08);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.08);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.08);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.08);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.08);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.16);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.08);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.16);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.08);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.08);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.16);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.16);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.16);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.24);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.16);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.24);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.24);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.16);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.24);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.16);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.24);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.24);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.24);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.24);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.24);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.33);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.24);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.33);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.33);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.24);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.33);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.24);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.33);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.33);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.33);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.33);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.33);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.4);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.33);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.4);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.4);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.33);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.4);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.33);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.4);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.4);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.4);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.4);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.4);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.48);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.4);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.48);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.48);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.4);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.48);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.4);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.48);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.48);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.48);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.48);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.48);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.575);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.48);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.575);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.575);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.48);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.575);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.48);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.575);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.575);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.575);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.575);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.575);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.625);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.575);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.625);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  //
  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.625);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.575);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.625);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.575);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  binx  = CorrectionMatrix->GetXaxis()->FindBin(0.625);
  biny  = CorrectionMatrix->GetYaxis()->FindBin(0.625);
  binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.625);
  biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.625);
  bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));

  // binx  = CorrectionMatrix->GetXaxis()->FindBin(0.625);
  // biny  = CorrectionMatrix->GetYaxis()->FindBin(0.625);
  // binx2 = CorrectionMatrix->GetXaxis()->FindBin(-0.625);
  // biny2 = CorrectionMatrix->GetYaxis()->FindBin(0.625);
  // bin   = CorrectionMatrix->GetBin(binx,biny,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  // bin2  = CorrectionMatrix->GetBin(binx2,biny2,0); //see doc of TH1::GetBin h2->SetBinContent(bin,c); //where bin is the linearized bin number
  // CorrectionMatrix->SetBinContent(bin, CorrectionMatrix->GetBinContent(bin2));







  new TCanvas;
  CorrectionMatrix->Draw("colZ text");


  new TCanvas;
  fBinMigration2->Add(CorrectionMatrix, -1);
  fBinMigration2->Draw("surf");




  new TCanvas;
  TH1F* fCosThetaReconstructed = (TH1F*)  fBinMigration2->ProjectionY("fCosThetaReconstructed");
  fCosThetaReconstructed->Draw();


  TFile* file = new TFile("SavedCosThetaLongitudinalAfterBkg.root", "recreate");
  fCosThetaReconstructed->Write();
  file->Close();
























  // TH2F *fBinMigrationPhi = (TH2F*)listingsMC->FindObject("fBinMigrationForPhiHelicityH");
  // fBinMigrationPhi->Sumw2();
  // TCanvas* ZNAEnergy2 = new TCanvas( "BinMigration2", "BinMigration2", 900, 800 );
  // gPad->SetMargin(0.13,0.01,0.12,0.01);
  // // gPad ->SetLogy();
  // gStyle->SetOptFit(0);
  // gStyle->SetOptStat(0);
  //
  //
  //
  // new TCanvas;
  // fBinMigrationPhi->Rebin2D(10, 10);
  // fBinMigrationPhi->Draw("colZ");
  // gPad->SetGrid();
  // gPad->SaveAs("pngResults/BinMigration2D.png", "RECREATE");
  //
  //
  //
  //
  //
  // fBinMigrationPhi->Rebin2D(2, 2);
  // fBinMigrationPhi->Rebin2D(2, 2);
  //
  //
  //
  // // fBinMigrationPhi->Rebin2D(2, 2);
  //
  //
  // new TCanvas;
  // gStyle->SetOptStat("nemr");
  // gPad->SetMargin(0.13,0.1,0.12,0.1);
  // // fBinMigrationPhi->Rebin2D(10, 10);
  // fBinMigrationPhi->GetXaxis()->SetTitle("cos(#theta_{generated})");
  // fBinMigrationPhi->GetYaxis()->SetTitle("cos(#theta_{reconstructed})");
  // fBinMigrationPhi->GetXaxis()->SetTitleOffset(1.25);
  // // fBinMigrationPhi->GetYaxis()->SetTitleOffset(1.25);
  // fBinMigrationPhi->GetYaxis()->SetTitleOffset(1.45);
  // fBinMigrationPhi->GetXaxis()->SetTitleSize(0.045);
  // fBinMigrationPhi->GetYaxis()->SetTitleSize(0.045);
  // fBinMigrationPhi->GetXaxis()->SetLabelSize(0.045);
  // fBinMigrationPhi->GetYaxis()->SetLabelSize(0.045);
  // fBinMigrationPhi->GetXaxis()->SetTitleFont(42);
  // fBinMigrationPhi->GetYaxis()->SetTitleFont(42);
  // fBinMigrationPhi->GetXaxis()->SetLabelFont(42);
  // fBinMigrationPhi->GetYaxis()->SetLabelFont(42);
  // // fBinMigrationPhi->GetXaxis()->SetNdivisions(408);
  // // fBinMigrationPhi->GetYaxis()->SetRangeUser(-1, 1);
  // // fBinMigrationPhi->GetXaxis()->SetRangeUser(-1, 1);
  // // fBinMigrationPhi->Draw("text colZ");
  // fBinMigrationPhi->Draw("text colZ");
  //
  //
  //
  //
  //
  //
  //
  //
  // TH2F *fBinMigrationPhi2 = (TH2F*)listingsMC2->FindObject("fBinMigrationForPhiHelicityH");
  // fBinMigrationPhi2->Sumw2();
  // fBinMigrationPhi2->Rebin2D(10, 10);
  // fBinMigrationPhi2->Rebin2D(2, 2);
  // fBinMigrationPhi2->Rebin2D(2, 2);
  //
  //
  // new TCanvas;
  // fBinMigrationPhi2->Draw("colZ text");
  // new TCanvas;
  // fBinMigrationPhi2->Draw("surf");
  //
  //
  //
  // new TCanvas;
  // fBinMigrationPhi2->Add(fBinMigrationPhi, -1);
  // fBinMigrationPhi2->Draw("surf");

}
