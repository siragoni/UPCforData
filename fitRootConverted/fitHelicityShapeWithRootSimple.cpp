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
#include "TStyle.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"


//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void fitHelicityShape(const char* AnalysisName, const char* MonteCarloName){
  /* - There are three cases for the selectionFlag:
     - 1) = 0 ; this implies the traditional pt-integrated plot;
     - 2) = 1 ; this is instead the coherent component;
     - 3) = 2 ; this is the incoherent component;
     - 4) = 3 ; ******************* ;
     -
   */
  TFile* dataList = new TFile(AnalysisName);
  TDirectory* dirData = dataList->GetDirectory("MyTask");
  TFile* mcList = new TFile(MonteCarloName);
  TDirectory* dirMC = mcList->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listingsData;
  dirData->GetObject("MyOutputContainer", listingsData);
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
  TH1F *fAngularDistribOfPositiveMuonRestFrameJPsiH = (TH1F*)listingsData->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(4);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw();
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->SaveAs("rawDistr.png", "RECREATE");


  TH1F *fMCAngularDistribOfPositiveMuonRestFrameJPsiH = (TH1F*)listingsMC->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(4);
  fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Draw();
  // fMCAngularDistribOfPositiveMuonRestFrameJPsiH->SaveAs("mcReconDistr.png", "RECREATE");


  TH1F *fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH = (TH1F*)listingsMC->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Sumw2();
  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Rebin(4);
  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Draw();
  // fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->SaveAs("mcGenDistr.png", "RECREATE");


  TH1F *acceptance = (TH1F*)fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Clone("acceptance");
  acceptance->Sumw2();
  acceptance->Divide(fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH);
  acceptance->SetMarkerStyle(21);
  acceptance->SetLineColor(kBlue+1);
  acceptance->SetLineWidth(2);
  acceptance->Draw("ep");
  // acceptance->SaveAs("acceptance.png", "RECREATE");


  // TH1F *fCorrectedShape = (TH1F*)fAngularDistribOfPositiveMuonRestFrameJPsiH->Clone("fCorrectedShape");
  // fCorrectedShape->Sumw2();
  // fCorrectedShape->Divide(acceptance);
  // fCorrectedShape->Rebin(2);
  // fCorrectedShape->SetMarkerStyle(21);
  // fCorrectedShape->SetLineColor(kBlue+1);
  // fCorrectedShape->SetLineWidth(2);
  // fCorrectedShape->Draw("ep");
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Divide(acceptance);
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(5);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->SetMarkerStyle(21);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->SetLineColor(kBlue+1);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->SetLineWidth(2);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("ep");
  TH1F *fCorrectedShape = fAngularDistribOfPositiveMuonRestFrameJPsiH;

  /* - Do the parabolic fit with Root only.
     - Retrieve the parameters later for Roofit plotting.
     -
  */
  TF1* ParabolicFit = new TF1("ParabolicFit","[0]*(1+[1]*x*x)",-1., 1.);
  ParabolicFit->SetNpx(1000);
  ParabolicFit->SetParameter(1, 1);
  ParabolicFit->SetParLimits(1, -3, +3);
  fCorrectedShape->Fit( ParabolicFit,"","", -0.5, 0.5 );
  Double_t normalization = ParabolicFit->Integral(-0.5, 0.5);
  cout << "normalization " << normalization << endl;

  fCorrectedShape->SetLineColor(kBlue);
  fCorrectedShape->SetLineStyle(kSolid);
  fCorrectedShape->SetLineWidth(3);
  fCorrectedShape->SetMarkerStyle(kFullCircle);
  fCorrectedShape->SetMarkerSize(1);
  fCorrectedShape->GetXaxis()->SetTitle("cos(#theta)");
  fCorrectedShape->GetYaxis()->SetTitle( Form( "Counts / (%.3f a.u.)",
                                               fCorrectedShape->GetXaxis()->GetBinWidth(1)
                                               )
                                              );
  fCorrectedShape->SetTitle("");
  TCanvas* ZNAEnergy = new TCanvas( "CosTheta", "CosTheta", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  /* - Beautifying is starting now.
     -
   */
  fCorrectedShape->GetXaxis()->SetTitleOffset(1.25);
  // fZNAEnergyAgainstEntriesH->GetYaxis()->SetTitleOffset(1.25);
  fCorrectedShape->GetYaxis()->SetTitleOffset(1.45);
  fCorrectedShape->GetXaxis()->SetTitleSize(0.045);
  fCorrectedShape->GetYaxis()->SetTitleSize(0.045);
  fCorrectedShape->GetXaxis()->SetLabelSize(0.045);
  fCorrectedShape->GetYaxis()->SetLabelSize(0.045);
  fCorrectedShape->GetXaxis()->SetTitleFont(42);
  fCorrectedShape->GetYaxis()->SetTitleFont(42);
  fCorrectedShape->GetXaxis()->SetLabelFont(42);
  fCorrectedShape->GetYaxis()->SetLabelFont(42);
  fCorrectedShape->GetXaxis()->SetNdivisions(408);
  fCorrectedShape->GetYaxis()->SetRangeUser(5, fCorrectedShape->GetMaximum()*10.);
  // gPad ->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fCorrectedShape->Draw("SAME");
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
  latex->DrawLatex(0.55,0.72,Form("#lambda_{#theta} = %.2f #pm %.2f",
                                  ParabolicFit->GetParameter(1),
                                  ParabolicFit->GetParError(1)
                                  )
                                );
  latex->DrawLatex(0.55,0.18,Form( "      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                   ParabolicFit->GetChisquare(),
                                   ParabolicFit->GetNDF(),
                                   ParabolicFit->GetChisquare()/ParabolicFit->GetNDF()
                                  )
                                 );


  gPad->SaveAs("pngResults/CosThetaHeRooFit.png", "RECREATE");

  new TCanvas;
  acceptance->Draw("ep");
  gPad->SaveAs("pngResults/CosThetaHeACCxEFF.png", "RECREATE");

}
