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
/* - Histograms to be used for the fit.
 * - What happens is that we will interpolate the many points together...
 * -
 */
TH1F*     fThetaMC;
TH1F*     fTheta;
TH1F*     fPhiMC;
TH1F*     fPhi;
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void ThetaPhiMaps(){

  TFile* fileMC;
  fileMC = new TFile("AnalysisResultsLHC18l7_pol.root");
  TDirectory* dirMC;
  dirMC = fileMC->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listingsMC;
  dirMC->GetObject("MyOutputContainer", listingsMC);
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
  TFile* file;
  file = new TFile("AnalysisResultsLHC18qr_10012020.root");
  TDirectory* dir;
  dir = file->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);

  fThetaMC = (TH1F*)listingsMC->FindObject("fThetaMuonH");
  fPhiMC   = (TH1F*)listingsMC->FindObject("fPhiMuonH");
  fTheta   = (TH1F*)listings  ->FindObject("fThetaMuonH");
  fPhi     = (TH1F*)listings  ->FindObject("fPhiMuonH");


  fThetaMC->Rebin(4);
  fTheta  ->Rebin(4);
  fPhiMC  ->Rebin(200);
  fPhi    ->Rebin(200);

  fThetaMC->Sumw2();
  fTheta  ->Sumw2();
  fPhiMC  ->Sumw2();
  fPhi    ->Sumw2();

  Double_t Integral_fThetaMC = fThetaMC->Integral();
  Double_t Integral_fTheta   = fTheta  ->Integral();
  Double_t Integral_fPhiMC   = fPhiMC  ->Integral();
  Double_t Integral_fPhi     = fPhi    ->Integral();
  fThetaMC->Scale( 1/Integral_fThetaMC );
  fTheta  ->Scale( 1/Integral_fTheta   );
  fPhiMC  ->Scale( 1/Integral_fPhiMC   );
  fPhi    ->Scale( 1/Integral_fPhi     );


  fThetaMC-> SetLineColor(kRed);
  fPhiMC  -> SetLineColor(kRed);
  fTheta  -> SetLineColor(kBlue);
  fPhi    -> SetLineColor(kBlue);
  fThetaMC-> SetLineWidth(3);
  fPhiMC  -> SetLineWidth(3);
  fTheta  -> SetLineWidth(3);
  fPhi    -> SetLineWidth(3);



  // TLegend* l = new TLegend(0.2,0.55,0.5,0.85);
  // l->SetMargin(0.1);
  // l->SetBorderSize(0);
  // l->AddEntry(  fCohJpsiToMu,        "fCohJpsiToMu"       );
  // l->AddEntry(  fCohPsi2sToMu,       "fCohPsi2sToMu"      );
  // // l->AddEntry(  fCohPsi2sToMuPi,     "fCohPsi2sToMuPi"    );
  // l->AddEntry(  fIncohJpsiToMu,      "fIncohJpsiToMu"     );
  // l->AddEntry(  fIncohPsi2sToMu,     "fIncohPsi2sToMu"    );
  // // l->AddEntry(  fIncohPsi2sToMuPi,   "fIncohPsi2sToMuPi"  );
  // l->AddEntry(  fTwoGammaToMuHigh,   "fTwoGammaToMuHigh"  );
  // l->AddEntry(  fTwoGammaToMuMedium, "fTwoGammaToMuMedium");
  // l->Draw("same");

  // gPad->SaveAs("pngResults/Coh0N0Nleg.png",      "RECREATE");







  // ------------------------------------------------------------------------------------------------------------------------------------
  // Pt plot
  //---Create pt canvas
  TCanvas *cPt = new TCanvas("cPt","cPt",800,800);
  TPad *cPt1 = new TPad("cPt1", "cPt1",0.001,0.001,0.999,0.999);
  // ---Pad option
  // cPt1->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  cPt1->SetLeftMargin(0.14);
  cPt1->SetRightMargin(0.01);
  cPt1->SetBottomMargin(0.12);
  cPt1->Draw();
  cPt1->cd();

  fThetaMC->SetTitle("");
  // fThetaMC->GetXaxis()->SetRangeUser(2.96, 3.12);
  fThetaMC->GetXaxis()->SetRangeUser(2.976, 3.1051);
  fThetaMC->GetYaxis()->SetRangeUser(0, 1.6*fThetaMC->GetMaximum());
  fThetaMC->GetXaxis()->SetTitle("Single muon #theta [rad]");
  fThetaMC->GetYaxis()->SetTitleOffset(2.2);
  fThetaMC->GetXaxis()->SetTitleOffset(1.1);
  // fThetaMC->GetXaxis()->SetTitleSize(42);
  // fThetaMC->SetTitleSize(22, "x");
  fThetaMC->GetYaxis()->SetTitle("Counts/Total number of counts [a.u.]");
  fThetaMC->Draw();
  fTheta  ->Draw("same");



  // TLatex * text_pt1 = new TLatex (0.1,0.7,"Preliminary");
  // text_pt1->SetTextSize(0.035);
  // text_pt1->Draw();
  // TLatex * text_pt2 = new TLatex (0.85,0.7,"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
  // text_pt2->SetTextSize(0.035);
  // text_pt2->Draw();
  // TLatex * text_pt3 = new TLatex (0.1,0.7,"0N0N, -4.0 < y < -2.5");
  // text_pt3->SetTextSize(0.035);
  // text_pt3->Draw();
  //
  // TLegend *leg_pt = new TLegend(0.6,0.45,0.95,0.79);
  // leg_pt->SetFillStyle(0);
  // leg_pt->SetBorderSize(0);
  // leg_pt->SetTextSize(0.04);
  // leg_pt->AddEntry("fThetaMC","ALICE Simulation", "L");
  // leg_pt->AddEntry("fTheta"  ,  "ALICE Data",     "L");
  // leg_pt->Draw();


  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, single muons");
  latex->DrawLatex(0.55,0.78,"2.9 < M_{#mu#mu} < 3.2 GeV/#it{c^{2}}");
  latex->DrawLatex(0.55,0.72,"#it{p}_{T}^{#mu#mu} < 0.2 GeV/#it{c}");
  latex->DrawLatex(0.55,0.66,Form("%.1f < y < %.1f",-4.0,-2.5));

  TLegend* l = new TLegend(0.15,0.45,0.55,0.65);
  // TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
  l->SetMargin(0.1);
  l->SetBorderSize(0);
  l->AddEntry(  fThetaMC, "ALICE Simulation");
  l->AddEntry(  fTheta  , "ALICE Data"      );
  l->Draw();



  gPad->SaveAs("pngResults/ThetaSingleMuon.png",  "RECREATE");

















  // ------------------------------------------------------------------------------------------------------------------------------------
  // Pt plot
  //---Create pt canvas
  TCanvas *c2 = new TCanvas("cPt2","cPt2",800,800);
  TPad *cPt2 = new TPad("cPt2", "cPt2",0.001,0.001,0.999,0.999);
  // ---Pad option
  // cPt1->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  cPt2->SetLeftMargin(0.14);
  cPt2->SetRightMargin(0.01);
  cPt2->SetBottomMargin(0.12);
  cPt2->Draw();
  cPt2->cd();

  fPhiMC->SetTitle("");
  fPhiMC->GetXaxis()->SetRangeUser(0, 3.14*2.);
  fPhiMC->GetYaxis()->SetRangeUser(0, 2.*fPhiMC->GetMaximum());
  fPhiMC->GetXaxis()->SetTitle("Single muon #varphi [rad]");
  fPhiMC->GetYaxis()->SetTitleOffset(1.6);
  fPhiMC->GetXaxis()->SetTitleOffset(1.1);
  fPhiMC->GetYaxis()->SetTitle("Counts/Total number of counts [a.u.]");
  fPhiMC->Draw();
  fPhi  ->Draw("same");



  // TLatex * text_pt1 = new TLatex (0.1,0.7,"Preliminary");
  // text_pt1->SetTextSize(0.035);
  // text_pt1->Draw();
  // TLatex * text_pt2 = new TLatex (0.85,0.7,"#bf{ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}");
  // text_pt2->SetTextSize(0.035);
  // text_pt2->Draw();
  // TLatex * text_pt3 = new TLatex (0.1,0.7,"0N0N, -4.0 < y < -2.5");
  // text_pt3->SetTextSize(0.035);
  // text_pt3->Draw();
  //
  // TLegend *leg_pt = new TLegend(0.6,0.45,0.95,0.79);
  // leg_pt->SetFillStyle(0);
  // leg_pt->SetBorderSize(0);
  // leg_pt->SetTextSize(0.04);
  // leg_pt->AddEntry("fThetaMC","ALICE Simulation", "L");
  // leg_pt->AddEntry("fTheta"  ,  "ALICE Data",     "L");
  // leg_pt->Draw();


  TLatex* latex2 = new TLatex();
  latex2->SetTextSize(0.05);
  latex2->SetTextFont(42);
  latex2->SetTextAlign(11);
  latex2->SetNDC();
  latex2->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2->SetTextSize(0.045);
  // latex2->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex2->DrawLatex(0.55,0.84,"UPC, single muons");
  latex2->DrawLatex(0.55,0.78,"2.9 < M_{#mu#mu} < 3.2 GeV/#it{c^{2}}");
  latex2->DrawLatex(0.55,0.72,"#it{p}_{T}^{#mu#mu} < 0.2 GeV/#it{c}");
  latex2->DrawLatex(0.55,0.66,Form("%.1f < y < %.1f",-4.0,-2.5));

  TLegend* l2 = new TLegend(0.15,0.85,0.45,0.75);
  // TLegend* l2 = new TLegend(0.45,0.55,0.98,0.85);
  l2->SetMargin(0.1);
  l2->SetBorderSize(0);
  l2->AddEntry(  fPhiMC, "ALICE Simulation");
  l2->AddEntry(  fPhi  , "ALICE Data"      );
  l2->Draw();



  gPad->SaveAs("pngResults/PhiSingleMuon.png",  "RECREATE");



















  TCanvas *c3 = new TCanvas("c3","multipads",1200,700);
  gStyle->SetOptStat(0);
  c3->Divide(2,1,0,0);
  TH2F h1("h1","test1",10,0,1,20,0,20);
  TH2F h2("h2","test2",10,0,1,20,0,100);
  TH2F h3("h3","test3",10,0,1,20,-1,1);
  TH2F h4("h4","test4",10,0,1,20,0,1000);

  c3->cd(1);
  gPad->SetTickx(2);
  gPad->SetLeftMargin(0.2);
  // gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.12);
  // fThetaMC->GetYaxis()->SetRangeUser(0, 2.*fPhiMC->GetMaximum());
  fThetaMC->GetYaxis()->SetRangeUser(0, 0.0195);
  fThetaMC->Draw("L HIST");
  fTheta  ->Draw("same");

  TLatex* latex3 = new TLatex();
  latex3->SetTextSize(0.045);
  latex3->SetTextFont(42);
  latex3->SetTextAlign(11);
  latex3->SetNDC();
  latex3->DrawLatex(0.25,0.94,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3->SetTextSize(0.036);
  // latex3->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  // latex3->DrawLatex(0.55,0.84,"UPC, single muons");
  // latex3->DrawLatex(0.55,0.78,"2.9 < M_{#mu#mu} < 3.2 GeV/#it{c^{2}}");
  // latex3->DrawLatex(0.55,0.72,"#it{p}_{T}^{#mu#mu} < 0.2 GeV/#it{c}");
  // latex3->DrawLatex(0.55,0.66,Form("%.1f < y < %.1f",-4.0,-2.5));
  latex3->DrawLatex(0.25,0.84,"UPC, single muons");
  latex3->DrawLatex(0.25,0.78,"2.9 < M_{#mu#mu} < 3.2 GeV/#it{c^{2}}");
  latex3->DrawLatex(0.25,0.72,"#it{p}_{T}^{#mu#mu} < 0.2 GeV/#it{c}");
  latex3->DrawLatex(0.25,0.66,Form("%.1f < y < %.1f",-4.0,-2.5));

  // TLegend* l3 = new TLegend(0.27,0.55,0.5,0.75);
  // // TLegend* l3 = new TLegend(0.45,0.55,0.98,0.85);
  // l3->SetMargin(0.1);
  // l3->SetBorderSize(0);
  // l3->AddEntry(  fThetaMC, "ALICE Simulation");
  // l3->AddEntry(  fTheta  , "ALICE Data"      );
  // l3->Draw();





  c3->cd(2);
  gPad->SetTickx(2);
  // gPad->SetTicky(2);
  // gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  // gPad->SetBottomMargin(0.12);
  fPhiMC->GetYaxis()->SetRangeUser(0, 0.0195);
  fPhiMC->GetXaxis()->SetRangeUser(-0.05, 3.14*2.);
  fPhiMC->GetYaxis()->SetLabelOffset(0.01);
  // fPhiMC->GetXaxis()->SetTitleSize(42);
  fPhiMC->SetTitleSize(1.1*fThetaMC->GetXaxis()->GetTitleSize(), "x");
  fPhiMC->GetXaxis()->SetTitleOffset(0.95);
  fPhiMC->Draw("L HIST");
  fPhi  ->Draw("same");



  // TLegend* l4 = new TLegend(0.27,0.55,0.75,0.75);
  TLegend* l4 = new TLegend(0.2,0.15,0.75,0.28);
  // TLegend* l4 = new TLegend(0.15,0.85,0.45,0.75);
  // TLegend* l4 = new TLegend(0.45,0.55,0.98,0.85);
  gStyle->SetLegendFont(42);
  l4->SetMargin(0.1);
  l4->SetBorderSize(0);
  l4->AddEntry(  fPhiMC, "ALICE Simulation");
  l4->AddEntry(  fPhi  , "ALICE Data"      );
  l4->Draw();

  c3->SaveAs("pngResults/ThetaPhiSingleMuon.png",  "RECREATE");
  c3->SaveAs("pngResults/ThetaPhiSingleMuon.pdf",  "RECREATE");


}
