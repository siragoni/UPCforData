#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <vector>



//_____________________________________________________________________________
/* - Histograms to be used for the fit.
 * - What happens is that we will interpolate the many points together...
 * -
 */
TH1F* fCohJpsiToMu;
TH1F* fCohPsi2sToMu;
TH1F* fCohPsi2sToMuPi;
TH1F* fIncohJpsiToMu;
TH1F* fIncohPsi2sToMu;
TH1F* fIncohPsi2sToMuPi;
TH1F* fTwoGammaToMuMedium;
TH1F* fTwoGammaToMuHigh;
TH1F* fHighPtTail;
//_____________________________________________________________________________
/* - Fit function for the Pt plots.
 * - I am using simple ROOT to make gaussian fits to the plot.
 */
Double_t fPtDistr(Double_t* x,Double_t* par)
{
  /* - Par 0, 1, 2:   coherent.
     - Par 3, 4, 5:   incoherent.
     - Par 6      :   gamma+gamma.
     -
   */
  Double_t val = 0;
  val += par[0]* ( fCohJpsiToMu       ->Interpolate(x[0]) );
  // val += par[1]* ( fCohPsi2sToMu      ->Interpolate(x[0]) );
  val += par[2]* ( fCohPsi2sToMuPi    ->Interpolate(x[0]) );
  val += par[3]* ( fIncohJpsiToMu     ->Interpolate(x[0]) );
  // val += par[4]* ( fIncohPsi2sToMu    ->Interpolate(x[0]) );
  val += par[5]* ( fIncohPsi2sToMuPi  ->Interpolate(x[0]) );
  val += par[6]* ( fTwoGammaToMuMedium->Interpolate(x[0]) );
  val += par[7]* ( fHighPtTail        ->Interpolate(x[0]) );

  return val;
}
//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void fitPtDistr(const char* AnalysisName){
  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TFile* fileMC[8];
  fileMC[0] = new TFile("MCtrainResults/2019-05-17/kCohJpsiToMu/AnalysisResults.root");
  fileMC[1] = new TFile("MCtrainResults/2019-05-17/kCohPsi2sToMu/AnalysisResults.root");
  fileMC[2] = new TFile("MCtrainResults/2019-05-17/kCohPsi2sToMuPi/AnalysisResults.root");
  fileMC[3] = new TFile("MCtrainResults/2019-05-17/kIncohJpsiToMu/AnalysisResults.root");
  fileMC[4] = new TFile("MCtrainResults/2019-05-17/kIncohPsi2sToMu/AnalysisResults.root");
  fileMC[5] = new TFile("MCtrainResults/2019-05-17/kIncohPsi2sToMuPi/AnalysisResults.root");
  fileMC[6] = new TFile("MCtrainResults/2019-05-17/kTwoGammaToMuHigh/AnalysisResults.root");
  fileMC[7] = new TFile("MCtrainResults/2019-05-17/kTwoGammaToMuMedium/AnalysisResults.root");
  TDirectory* dirMC[8];
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory] = fileMC[iDirectory]->GetDirectory("MyTask");
  }
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);
  TList* listingsMC[8];
  for(Int_t iDirectory = 0; iDirectory < 8; iDirectory++) {
    dirMC[iDirectory]->GetObject("MyOutputContainer", listingsMC[iDirectory]);
  }
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
  fCohJpsiToMu        = (TH1F*)listingsMC[0]->FindObject("fDimuonPtDistributionH");
  fCohPsi2sToMu       = (TH1F*)listingsMC[1]->FindObject("fDimuonPtDistributionH");
  fCohPsi2sToMuPi     = (TH1F*)listingsMC[2]->FindObject("fDimuonPtDistributionH");
  fIncohJpsiToMu      = (TH1F*)listingsMC[3]->FindObject("fDimuonPtDistributionH");
  fIncohPsi2sToMu     = (TH1F*)listingsMC[4]->FindObject("fDimuonPtDistributionH");
  fIncohPsi2sToMuPi   = (TH1F*)listingsMC[5]->FindObject("fDimuonPtDistributionH");
  fTwoGammaToMuMedium = (TH1F*)listingsMC[6]->FindObject("fDimuonPtDistributionH");
  fTwoGammaToMuHigh   = (TH1F*)listingsMC[7]->FindObject("fDimuonPtDistributionH");
  /* - Rebin
     -
   */
  fCohJpsiToMu        -> Rebin(4);
  fCohPsi2sToMu       -> Rebin(4);
  fCohPsi2sToMuPi     -> Rebin(4);
  fIncohJpsiToMu      -> Rebin(4);
  fIncohPsi2sToMu     -> Rebin(4);
  fIncohPsi2sToMuPi   -> Rebin(4);
  fTwoGammaToMuMedium -> Rebin(4);
  fTwoGammaToMuHigh   -> Rebin(4);
  /* - Firstly we normalize the histograms.
     - Remember to always Sumw2()!!
     -
   */
  fCohJpsiToMu        -> Sumw2();
  fCohPsi2sToMu       -> Sumw2();
  fCohPsi2sToMuPi     -> Sumw2();
  fIncohJpsiToMu      -> Sumw2();
  fIncohPsi2sToMu     -> Sumw2();
  fIncohPsi2sToMuPi   -> Sumw2();
  fTwoGammaToMuMedium -> Sumw2();
  fTwoGammaToMuHigh   -> Sumw2();
  // Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral(0, 20);
  // Double_t Integral_fCohPsi2sToMu       = fCohPsi2sToMu       -> Integral(0, 20);
  // Double_t Integral_fCohPsi2sToMuPi     = fCohPsi2sToMuPi     -> Integral(0, 20);
  // Double_t Integral_fIncohJpsiToMu      = fIncohJpsiToMu      -> Integral(0, 20);
  // Double_t Integral_fIncohPsi2sToMu     = fIncohPsi2sToMu     -> Integral(0, 20);
  // Double_t Integral_fIncohPsi2sToMuPi   = fIncohPsi2sToMuPi   -> Integral(0, 20);
  // Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral(0, 20);
  // Double_t Integral_fTwoGammaToMuHigh   = fTwoGammaToMuHigh   -> Integral(0, 20);
  Double_t Integral_fCohJpsiToMu        = fCohJpsiToMu        -> Integral();
  Double_t Integral_fCohPsi2sToMu       = fCohPsi2sToMu       -> Integral();
  Double_t Integral_fCohPsi2sToMuPi     = fCohPsi2sToMuPi     -> Integral();
  Double_t Integral_fIncohJpsiToMu      = fIncohJpsiToMu      -> Integral();
  Double_t Integral_fIncohPsi2sToMu     = fIncohPsi2sToMu     -> Integral();
  Double_t Integral_fIncohPsi2sToMuPi   = fIncohPsi2sToMuPi   -> Integral();
  Double_t Integral_fTwoGammaToMuMedium = fTwoGammaToMuMedium -> Integral();
  Double_t Integral_fTwoGammaToMuHigh   = fTwoGammaToMuHigh   -> Integral();
  fCohJpsiToMu        -> Scale( 1/Integral_fCohJpsiToMu        );
  fCohPsi2sToMu       -> Scale( 1/Integral_fCohPsi2sToMu       );
  fCohPsi2sToMuPi     -> Scale( 1/Integral_fCohPsi2sToMuPi     );
  fIncohJpsiToMu      -> Scale( 1/Integral_fIncohJpsiToMu      );
  fIncohPsi2sToMu     -> Scale( 1/Integral_fIncohPsi2sToMu     );
  fIncohPsi2sToMuPi   -> Scale( 1/Integral_fIncohPsi2sToMuPi   );
  fTwoGammaToMuMedium -> Scale( 1/Integral_fTwoGammaToMuMedium );
  fTwoGammaToMuHigh   -> Scale( 1/Integral_fTwoGammaToMuHigh   );

  /* - High Pt-tail, with HERA's data.
     -
   */
  TF1* fModelForHighPtTail = new TF1("fModelForHighPtTail","[0]*x*(1+[1]/[2]*x*x)^(-[2])",0,4);
  fModelForHighPtTail->SetParameter(0,1);
//  fModelForHighPtTail->SetParameter(1,debug==4 ? 1.25 : 1.);
//  fModelForHighPtTail->SetParameter(2,debug==4 ? 6.1 : 1.);
  fModelForHighPtTail->SetParameter(1, 1.6/*1.79*/);
  fModelForHighPtTail->SetParameter(2, 3.58);
  fModelForHighPtTail->SetNpx( fCohJpsiToMu->GetNbinsX() );
  fHighPtTail = (TH1F*) fModelForHighPtTail->GetHistogram()->Clone("fHighPtTail");
  for (Int_t ibin=1; ibin<=fHighPtTail->GetNbinsX(); ibin++) {
    fHighPtTail->SetBinError(ibin,0);
  }



  TH1F *fDimuonPtDistributionDataH = (TH1F*)listings->FindObject("fDimuonPtDistributionH");
  fDimuonPtDistributionDataH->Rebin(4);
  fDimuonPtDistributionDataH->Draw("PE");


  TF1* FitPtDistr = new TF1(  "FitPtDistr",
                              fPtDistr,
                              0, 3,
                              /*7*/8
                              );
  FitPtDistr->SetNpx(1000);

  /* - Setting parameters for the fit.
     -
   */
  // Coherent J/Psi
  // FitPtDistr->SetParameter(0, 38000);
  // FitPtDistr->SetParameter(0, 20144);
  // FitPtDistr->SetParLimits(0, FitPtDistr->GetParameter(0)*0.95, FitPtDistr->GetParameter(0)*1.05);

  // FEED-DOWN COHERENT
  // Double_t kFeedDownCoherent                       = 1;     // neutral element
  // Double_t kCrossSectionCoherentPsiPrimeToJPsiPiPi = 0.818; // millibarn
  // Double_t kCrossSectionCoherentJPsiDirect         = 4.370; // millibarn
  // Double_t AccTimesEffCoherentPsiPrimeToJPsiPiPi   = 0.15840;
  // Double_t AccTimesEffCoherentJPsiDirect           = 0.12251;
  // Double_t kBrachingRationCoherentPsiPrime         = 0.609;
  // kFeedDownCoherent *= kCrossSectionCoherentPsiPrimeToJPsiPiPi;
  // kFeedDownCoherent *= kBrachingRationCoherentPsiPrime;
  // kFeedDownCoherent *= AccTimesEffCoherentPsiPrimeToJPsiPiPi;
  // kFeedDownCoherent =  kFeedDownCoherent / (kCrossSectionCoherentJPsiDirect*AccTimesEffCoherentJPsiDirect);
  // FitPtDistr->SetParameter(2, kFeedDownCoherent);
  // FitPtDistr->SetParLimits(2, FitPtDistr->GetParameter(2)*0.95, FitPtDistr->GetParameter(2)*1.05);
  FitPtDistr->FixParameter(1, 0);
  FitPtDistr->FixParameter(4, 0);

  // FitPtDistr->SetParameter(0, 4.7610e+04);
  // FitPtDistr->SetParLimits(0, FitPtDistr->GetParameter(0)*0.8, FitPtDistr->GetParameter(0)*1.2);


  Double_t kFeedDownCoherent = 1705;     // neutral element
  Double_t kError            = 0.175;
  FitPtDistr->SetParameter(2, kFeedDownCoherent);
  FitPtDistr->SetParLimits(2, FitPtDistr->GetParameter(2)*(1-kError), FitPtDistr->GetParameter(2)*(1+kError));

  // FEED-DOWN COHERENT
  Double_t kFeedDownIncoherent = 1705;     // neutral element
  Double_t kErrorInc           = 0.175;
  FitPtDistr->SetParameter(5, kFeedDownIncoherent);
  FitPtDistr->SetParLimits(5, FitPtDistr->GetParameter(5)*(1-kErrorInc), FitPtDistr->GetParameter(5)*(1+kErrorInc));

  // Gamma+Gamma Medium
  // FitPtDistr->SetParameter(6, 4800);
  FitPtDistr->SetParameter(6, 6531);
  // FitPtDistr->FixParameter(6, 0);
  FitPtDistr->SetParLimits(6, FitPtDistr->GetParameter(6)*0.95, FitPtDistr->GetParameter(6)*1.05);

  // Unknown component
  FitPtDistr->SetParameter(7, 743);
  FitPtDistr->SetParLimits(7, 720, 745);
  // FitPtDistr->SetParLimits(7, 740, 745);

  fDimuonPtDistributionDataH->SetLineColor(kBlue);
  fDimuonPtDistributionDataH->SetLineStyle(kSolid);
  fDimuonPtDistributionDataH->SetLineWidth(3);
  fDimuonPtDistributionDataH->SetMarkerStyle(kFullCircle);
  fDimuonPtDistributionDataH->SetMarkerColor(kBlue);
  fDimuonPtDistributionDataH->SetMarkerSize(1);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  fDimuonPtDistributionDataH->GetYaxis()->SetTitle( Form( "Counts / (%.3f GeV/#it{c})",
                                                          fDimuonPtDistributionDataH->GetXaxis()->GetBinWidth(1)
                                                        )
                                                    );
  fDimuonPtDistributionDataH->SetTitle("");
  fDimuonPtDistributionDataH->Fit( FitPtDistr,"","", 0., 3. );
  TCanvas* PtDistrCanvas = new TCanvas( "PtDistrCanvas", "PtDistrCanvas", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  gPad->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  /* - Beautifying is starting now.
     -
   */
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleOffset(1.25);
  // fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.25);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleOffset(1.45);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleSize(0.045);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleSize(0.045);
  fDimuonPtDistributionDataH->GetXaxis()->SetLabelSize(0.045);
  fDimuonPtDistributionDataH->GetYaxis()->SetLabelSize(0.045);
  fDimuonPtDistributionDataH->GetXaxis()->SetTitleFont(42);
  fDimuonPtDistributionDataH->GetYaxis()->SetTitleFont(42);
  fDimuonPtDistributionDataH->GetXaxis()->SetLabelFont(42);
  fDimuonPtDistributionDataH->GetYaxis()->SetLabelFont(42);
  // fDimuonPtDistributionDataH->GetXaxis()->SetNdivisions(408);
  fDimuonPtDistributionDataH->GetYaxis()->SetRangeUser(10, fDimuonPtDistributionDataH->GetMaximum()*10.);
  fDimuonPtDistributionDataH->GetXaxis()->SetRangeUser(0, 3);
  gPad ->SetLogy();
  fDimuonPtDistributionDataH->Draw("PEsame");
  fCohJpsiToMu        -> SetLineColor(kRed);
  fCohPsi2sToMu       -> SetLineColor(kMagenta);
  fCohPsi2sToMuPi     -> SetLineColor(kYellow+1);
  fIncohJpsiToMu      -> SetLineColor(kCyan);
  fIncohPsi2sToMu     -> SetLineColor(kYellow);
  fIncohPsi2sToMuPi   -> SetLineColor(kBlue+2);
  fTwoGammaToMuMedium -> SetLineColor(kGreen);
  fTwoGammaToMuHigh   -> SetLineColor(kBlue+3);
  fHighPtTail         -> SetLineColor(kGreen+1);
  fCohJpsiToMu        -> SetLineWidth(3);
  fCohPsi2sToMu       -> SetLineWidth(3);
  fCohPsi2sToMuPi     -> SetLineWidth(3);
  fIncohJpsiToMu      -> SetLineWidth(3);
  fIncohPsi2sToMu     -> SetLineWidth(3);
  fIncohPsi2sToMuPi   -> SetLineWidth(3);
  fTwoGammaToMuMedium -> SetLineWidth(3);
  fTwoGammaToMuHigh   -> SetLineWidth(3);
  fHighPtTail         -> SetLineWidth(3);
  TH1F* fCohJpsiToMuC        = (TH1F*) fCohJpsiToMu        -> Clone("fCohJpsiToMuC");
  TH1F* fCohPsi2sToMuC       = (TH1F*) fCohPsi2sToMu       -> Clone("fCohPsi2sToMuC");
  TH1F* fCohPsi2sToMuPiC     = (TH1F*) fCohPsi2sToMuPi     -> Clone("fCohPsi2sToMuPiC");
  TH1F* fIncohJpsiToMuC      = (TH1F*) fIncohJpsiToMu      -> Clone("fIncohJpsiToMuC");
  TH1F* fIncohPsi2sToMuC     = (TH1F*) fIncohPsi2sToMu     -> Clone("fIncohPsi2sToMuC");
  TH1F* fIncohPsi2sToMuPiC   = (TH1F*) fIncohPsi2sToMuPi   -> Clone("fIncohPsi2sToMuPiC");
  TH1F* fTwoGammaToMuMediumC = (TH1F*) fTwoGammaToMuMedium -> Clone("fTwoGammaToMuMediumC");
  TH1F* fTwoGammaToMuHighC   = (TH1F*) fTwoGammaToMuHigh   -> Clone("fTwoGammaToMuHighC");
  TH1F* fHighPtTailC         = (TH1F*) fHighPtTail         -> Clone("fHighPtTailC");
  fCohJpsiToMuC        -> Scale( FitPtDistr->GetParameter(0) );
  fCohPsi2sToMuC       -> Scale( FitPtDistr->GetParameter(1) );
  fCohPsi2sToMuPiC     -> Scale( FitPtDistr->GetParameter(2) );
  fIncohJpsiToMuC      -> Scale( FitPtDistr->GetParameter(3) );
  fIncohPsi2sToMuC     -> Scale( FitPtDistr->GetParameter(4) );
  fIncohPsi2sToMuPiC   -> Scale( FitPtDistr->GetParameter(5) );
  fTwoGammaToMuMediumC -> Scale( FitPtDistr->GetParameter(6) );
  fTwoGammaToMuHighC   -> Scale( FitPtDistr->GetParameter(6) );
  fHighPtTailC         -> Scale( FitPtDistr->GetParameter(7) );
  fCohJpsiToMuC        -> Draw("same");
  // fCohPsi2sToMuC       -> Draw("same");
  fCohPsi2sToMuPiC     -> Draw("same");
  fIncohJpsiToMuC      -> Draw("same");
  // fIncohPsi2sToMuC     -> Draw("same");
  fIncohPsi2sToMuPiC   -> Draw("same");
  fTwoGammaToMuMediumC -> Draw("same");
  // fTwoGammaToMuHighC   -> Draw("same");
  fHighPtTailC         -> Draw("Esame");
  // fCohJpsiToMu        -> Scale( FitPtDistr->GetParameter(0) );
  // fCohPsi2sToMu       -> Scale( FitPtDistr->GetParameter(1) );
  // fCohPsi2sToMuPi     -> Scale( FitPtDistr->GetParameter(2) );
  // fIncohJpsiToMu      -> Scale( FitPtDistr->GetParameter(3) );
  // fIncohPsi2sToMu     -> Scale( FitPtDistr->GetParameter(4) );
  // fIncohPsi2sToMuPi   -> Scale( FitPtDistr->GetParameter(5) );
  // fTwoGammaToMuMedium -> Scale( FitPtDistr->GetParameter(6) );
  // fTwoGammaToMuHigh   -> Scale( FitPtDistr->GetParameter(0) );
  // fCohJpsiToMu        -> Draw("same");
  // fCohPsi2sToMu       -> Draw("same");
  // fCohPsi2sToMuPi     -> Draw("same");
  // fIncohJpsiToMu      -> Draw("same");
  // fIncohPsi2sToMu     -> Draw("same");
  // fIncohPsi2sToMuPi   -> Draw("same");
  // fTwoGammaToMuMedium -> Draw("same");
  // fTwoGammaToMuHigh   -> Draw("same");


  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  // latex->DrawLatex(0.55,0.84,"UPC, LHC18q+LHC18r data");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  // latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));
  /* - Chi square computation.
     -
   */
  // latex->DrawLatex(0.55,0.68,Form( "  #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
  //                                  FitPtDistr->GetChisquare(),
  //                                  FitPtDistr->GetNDF(),
  //                                  FitPtDistr->GetChisquare()/FitPtDistr->GetNDF()
  //                                 )
  //                                );

  TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
  l->SetMargin(0.1);
  l->SetBorderSize(0);
  l->AddEntry(  fDimuonPtDistributionDataH, "ALICE data 2018");
  l->AddEntry(  FitPtDistr, Form( "Fit: #chi^{2}/NDF = %.2f / %.2d = %.2f  ",
                                   FitPtDistr->GetChisquare(),
                                   FitPtDistr->GetNDF(),
                                   FitPtDistr->GetChisquare()/FitPtDistr->GetNDF()
                                   )
                                  );
  l->AddEntry(  fCohJpsiToMuC,        "Coherent   J/#psi");
  l->AddEntry(  fIncohJpsiToMuC,      "Incoherent J/#psi");
  l->AddEntry(  fHighPtTailC,         "Incoherent dissociative J/#psi");
  l->AddEntry(  fCohPsi2sToMuPiC,     "Coherent   #psi(2S) feeddown");
  l->AddEntry(  fIncohPsi2sToMuPiC,   "Incoherent #psi(2S) feeddown");
  l->AddEntry(  fTwoGammaToMuMediumC, "Continuum  #gamma#gamma #rightarrow #mu#mu");
  l->Draw();

  gPad->SaveAs("pngResults/fitPtDistr.png", "RECREATE");

}
