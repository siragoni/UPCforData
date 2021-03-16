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

//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void PolarisationAxEratio(){

  TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);

  TFile* fileList2 = new TFile("AnalysisResultsLHC18l7_coh_nophotonsPol.root");
  TDirectory* dir2 = fileList2->GetDirectory("MyTask");
  TList* listings2;
  dir2->GetObject("MyOutputContainer", listings2);

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
  TH1F* fReconCosThetaH  = (TH1F*)listings->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH  = (TH1F*)listings->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fReconCosThetaH2 = (TH1F*)listings2->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH2 = (TH1F*)listings2->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");

  // fReconCosThetaH->Rebin(4);
  // fGenerCosThetaH->Rebin(4);
  fReconCosThetaH ->Sumw2();
  fGenerCosThetaH ->Sumw2();
  fReconCosThetaH2->Sumw2();
  fGenerCosThetaH2->Sumw2();
  // new TCanvas;
  // fGenerCosThetaH->Draw();
  TH1F* ReconTheta  = (TH1F*) fReconCosThetaH->Clone("ReconTheta");
  TH1F* ReconTheta2 = (TH1F*) fReconCosThetaH2->Clone("ReconTheta");

  TH1F* fReconPhiH = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerPhiH = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsH");
  fReconPhiH->Sumw2();
  fGenerPhiH->Sumw2();

  TH1F* fReconTildePhiH = (TH1F*)listings->FindObject("fTildePhiHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerTildePhiH = (TH1F*)listings->FindObject("fMCTildePhiHelicityFrameTwentyfiveBinsH");
  fReconTildePhiH->Sumw2();
  fGenerTildePhiH->Sumw2();





  TCanvas* AcceptanceCanvasCosTheta = new TCanvas("AcceptanceCanvasCosTheta","AcceptanceCanvasCosTheta",900,800);
  TH1F* acceptanceCosTheta = (TH1F*) fReconCosThetaH->Clone("acceptanceCosTheta");
  acceptanceCosTheta->Divide(fGenerCosThetaH);
  acceptanceCosTheta->SetLineWidth(2);
  acceptanceCosTheta->SetLineColor(kRed);
  acceptanceCosTheta->Draw("ep");

  TH1F* acceptanceCosTheta2 = (TH1F*) fReconCosThetaH2->Clone("acceptanceCosTheta2");
  acceptanceCosTheta2->Divide(fGenerCosThetaH2);
  acceptanceCosTheta2->SetLineWidth(2);
  acceptanceCosTheta2->SetLineColor(kMagenta);
  acceptanceCosTheta2->Draw("epsame");

  TCanvas* RatioWithWithoutPhotons = new TCanvas("RatioWithWithoutPhotons","RatioWithWithoutPhotons",900,800);
  TH1F* acceptanceCosTheta2a  = (TH1F*) acceptanceCosTheta2->Clone("acceptanceCosTheta2a");
  acceptanceCosTheta2a->Divide(acceptanceCosTheta);
  acceptanceCosTheta2a->SetLineWidth(2);
  acceptanceCosTheta2a->SetLineColor(kMagenta);
  acceptanceCosTheta2a->GetYaxis()->SetRangeUser(0.935, 1.03);
  acceptanceCosTheta2a->Draw("ep");

  // TCanvas* AcceptanceCanvasPhi = new TCanvas("AcceptanceCanvasPhi","AcceptanceCanvasPhi",900,800);
  // TH1F* acceptancePhi = (TH1F*) fReconPhiH->Clone("acceptancePhi");
  // acceptancePhi->Divide(fGenerPhiH);
  // acceptancePhi->Draw("ep");
  //
  // TCanvas* CorrCanvasPhi = new TCanvas("CorrCanvasPhi","CorrCanvasPhi",900,800);
  // RawPhiH->Divide(acceptancePhi);
  // RawPhiH->Draw("ep");
  //
  // TCanvas* AcceptanceCanvasTildePhi = new TCanvas("AcceptanceCanvasTildePhi","AcceptanceCanvasTildePhi",900,800);
  // TH1F* acceptanceTildePhi = (TH1F*) fReconTildePhiH->Clone("acceptanceTildePhi");
  // acceptanceTildePhi->Divide(fGenerTildePhiH);
  // acceptanceTildePhi->Draw("ep");
  //
  // TCanvas* CorrCanvasTildePhi = new TCanvas("CorrCanvasTildePhi","CorrCanvasTildePhi",900,800);
  // RawTildePhiH->Divide(acceptanceTildePhi);
  // RawTildePhiH->Draw("ep");



}
