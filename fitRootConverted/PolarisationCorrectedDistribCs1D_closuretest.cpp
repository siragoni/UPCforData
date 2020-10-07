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
void PolarisationCorrectedDistribCs1D( Int_t RangeSelectionMode = 0 ){

  TFile* fileListFullT   = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root"); // reweighting
  TFile* fileListShortT  = new TFile("AnalysisResultsLHC18l7Michal.root"); // reweighting
  TFile* fileListShortL  = new TFile("AnalysisResultsMichal.root"); // trial longitudinal
  TDirectory* dirFullT   = fileListFullT ->GetDirectory("MyTask");
  TDirectory* dirShortT  = fileListShortT->GetDirectory("MyTask");
  TDirectory* dirShortL  = fileListShortL->GetDirectory("MyTask");
  TList* listingsFullT;
  TList* listingsShortT;
  TList* listingsShortL;
  dirFullT ->GetObject("MyOutputContainer", listingsFullT );
  dirShortT->GetObject("MyOutputContainer", listingsShortT);
  dirShortL->GetObject("MyOutputContainer", listingsShortL);
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

  //_____________________________________________________________________________
  // Full Transverse
  //
  TH1F* fReconCosThetaFullTH = (TH1F*)listingsFullT->FindObject("fCosThetaCsFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaFullTH = (TH1F*)listingsFullT->FindObject("fMCCosThetaCsFrameTwentyfiveBinsH");
  fReconCosThetaFullTH->Sumw2();
  fGenerCosThetaFullTH->Sumw2();
  TH1F* ReconThetaFullT = (TH1F*) fReconCosThetaFullTH->Clone("ReconThetaFullT");
  TH1F* DataThetaFullT0 = (TH1F*) fReconCosThetaFullTH->Clone("DataThetaFullT0");
  TH1F* DataThetaFullT1 = (TH1F*) fReconCosThetaFullTH->Clone("DataThetaFullT1");
  TH1F* DataThetaFullT2 = (TH1F*) fReconCosThetaFullTH->Clone("DataThetaFullT2");

  TH1F* fReconPhiFullTH = (TH1F*)listingsFullT->FindObject("fPhiCsFrameTwentyfiveBinsH");
  TH1F* fGenerPhiFullTH = (TH1F*)listingsFullT->FindObject("fMCPhiCsFrameTwentyfiveBinsH");
  fReconPhiFullTH->Sumw2();
  fGenerPhiFullTH->Sumw2();
  TH1F* ReconPhiFullT = (TH1F*) fReconPhiFullTH->Clone("ReconPhiFullT");
  TH1F* DataPhiFullT0 = (TH1F*) fReconPhiFullTH->Clone("DataPhiFullT0");
  TH1F* DataPhiFullT1 = (TH1F*) fReconPhiFullTH->Clone("DataPhiFullT1");
  TH1F* DataPhiFullT2 = (TH1F*) fReconPhiFullTH->Clone("DataPhiFullT2");

  TH1F* fReconTildePhiFullTH = (TH1F*)listingsFullT->FindObject("fTildePhiCsFrameTwentyfiveBinsH");
  TH1F* fGenerTildePhiFullTH = (TH1F*)listingsFullT->FindObject("fMCTildePhiCsFrameTwentyfiveBinsH");
  fReconTildePhiFullTH->Sumw2();
  fGenerTildePhiFullTH->Sumw2();
  TH1F* ReconTildePhiFullT = (TH1F*) fReconTildePhiFullTH->Clone("ReconTildePhiFullT");
  TH1F* DataTildePhiFullT0 = (TH1F*) fReconTildePhiFullTH->Clone("DataTildePhiFullT0");
  TH1F* DataTildePhiFullT1 = (TH1F*) fReconTildePhiFullTH->Clone("DataTildePhiFullT1");
  TH1F* DataTildePhiFullT2 = (TH1F*) fReconTildePhiFullTH->Clone("DataTildePhiFullT2");





  //_____________________________________________________________________________
  // Short Transverse
  //
  TH1F* fReconCosThetaShortTH = (TH1F*)listingsShortT->FindObject("fCosThetaCsFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaShortTH = (TH1F*)listingsShortT->FindObject("fMCCosThetaCsFrameTwentyfiveBinsH");
  fReconCosThetaShortTH->Sumw2();
  fGenerCosThetaShortTH->Sumw2();
  TH1F* ReconThetaShortT  = (TH1F*) fReconCosThetaShortTH->Clone("ReconThetaShortT");
  TH1F* DataThetaShortT0  = (TH1F*) fReconCosThetaShortTH->Clone("DataThetaShortT0");
  TH1F* DataThetaShortT1  = (TH1F*) fReconCosThetaShortTH->Clone("DataThetaShortT1");
  TH1F* DataThetaShortT2  = (TH1F*) fReconCosThetaShortTH->Clone("DataThetaShortT2");

  TH1F* fReconPhiShortTH = (TH1F*)listingsShortT->FindObject("fPhiCsFrameTwentyfiveBinsH");
  TH1F* fGenerPhiShortTH = (TH1F*)listingsShortT->FindObject("fMCPhiCsFrameTwentyfiveBinsH");
  fReconPhiShortTH->Sumw2();
  fGenerPhiShortTH->Sumw2();
  TH1F* ReconPhiShortT  = (TH1F*) fReconPhiShortTH->Clone("ReconPhiShortT");
  TH1F* DataPhiShortT0  = (TH1F*) fReconPhiShortTH->Clone("DataPhiShortT0");
  TH1F* DataPhiShortT1  = (TH1F*) fReconPhiShortTH->Clone("DataPhiShortT1");
  TH1F* DataPhiShortT2  = (TH1F*) fReconPhiShortTH->Clone("DataPhiShortT2");

  TH1F* fReconTildePhiShortTH = (TH1F*)listingsShortT->FindObject("fTildePhiCsFrameTwentyfiveBinsH");
  TH1F* fGenerTildePhiShortTH = (TH1F*)listingsShortT->FindObject("fMCTildePhiCsFrameTwentyfiveBinsH");
  fReconTildePhiShortTH->Sumw2();
  fGenerTildePhiShortTH->Sumw2();
  TH1F* ReconTildePhiShortT  = (TH1F*) fReconTildePhiShortTH->Clone("ReconTildePhiShortT");
  TH1F* DataTildePhiShortT0  = (TH1F*) fReconTildePhiShortTH->Clone("DataTildePhiShortT0");
  TH1F* DataTildePhiShortT1  = (TH1F*) fReconTildePhiShortTH->Clone("DataTildePhiShortT1");
  TH1F* DataTildePhiShortT2  = (TH1F*) fReconTildePhiShortTH->Clone("DataTildePhiShortT2");







  //_____________________________________________________________________________
  // Short Longitudinal
  //
  TH1F* fReconCosThetaShortLH = (TH1F*)listingsShortL->FindObject("fCosThetaCsFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaShortLH = (TH1F*)listingsShortL->FindObject("fMCCosThetaCsFrameTwentyfiveBinsH");
  fReconCosThetaShortLH->Sumw2();
  fGenerCosThetaShortLH->Sumw2();
  TH1F* ReconThetaShortL = (TH1F*) fReconCosThetaShortLH->Clone("ReconThetaShortL");
  TH1F* DataThetaShortL0 = (TH1F*) fReconCosThetaShortLH->Clone("DataThetaShortL0");
  TH1F* DataThetaShortL1 = (TH1F*) fReconCosThetaShortLH->Clone("DataThetaShortL1");
  TH1F* DataThetaShortL2 = (TH1F*) fReconCosThetaShortLH->Clone("DataThetaShortL2");

  TH1F* fReconPhiShortLH = (TH1F*)listingsShortL->FindObject("fPhiCsFrameTwentyfiveBinsH");
  TH1F* fGenerPhiShortLH = (TH1F*)listingsShortL->FindObject("fMCPhiCsFrameTwentyfiveBinsH");
  fReconPhiShortLH->Sumw2();
  fGenerPhiShortLH->Sumw2();
  TH1F* ReconPhiShortL = (TH1F*) fReconPhiShortLH->Clone("ReconPhiShortL");
  TH1F* DataPhiShortL0 = (TH1F*) fReconPhiShortLH->Clone("DataPhiShortL0");
  TH1F* DataPhiShortL1 = (TH1F*) fReconPhiShortLH->Clone("DataPhiShortL1");
  TH1F* DataPhiShortL2 = (TH1F*) fReconPhiShortLH->Clone("DataPhiShortL2");

  TH1F* fReconTildePhiShortLH = (TH1F*)listingsShortL->FindObject("fTildePhiCsFrameTwentyfiveBinsH");
  TH1F* fGenerTildePhiShortLH = (TH1F*)listingsShortL->FindObject("fMCTildePhiCsFrameTwentyfiveBinsH");
  fReconTildePhiShortLH->Sumw2();
  fGenerTildePhiShortLH->Sumw2();
  TH1F* ReconTildePhiShortL = (TH1F*) fReconTildePhiShortLH->Clone("ReconTildePhiShortL");
  TH1F* DataTildePhiShortL0 = (TH1F*) fReconTildePhiShortLH->Clone("DataTildePhiShortL0");
  TH1F* DataTildePhiShortL1 = (TH1F*) fReconTildePhiShortLH->Clone("DataTildePhiShortL1");
  TH1F* DataTildePhiShortL2 = (TH1F*) fReconTildePhiShortLH->Clone("DataTildePhiShortL2");



  // TDatime d;
  //
  // TFile* fileDataRawCosTheta = 0x0;
  // TFile* fileDataRawPhi      = 0x0;
  // TFile* fileDataRawTildePhi = 0x0;
  // fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaCS/CosThetaCSFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiCS/PhiCsFrame.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiCS/TildePhiCsFrame.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  //
  //
  // TH1F* CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  // CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  // TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");
  //
  // TH1F* PhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  // PhiAfterSignalExtractionErrorsRawH->Sumw2();
  // TH1F* RawPhiH = (TH1F*) PhiAfterSignalExtractionErrorsRawH->Clone("RawPhiH");
  //
  // TH1F* TildePhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
  // TildePhiAfterSignalExtractionErrorsRawH->Sumw2();
  // TH1F* RawTildePhiH = (TH1F*) TildePhiAfterSignalExtractionErrorsRawH->Clone("RawTildePhiH");






  //_____________________________________________________________________________
  // Full Transverse
  //
  TCanvas* AcceptanceCanvasCosThetaFullT = new TCanvas("AcceptanceCanvasCosThetaFullT","AcceptanceCanvasCosThetaFullT",900,800);
  TH1F* acceptanceCosThetaFullT = (TH1F*) fReconCosThetaFullTH->Clone("acceptanceCosThetaFullT");
  // acceptanceCosThetaFullT->Divide(fGenerCosThetaFullTH);
  acceptanceCosThetaFullT->Divide(acceptanceCosThetaFullT, fGenerCosThetaFullTH, 1, 1, "b");
  acceptanceCosThetaFullT->Draw("ep");

  TCanvas* CorrCanvasCosThetaFullT = new TCanvas("CorrCanvasCosThetaFullT","CorrCanvasCosThetaFullT",900,800);
  DataThetaFullT0->Divide(acceptanceCosThetaFullT);
  DataThetaFullT0->Draw("ep");
  DataThetaShortT0->Divide(acceptanceCosThetaFullT);
  DataThetaShortT0->Draw("ep");
  DataThetaShortL0->Divide(acceptanceCosThetaFullT);
  DataThetaShortL0->Draw("ep");

  TCanvas* AcceptanceCanvasPhiFullT = new TCanvas("AcceptanceCanvasPhiFullT","AcceptanceCanvasPhiFullT",900,800);
  TH1F* acceptancePhiFullT = (TH1F*) fReconPhiFullTH->Clone("acceptancePhiFullT");
  // acceptancePhiFullT->Divide(fGenerPhiFullTH);
  acceptancePhiFullT->Divide(acceptancePhiFullT, fGenerPhiFullTH, 1, 1, "b");
  acceptancePhiFullT->Draw("ep");

  TCanvas* CorrCanvasPhiFullT = new TCanvas("CorrCanvasPhiFullT","CorrCanvasPhiFullT",900,800);
  DataPhiFullT0->Divide(acceptancePhiFullT);
  DataPhiFullT0->Draw("ep");
  DataPhiShortT0->Divide(acceptancePhiFullT);
  DataPhiShortT0->Draw("ep");
  DataPhiShortL0->Divide(acceptancePhiFullT);
  DataPhiShortL0->Draw("ep");

  TCanvas* AcceptanceCanvasTildePhiFullT = new TCanvas("AcceptanceCanvasTildePhiFullT","AcceptanceCanvasTildePhiFullT",900,800);
  TH1F* acceptanceTildePhiFullT = (TH1F*) fReconTildePhiFullTH->Clone("acceptanceTildePhiFullT");
  // acceptanceTildePhiFullT->Divide(fGenerTildePhiFullTH);
  acceptanceTildePhiFullT->Divide(acceptanceTildePhiFullT, fGenerTildePhiFullTH, 1, 1, "b");
  acceptanceTildePhiFullT->Draw("ep");

  TCanvas* CorrCanvasTildePhiFullT = new TCanvas("CorrCanvasTildePhiFullT","CorrCanvasTildePhiFullT",900,800);
  DataTildePhiFullT0->Divide(acceptanceTildePhiFullT);
  DataTildePhiFullT0->Draw("ep");
  DataTildePhiShortT0->Divide(acceptanceTildePhiFullT);
  DataTildePhiShortT0->Draw("ep");
  DataTildePhiShortL0->Divide(acceptanceTildePhiFullT);
  DataTildePhiShortL0->Draw("ep");







  //_____________________________________________________________________________
  // Short Transverse
  //
  TCanvas* AcceptanceCanvasCosThetaShortT = new TCanvas("AcceptanceCanvasCosThetaShortT","AcceptanceCanvasCosThetaShortT",900,800);
  TH1F* acceptanceCosThetaShortT = (TH1F*) fReconCosThetaShortTH->Clone("acceptanceCosThetaShortT");
  // acceptanceCosThetaShortT->Divide(fGenerCosThetaShortTH);
  acceptanceCosThetaShortT->Divide(acceptanceCosThetaShortT, fGenerCosThetaShortTH, 1, 1, "b");
  acceptanceCosThetaShortT->Draw("ep");

  TCanvas* CorrCanvasCosThetaShortT = new TCanvas("CorrCanvasCosThetaShortT","CorrCanvasCosThetaShortT",900,800);
  DataThetaFullT1->Divide(acceptanceCosThetaShortT);
  DataThetaFullT1->Draw("ep");
  DataThetaShortT1->Divide(acceptanceCosThetaShortT);
  DataThetaShortT1->Draw("ep");
  DataThetaShortL1->Divide(acceptanceCosThetaShortT);
  DataThetaShortL1->Draw("ep");

  TCanvas* AcceptanceCanvasPhiShortT = new TCanvas("AcceptanceCanvasPhiShortT","AcceptanceCanvasPhiShortT",900,800);
  TH1F* acceptancePhiShortT = (TH1F*) fReconPhiShortTH->Clone("acceptancePhiShortT");
  // acceptancePhiShortT->Divide(fGenerPhiShortTH);
  acceptancePhiShortT->Divide(acceptancePhiShortT, fGenerPhiShortTH, 1, 1, "b");
  acceptancePhiShortT->Draw("ep");

  TCanvas* CorrCanvasPhiShortT = new TCanvas("CorrCanvasPhiShortT","CorrCanvasPhiShortT",900,800);
  DataPhiFullT1->Divide(acceptancePhiShortT);
  DataPhiFullT1->Draw("ep");
  DataPhiShortT1->Divide(acceptancePhiShortT);
  DataPhiShortT1->Draw("ep");
  DataPhiShortL1->Divide(acceptancePhiShortT);
  DataPhiShortL1->Draw("ep");

  TCanvas* AcceptanceCanvasTildePhiShortT = new TCanvas("AcceptanceCanvasTildePhiShortT","AcceptanceCanvasTildePhiShortT",900,800);
  TH1F* acceptanceTildePhiShortT = (TH1F*) fReconTildePhiShortTH->Clone("acceptanceTildePhiShortT");
  // acceptanceTildePhiShortT->Divide(fGenerTildePhiShortTH);
  acceptanceTildePhiShortT->Divide(acceptanceTildePhiShortT, fGenerTildePhiShortTH, 1, 1, "b");
  acceptanceTildePhiShortT->Draw("ep");

  TCanvas* CorrCanvasTildePhiShortT = new TCanvas("CorrCanvasTildePhiShortT","CorrCanvasTildePhiShortT",900,800);
  DataTildePhiFullT1->Divide(acceptanceTildePhiShortT);
  DataTildePhiFullT1->Draw("ep");
  DataTildePhiShortT1->Divide(acceptanceTildePhiShortT);
  DataTildePhiShortT1->Draw("ep");
  DataTildePhiShortL1->Divide(acceptanceTildePhiShortT);
  DataTildePhiShortL1->Draw("ep");








  //_____________________________________________________________________________
  // Short Longitudinal
  //
  TCanvas* AcceptanceCanvasCosThetaShortL = new TCanvas("AcceptanceCanvasCosThetaShortL","AcceptanceCanvasCosThetaShortL",900,800);
  TH1F* acceptanceCosThetaShortL = (TH1F*) fReconCosThetaShortLH->Clone("acceptanceCosThetaShortL");
  // acceptanceCosThetaShortL->Divide(fGenerCosThetaShortLH);
  acceptanceCosThetaShortL->Divide(acceptanceCosThetaShortL, fGenerCosThetaShortLH, 1, 1, "b");
  acceptanceCosThetaShortL->Draw("ep");

  TCanvas* CorrCanvasCosThetaShortL = new TCanvas("CorrCanvasCosThetaShortL","CorrCanvasCosThetaShortL",900,800);
  DataThetaFullT2->Divide(acceptanceCosThetaShortL);
  DataThetaFullT2->Draw("ep");
  DataThetaShortT2->Divide(acceptanceCosThetaShortL);
  DataThetaShortT2->Draw("ep");
  DataThetaShortL2->Divide(acceptanceCosThetaShortL);
  DataThetaShortL2->Draw("ep");

  TCanvas* AcceptanceCanvasPhiShortL = new TCanvas("AcceptanceCanvasPhiShortL","AcceptanceCanvasPhiShortL",900,800);
  TH1F* acceptancePhiShortL = (TH1F*) fReconPhiShortLH->Clone("acceptancePhiShortL");
  // acceptancePhiShortL->Divide(fGenerPhiShortLH);
  acceptancePhiShortL->Divide(acceptancePhiShortL, fGenerPhiShortLH, 1, 1, "b");
  acceptancePhiShortL->Draw("ep");

  TCanvas* CorrCanvasPhiShortL = new TCanvas("CorrCanvasPhiShortL","CorrCanvasPhiShortL",900,800);
  DataPhiFullT2->Divide(acceptancePhiShortL);
  DataPhiFullT2->Draw("ep");
  DataPhiShortT2->Divide(acceptancePhiShortL);
  DataPhiShortT2->Draw("ep");
  DataPhiShortL2->Divide(acceptancePhiShortL);
  DataPhiShortL2->Draw("ep");

  TCanvas* AcceptanceCanvasTildePhiShortL = new TCanvas("AcceptanceCanvasTildePhiShortL","AcceptanceCanvasTildePhiShortL",900,800);
  TH1F* acceptanceTildePhiShortL = (TH1F*) fReconTildePhiShortLH->Clone("acceptanceTildePhiShortL");
  // acceptanceTildePhiShortL->Divide(fGenerTildePhiShortLH);
  acceptanceTildePhiShortL->Divide(acceptanceTildePhiShortL, fGenerTildePhiShortLH, 1, 1, "b");
  acceptanceTildePhiShortL->Draw("ep");

  TCanvas* CorrCanvasTildePhiShortL = new TCanvas("CorrCanvasTildePhiShortL","CorrCanvasTildePhiShortL",900,800);
  DataTildePhiFullT2->Divide(acceptanceTildePhiShortL);
  DataTildePhiFullT2->Draw("ep");
  DataTildePhiShortT2->Divide(acceptanceTildePhiShortL);
  DataTildePhiShortT2->Draw("ep");
  DataTildePhiShortL2->Divide(acceptanceTildePhiShortL);
  DataTildePhiShortL2->Draw("ep");





  // TH1F* CorrCosThetaH = (TH1F*) RawCosThetaH->Clone("CorrCosThetaH");
  // TH1F* CorrPhiH      = (TH1F*) RawPhiH     ->Clone("CorrPhiH");
  // TH1F* CorrTildePhiH = (TH1F*) RawTildePhiH->Clone("CorrTildePhiH");

  // CLOSURE TEST 2
  TH1F* CorrectedDistribCosTheta[9]  = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
  CorrectedDistribCosTheta[0]  = (TH1F*) DataThetaFullT0 ->Clone("CorrectedDistribCosTheta_0");
  CorrectedDistribCosTheta[1]  = (TH1F*) DataThetaShortT0->Clone("CorrectedDistribCosTheta_1");
  CorrectedDistribCosTheta[2]  = (TH1F*) DataThetaShortL0->Clone("CorrectedDistribCosTheta_2");
  CorrectedDistribCosTheta[3]  = (TH1F*) DataThetaFullT1 ->Clone("CorrectedDistribCosTheta_3");
  CorrectedDistribCosTheta[4]  = (TH1F*) DataThetaShortT1->Clone("CorrectedDistribCosTheta_4");
  CorrectedDistribCosTheta[5]  = (TH1F*) DataThetaShortL1->Clone("CorrectedDistribCosTheta_5");
  CorrectedDistribCosTheta[6]  = (TH1F*) DataThetaFullT2 ->Clone("CorrectedDistribCosTheta_6");
  CorrectedDistribCosTheta[7]  = (TH1F*) DataThetaShortT2->Clone("CorrectedDistribCosTheta_7");
  CorrectedDistribCosTheta[8]  = (TH1F*) DataThetaShortL2->Clone("CorrectedDistribCosTheta_8");
  TH1F* CorrectedDistribCosTheta2[9] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
  CorrectedDistribCosTheta2[0] = (TH1F*) DataThetaFullT0 ->Clone("CorrectedDistribCosTheta2_0");
  CorrectedDistribCosTheta2[1] = (TH1F*) DataThetaShortT0->Clone("CorrectedDistribCosTheta2_1");
  CorrectedDistribCosTheta2[2] = (TH1F*) DataThetaShortL0->Clone("CorrectedDistribCosTheta2_2");
  CorrectedDistribCosTheta2[3] = (TH1F*) DataThetaFullT1 ->Clone("CorrectedDistribCosTheta2_3");
  CorrectedDistribCosTheta2[4] = (TH1F*) DataThetaShortT1->Clone("CorrectedDistribCosTheta2_4");
  CorrectedDistribCosTheta2[5] = (TH1F*) DataThetaShortL1->Clone("CorrectedDistribCosTheta2_5");
  CorrectedDistribCosTheta2[6] = (TH1F*) DataThetaFullT2 ->Clone("CorrectedDistribCosTheta2_6");
  CorrectedDistribCosTheta2[7] = (TH1F*) DataThetaShortT2->Clone("CorrectedDistribCosTheta2_7");
  CorrectedDistribCosTheta2[8] = (TH1F*) DataThetaShortL2->Clone("CorrectedDistribCosTheta2_8");


  TH1F* CorrectedDistribPhi[9]  = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
  CorrectedDistribPhi[0]  = (TH1F*) DataPhiFullT0 ->Clone("CorrectedDistribPhi_0");
  CorrectedDistribPhi[1]  = (TH1F*) DataPhiShortT0->Clone("CorrectedDistribPhi_1");
  CorrectedDistribPhi[2]  = (TH1F*) DataPhiShortL0->Clone("CorrectedDistribPhi_2");
  CorrectedDistribPhi[3]  = (TH1F*) DataPhiFullT1 ->Clone("CorrectedDistribPhi_3");
  CorrectedDistribPhi[4]  = (TH1F*) DataPhiShortT1->Clone("CorrectedDistribPhi_4");
  CorrectedDistribPhi[5]  = (TH1F*) DataPhiShortL1->Clone("CorrectedDistribPhi_5");
  CorrectedDistribPhi[6]  = (TH1F*) DataPhiFullT2 ->Clone("CorrectedDistribPhi_6");
  CorrectedDistribPhi[7]  = (TH1F*) DataPhiShortT2->Clone("CorrectedDistribPhi_7");
  CorrectedDistribPhi[8]  = (TH1F*) DataPhiShortL2->Clone("CorrectedDistribPhi_8");
  TH1F* CorrectedDistribPhi2[9] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
  CorrectedDistribPhi2[0] = (TH1F*) DataPhiFullT0 ->Clone("CorrectedDistribPhi2_0");
  CorrectedDistribPhi2[1] = (TH1F*) DataPhiShortT0->Clone("CorrectedDistribPhi2_1");
  CorrectedDistribPhi2[2] = (TH1F*) DataPhiShortL0->Clone("CorrectedDistribPhi2_2");
  CorrectedDistribPhi2[3] = (TH1F*) DataPhiFullT1 ->Clone("CorrectedDistribPhi2_3");
  CorrectedDistribPhi2[4] = (TH1F*) DataPhiShortT1->Clone("CorrectedDistribPhi2_4");
  CorrectedDistribPhi2[5] = (TH1F*) DataPhiShortL1->Clone("CorrectedDistribPhi2_5");
  CorrectedDistribPhi2[6] = (TH1F*) DataPhiFullT2 ->Clone("CorrectedDistribPhi2_6");
  CorrectedDistribPhi2[7] = (TH1F*) DataPhiShortT2->Clone("CorrectedDistribPhi2_7");
  CorrectedDistribPhi2[8] = (TH1F*) DataPhiShortL2->Clone("CorrectedDistribPhi2_8");



  TH1F* CorrectedDistribTildePhi[9] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
  CorrectedDistribTildePhi[0] = (TH1F*) DataTildePhiFullT0 ->Clone("CorrectedDistribTildePhi_0");
  CorrectedDistribTildePhi[1] = (TH1F*) DataTildePhiShortT0->Clone("CorrectedDistribTildePhi_1");
  CorrectedDistribTildePhi[2] = (TH1F*) DataTildePhiShortL0->Clone("CorrectedDistribTildePhi_2");
  CorrectedDistribTildePhi[3] = (TH1F*) DataTildePhiFullT1 ->Clone("CorrectedDistribTildePhi_3");
  CorrectedDistribTildePhi[4] = (TH1F*) DataTildePhiShortT1->Clone("CorrectedDistribTildePhi_4");
  CorrectedDistribTildePhi[5] = (TH1F*) DataTildePhiShortL1->Clone("CorrectedDistribTildePhi_5");
  CorrectedDistribTildePhi[6] = (TH1F*) DataTildePhiFullT2 ->Clone("CorrectedDistribTildePhi_6");
  CorrectedDistribTildePhi[7] = (TH1F*) DataTildePhiShortT2->Clone("CorrectedDistribTildePhi_7");
  CorrectedDistribTildePhi[8] = (TH1F*) DataTildePhiShortL2->Clone("CorrectedDistribTildePhi_8");
  TH1F* CorrectedDistribTildePhi2[9] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
  CorrectedDistribTildePhi2[0] = (TH1F*) DataTildePhiFullT0 ->Clone("CorrectedDistribTildePhi2_0");
  CorrectedDistribTildePhi2[1] = (TH1F*) DataTildePhiShortT0->Clone("CorrectedDistribTildePhi2_1");
  CorrectedDistribTildePhi2[2] = (TH1F*) DataTildePhiShortL0->Clone("CorrectedDistribTildePhi2_2");
  CorrectedDistribTildePhi2[3] = (TH1F*) DataTildePhiFullT1 ->Clone("CorrectedDistribTildePhi2_3");
  CorrectedDistribTildePhi2[4] = (TH1F*) DataTildePhiShortT1->Clone("CorrectedDistribTildePhi2_4");
  CorrectedDistribTildePhi2[5] = (TH1F*) DataTildePhiShortL1->Clone("CorrectedDistribTildePhi2_5");
  CorrectedDistribTildePhi2[6] = (TH1F*) DataTildePhiFullT2 ->Clone("CorrectedDistribTildePhi2_6");
  CorrectedDistribTildePhi2[7] = (TH1F*) DataTildePhiShortT2->Clone("CorrectedDistribTildePhi2_7");
  CorrectedDistribTildePhi2[8] = (TH1F*) DataTildePhiShortL2->Clone("CorrectedDistribTildePhi2_8");

  for( Int_t i = 0; i < 9; i++ ){
    CorrectedDistribCosTheta[i] ->SetLineWidth(3);
    CorrectedDistribPhi[i]      ->SetLineWidth(3);
    CorrectedDistribTildePhi[i] ->SetLineWidth(3);
    CorrectedDistribCosTheta2[i]->SetLineWidth(3);
    CorrectedDistribPhi2[i]     ->SetLineWidth(3);
    CorrectedDistribTildePhi2[i]->SetLineWidth(3);
  }


  Int_t binThetaLow  = CorrectedDistribCosTheta[1] ->GetXaxis()->FindBin(-0.5);
  Int_t binThetaHigh = CorrectedDistribCosTheta[1] ->GetXaxis()->FindBin(0.5);


  fGenerCosThetaFullTH ->SetLineWidth(3);
  fGenerCosThetaShortTH->SetLineWidth(3);
  fGenerCosThetaShortLH->SetLineWidth(3);
  fGenerCosThetaFullTH ->SetLineColor(kRed);
  fGenerCosThetaShortTH->SetLineColor(kRed);
  fGenerCosThetaShortLH->SetLineColor(kRed);

  fGenerPhiFullTH ->SetLineWidth(3);
  fGenerPhiShortTH->SetLineWidth(3);
  fGenerPhiShortLH->SetLineWidth(3);
  fGenerPhiFullTH ->SetLineColor(kRed);
  fGenerPhiShortTH->SetLineColor(kRed);
  fGenerPhiShortLH->SetLineColor(kRed);

  fGenerTildePhiFullTH ->SetLineWidth(3);
  fGenerTildePhiShortTH->SetLineWidth(3);
  fGenerTildePhiShortLH->SetLineWidth(3);
  fGenerTildePhiFullTH ->SetLineColor(kRed);
  fGenerTildePhiShortTH->SetLineColor(kRed);
  fGenerTildePhiShortLH->SetLineColor(kRed);


  // Full T   vs Full T
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 800);
  c1->Divide(3,2);
  c1->cd(1);
  CorrectedDistribCosTheta[0] ->Draw("ep");
  TH1F* fGenerCosThetaFullTHc1 = (TH1F*) fGenerCosThetaFullTH->Clone("fGenerCosThetaFullTHc1");
  fGenerCosThetaFullTHc1      ->Draw("same ep");
  c1->cd(2);
  CorrectedDistribPhi[0]->GetYaxis()->SetRangeUser(CorrectedDistribPhi[0]->GetMaximum()*0.6, CorrectedDistribPhi[0]->GetMaximum()*1.5);
  CorrectedDistribPhi[0]      ->Draw("ep");
  TH1F* fGenerPhiFullTHc1 = (TH1F*) fGenerPhiFullTH->Clone("fGenerPhiFullTHc1");
  fGenerPhiFullTHc1           ->Draw("same ep");
  c1->cd(3);
  CorrectedDistribTildePhi[0]->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi[0]->GetMaximum()*0.6, CorrectedDistribTildePhi[0]->GetMaximum()*1.5);
  CorrectedDistribTildePhi[0] ->Draw("ep");
  TH1F* fGenerTildePhiFullTHc1 = (TH1F*) fGenerTildePhiFullTH->Clone("fGenerTildePhiFullTHc1");
  fGenerTildePhiFullTH        ->Draw("same ep");
  c1->cd(4);
  CorrectedDistribCosTheta2[0]->Divide(fGenerCosThetaFullTH);
  CorrectedDistribCosTheta2[0]->Draw("ep");
  c1->cd(5);
  CorrectedDistribPhi2[0]     ->Divide(fGenerPhiFullTH);
  CorrectedDistribPhi2[0]     ->GetYaxis()->SetRangeUser(CorrectedDistribPhi2[0]->GetMaximum()*0.6, CorrectedDistribPhi2[0]->GetMaximum()*1.5);
  CorrectedDistribPhi2[0]     ->Draw("ep");
  c1->cd(6);
  CorrectedDistribTildePhi2[0]->Divide(fGenerTildePhiFullTH);
  CorrectedDistribTildePhi2[0]->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi2[0]->GetMaximum()*0.6, CorrectedDistribTildePhi2[0]->GetMaximum()*1.5);
  CorrectedDistribTildePhi2[0]->Draw("ep");






  // Small T vs Small T
  TCanvas *c2 = new TCanvas("c2", "c2", 900, 800);
  c2->Divide(3,2);
  c2->cd(1);
  CorrectedDistribCosTheta[4] ->Draw("ep");
  TH1F* fGenerCosThetaShortTHc2 = (TH1F*) fGenerCosThetaShortTH->Clone("fGenerCosThetaShortTHc2");
  fGenerCosThetaShortTHc2     ->Draw("same ep");
  c2->cd(2);
  CorrectedDistribPhi[4]      ->GetYaxis()->SetRangeUser(CorrectedDistribPhi[4]      ->GetMaximum()*0.6, CorrectedDistribPhi[4]      ->GetMaximum()*1.5);
  CorrectedDistribPhi[4]      ->Draw("ep");
  TH1F* fGenerPhiShortTHc2 = (TH1F*) fGenerPhiShortTH->Clone("fGenerPhiShortTHc2");
  fGenerPhiShortTHc2          ->Draw("same ep");
  c2->cd(3);
  CorrectedDistribTildePhi[4] ->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi[4] ->GetMaximum()*0.6, CorrectedDistribTildePhi[4] ->GetMaximum()*1.5);
  CorrectedDistribTildePhi[4] ->Draw("ep");
  TH1F* fGenerTildePhiShortTHc2 = (TH1F*) fGenerTildePhiShortTH->Clone("fGenerTildePhiShortTHc2");
  fGenerTildePhiShortTHc2     ->Draw("same ep");
  c2->cd(4);
  CorrectedDistribCosTheta2[4]->Divide(fGenerCosThetaShortTH);
  CorrectedDistribCosTheta2[4]->GetYaxis()->SetRangeUser(CorrectedDistribCosTheta2[4]->GetMaximum()*0.6, CorrectedDistribCosTheta2[4]->GetMaximum()*1.5);
  CorrectedDistribCosTheta2[4]->Draw("ep");
  c2->cd(5);
  CorrectedDistribPhi2[4]     ->Divide(fGenerPhiShortTH);
  CorrectedDistribPhi2[4]     ->GetYaxis()->SetRangeUser(CorrectedDistribPhi2[4]     ->GetMaximum()*0.6, CorrectedDistribPhi2[4]     ->GetMaximum()*1.5);
  CorrectedDistribPhi2[4]     ->Draw("ep");
  c2->cd(6);
  CorrectedDistribTildePhi2[4]->Divide(fGenerTildePhiShortTH);
  CorrectedDistribTildePhi2[4]->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi2[4]->GetMaximum()*0.6, CorrectedDistribTildePhi2[4]->GetMaximum()*1.5);
  CorrectedDistribTildePhi2[4]->Draw("ep");







  // Full T vs    Small T / Full T
  TCanvas *c3 = new TCanvas("c3", "c3", 900, 800);
  c3->Divide(3,2);
  // c3->cd(1);
  // CorrectedDistribCosTheta[1] ->Draw("ep");
  // fGenerCosThetaFullTH        ->Draw("same ep");
  // c3->cd(2);
  // CorrectedDistribPhi[1]      ->Draw("ep");
  // fGenerPhiFullTH             ->Draw("same ep");
  // c3->cd(3);
  // CorrectedDistribTildePhi[1] ->Draw("ep");
  // fGenerTildePhiFullTH        ->Draw("same ep");
  TH1F* fGenerCosThetaFullTHc3 = (TH1F*) fGenerCosThetaFullTH->Clone("fGenerCosThetaFullTHc3");
  TH1F* fGenerPhiFullTHc3      = (TH1F*) fGenerPhiFullTH     ->Clone("fGenerPhiFullTHc3");
  TH1F* fGenerTildePhiFullTHc3 = (TH1F*) fGenerTildePhiFullTH->Clone("fGenerTildePhiFullTHc3");
  c3->cd(4);
  CorrectedDistribCosTheta2[1]->Divide(fGenerCosThetaFullTH);
  CorrectedDistribCosTheta2[1]->Scale( 2./( CorrectedDistribCosTheta2[1]->GetMaximum()+CorrectedDistribCosTheta2[1]->GetMinimum()  ) );
  CorrectedDistribCosTheta2[1]->Scale( 0.7 );
  // CorrectedDistribCosTheta2[1]->Scale( 1./( CorrectedDistribCosTheta2[1]->Integral() ) );
  // CorrectedDistribCosTheta2[1]->Scale( 1.2/CorrectedDistribCosTheta2[1]->GetXaxis()->GetBinWidth(1) );
  // CorrectedDistribCosTheta2[1]->Scale( 1./( fGenerCosThetaFullTHc3->Integral(binThetaLow, binThetaHigh) ) );
  // CorrectedDistribCosTheta2[1]->Scale( CorrectedDistribCosTheta[1]->Integral(binThetaLow, binThetaHigh) );
  CorrectedDistribCosTheta2[1]->Draw("ep");
  c3->cd(5);
  CorrectedDistribPhi2[1]     ->Divide(fGenerPhiFullTH);
  CorrectedDistribPhi2[1]     ->Scale( 2./( CorrectedDistribPhi2[1]->GetMaximum()+CorrectedDistribPhi2[1]->GetMinimum()  ) );
  // CorrectedDistribPhi2[1]     ->Scale( 1./( CorrectedDistribPhi2[1]     ->Integral() ) );
  // CorrectedDistribPhi2[1]     ->Scale( 100.*CorrectedDistribPhi2[1]     ->GetXaxis()->GetBinWidth(1) );
  // CorrectedDistribPhi2[1]     ->Scale( 1./( fGenerPhiFullTHc3->Integral() ) );
  // CorrectedDistribPhi2[1]     ->Scale( CorrectedDistribPhi[1]->Integral() );
  CorrectedDistribPhi2[1]     ->GetYaxis()->SetRangeUser(CorrectedDistribPhi2[1]     ->GetMaximum()*0.6, CorrectedDistribPhi2[1]     ->GetMaximum()*1.5);
  CorrectedDistribPhi2[1]     ->Draw("ep");
  c3->cd(6);
  CorrectedDistribTildePhi2[1]->Divide(fGenerTildePhiFullTH);
  CorrectedDistribTildePhi2[1]->Scale( 2./( CorrectedDistribTildePhi2[1]->GetMaximum()+CorrectedDistribTildePhi2[1]->GetMinimum()  ) );
  // CorrectedDistribTildePhi2[1]->Scale( 1./( fGenerTildePhiFullTHc3->Integral() ) );
  // CorrectedDistribTildePhi2[1]->Scale( 10.*CorrectedDistribTildePhi2[1]->GetXaxis()->GetBinWidth(1) );
  // CorrectedDistribTildePhi2[1]->Scale( 1./( fGenerTildePhiFullTHc3->Integral() ) );
  // CorrectedDistribTildePhi2[1]->Scale( CorrectedDistribTildePhi[1]->Integral() );
  CorrectedDistribTildePhi2[1]->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi2[1]->GetMaximum()*0.6, CorrectedDistribTildePhi2[1]->GetMaximum()*1.5);
  CorrectedDistribTildePhi2[1]->Draw("ep");
  c3->cd(1);
  // CorrectedDistribCosTheta[1] ->GetYaxis()->SetRangeUser(CorrectedDistribCosTheta[1] ->GetMaximum()*0.6, CorrectedDistribCosTheta[1] ->GetMaximum()*1.5);
  // CorrectedDistribCosTheta[1] ->Draw("ep");
  // TH1F* fGenerCosThetaFullTHc3 = (TH1F*) fGenerCosThetaFullTH->Clone("fGenerCosThetaFullTHc3");
  fGenerCosThetaFullTHc3      ->Scale( 1./( fGenerCosThetaFullTHc3->Integral(binThetaLow, binThetaHigh) ) );
  fGenerCosThetaFullTHc3      ->Scale( CorrectedDistribCosTheta[1]->Integral(binThetaLow, binThetaHigh) );
  CorrectedDistribCosTheta[1] ->GetYaxis()->SetRangeUser(CorrectedDistribCosTheta[1] ->GetMaximum()*0.6, CorrectedDistribCosTheta[1] ->GetMaximum()*1.5);
  CorrectedDistribCosTheta[1] ->Draw("ep");
  fGenerCosThetaFullTHc3      ->Draw("same ep");
  c3->cd(2);
  CorrectedDistribPhi[1]      ->GetYaxis()->SetRangeUser(CorrectedDistribPhi[1]      ->GetMaximum()*0.6, CorrectedDistribPhi[1]      ->GetMaximum()*1.5);
  CorrectedDistribPhi[1]      ->Draw("ep");
  // TH1F* fGenerPhiFullTHc3 = (TH1F*) fGenerPhiFullTH->Clone("fGenerPhiFullTHc3");
  fGenerPhiFullTHc3           ->Scale( 1./( fGenerPhiFullTHc3->Integral() ) );
  fGenerPhiFullTHc3           ->Scale( CorrectedDistribPhi[1]->Integral() );
  fGenerPhiFullTHc3           ->Draw("same ep");
  c3->cd(3);
  CorrectedDistribTildePhi[1] ->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi[1] ->GetMaximum()*0.6, CorrectedDistribTildePhi[1] ->GetMaximum()*1.5);
  CorrectedDistribTildePhi[1] ->Draw("ep");
  // TH1F* fGenerTildePhiFullTHc3 = (TH1F*) fGenerTildePhiFullTH->Clone("fGenerTildePhiFullTHc3");
  fGenerTildePhiFullTHc3      ->Scale( 1./( fGenerTildePhiFullTHc3->Integral() ) );
  fGenerTildePhiFullTHc3      ->Scale( CorrectedDistribTildePhi[1]->Integral() );
  fGenerTildePhiFullTHc3      ->Draw("same ep");






  // Small L vs Small L
  TCanvas *c4 = new TCanvas("c4", "c4", 900, 800);
  c4->Divide(3,2);
  c4->cd(1);
  CorrectedDistribCosTheta[8] ->Draw("ep");
  TH1F* fGenerCosThetaShortLHc4 = (TH1F*) fGenerCosThetaShortLH->Clone("fGenerCosThetaShortLHc4");
  fGenerCosThetaShortLHc4     ->Draw("same ep");
  c4->cd(2);
  CorrectedDistribPhi[8]      ->GetYaxis()->SetRangeUser(CorrectedDistribPhi[8]      ->GetMaximum()*0.6, CorrectedDistribPhi[8]      ->GetMaximum()*1.5);
  CorrectedDistribPhi[8]      ->Draw("ep");
  TH1F* fGenerPhiShortLHc4 = (TH1F*) fGenerPhiShortLH->Clone("fGenerPhiShortLHc4");
  fGenerPhiShortLHc4          ->Draw("same ep");
  fGenerPhiShortLHc4          ->Draw("same ep");
  c4->cd(3);
  CorrectedDistribTildePhi[8] ->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi[8] ->GetMaximum()*0.6, CorrectedDistribTildePhi[8] ->GetMaximum()*1.5);
  CorrectedDistribTildePhi[8] ->Draw("ep");
  TH1F* fGenerTildePhiShortLHc4 = (TH1F*) fGenerTildePhiShortLH->Clone("fGenerTildePhiShortLHc4");
  fGenerTildePhiShortLHc4     ->Draw("same ep");
  fGenerTildePhiShortLHc4     ->Draw("same ep");
  c4->cd(4);
  CorrectedDistribCosTheta2[8]->Divide(fGenerCosThetaShortLH);
  CorrectedDistribCosTheta2[8]->Draw("ep");
  c4->cd(5);
  CorrectedDistribPhi2[8]     ->Divide(fGenerPhiShortLH);
  CorrectedDistribPhi2[8]     ->GetYaxis()->SetRangeUser(CorrectedDistribPhi2[8]     ->GetMaximum()*0.6, CorrectedDistribPhi2[8]     ->GetMaximum()*1.5);
  CorrectedDistribPhi2[8]     ->Draw("ep");
  c4->cd(6);
  CorrectedDistribTildePhi2[8]->Divide(fGenerTildePhiShortLH);
  CorrectedDistribTildePhi2[8]->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi2[8]->GetMaximum()*0.6, CorrectedDistribTildePhi2[8]->GetMaximum()*1.5);
  CorrectedDistribTildePhi2[8]->Draw("ep");





  // Small L vs    Small L / Full T
  TCanvas *c5 = new TCanvas("c5", "c5", 900, 800);
  c5->Divide(3,2);
  // c5->cd(1);
  // CorrectedDistribCosTheta[2] ->Draw("ep");
  // fGenerCosThetaShortLH       ->Draw("same ep");
  // c5->cd(2);
  // CorrectedDistribPhi[2]      ->Draw("ep");
  // fGenerPhiShortLH            ->Draw("same ep");
  // c5->cd(3);
  // CorrectedDistribTildePhi[2] ->Draw("ep");
  // fGenerTildePhiShortLH       ->Draw("same ep");
  c5->cd(4);
  CorrectedDistribCosTheta2[2]->Divide(fGenerCosThetaShortLH);
  CorrectedDistribCosTheta2[2]->Draw("ep");
  c5->cd(5);
  CorrectedDistribPhi2[2]     ->Divide(fGenerPhiShortLH);
  CorrectedDistribPhi2[2]     ->Scale( 2./( CorrectedDistribPhi2[2]     ->GetMaximum()+CorrectedDistribPhi2[2]     ->GetMinimum()  ) );
  CorrectedDistribPhi2[2]     ->GetYaxis()->SetRangeUser(CorrectedDistribPhi2[2]     ->GetMaximum()*0.6, CorrectedDistribPhi2[2]     ->GetMaximum()*1.5);
  CorrectedDistribPhi2[2]     ->Draw("ep");
  c5->cd(6);
  CorrectedDistribTildePhi2[2]->Divide(fGenerTildePhiShortLH);
  CorrectedDistribTildePhi2[2]->Scale( 2./( CorrectedDistribTildePhi2[2]->GetMaximum()+CorrectedDistribTildePhi2[2]->GetMinimum()  ) );
  CorrectedDistribTildePhi2[2]->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi2[2]->GetMaximum()*0.6, CorrectedDistribTildePhi2[2]->GetMaximum()*1.5);
  CorrectedDistribTildePhi2[2]->Draw("ep");
  c5->cd(1);
  CorrectedDistribCosTheta[2] ->Draw("ep");
  TH1F* fGenerCosThetaShortLHc5 = (TH1F*) fGenerCosThetaShortLH->Clone("fGenerCosThetaShortLHc5");
  fGenerCosThetaShortLHc5     ->Draw("same ep");
  c5->cd(2);
  CorrectedDistribPhi[2]      ->GetYaxis()->SetRangeUser(CorrectedDistribPhi[2]      ->GetMaximum()*0.6, CorrectedDistribPhi[2]      ->GetMaximum()*1.5);
  CorrectedDistribPhi[2]      ->Draw("ep");
  TH1F* fGenerPhiShortLHc5 = (TH1F*) fGenerPhiShortLH->Clone("fGenerPhiShortLHc5");
  fGenerPhiShortLHc5          ->Scale( 1./( fGenerPhiShortLHc5->Integral() ) );
  fGenerPhiShortLHc5          ->Scale( CorrectedDistribPhi[2]->Integral() );
  fGenerPhiShortLHc5          ->Draw("same ep");
  c5->cd(3);
  CorrectedDistribTildePhi[2] ->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi[2] ->GetMaximum()*0.6, CorrectedDistribTildePhi[2] ->GetMaximum()*1.5);
  CorrectedDistribTildePhi[2] ->Draw("ep");
  TH1F* fGenerTildePhiShortLHc5 = (TH1F*) fGenerTildePhiShortLH->Clone("fGenerTildePhiShortLHc5");
  fGenerTildePhiShortLHc5     ->Scale( 1./( fGenerTildePhiShortLHc5->Integral() ) );
  fGenerTildePhiShortLHc5     ->Scale( CorrectedDistribTildePhi[2]->Integral() );
  fGenerTildePhiShortLHc5     ->Draw("same ep");





  // Small T vs Small T / Small L
  TCanvas *c6 = new TCanvas("c6", "c6", 900, 800);
  c6->Divide(3,2);
  // c6->cd(1);
  // CorrectedDistribCosTheta[7] ->Draw("ep");
  // fGenerCosThetaShortTH       ->Draw("same ep");
  // c6->cd(2);
  // CorrectedDistribPhi[7]      ->Draw("ep");
  // fGenerPhiShortTH            ->Draw("same ep");
  // c6->cd(3);
  // CorrectedDistribTildePhi[7] ->Draw("ep");
  // fGenerTildePhiShortTH       ->Draw("same ep");
  c6->cd(4);
  CorrectedDistribCosTheta2[7]->Divide(fGenerCosThetaShortTH);
  CorrectedDistribCosTheta2[7]->Draw("ep");
  c6->cd(5);
  CorrectedDistribPhi2[7]     ->Divide(fGenerPhiShortTH);
  CorrectedDistribPhi2[7]     ->Scale( 2./( CorrectedDistribPhi2[7]     ->GetMaximum()+CorrectedDistribPhi2[7]     ->GetMinimum()  ) );
  CorrectedDistribPhi2[7]     ->GetYaxis()->SetRangeUser(CorrectedDistribPhi2[7]     ->GetMaximum()*0.6, CorrectedDistribPhi2[7]     ->GetMaximum()*1.5);
  CorrectedDistribPhi2[7]     ->Draw("ep");
  c6->cd(6);
  CorrectedDistribTildePhi2[7]->Divide(fGenerTildePhiShortTH);
  CorrectedDistribTildePhi2[7]->Scale( 2./( CorrectedDistribTildePhi2[7]->GetMaximum()+CorrectedDistribTildePhi2[7]->GetMinimum()  ) );
  CorrectedDistribTildePhi2[7]->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi2[7]->GetMaximum()*0.6, CorrectedDistribTildePhi2[7]->GetMaximum()*1.5);
  CorrectedDistribTildePhi2[7]->Draw("ep");
  c6->cd(1);
  CorrectedDistribCosTheta[7] ->Draw("ep");
  TH1F* fGenerCosThetaShortTHc6 = (TH1F*) fGenerCosThetaShortTH->Clone("fGenerCosThetaShortTHc6");
  fGenerCosThetaShortTHc6       ->Draw("same ep");
  c6->cd(2);
  CorrectedDistribPhi[7]      ->GetYaxis()->SetRangeUser(CorrectedDistribPhi[7]      ->GetMaximum()*0.6, CorrectedDistribPhi[7]      ->GetMaximum()*1.5);
  CorrectedDistribPhi[7]      ->Draw("ep");
  TH1F* fGenerPhiShortTHc6 = (TH1F*) fGenerPhiShortTH->Clone("fGenerPhiShortTHc6");
  fGenerPhiShortTHc6          ->Scale( 1./( fGenerPhiShortTHc6->Integral() ) );
  fGenerPhiShortTHc6          ->Scale( CorrectedDistribPhi[7]->Integral() );
  fGenerPhiShortTHc6          ->Draw("same ep");
  c6->cd(3);
  CorrectedDistribTildePhi[7] ->GetYaxis()->SetRangeUser(CorrectedDistribTildePhi[7] ->GetMaximum()*0.6, CorrectedDistribTildePhi[7] ->GetMaximum()*1.5);
  CorrectedDistribTildePhi[7] ->Draw("ep");
  TH1F* fGenerTildePhiShortTHc6 = (TH1F*) fGenerTildePhiShortTH->Clone("fGenerTildePhiShortTHc6");
  fGenerTildePhiShortTHc6     ->Scale( 1./( fGenerTildePhiShortTHc6->Integral() ) );
  fGenerTildePhiShortTHc6     ->Scale( CorrectedDistribTildePhi[7]->Integral() );
  fGenerTildePhiShortTHc6     ->Draw("same ep");






  TFile f("pngResults/PolarisationCorrectedCs1D_closuretest.root", "recreate");
  acceptanceCosThetaFullT ->Write();
  acceptancePhiFullT      ->Write();
  acceptanceTildePhiFullT ->Write();
  DataThetaFullT0         ->Write();
  DataThetaShortT0        ->Write();
  DataThetaShortL0        ->Write();
  DataPhiFullT0           ->Write();
  DataPhiShortT0          ->Write();
  DataPhiShortL0          ->Write();
  DataTildePhiFullT0      ->Write();
  DataTildePhiShortT0     ->Write();
  DataTildePhiShortL0     ->Write();
  acceptanceCosThetaShortT->Write();
  acceptancePhiShortT     ->Write();
  acceptanceTildePhiShortT->Write();
  DataThetaFullT1         ->Write();
  DataThetaShortT1        ->Write();
  DataThetaShortL1        ->Write();
  DataPhiFullT1           ->Write();
  DataPhiShortT1          ->Write();
  DataPhiShortL1          ->Write();
  DataTildePhiFullT1      ->Write();
  DataTildePhiShortT1     ->Write();
  DataTildePhiShortL1     ->Write();
  acceptanceCosThetaShortL->Write();
  acceptancePhiShortL     ->Write();
  acceptanceTildePhiShortL->Write();
  DataThetaFullT2         ->Write();
  DataThetaShortT2        ->Write();
  DataThetaShortL2        ->Write();
  DataPhiFullT2           ->Write();
  DataPhiShortT2          ->Write();
  DataPhiShortL2          ->Write();
  DataTildePhiFullT2      ->Write();
  DataTildePhiShortT2     ->Write();
  DataTildePhiShortL2     ->Write();
  f.Close();
}
