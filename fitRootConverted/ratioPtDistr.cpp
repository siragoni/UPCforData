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
/* - Fit function for the ZNC plots.
 * -
 */
void ratioPtDistr(const char* AnalysisName, const char* AnalysisName2){
  TFile* fileList  = new TFile(AnalysisName);
  TDirectory* dir  = fileList ->GetDirectory("MyTask");
  TFile* fileList2 = new TFile(AnalysisName2);
  TDirectory* dir2 = fileList2->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listings;
  dir ->GetObject("MyOutputContainer", listings);
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

  TCanvas* RatioCanvas = new TCanvas("RatioCanvas", "RatioCanvas", 1000, 200);

  TH1F* fDimuonPtDistributionDataH  = (TH1F*)listings ->FindObject("fDimuonPtDistributionH");
  TH1F* fDimuonPtDistributionDataH2 = (TH1F*)listings2->FindObject("fDimuonPtDistributionH");
  TCanvas* PlottingSame = new TCanvas("PlottingSame", "PlottingSame", 1000, 900);
  fDimuonPtDistributionDataH ->SetLineColor(kRed);
  fDimuonPtDistributionDataH ->Draw();
  fDimuonPtDistributionDataH2->Draw("same");

  PlottingSame->SaveAs("pngResults/PlottingSame.png", "RECREATE");

  // fDimuonPtDistributionDataH->Rebin(5);
  // fDimuonPtDistributionDataH->Sumw2();
  //
  // fDimuonPtDistributionDataH2->Rebin(5);
  // fDimuonPtDistributionDataH2->Sumw2();
  //
  // TH1F* RatioH = (TH1F*)fDimuonPtDistributionDataH2->Clone("RatioH");
  // // RatioH->Rebin(5);
  // RatioH->Sumw2();
  // RatioH->Divide(fDimuonPtDistributionDataH);
  //
  //
  // RatioH->SetMarkerStyle(21);
  // RatioH->SetLineColor(kBlue);
  // RatioH->SetLineWidth(2);
  // RatioH->Draw("ep");
  // gPad->SaveAs("pngResults/RatioPtDistr.png", "RECREATE");



}
