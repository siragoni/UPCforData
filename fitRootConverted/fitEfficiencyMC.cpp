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
/* - Computes the efficiency of the MC on a
 * - run-by-run basis.
 * -
 */
void fitEfficiencyMC(){

  // TFile* fileList = new TFile("../MyUPC-MonteCarlo/AnalysisResults.root");
  // TFile* fileList = new TFile("AnalysisResultsMC06062019.root");
  // TFile* fileList = new TFile("MCtrainResults/2019-06-16-LHC18qr/kCohJpsiToMu/AnalysisResults.root");
  TFile* fileList = new TFile("MCtrainResults/2019-06-16-LHC15o/kCohJpsiToMu/AnalysisResults.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);
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
  TH1F* fEfficiencyPerRunH   = (TH1F*)listings->FindObject("fEfficiencyPerRunH");
  TH1F* fMCEfficiencyPerRunH = (TH1F*)listings->FindObject("fMCEfficiencyPerRunH");
  fEfficiencyPerRunH  ->Sumw2();
  fMCEfficiencyPerRunH->Sumw2();

  TCanvas* EffCanvas = new TCanvas("EffCanvas","EffCanvas",900,800);
  TH1F* RealEfficiency = (TH1F*) fEfficiencyPerRunH->Clone("RealEfficiency");
  TH1F* MCEfficiency = (TH1F*) fEfficiencyPerRunH->Clone("RealEfficiency");
  RealEfficiency->Divide(fMCEfficiencyPerRunH);
  RealEfficiency->Draw("ep");




  TFile f("pngResults/efficiency.root", "recreate");
  RealEfficiency->Write();
  f.Close();
}
