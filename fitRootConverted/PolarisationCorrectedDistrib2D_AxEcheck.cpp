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
void AxEcheck(){

  // TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  TFile* fileList = new TFile("AnalysisResultsMC_flat2D.root");
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
  // TH2F* fReconH = (TH2F*)listings->FindObject("fCosThetaAndPhiHelicityFrameMyBinningH");
  // TH2F* fGenerH = (TH2F*)listings->FindObject("fMCCosThetaAndPhiHelicityFrameMyBinningH");
  TH2F* fReconHhe = (TH2F*)listings->FindObject("fCosThetaAndPhiHelicityFrameMyBinningReweightingH");
  TH2F* fGenerHhe = (TH2F*)listings->FindObject("fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH");
  fReconHhe->Sumw2();
  fGenerHhe->Sumw2();
  TH2F* fReconHcs = (TH2F*)listings->FindObject("fCosThetaAndPhiCsFrameMyBinningReweightingH");
  TH2F* fGenerHcs = (TH2F*)listings->FindObject("fMCCosThetaAndPhiCsFrameMyBinningReweightingH");
  fReconHcs->Sumw2();
  fGenerHcs->Sumw2();

  TH2F* AxEhe = (TH2F*) fReconHhe->Clone("AxEhe");
  AxEhe->Divide(fGenerHhe);
  new TCanvas;
  AxEhe->Draw("colZ");


  TH2F* AxEcs = (TH2F*) fReconHcs->Clone("AxEcs");
  AxEcs->Divide(fGenerHcs);
  new TCanvas;
  AxEcs->Draw("colZ");

}
