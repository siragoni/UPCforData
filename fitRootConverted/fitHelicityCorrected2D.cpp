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
/* - Fit function for the ZNC plots.
 * -
 */
void fitHelicityCorrected2D(){

  // TFile* fileList = new TFile("AnalysisResultsMC30052019.root");
  TFile* fileList = new TFile("MCtrainResults/2019-06-16-LHC18qr/kCohJpsiToMu/AnalysisResults.root");
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
  TH2F* fReconH = (TH2F*)listings->FindObject("fInvariantMassDistributionForSignalExtractionHelicityFrameH");
  TH2F* fGenerH = (TH2F*)listings->FindObject("fMCInvariantMassDistributionForSignalExtractionHelicityFrameH");
  fReconH->Sumw2();
  fGenerH->Sumw2();

  // gSystem->cd("pngResults/2019-05-30-18qr15o-NoSPD/SignalExtraction/");

  // TFile* fileDataRaw = new TFile("pngResults/2019-05-30-18qr15o-NoSPD/SignalExtraction/TH2signalEX.root");
  TFile* fileDataRaw = new TFile("pngResults/2019-06-19/SignalExtraction/TH2signalEX.root");
  TH2F* helicity2DafterSignalExtractionErrors = (TH2F*)fileDataRaw->Get("helicity2DafterSignalExtractionErrors");
  helicity2DafterSignalExtractionErrors->Sumw2();
  TH2F* RawH = (TH2F*) helicity2DafterSignalExtractionErrors->Clone("RawH");


  TCanvas* AcceptanceCanvas = new TCanvas("AcceptanceCanvas","AcceptanceCanvas",900,800);
  TH2F* acceptance = (TH2F*) fReconH->Clone("acceptance");
  acceptance->Divide(fGenerH);
  acceptance->Draw("ep");

  TCanvas* CorrCanvas = new TCanvas("CorrCanvas","CorrCanvas",900,800);
  RawH->Divide(acceptance);
  RawH->Draw("colZ");


  TFile f("pngResults/TH2corr.root", "new");
  acceptance->Write();
  RawH      ->Write();
  f.Close();
}
