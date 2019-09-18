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
void PolarisationCorrectedDistribCs2D(){

  TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
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
  TH2F* fReconH = (TH2F*)listings->FindObject("fCosThetaAndPhiHelicityFrameMyBinningH");
  TH2F* fGenerH = (TH2F*)listings->FindObject("fMCCosThetaAndPhiHelicityFrameMyBinningH");
  // TH2F* fReconH = (TH2F*)listings->FindObject("fCosThetaAndPhiCsFrameMyBinningH");
  // TH2F* fGenerH = (TH2F*)listings->FindObject("fMCCosThetaAndPhiCsFrameMyBinningH");
  fReconH->Sumw2();
  fGenerH->Sumw2();

  TFile* fileDataRaw = new TFile("pngResults/2019-09-18/2DCS/Polarisation2DCs.root");
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
  TCanvas* CorrCanvas2 = new TCanvas("CorrCanvas2","CorrCanvas2",900,800);
  RawH->Draw("ep");

  // TH1F* AccErrors = new TH1F("AccErrors", "AccErrors", 40, -1, 1);
  // for( Int_t iLoop = 1; iLoop <= 40; iLoop++ ){
  //   Double_t AcceptanceEntries = fReconH   ->GetBinContent(iLoop);
  //   Double_t AcceptanceValue   = acceptance->GetBinContent(iLoop);
  //   if( AcceptanceEntries >= 1  ){
  //     AccErrors->SetBinContent( iLoop, AcceptanceValue/TMath::Sqrt(AcceptanceEntries) );
  //   }
  // }

  TFile f("pngResults/PolarisationCorrectedCs2D.root", "recreate");
  acceptance->Write();
  RawH      ->Write();
  // AccErrors ->Write();
  f.Close();
}
