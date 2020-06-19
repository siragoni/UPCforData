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
void PolarisationCorrectedDistribCs2D( Int_t RangeSelectionMode = 0 ){

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
  TH2F* fReconH = (TH2F*)listings->FindObject("fCosThetaAndPhiCsFrameMyBinningReweightingH");
  TH2F* fGenerH = (TH2F*)listings->FindObject("fMCCosThetaAndPhiCsFrameMyBinningReweightingH");
  fReconH->Rebin2D(4,4);
  fGenerH->Rebin2D(4,4);
  fReconH->Sumw2();
  fGenerH->Sumw2();
  new TCanvas;
  fReconH->Draw("colZ");
  new TCanvas;
  fGenerH->Draw("colZ");
  new TCanvas;


  TDatime d;
  // TFile* fileDataRaw = new TFile(Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* fileDataRaw = 0x0;
  if (        RangeSelectionMode == 0 ) {
    fileDataRaw = new TFile( Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 1 ) {
    fileDataRaw = new TFile( Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs_1.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 2 ) {
    fileDataRaw = new TFile( Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 3 ) {
    fileDataRaw = new TFile( Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 4 ) {
    fileDataRaw = new TFile( Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( RangeSelectionMode == 5 ) {
    fileDataRaw = new TFile( Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else {
    fileDataRaw = new TFile( Form("pngResults/%d-%2.2d-%2.2d/2DCS/Polarisation2DCs.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  }




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

  // TFile f("pngResults/PolarisationCorrectedCs2D.root", "recreate");
  // acceptance->Write();
  // RawH      ->Write();
  // // AccErrors ->Write();
  // f.Close();
  if        ( RangeSelectionMode == 0 ) {
    TFile f("pngResults/PolarisationCorrectedCs2D_flat.root", "recreate");
    acceptance->Write();
    RawH      ->Write();
    // AccErrors ->Write();
    f.Close();
  } else if ( RangeSelectionMode == 1 ) {
    TFile f("pngResults/PolarisationCorrectedCs2D_1.root", "recreate");
    acceptance->Write();
    RawH      ->Write();
    // AccErrors ->Write();
    f.Close();
  } else if ( RangeSelectionMode == 2 ) {
    TFile f("pngResults/PolarisationCorrectedCs2D_2.root", "recreate");
    acceptance->Write();
    RawH      ->Write();
    // AccErrors ->Write();
    f.Close();
  } else if ( RangeSelectionMode == 3 ) {
    TFile f("pngResults/PolarisationCorrectedCs2D_3.root", "recreate");
    acceptance->Write();
    RawH      ->Write();
    // AccErrors ->Write();
    f.Close();
  } else if ( RangeSelectionMode == 4 ) {
    TFile f("pngResults/PolarisationCorrectedCs2D_4.root", "recreate");
    acceptance->Write();
    RawH      ->Write();
    // AccErrors ->Write();
    f.Close();
  } else if ( RangeSelectionMode == 5 ) {
    TFile f("pngResults/PolarisationCorrectedCs2D_5.root", "recreate");
    acceptance->Write();
    RawH      ->Write();
    // AccErrors ->Write();
    f.Close();
  } else {
    TFile f("pngResults/PolarisationCorrectedCs2D.root", "recreate");
    acceptance->Write();
    RawH      ->Write();
    // AccErrors ->Write();
    f.Close();
  }

}
