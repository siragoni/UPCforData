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
void PolarisationCorrectedDistribHe1D( Int_t selectionFlag = 0 ){

  TDatime d;
  TFile *file1D_T = new TFile(Form("pngResults/%d-%2.2d-%2.2d/Closure/PolarisationCorrectedHe1Dv2_RecTAxET.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile *file1D_L = new TFile(Form("pngResults/%d-%2.2d-%2.2d/Closure/PolarisationCorrectedHe1Dv2_RecTAxEL.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );

  TH1F* acceptanceCosTheta_T = (TH1F*) file1D_T->Get("acceptanceCosTheta");
  TH1F* acceptancePhi_T      = (TH1F*) file1D_T->Get("acceptancePhi");
  TH1F* acceptanceTildePhi_T = (TH1F*) file1D_T->Get("acceptanceTildePhi");

  TH1F* acceptanceCosTheta_L = (TH1F*) file1D_L->Get("acceptanceCosTheta");
  TH1F* acceptancePhi_L      = (TH1F*) file1D_L->Get("acceptancePhi");
  TH1F* acceptanceTildePhi_L = (TH1F*) file1D_L->Get("acceptanceTildePhi");


  TH1F* differenceCosTheta = (TH1F*) acceptanceCosTheta_T->Clone("differenceCosTheta");
  TH1F* differencePhi      = (TH1F*) acceptancePhi_T     ->Clone("differencePhi");
  TH1F* differenceTildePhi = (TH1F*) acceptanceTildePhi_T->Clone("differenceTildePhi");


  differenceCosTheta->Add(acceptanceCosTheta_L, -1.);
  differencePhi     ->Add(acceptancePhi_L,      -1.);
  differenceTildePhi->Add(acceptanceTildePhi_L, -1.);


  differenceCosTheta->Sumw2();
  differencePhi     ->Sumw2();
  differenceTildePhi->Sumw2();



  differenceCosTheta->Divide(acceptanceCosTheta_T);
  differencePhi     ->Divide(acceptancePhi_T);
  differenceTildePhi->Divide(acceptanceTildePhi_T);

  new TCanvas;
  differenceCosTheta->Draw();
  new TCanvas;
  differencePhi     ->Draw();
  new TCanvas;
  differenceTildePhi->Draw();








  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_19032021.root");
  // TFile* fileList2 = new TFile("AnalysisResultsLHC18l7_long_19032021.root"); // trial longitudinal
  TFile* fileList = new TFile("AnalysisResultsLHC18l7_trans_30032021.root");
  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_nonpol_rapidity_30032021.root");
  TFile* fileList3 = new TFile("AnalysisResultsLHC18l7_long_rapidity_30032021.root"); // trial longitudinal
  TFile* fileList2 = new TFile("AnalysisResultsLHC18l7_long_faultyevents_31032021.root"); // trial longitudinal
  // TFile* fileList3 = new TFile("AnalysisResultsLHC18l7_long_niceevents_31032021.root"); // trial longitudinal
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TDirectory* dir2 = fileList2->GetDirectory("MyTask");
  TDirectory* dir3 = fileList3->GetDirectory("MyTask");
  TList* listings;
  TList* listings2;
  TList* listings3;
  dir->GetObject("MyOutputContainer", listings);
  dir2->GetObject("MyOutputContainer", listings2);
  dir3->GetObject("MyOutputContainer", listings3);


  TH2F* fRecon2DH_T = (TH2F*)listings->FindObject("fCosThetaAndPhiHelicityFrameMyBinningH");
  TH2F* fGener2DH_T = (TH2F*)listings->FindObject("fMCCosThetaAndPhiHelicityFrameMyBinningH");
  fRecon2DH_T->Sumw2();
  fGener2DH_T->Sumw2();


  TH2F* fRecon2DH_L = (TH2F*)listings2->FindObject("fCosThetaAndPhiHelicityFrameMyBinningH");
  TH2F* fGener2DH_L = (TH2F*)listings2->FindObject("fMCCosThetaAndPhiHelicityFrameMyBinningH");
  fRecon2DH_L->Sumw2();
  fGener2DH_L->Sumw2();
  TH1F* fReconCosThetaH = (TH1F*)listings2->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH = (TH1F*)listings2->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  fReconCosThetaH->Sumw2();
  fGenerCosThetaH->Sumw2();
  TH1F* fReconCosThetaH_T = (TH1F*)listings->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH_T = (TH1F*)listings->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  fReconCosThetaH_T->Sumw2();
  fGenerCosThetaH_T->Sumw2();
  TH1F* fReconPhiH_T = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerPhiH_T = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsH");
  fReconPhiH_T->Sumw2();
  fGenerPhiH_T->Sumw2();
  // TH1F* fReconCosThetaH_nice = (TH1F*)listings3->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TFile* SavedFile = new TFile("SavedCosThetaLongitudinalAfterBkg.root");
  TH1F* fReconCosThetaH_nice = (TH1F*)SavedFile->Get("fCosThetaReconstructed");
  TH1F* fGenerCosThetaH_nice = (TH1F*)listings3->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  fReconCosThetaH_nice->Sumw2();
  fGenerCosThetaH_nice->Sumw2();
  TH1F* fGenerCosThetaH_nice_2 = (TH1F*)fGenerCosThetaH_nice->Clone("fGenerCosThetaH_nice_2");
  TH1F* fReconCosThetaH_2 = (TH1F*)fReconCosThetaH->Clone("fReconCosThetaH_2");
  TH1F* fReconCosThetaH_3 = (TH1F*)fReconCosThetaH->Clone("fReconCosThetaH_3");
  TH1F* fReconCosThetaH_4 = (TH1F*)fReconCosThetaH->Clone("fReconCosThetaH_4");
  fReconCosThetaH_2->Divide(fGenerCosThetaH);

  fReconCosThetaH     ->Divide(fReconCosThetaH_2);
  fReconCosThetaH_T   ->Divide(fGenerCosThetaH_T);
  fReconCosThetaH_3   ->Divide(fReconCosThetaH_T);
  fReconCosThetaH_nice->Divide(fGenerCosThetaH_nice);
  fReconCosThetaH_4   ->Divide(fReconCosThetaH_nice);


  TH2F* AxE_T = (TH2F*) fRecon2DH_T->Clone("AxE_T");
  TH2F* AxE_L = (TH2F*) fRecon2DH_L->Clone("AxE_L");
  AxE_T->Divide(fGener2DH_T);
  AxE_L->Divide(fGener2DH_L);

  new TCanvas;
  AxE_T->Draw("ColZ text");
  new TCanvas;
  AxE_L->Draw("ColZ text");


  TH2F* diffAxE = (TH2F*) AxE_T->Clone("diffAxE");
  diffAxE->Add(AxE_L, -1.);
  diffAxE->Divide(AxE_T);

  new TCanvas;
  diffAxE->Draw("ColZ text");


  Int_t nXbins = diffAxE->GetNbinsX();
  Int_t nYbins = diffAxE->GetNbinsY();


  TH1F* diffAxE_1D = new TH1F("diffAxE_1D","diffAxE_1D", nXbins*nYbins, 0., nXbins*nYbins );
  Int_t ibin = 1;
  for( Int_t ix = 1; ix < nXbins+1; ix++){
    for( Int_t iy = 1; iy < nYbins+1; iy++){
      diffAxE_1D->SetBinContent(ibin, diffAxE->GetBinContent(diffAxE->GetBin(ix, iy)));
      ibin++;
    }
  }
  new TCanvas;
  diffAxE_1D->Draw();



  new TCanvas;
  gPad->SetLogy();
  acceptanceCosTheta_T->Draw();
  acceptanceCosTheta_L->SetLineColor(kRed);
  acceptanceCosTheta_L->Draw("same");

  new TCanvas;
  gPad->SetLogy();
  acceptancePhi_T->Draw();
  acceptancePhi_L->SetLineColor(kRed);
  acceptancePhi_L->Draw("same");

  new TCanvas;
  gPad->SetLogy();
  acceptanceTildePhi_T->Draw();
  acceptanceTildePhi_L->SetLineColor(kRed);
  acceptanceTildePhi_L->Draw("same");


  new TCanvas;
  TH2F* ratioAxE = (TH2F*) AxE_T->Clone("diffAxE");
  ratioAxE->Divide(AxE_L);
  ratioAxE->Draw("ColZ text");

  new TCanvas;
  fReconCosThetaH     ->SetLineColor(kRed);
  fReconCosThetaH     ->SetLineWidth(3);
  fReconCosThetaH     ->Draw();
  new TCanvas;
  fReconCosThetaH_3   ->SetLineColor(kMagenta);
  fReconCosThetaH_3   ->SetLineWidth(3);
  fReconCosThetaH_3   ->Draw();
  new TCanvas;
  fReconCosThetaH_4   ->SetLineWidth(3);
  fReconCosThetaH_4   ->Draw();
  new TCanvas;
  fReconCosThetaH     ->Draw();
  fReconCosThetaH_3   ->Draw("same");
  fReconCosThetaH_4   ->Draw("same");
  // fGenerCosThetaH_nice_2->SetLineColor(kBlack);
  // fGenerCosThetaH_nice_2->SetLineWidth(3);
  // fGenerCosThetaH_nice_2->Draw("same");


  // new TCanvas;
  // TH1F* ratioAxE_1D_TvsL = (TH1F*) fReconCosThetaH_T->Clone("ratioAxE_1D_TvsL");
  // ratioAxE_1D_TvsL->Divide(fReconCosThetaH_nice);
  // ratioAxE_1D_TvsL->Draw();




  new TCanvas;
  fReconPhiH_T->Divide(fGenerPhiH_T);
  TH1F* ratioAxE_1D_TvsL_onlyPhi = (TH1F*) fReconPhiH_T->Clone("ratioAxE_1D_TvsL_onlyPhi");
  TFile* SavedFilePhi = new TFile("SavedPhiLongitudinalAfterBkg.root");
  TH1F* fReconPhiH_nice = (TH1F*)SavedFilePhi->Get("fPhiReconstructed");
  TH1F* fGenerPhiH_nice = (TH1F*)listings3->FindObject("fMCPhiHelicityFrameTwentyfiveBinsH");
  fReconPhiH_nice->Sumw2();
  fGenerPhiH_nice->Sumw2();
  fReconPhiH_nice->Divide(fGenerPhiH_nice);
  ratioAxE_1D_TvsL_onlyPhi->Divide(fReconPhiH_nice);
  ratioAxE_1D_TvsL_onlyPhi->Draw();

}






//_____________________________________________________________________________
/* - Final closure.
 * -
 */
void PolarisationCorrectedDistribHe1D_2( Int_t selectionFlag = 0 ){

  TDatime d;
  TFile *file1D_T = new TFile(Form("pngResults/%d-%2.2d-%2.2d/Closure/PolarisationCorrectedHe1Dv2_RecTAxET.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile *file1D_L = new TFile(Form("pngResults/%d-%2.2d-%2.2d/Closure/PolarisationCorrectedHe1Dv2_RecTAxEL.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );

  TH1F* acceptanceCosTheta_T = (TH1F*) file1D_T->Get("acceptanceCosTheta");
  TH1F* acceptancePhi_T      = (TH1F*) file1D_T->Get("acceptancePhi");
  TH1F* acceptanceTildePhi_T = (TH1F*) file1D_T->Get("acceptanceTildePhi");

  TH1F* acceptanceCosTheta_L = (TH1F*) file1D_L->Get("acceptanceCosTheta");
  TH1F* acceptancePhi_L      = (TH1F*) file1D_L->Get("acceptancePhi");
  TH1F* acceptanceTildePhi_L = (TH1F*) file1D_L->Get("acceptanceTildePhi");


  // TH1F* differenceCosTheta = (TH1F*) acceptanceCosTheta_T->Clone("differenceCosTheta");
  // TH1F* differencePhi      = (TH1F*) acceptancePhi_T     ->Clone("differencePhi");
  // TH1F* differenceTildePhi = (TH1F*) acceptanceTildePhi_T->Clone("differenceTildePhi");
  //
  //
  // differenceCosTheta->Add(acceptanceCosTheta_L, -1.);
  // differencePhi     ->Add(acceptancePhi_L,      -1.);
  // differenceTildePhi->Add(acceptanceTildePhi_L, -1.);
  //
  //
  // differenceCosTheta->Sumw2();
  // differencePhi     ->Sumw2();
  // differenceTildePhi->Sumw2();
  //
  //
  //
  // differenceCosTheta->Divide(acceptanceCosTheta_T);
  // differencePhi     ->Divide(acceptancePhi_T);
  // differenceTildePhi->Divide(acceptanceTildePhi_T);
  //
  // new TCanvas;
  // differenceCosTheta->Draw();
  // new TCanvas;
  // differencePhi     ->Draw();
  // new TCanvas;
  // differenceTildePhi->Draw();








  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_19032021.root");
  // TFile* fileList2 = new TFile("AnalysisResultsLHC18l7_long_19032021.root"); // trial longitudinal
  TFile* fileList = new TFile("AnalysisResultsLHC18l7_trans_30032021.root");
  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_nonpol_rapidity_30032021.root");
  // TFile* fileList3 = new TFile("AnalysisResultsLHC18l7_long_rapidity_30032021.root"); // trial longitudinal
  TFile* fileList2 = new TFile("AnalysisResultsLHC18l7_long_faultyevents_31032021.root"); // trial longitudinal
  TFile* fileList3 = new TFile("AnalysisResultsLHC18l7_long_niceevents_31032021.root"); // trial longitudinal
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TDirectory* dir2 = fileList2->GetDirectory("MyTask");
  TDirectory* dir3 = fileList3->GetDirectory("MyTask");
  TList* listings;
  TList* listings2;
  TList* listings3;
  dir->GetObject("MyOutputContainer", listings);
  dir2->GetObject("MyOutputContainer", listings2);
  dir3->GetObject("MyOutputContainer", listings3);


  // TH2F* fRecon2DH_T = (TH2F*)listings->FindObject("fCosThetaAndPhiHelicityFrameMyBinningH");
  // TH2F* fGener2DH_T = (TH2F*)listings->FindObject("fMCCosThetaAndPhiHelicityFrameMyBinningH");
  // fRecon2DH_T->Sumw2();
  // fGener2DH_T->Sumw2();
  //
  //
  // TH2F* fRecon2DH_L = (TH2F*)listings2->FindObject("fCosThetaAndPhiHelicityFrameMyBinningH");
  // TH2F* fGener2DH_L = (TH2F*)listings2->FindObject("fMCCosThetaAndPhiHelicityFrameMyBinningH");
  // fRecon2DH_L->Sumw2();
  // fGener2DH_L->Sumw2();
  // TH1F* fReconCosThetaH = (TH1F*)listings2->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  // TH1F* fGenerCosThetaH = (TH1F*)listings2->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fReconCosThetaH = (TH1F*)listings3->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH = (TH1F*)listings3->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  fReconCosThetaH->Sumw2();
  fGenerCosThetaH->Sumw2();
  TH1F* fReconCosThetaH_T = (TH1F*)listings->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TH1F* fGenerCosThetaH_T = (TH1F*)listings->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  fReconCosThetaH_T->Sumw2();
  fGenerCosThetaH_T->Sumw2();
  // TH1F* fReconPhiH_T = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsH");
  // TH1F* fGenerPhiH_T = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsH");
  // fReconPhiH_T->Sumw2();
  // fGenerPhiH_T->Sumw2();
  // TH1F* fReconCosThetaH_nice = (TH1F*)listings3->FindObject("fCosThetaHelicityFrameTwentyfiveBinsH");
  TFile* SavedFile = new TFile("SavedCosThetaLongitudinalAfterBkg.root");
  TH1F* fReconCosThetaH_nice = (TH1F*)SavedFile->Get("fCosThetaReconstructed");
  TH1F* fGenerCosThetaH_nice = (TH1F*)listings3->FindObject("fMCCosThetaHelicityFrameTwentyfiveBinsH");
  fReconCosThetaH_nice->Sumw2();
  fGenerCosThetaH_nice->Sumw2();
  TH1F* fGenerCosThetaH_nice_2 = (TH1F*)fGenerCosThetaH_nice->Clone("fGenerCosThetaH_nice_2");
  // TH1F* fReconCosThetaH_2 = (TH1F*)fReconCosThetaH->Clone("fReconCosThetaH_2");
  // TH1F* fReconCosThetaH_3 = (TH1F*)fReconCosThetaH->Clone("fReconCosThetaH_3");
  // TH1F* fReconCosThetaH_4 = (TH1F*)fReconCosThetaH->Clone("fReconCosThetaH_4");
  TH1F* fReconCosThetaH_2 = (TH1F*)fReconCosThetaH_nice->Clone("fReconCosThetaH_2");
  TH1F* fReconCosThetaH_3 = (TH1F*)fReconCosThetaH_nice->Clone("fReconCosThetaH_3");
  TH1F* fReconCosThetaH_4 = (TH1F*)fReconCosThetaH_nice->Clone("fReconCosThetaH_4");
  fReconCosThetaH_2->Divide(fGenerCosThetaH);

  fReconCosThetaH     ->Divide(fReconCosThetaH_2);
  fReconCosThetaH_T   ->Divide(fGenerCosThetaH_T);
  fReconCosThetaH_3   ->Divide(fReconCosThetaH_T);
  fReconCosThetaH_nice->Divide(fGenerCosThetaH_nice);
  fReconCosThetaH_4   ->Divide(fReconCosThetaH_nice);


  // TH2F* AxE_T = (TH2F*) fRecon2DH_T->Clone("AxE_T");
  // TH2F* AxE_L = (TH2F*) fRecon2DH_L->Clone("AxE_L");
  // AxE_T->Divide(fGener2DH_T);
  // AxE_L->Divide(fGener2DH_L);

  // new TCanvas;
  // AxE_T->Draw("ColZ text");
  // new TCanvas;
  // AxE_L->Draw("ColZ text");


  // TH2F* diffAxE = (TH2F*) AxE_T->Clone("diffAxE");
  // diffAxE->Add(AxE_L, -1.);
  // diffAxE->Divide(AxE_T);
  //
  // new TCanvas;
  // diffAxE->Draw("ColZ text");
  //
  //
  // Int_t nXbins = diffAxE->GetNbinsX();
  // Int_t nYbins = diffAxE->GetNbinsY();
  //
  //
  // TH1F* diffAxE_1D = new TH1F("diffAxE_1D","diffAxE_1D", nXbins*nYbins, 0., nXbins*nYbins );
  // Int_t ibin = 1;
  // for( Int_t ix = 1; ix < nXbins+1; ix++){
  //   for( Int_t iy = 1; iy < nYbins+1; iy++){
  //     diffAxE_1D->SetBinContent(ibin, diffAxE->GetBinContent(diffAxE->GetBin(ix, iy)));
  //     ibin++;
  //   }
  // }
  // new TCanvas;
  // diffAxE_1D->Draw();



  // new TCanvas;
  // gPad->SetLogy();
  // acceptanceCosTheta_T->Draw();
  // acceptanceCosTheta_L->SetLineColor(kRed);
  // acceptanceCosTheta_L->Draw("same");
  //
  // new TCanvas;
  // gPad->SetLogy();
  // acceptancePhi_T->Draw();
  // acceptancePhi_L->SetLineColor(kRed);
  // acceptancePhi_L->Draw("same");
  //
  // new TCanvas;
  // gPad->SetLogy();
  // acceptanceTildePhi_T->Draw();
  // acceptanceTildePhi_L->SetLineColor(kRed);
  // acceptanceTildePhi_L->Draw("same");


  // new TCanvas;
  // TH2F* ratioAxE = (TH2F*) AxE_T->Clone("diffAxE");
  // ratioAxE->Divide(AxE_L);
  // ratioAxE->Draw("ColZ text");

  new TCanvas;
  fReconCosThetaH     ->SetLineColor(kRed);
  fReconCosThetaH     ->SetLineWidth(3);
  fReconCosThetaH     ->Draw();
  new TCanvas;
  fReconCosThetaH_3   ->SetLineColor(kMagenta);
  fReconCosThetaH_3   ->SetLineWidth(3);
  fReconCosThetaH_3   ->Draw();
  new TCanvas;
  fReconCosThetaH_4   ->SetLineWidth(3);
  fReconCosThetaH_4   ->Draw();
  new TCanvas;
  fReconCosThetaH     ->Draw();
  fReconCosThetaH_3   ->Draw("same");
  fReconCosThetaH_4   ->Draw("same");
  fGenerCosThetaH_nice_2->SetLineColor(kBlack);
  fGenerCosThetaH_nice_2->SetLineWidth(3);
  fGenerCosThetaH_nice_2->Draw("same");


  // new TCanvas;
  // TH1F* ratioAxE_1D_TvsL = (TH1F*) fReconCosThetaH_T->Clone("ratioAxE_1D_TvsL");
  // ratioAxE_1D_TvsL->Divide(fReconCosThetaH_nice);
  // ratioAxE_1D_TvsL->Draw();




  // new TCanvas;
  // fReconPhiH_T->Divide(fGenerPhiH_T);
  // TH1F* ratioAxE_1D_TvsL_onlyPhi = (TH1F*) fReconPhiH_T->Clone("ratioAxE_1D_TvsL_onlyPhi");
  // TFile* SavedFilePhi = new TFile("SavedPhiLongitudinalAfterBkg.root");
  // TH1F* fReconPhiH_nice = (TH1F*)SavedFilePhi->Get("fPhiReconstructed");
  // TH1F* fGenerPhiH_nice = (TH1F*)listings3->FindObject("fMCPhiHelicityFrameTwentyfiveBinsH");
  // fReconPhiH_nice->Sumw2();
  // fGenerPhiH_nice->Sumw2();
  // fReconPhiH_nice->Divide(fGenerPhiH_nice);
  // ratioAxE_1D_TvsL_onlyPhi->Divide(fReconPhiH_nice);
  // ratioAxE_1D_TvsL_onlyPhi->Draw();

}
