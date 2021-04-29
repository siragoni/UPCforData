#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TStyle.h"
using namespace std;
#include <vector>


//_____________________________________________________________________________
/* -
 * -
 */
void Ratio(){

  TDatime d;
  // TFile* file1D  = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2_AxE.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* file1D_ = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2.root",       d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* file1D  = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedCs1Dv2_AxE.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* file1D_ = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedCs1Dv2.root",       d.GetYear(), d.GetMonth(), d.GetDay() ) );

  // TH1F* CorrectedCosTheta  = (TH1F*) file1D ->Get("CorrCosThetaH");
  // TH1F* CorrectedCosTheta_ = (TH1F*) file1D_->Get("CorrCosThetaH");
  TH1F* CorrectedCosTheta  = (TH1F*) file1D ->Get("acceptanceCosTheta");
  TH1F* CorrectedCosTheta_ = (TH1F*) file1D_->Get("acceptanceCosTheta");
  CorrectedCosTheta ->Sumw2();
  CorrectedCosTheta_->Sumw2();

  // CorrectedCosTheta_->Divide(CorrectedCosTheta);
  new TCanvas;
  CorrectedCosTheta_->Draw();



  Double_t IntegralCosTheta  = 0.;
  Double_t IntegralCosTheta_ = 0.;
  Double_t ErrorCosTheta     = 0.;
  Double_t ErrorCosTheta_    = 0.;


  for (size_t i = 1; i < 25; i++) {

    if ( CorrectedCosTheta->GetBinContent(i) > 0.0001 ) {
      IntegralCosTheta += CorrectedCosTheta->GetBinContent(i);
      ErrorCosTheta    += CorrectedCosTheta->GetBinError(i);
    }
    if ( CorrectedCosTheta_->GetBinContent(i) > 0.0001 ) {
      IntegralCosTheta_ += CorrectedCosTheta_->GetBinContent(i);
      ErrorCosTheta_    += CorrectedCosTheta_->GetBinError(i);
    }


  }



  cout << "IntegralCosTheta  = " << IntegralCosTheta  << " +/- " << ErrorCosTheta  << endl;
  cout << "IntegralCosTheta_ = " << IntegralCosTheta_ << " +/- " << ErrorCosTheta_ << endl;
  // cout << "Ratio             = " << IntegralCosTheta_/IntegralCosTheta  << endl;
  cout << "Ratio             = " << IntegralCosTheta/IntegralCosTheta_  << endl;

}
