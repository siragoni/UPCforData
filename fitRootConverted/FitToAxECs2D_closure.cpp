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
Double_t CosThetaFit( Double_t CosThetaToBeWeighted )
{
  Double_t val = 0;
  // Double_t par[13] = { 1.84270e-04, 1.35000e-02, 2.52572e+00, 4.49481e-01, 7.71374e-01,
  //                      4.10537e-01,-1.91331e-03,-4.33903e+00, 4.44258e-01,-7.64919e-01,
  //                      4.32631e-04, 1.48228e-02, 2.75000e+00 };
  Double_t par[13] = {-1.90267e-03, 1.64993e-02, 2.75000e+00, 4.63378e-01, 7.77394e-01,
                       4.37958e-01,-5.29215e-02,-4.89512e+00, 4.53152e-01,-7.93050e-01,
                      -5.20841e-04, 1.35817e-02, 2.74615e+00 };
  if (        CosThetaToBeWeighted < -0.550 ) {
    val = par[0] + par[1] * ( CosThetaToBeWeighted + 0.650 ) + par[2] * ( CosThetaToBeWeighted + 0.650 ) * ( CosThetaToBeWeighted + 0.650 );
  } else if ( CosThetaToBeWeighted < -0.125 ) {
    val = par[3] + par[4] * CosThetaToBeWeighted;
  } else if ( CosThetaToBeWeighted <  0.125 ) {
    val = par[5] + par[6] * CosThetaToBeWeighted + par[7] * CosThetaToBeWeighted * CosThetaToBeWeighted;
  } else if ( CosThetaToBeWeighted <  0.550 ) {
    val = par[8] + par[9] * CosThetaToBeWeighted;
  } else {
    val = par[10] + par[11] * ( CosThetaToBeWeighted - 0.650 ) + par[12] * ( CosThetaToBeWeighted - 0.650 ) * ( CosThetaToBeWeighted - 0.650 );
  }
  return val;

}
//_____________________________________________________________________________
Double_t fitf(Double_t *x, Double_t *par)
{
  // Double_t arg = 0;
  // if (par[2] != 0) arg = (x[0] - par[1])/par[2];

  Double_t fitval = par[0]*CosThetaFit(x[0])*(par[1]+par[2]*TMath::Cos(2*x[1]));
  return fitval;



}//_____________________________________________________________________________
/* - Fit function for the 2D AxE.
 * -
 */
void FitToAxECs(){

  // TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  TFile* fileList = new TFile("AnalysisResultsCohJpsiLocal.root");
  // TFile* fileList = new TFile("MCtrainResults/2020-06-25/kCohJpsiToMu/AnalysisResults.root");
  // TFile* fileList = new TFile("AnalysisResults_flatpolarisation_true_lumi.root");
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
  TH2F* fReconH = (TH2F*)listings->FindObject("fCosThetaAndPhiCsFrameMyBinningH");
  TH2F* fGenerH = (TH2F*)listings->FindObject("fMCCosThetaAndPhiCsFrameMyBinningH");
  fReconH->Sumw2();
  fGenerH->Sumw2();








  TCanvas* AcceptanceCanvas = new TCanvas("AcceptanceCanvas","AcceptanceCanvas",900,800);
  TH2F* acceptance = (TH2F*) fReconH->Clone("acceptance");
  acceptance->Divide(fGenerH);
  acceptance->Draw("ep");




  TF2 *fitting = new TF2("fitting",fitf,-0.65,0.65,-3.14,3.14,3);
  fitting->SetParameter(0, 1);
  fitting->SetParLimits(0, 0.0001, 100000);
  fitting->SetParameter(1, 0.0003);
  fitting->SetParLimits(1, 0.0000001, 0.0007);
  fitting->SetParameter(2, 1);
  fitting->SetParLimits(2, 0.0000000001, 100000);

  acceptance->Fit("fitting");
  Double_t chi = fitting->GetChisquare()/fitting->GetNDF();
  cout << chi << endl;
}
