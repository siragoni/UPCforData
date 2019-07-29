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
/* - Fit function for acceptance.
 * - I am using simple ROOT to make sigmoid fits to the plot.
 */
Double_t fDoubleSigmoid(Double_t* x,Double_t* par)
{
  /* - Par 0: normalisation
   * - Par 1: shift of each sigmoid.
   * -
   */
  Double_t val = 0;
  val = -1 * par[0] / ( 1 + TMath::Exp( -par[2] * TMath::Abs( x[0] ) + par[1] ) );
  return val;
}
//_____________________________________________________________________________
/* - Fit function for acceptance.
 * - I am using simple ROOT to make sigmoid fits to the plot.
 */
Double_t fTripleParabolic(Double_t* x,Double_t* par)
{
  /* - Par 0: normalisation
   * - Par 1: shift of each sigmoid.
   * -
   */
  Double_t val = 0;
  if (        x[0] < -0.550 ) {
    val = par[0] + par[1] * ( x[0] + 0.650 ) + par[2] * ( x[0] + 0.650 ) * ( x[0] + 0.650 );
    // if (  x[0]     > -0.560 ) {
    //   val = par[3] + par[4] * x[0];
    // }
  } else if ( x[0] < -0.125 ) {
    // if (  x[0]     < -0.530 ){
    //   val = par[0] + par[1] * ( x[0] + 0.650 ) + par[2] * ( x[0] + 0.650 ) * ( x[0] + 0.650 );
    // }
    val = par[3] + par[4] * x[0];
  } else if ( x[0] <  0.125 ) {
    val = par[5] + par[6] * x[0] + par[7] * x[0] * x[0];
  } else if ( x[0] <  0.550 ) {
    val = par[8] + par[9] * x[0];
  } else {
    val = par[10] + par[11] * ( x[0] - 0.650 ) + par[12] * ( x[0] - 0.650 ) * ( x[0] - 0.650 );
  }
  return val;
}
//_____________________________________________________________________________
/* - Fit function for acceptance.
 * - I am using simple ROOT to make sigmoid fits to the plot.
 */
Double_t fTripleParabolicToCheck(Double_t* x,Double_t* par)
{
  /* - Par 0: normalisation
   * - Par 1: shift of each sigmoid.
   * -
   */
  Double_t val = 0;
  if (        x[0] < -0.550 ) {
    val = par[0] + par[1] * ( x[0] + 0.650 ) + par[2] * ( x[0] + 0.650 ) * ( x[0] + 0.650 );
    // if (  x[0]     > -0.560 ) {
    //   val = par[3] + par[4] * x[0];
    // }
  } else if ( x[0] < -0.125 ) {
    // if (  x[0]     < -0.530 ){
    //   val = par[0] + par[1] * ( x[0] + 0.650 ) + par[2] * ( x[0] + 0.650 ) * ( x[0] + 0.650 );
    // }
    val = par[3] + par[4] * x[0];
  } else if ( x[0] <  0.125 ) {
    val = par[5] + par[6] * x[0] + par[7] * x[0] * x[0];
  } else if ( x[0] <  0.550 ) {
    val = par[8] + par[9] * x[0];
  } else {
    val = par[10] + par[11] * ( x[0] - 0.650 ) + par[12] * ( x[0] - 0.650 ) * ( x[0] - 0.650 );
  }
  return val;
}
//_____________________________________________________________________________
/* - Fit function for the acceptance.
 * - This way I can reweigh the data directly
 * - to account for the acceptance.
 * - Hence the acceptance binning effect
 * - would be less important...
 * -
 */
void ManyBinsAcceptanceForFit(){

  // TFile* fileList = new TFile("AnalysisResultsMC08072019.root");
  // TFile* fileList = new TFile("MCtrainResults/2019-06-08/kCohJpsiToMu/AnalysisResults.root");
  TFile* fileList = new TFile("MCtrainResults/2019-06-24/kCohJpsiToMu/AnalysisResults.root");
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
  TH1F* fReconCosThetaH = (TH1F*)listings->FindObject("fAngularDistribOfPositiveMuonRestFrameJPsiH");
  TH1F* fGenerCosThetaH = (TH1F*)listings->FindObject("fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH");
  fReconCosThetaH->Sumw2();
  fGenerCosThetaH->Sumw2();
  TH1F* ReconTheta = (TH1F*) fReconCosThetaH->Clone("ReconTheta");




  TCanvas* AcceptanceCanvasCosTheta = new TCanvas("AcceptanceCanvasCosTheta","AcceptanceCanvasCosTheta",900,800);
  TH1F* acceptanceCosTheta = (TH1F*) fReconCosThetaH->Clone("acceptanceCosTheta");
  acceptanceCosTheta->Divide(fGenerCosThetaH);
  acceptanceCosTheta->Draw("ep");


  TF1* fFitAcceptance = new TF1( "fFitAcceptance", fDoubleSigmoid, -0.5, 0.5, 3 );
  fFitAcceptance->SetParameter(0, 1.   );
  fFitAcceptance->SetParameter(1, 0.325);
  fFitAcceptance->SetParameter(2, 5);
  fFitAcceptance->SetNpx( acceptanceCosTheta->GetNbinsX() );

  TF1* fFitAcceptanceParab = new TF1( "fFitAcceptanceParab", fTripleParabolic, -0.65, 0.65, 13 );
  fFitAcceptanceParab->SetParameter(0, 0.    );
  fFitAcceptanceParab->SetParameter(1, 0.015 );
  fFitAcceptanceParab->SetParameter(2, 2.5   );
  fFitAcceptanceParab->SetParameter(3, 0.45  );
  fFitAcceptanceParab->SetParameter(4, 0.8   );
  fFitAcceptanceParab->SetParameter(5, 0.41  );
  // fFitAcceptanceParab->SetParameter(6, 0.01 );
  fFitAcceptanceParab->SetParameter(6, 0.    );
  // fFitAcceptanceParab->SetParameter(7,-2.4  );
  fFitAcceptanceParab->SetParameter(7,-4.3   );
  fFitAcceptanceParab->SetParameter(8, 0.44  );
  fFitAcceptanceParab->SetParameter(9,-0.76  );
  fFitAcceptanceParab->SetParameter(10, 0.   );
  fFitAcceptanceParab->SetParameter(11, 0.015);
  fFitAcceptanceParab->SetParameter(12, 2.5  );

  fFitAcceptanceParab->SetParLimits(0,  fFitAcceptanceParab->GetParameter(0) * 0.9,  fFitAcceptanceParab->GetParameter(0) * 1.1 );
  fFitAcceptanceParab->SetParLimits(1,  fFitAcceptanceParab->GetParameter(1) * 0.9,  fFitAcceptanceParab->GetParameter(1) * 1.1 );
  fFitAcceptanceParab->SetParLimits(2,  fFitAcceptanceParab->GetParameter(2) * 0.9,  fFitAcceptanceParab->GetParameter(2) * 1.1 );
  fFitAcceptanceParab->SetParLimits(3,  fFitAcceptanceParab->GetParameter(3) * 0.9,  fFitAcceptanceParab->GetParameter(3) * 1.1 );
  fFitAcceptanceParab->SetParLimits(4,  fFitAcceptanceParab->GetParameter(4) * 0.9,  fFitAcceptanceParab->GetParameter(4) * 1.1 );
  fFitAcceptanceParab->SetParLimits(5,  fFitAcceptanceParab->GetParameter(5) * 0.9,  fFitAcceptanceParab->GetParameter(5) * 1.1 );
  fFitAcceptanceParab->SetParLimits(6,  fFitAcceptanceParab->GetParameter(6) * 0.,  fFitAcceptanceParab->GetParameter(6) * 10 );
  fFitAcceptanceParab->SetParLimits(7,  fFitAcceptanceParab->GetParameter(7) * 2.0,  0 );
  fFitAcceptanceParab->SetParLimits(8,  fFitAcceptanceParab->GetParameter(8) * 0.9,  fFitAcceptanceParab->GetParameter(8) * 1.1 );
  fFitAcceptanceParab->SetParLimits(9,  fFitAcceptanceParab->GetParameter(9) * 1.1,  fFitAcceptanceParab->GetParameter(9) * 0.9 );
  fFitAcceptanceParab->SetParLimits(10, fFitAcceptanceParab->GetParameter(10) * 0.9, fFitAcceptanceParab->GetParameter(10) * 1.1 );
  fFitAcceptanceParab->SetParLimits(11, fFitAcceptanceParab->GetParameter(11) * 0.9, fFitAcceptanceParab->GetParameter(11) * 1.1 );
  fFitAcceptanceParab->SetParLimits(12, fFitAcceptanceParab->GetParameter(12) * 0.9, fFitAcceptanceParab->GetParameter(12) * 1.1 );

  fFitAcceptanceParab->SetNpx( acceptanceCosTheta->GetNbinsX() );

  // acceptanceCosTheta->Fit( fFitAcceptance, "LR" );
  acceptanceCosTheta->Fit( fFitAcceptanceParab, "LR" );
  // acceptanceCosTheta->Fit("pol9");


  new TCanvas;
  TF1* CheckFunction = new TF1( "CheckFunction", fTripleParabolicToCheck, -0.65, 0.65, 13 );
  for ( Int_t iSetPar = 0; iSetPar < 13; iSetPar++ ) {
    CheckFunction->FixParameter( iSetPar, fFitAcceptanceParab->GetParameter(iSetPar) );
  }
  acceptanceCosTheta->Draw();
  CheckFunction->Draw("same");

  // TFile f("pngResults/TH1corrMyBinningSmall.root", "recreate");
  // acceptanceCosTheta->Write();
  // RawCosThetaH      ->Write();
  // f.Close();
}
