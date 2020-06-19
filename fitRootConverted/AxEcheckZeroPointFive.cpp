#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <math.h>
#include "TH2D.h"
#include "TF2.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"

#include "TDatime.h"

#include <vector>
#include <map>



//_____________________________________________________________________________
/* - Data needed for the global fit.
 * - They have to be global to be visible.
 */
// std::vector< std::pair< Double_t, Double_t > > coords;
std::vector< Double_t > coordsX;
std::vector< Double_t > coordsYpositiveCosTheta;
std::vector< Double_t > coordsYnegativeCosTheta;
std::vector< Double_t > valuespositiveCosTheta;
std::vector< Double_t > valuesnegativeCosTheta;
std::vector< Double_t > errorspositiveCosTheta;
std::vector< Double_t > errorsnegativeCosTheta;
std::vector< Double_t > coordsXmap;
std::vector< Double_t > coordsYmap;
std::vector< Double_t > valuesmap;
std::vector< Double_t > errorsmap;




void AxEcheck(){

  TDatime d;
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
  fReconH->Sumw2();
  fGenerH->Sumw2();
  TH2F*  AxEb   = (TH2F*) fReconH->Clone("AxEb");
  AxEb->Divide(fGenerH);
  AxEb->Sumw2();


  Int_t nBinsCosTheta2 = AxEb->GetNbinsX();
  Int_t nBinsPhi2      = AxEb->GetNbinsY();


  /// reset data structure
  coordsXmap = std::vector<Double_t>();
  coordsYmap = std::vector<Double_t>();
  valuesmap  = std::vector<Double_t>();
  errorsmap  = std::vector<Double_t>();
  for (Int_t ix = 1; ix <= nBinsCosTheta2; ++ix) {
    for (Int_t iy = 1; iy <= nBinsPhi2; ++iy) {
        coordsXmap.push_back( AxEb->GetXaxis()->GetBinCenter( ix ) );
        coordsYmap.push_back( AxEb->GetYaxis()->GetBinCenter( iy ) );
        valuesmap.push_back(  AxEb->GetBinContent( ix, iy )        );
        errorsmap.push_back(  AxEb->GetBinError(   ix, iy )        );
    }
  }



  Double_t CosThetaCenters[] = { -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5 };
  Double_t PhiCenters[] = { -3.14*19.5*0.05, -3.14*18.5*0.05, -3.14*17.5*0.05, -3.14*15*0.05,
                            -3.14*11*0.05,   -3.14*7.5*0.05,  -3.14*5*0.05,    -3.14*3*0.05,
                            -3.14*1.5*0.05,  -3.14*0.5*0.05,  +3.14*0.5*0.05,  +3.14*1.5*0.05,
                            +3.14*3*0.05,    +3.14*5*0.05,    +3.14*7.5*0.05,  +3.14*11*0.05,
                            +3.14*15*0.05,   +3.14*17.5*0.05, +3.14*18.5*0.05, +3.14*19.5*0.05 };


  Double_t MyVariableCosThetaBinning2[] = { -0.51, -0.35, -0.15, -0.05,
                                            0.05,  0.15,  0.35,  0.51 };
  Double_t MyVariablePhiBinning2[] = { -3.14*1,       -3.14*19*0.05, -3.14*18*0.05, -3.14*17*0.05,
                                      -3.14*13*0.05, -3.14*9*0.05,  -3.14*6*0.05,  -3.14*4*0.05,
                                      -3.14*2*0.05,  -3.14*1*0.05,   0,            +3.14*1*0.05,
                                      +3.14*2*0.05,  +3.14*4*0.05,  +3.14*6*0.05,  +3.14*9*0.05,
                                      +3.14*13*0.05, +3.14*17*0.05, +3.14*18*0.05, +3.14*19*0.05,
                                      +3.14*1 };

  TH2F* AxE = new TH2F( "AxE", "AxE", 7, MyVariableCosThetaBinning2, 20, MyVariablePhiBinning2 );
  Int_t BinX = 0;
  Int_t BinY = 0;
  Double_t PreviousError = 0;
  Double_t NewErrors[7][20];
  for ( Int_t ix = 0; ix < 7; ix++ ) {
    for ( Int_t iy = 0; iy < 20; iy++ ) {
      NewErrors[ix][iy] = 0;
    }
  }
  for ( Int_t i = 0; i < nBinsCosTheta2*nBinsPhi2; i++ ) {
    BinX          = 0;
    BinY          = 0;
    PreviousError = 0;
    AxE->Fill(coordsXmap[i], coordsYmap[i], valuesmap[i]);
    BinX          = AxE->GetXaxis()->FindBin(coordsXmap[i]);
    BinY          = AxE->GetYaxis()->FindBin(coordsYmap[i]);
    PreviousError = AxE->GetBinError( BinX, BinY );
    // AxE->SetBinError( BinX, BinY, TMath::Sqrt( PreviousError*PreviousError + errorsmap[i]*errorsmap[i] ) );
    AxE->SetBinError( BinX, BinY, 0. );
    NewErrors[BinX-1][BinY-1] += errorsmap[i]*errorsmap[i];
  }
  for ( Int_t ix = 0; ix < 7; ix++ ) {
    for ( Int_t iy = 0; iy < 20; iy++ ) {
      AxE->SetBinError( ix + 1, iy + 1, NewErrors[ix][iy] );
    }
  }
  AxE->Sumw2();
  new TCanvas;
  AxE->Draw("colZ");











  Int_t nBinsCosTheta  = AxE->GetNbinsX();
  Int_t nBinsPhi       = AxE->GetNbinsY();


  /// reset data structure
  coordsX                 = std::vector<Double_t>();
  coordsYnegativeCosTheta = std::vector<Double_t>();
  coordsYpositiveCosTheta = std::vector<Double_t>();
  valuesnegativeCosTheta  = std::vector<Double_t>();
  valuespositiveCosTheta  = std::vector<Double_t>();
  errorspositiveCosTheta  = std::vector<Double_t>();
  errorsnegativeCosTheta  = std::vector<Double_t>();
  /// fill data structure
  for (Int_t ix = 1; ix <= nBinsCosTheta; ++ix) {
    for (Int_t iy = 1; iy <= nBinsPhi; ++iy) {
      // coordsX.push_back( ((TAxis*) AxE->GetXaxis())->GetBinCenter( ix ) );
      if ( ix == 1 ) {
        coordsYnegativeCosTheta.push_back( AxE->GetYaxis()->GetBinCenter( iy ) );
        valuesnegativeCosTheta.push_back(  AxE->GetBinContent( ix, iy )        );
        errorsnegativeCosTheta.push_back(  AxE->GetBinError(   ix, iy )        );
      } else if ( ix == nBinsCosTheta ) {
        coordsYpositiveCosTheta.push_back( AxE->GetYaxis()->GetBinCenter( iy ) );
        valuespositiveCosTheta.push_back(  AxE->GetBinContent( ix, iy )        );
        errorspositiveCosTheta.push_back(  AxE->GetBinError(   ix, iy )        );
      }

    }
  }





  Double_t MyVariablePhiBinning[] = { -3.14*1,       -3.14*19*0.05, -3.14*18*0.05, -3.14*17*0.05,
                                      -3.14*13*0.05, -3.14*9*0.05,  -3.14*6*0.05,  -3.14*4*0.05,
                                      -3.14*2*0.05,  -3.14*1*0.05,   0,            +3.14*1*0.05,
                                      +3.14*2*0.05,  +3.14*4*0.05,  +3.14*6*0.05,  +3.14*9*0.05,
                                      +3.14*13*0.05, +3.14*17*0.05, +3.14*18*0.05, +3.14*19*0.05,
                                      +3.14*1 };

  TH1F* Negative = new TH1F( "Negative", "Negative", 20, MyVariablePhiBinning );
  TH1F* Positive = new TH1F( "Positive", "Positive", 20, MyVariablePhiBinning );

  for ( Int_t iBin = 1; iBin <= 20; iBin++ ) {

      Negative->SetBinContent( iBin, valuesnegativeCosTheta[iBin-1] );
      Positive->SetBinContent( iBin, valuespositiveCosTheta[iBin-1] );
      Negative->SetBinError(   iBin, errorsnegativeCosTheta[iBin-1] );
      Positive->SetBinError(   iBin, errorspositiveCosTheta[iBin-1] );
  }

  Negative->Sumw2();
  Positive->Sumw2();


  new TCanvas;
  Negative->SetLineColor(kRed);
  Negative->SetLineWidth(3);
  Negative->GetYaxis()->SetRangeUser(0., Negative->GetMaximum()*2.0);
  Negative->SetTitle( ";#phi; AxE" );
  Negative->Draw();
  Positive->SetLineColor(kBlue);
  Positive->SetLineWidth(3);
  Positive->Draw("same");
  gPad->BuildLegend();
  gPad->Modified();


  TH1F* PositiveOverNegative = (TH1F*) Positive->Clone("PositiveOverNegative");
  PositiveOverNegative->Divide(Negative);
  PositiveOverNegative->Sumw2();

  new TCanvas;
  PositiveOverNegative->Draw();
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  PositiveOverNegative->SetLineColor(kRed);
  PositiveOverNegative->SetLineWidth(3);
  PositiveOverNegative->GetXaxis()->SetTitleOffset(1.15);
  // PositiveOverNegative->GetYaxis()->SetTitleOffset(1.25);
  PositiveOverNegative->GetYaxis()->SetTitleOffset(1.);
  PositiveOverNegative->GetXaxis()->SetTitleSize(0.045);
  PositiveOverNegative->GetYaxis()->SetTitleSize(0.045);
  PositiveOverNegative->GetXaxis()->SetLabelSize(0.045);
  PositiveOverNegative->GetYaxis()->SetLabelSize(0.045);
  PositiveOverNegative->GetXaxis()->SetTitleFont(42);
  PositiveOverNegative->GetYaxis()->SetTitleFont(42);
  PositiveOverNegative->GetXaxis()->SetLabelFont(42);
  PositiveOverNegative->GetYaxis()->SetLabelFont(42);
  PositiveOverNegative->GetXaxis()->SetNdivisions(408);
  PositiveOverNegative->GetYaxis()->SetRangeUser(0., PositiveOverNegative->GetMaximum()*2.0);
  PositiveOverNegative->SetTitle( ";#phi; AxE ratio" );
  PositiveOverNegative->Draw("text");
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"Positive-Negative Ratio");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  // latex->DrawLatex(0.55,0.78,"Minuit 2D Fit");
  // latex->DrawLatex(0.55,0.70,Form("#lambda_{#theta} = %.3f #pm %.3f",     LambdaTheta,      LambdaThetaErr   ));
  // latex->DrawLatex(0.55,0.62,Form("#lambda_{#phi} = %.3f #pm %.3f",       LambdaPhi,        LambdaPhiErr     ));
  // latex->DrawLatex(0.55,0.54,Form("#lambda_{#theta#phi} = %.3f #pm %.3f", LambdaThetaPhi,   LambdaThetaPhiErr));
  // latex->DrawLatex(0.55,0.44,Form("#tilde{#chi} = %3.3f/%3.3d = %2.2f",      GlobalChi,        ndf, (Double_t)GlobalChi/(Double_t)ndf ));



}
