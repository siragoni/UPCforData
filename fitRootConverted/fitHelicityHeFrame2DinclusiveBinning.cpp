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
#include <vector>


#include "TH2.h"


//_____________________________________________________________________________
/* - Codign in the fit functions.
   - The fit is modelled as the sum of 3 different functions:
   - 1) CosTheta only
   - 2) Phi      only
   - 3) Mix of the two
   -
 */
Double_t CosTheta(Double_t *x, Double_t *par) {
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t returnValue     = 1 + par[0] * CosSquaredTheta;
  return   returnValue;
}
//______________________________________________
Double_t Phi(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi     = TMath::Cos( 2 * x[1] );
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t SinSquaredTheta = 1 - CosSquaredTheta;

  Double_t returnValue = par[0] * SinSquaredTheta * CosOfTwoPhi;
  return   returnValue;
}
//______________________________________________
Double_t Mix(Double_t *x, Double_t *par) {
  Double_t CosPhi               = TMath::Cos( x[1] );
  Double_t CosSquaredTheta      = x[0]*x[0];
  Double_t SinSquaredTheta      = 1 - CosSquaredTheta;
  Double_t SinSquaredOfTwoTheta = 4 * CosSquaredTheta * SinSquaredTheta;
  Double_t SinOfTwoTheta        = TMath::Sqrt(SinSquaredOfTwoTheta);
  Double_t returnValue          = par[0] * SinOfTwoTheta * CosPhi;
  return   returnValue;
}
//______________________________________________
Double_t helicity2D(Double_t *x, Double_t *par) {
   Double_t *lambdaTheta    = &par[0];
   Double_t *lambdaPhi      = &par[1];
   Double_t *lambdaThetaPhi = &par[2];
   Double_t sumOfTheSubFits = CosTheta( x, lambdaTheta) + Phi( x, lambdaPhi ) + Mix( x, lambdaThetaPhi );
   Double_t FinalResult     = par[3] * sumOfTheSubFits / ( 3 + par[0] );
   return   FinalResult;
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void fitHelicityShape(const char* AnalysisName, const char* MonteCarloName){
  /* - There are three cases for the selectionFlag:
     - 1) = 0 ; this implies the traditional pt-integrated plot;
     - 2) = 1 ; this is instead the coherent component;
     - 3) = 2 ; this is the incoherent component;
     - 4) = 3 ; ******************* ;
     -
   */
  TFile* dataList = new TFile(AnalysisName);
  TDirectory* dirData = dataList->GetDirectory("MyTask");
  TFile* mcList = new TFile(MonteCarloName);
  TDirectory* dirMC = mcList->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
  TList* listingsData;
  dirData->GetObject("MyOutputContainer", listingsData);
  TList* listingsMC;
  dirMC  ->GetObject("MyOutputContainer", listingsMC);
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
  new TCanvas;
  TH2F *fAngularDistribOfPositiveMuonRestFrameJPsiH = (TH2F*)listingsData->FindObject("fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH");
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin2D(2,2);
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(16);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("colZ");
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->SaveAs("rawDistr.png", "RECREATE");


  new TCanvas;
  TH2F *fMCAngularDistribOfPositiveMuonRestFrameJPsiH = (TH2F*)listingsMC->FindObject("fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH");
  fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  // fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin2D(2,2);
  // fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(16);
  fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("colZ");
  // fMCAngularDistribOfPositiveMuonRestFrameJPsiH->SaveAs("mcReconDistr.png", "RECREATE");


  new TCanvas;
  TH2F *fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH = (TH2F*)listingsMC->FindObject("fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH");
  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Sumw2();
  // fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Rebin2D(2,2);
  // fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Rebin(16);
  fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->Draw("colZ");
  // fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH->SaveAs("mcGenDistr.png", "RECREATE");


  new TCanvas;
  TH2F *acceptance = (TH2F*)fMCAngularDistribOfPositiveMuonRestFrameJPsiH->Clone("acceptance");
  acceptance->Sumw2();
  acceptance->Divide(fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH);
  acceptance->SetMarkerStyle(21);
  acceptance->SetLineColor(kBlue+1);
  acceptance->SetLineWidth(2);
  acceptance->Draw("ep");
  // acceptance->SaveAs("acceptance.png", "RECREATE");


  new TCanvas;
  // TH1F *fCorrectedShape = (TH1F*)fAngularDistribOfPositiveMuonRestFrameJPsiH->Clone("fCorrectedShape");
  // fCorrectedShape->Sumw2();
  // fCorrectedShape->Divide(acceptance);
  // fCorrectedShape->Rebin(2);
  // fCorrectedShape->SetMarkerStyle(21);
  // fCorrectedShape->SetLineColor(kBlue+1);
  // fCorrectedShape->SetLineWidth(2);
  // fCorrectedShape->Draw("ep");
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Sumw2();
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Divide(acceptance);
  // fAngularDistribOfPositiveMuonRestFrameJPsiH->Rebin(5);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->SetMarkerStyle(21);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->SetLineColor(kBlue+1);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->SetLineWidth(2);
  fAngularDistribOfPositiveMuonRestFrameJPsiH->Draw("ep");
  TH2F *fCorrectedShape = fAngularDistribOfPositiveMuonRestFrameJPsiH;

  /* - Do the parabolic fit with Root only.
     - Retrieve the parameters later for Roofit plotting.
     -
     - dN/dPhi = 1/2pi [ 1 - 4*rho_{1,-1}*Cos2Phi ]
  */
  // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-1, 1, -4, 4, 3);
  // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-0.5, 0.5, -3, 3, 3);
  /* - The global normalisation SHOULD be a parameter too!!!
     -
   */
  // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-1, 1, -4, 4, 3);
  // TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-0.9, 0.9, -3, 3, 3);
  /* - The global normalisation SHOULD be a parameter too!!!
    -
   */
  TF2* ParabolicFit = new TF2("ParabolicFit",helicity2D,-0.9, 0.9, -3, 3, 4);
  ParabolicFit->SetNpx(1000);
  ParabolicFit->SetParameter(0, 1);
  ParabolicFit->SetParameter(1, 0);
  ParabolicFit->SetParameter(2, 0);
  ParabolicFit->SetParameter(3, 1);
  ParabolicFit->SetParLimits(0, -2, +2);
  ParabolicFit->SetParLimits(1, -2, +2);
  ParabolicFit->SetParLimits(2, -2, +2);
  fCorrectedShape->Fit( "ParabolicFit" );
  // fCorrectedShape->Fit( ParabolicFit,"","", -0.5, 0.5, -2, 2 );

  fCorrectedShape->SetLineColor(kBlue);
  fCorrectedShape->SetLineStyle(kSolid);
  fCorrectedShape->SetLineWidth(3);
  fCorrectedShape->SetMarkerStyle(kFullCircle);
  fCorrectedShape->SetMarkerSize(1);
  fCorrectedShape->GetXaxis()->SetTitle("cos(#theta)");
  fCorrectedShape->GetYaxis()->SetTitle("#phi");
  // fCorrectedShape->GetYaxis()->SetTitle( Form( "Counts / (%.3f a.u.)",
  //                                              fCorrectedShape->GetXaxis()->GetBinWidth(1)
  //                                              )
  //                                             );
  fCorrectedShape->SetTitle("");
  new TCanvas;
  TCanvas* ZNAEnergy = new TCanvas( "CosTheta", "CosTheta", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  /* - Beautifying is starting now.
     -
   */
  fCorrectedShape->GetXaxis()->SetTitleOffset(1.25);
  fCorrectedShape->GetYaxis()->SetTitleOffset(1.25);
  // fCorrectedShape->GetYaxis()->SetTitleOffset(1.45);
  fCorrectedShape->GetXaxis()->SetTitleSize(0.045);
  fCorrectedShape->GetYaxis()->SetTitleSize(0.045);
  fCorrectedShape->GetXaxis()->SetLabelSize(0.045);
  fCorrectedShape->GetYaxis()->SetLabelSize(0.045);
  fCorrectedShape->GetXaxis()->SetTitleFont(42);
  fCorrectedShape->GetYaxis()->SetTitleFont(42);
  fCorrectedShape->GetXaxis()->SetLabelFont(42);
  fCorrectedShape->GetYaxis()->SetLabelFont(42);
  fCorrectedShape->GetXaxis()->SetNdivisions(408);
  fCorrectedShape->GetYaxis()->SetRangeUser(-3.1, 3.1);
  // fCorrectedShape->GetYaxis()->SetRangeUser(-2.5, 2.5);
  // fCorrectedShape->GetXaxis()->SetRangeUser(-1, 1);
  fCorrectedShape->GetXaxis()->SetRangeUser(-1, 1);
  // fCorrectedShape->GetYaxis()->SetRangeUser(5, fCorrectedShape->GetMaximum()*10.);
  // gPad ->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fCorrectedShape->Draw("colZ  SAME");
  // fCorrectedShape->Draw("surf3  ");
  ParabolicFit   ->Draw("cont1 SAME");
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, LHC18qr");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  latex->DrawLatex(0.55,0.78,Form("%.1f < y < %.1f",-4.0,-2.5));
  latex->DrawLatex(0.55,0.72,Form("#lambda_{#theta} = %.4f #pm %.4f",
                                  ParabolicFit->GetParameter(0),
                                  ParabolicFit->GetParError(0)
                                  )
                                );
  latex->DrawLatex(0.55,0.66,Form("#lambda_{#phi} = %.4f #pm %.4f",
                                  ParabolicFit->GetParameter(1),
                                  ParabolicFit->GetParError(1)
                                  )
                                );
  latex->DrawLatex(0.55,0.60,Form("#lambda_{#theta#phi} = %.4f #pm %.4f",
                                  ParabolicFit->GetParameter(2),
                                  ParabolicFit->GetParError(2)
                                  )
                                );
  latex->DrawLatex(0.55,0.18,Form( "      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                   ParabolicFit->GetChisquare(),
                                   ParabolicFit->GetNDF(),
                                   ParabolicFit->GetChisquare()/ParabolicFit->GetNDF()
                                  )
                                 );


  gPad->SaveAs("pngResults/fit2DinclusiveBinning.png",         "RECREATE");

  new TCanvas;
  acceptance->Draw("ep");
  gPad->SaveAs("pngResults/fit2DinclusiveBinning_ACCxEFF.png", "RECREATE");

}
