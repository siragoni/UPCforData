#include <TError.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include "TUnfoldDensity.h"
// #define VERBOSE_LCURVE_SCAN
using namespace std;
TRandom *rnd=0;
TH2 *gHistInvEMatrix;
TVirtualFitter *gFitter=0;
void chisquare_corr(Int_t &npar, Double_t * /*gin */, Double_t &f, Double_t *u, Int_t /*flag */) {
  //  Minimization function for H1s using a Chisquare method
  //  only one-dimensional histograms are supported
  //  Corelated errors are taken from an external inverse covariance matrix
  //  stored in a 2-dimensional histogram
  Double_t x;
  TH1 *hfit = (TH1*)gFitter->GetObjectFit();
  TF1 *f1   = (TF1*)gFitter->GetUserFunc();
  f1->InitArgs(&x,u);
  npar = f1->GetNpar();
  f = 0;
  Int_t npfit = 0;
  Int_t nPoints=hfit->GetNbinsX();
  Double_t *df=new Double_t[nPoints];
  for (Int_t i=0;i<nPoints;i++) {
    x     = hfit->GetBinCenter(i+1);
    TF1::RejectPoint(kFALSE);
    df[i] = f1->EvalPar(&x,u)-hfit->GetBinContent(i+1);
    if (TF1::RejectedPoint()) df[i]=0.0;
    else npfit++;
  }
  for (Int_t i=0;i<nPoints;i++) {
    for (Int_t j=0;j<nPoints;j++) {
      f += df[i]*df[j]*gHistInvEMatrix->GetBinContent(i+1,j+1);
    }
  }
  delete[] df;
  f1->SetNumberFitPoints(npfit);
}











int testUnfold1()
{
  // switch on histogram errors
  TH1::SetDefaultSumw2();
  // show fit result
  gStyle->SetOptFit(1111);

  // TFile* fileList = new TFile("AnalysisResultsLHC18l7_19032021.root");
  TFile* fileList = new TFile("AnalysisResultsLHC18l7_coherent_30112021.root");
  // TFile* fileList = new TFile("AnalysisResults.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);
  //============================================
  // generate MC distribution
  //
  TH1F *histMgenMC    = (TH1F*)listings->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  TH1F *histMdetMC    = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  TH2F *histMdetGenMC = (TH2F*)listings->FindObject("fPhiRecVsGenHelicityH");
  // generate data distribution
  //
  // TH1D *histMgenData = new TH1D("MgenData",";mass(gen)",nGen,xminGen,xmaxGen);
  TFile* fileDataRawPhi = new TFile("pngResults/2021-09-21/PhiHEv2/PhiHeFrameV2.root");
  // TH1F *histMdetData    = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  TH1F *histMdetData    = (TH1F*)listings->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");

  histMdetGenMC->Rebin2D(1,2);
  // histMdetGenMC->Rebin2D(1,3);
  Double_t integral_data = 0;
  for (size_t i = 0; i < 25; i++) {
    integral_data += histMdetData->GetBinContent(i+1);
  }
  histMdetGenMC->Scale(integral_data/histMdetGenMC->GetEntries());
  for (size_t i = 0; i < 25; i++) {
    for (size_t j = 0; j < 24; j++) {
      histMdetGenMC->SetBinContent(i+1, j+1, histMdetGenMC->GetBinContent(i+1, j+1)*integral_data/histMdetGenMC->GetEntries());
      histMdetGenMC->SetBinError(i+1, j+1, histMdetGenMC->GetBinError(i+1, j+1)*integral_data/histMdetGenMC->GetEntries());
    }
  }
  //=========================================================================
  // divide by bin withd to get density distributions
  // TH1D *histDensityGenData=new TH1D("DensityGenData",";mass(gen)",
  //                                   nGen,xminGen,xmaxGen);
  TH1F *histDensityGenMC=new TH1F("DensityGenMC",";mass(gen)",
                                    25, 0., 2.*TMath::Pi() );
  for(Int_t i=1;i<=25;i++) {
     // histDensityGenData->SetBinContent(i,histMgenData->GetBinContent(i)/
     //                                   histMgenData->GetBinWidth(i));
     histDensityGenMC->SetBinContent(i,histMgenMC->GetBinContent(i)/
                                       histMgenMC->GetBinWidth(i));
  }
  //=========================================================================
  // set up the unfolding
  // define migration matrix
  TUnfoldDensity unfold(histMdetGenMC,TUnfold::kHistMapOutputVert);
  // define input and bias scame
  // do not use the bias, because MC peak may be at the wrong place
  // watch out for error codes returned by the SetInput method
  // errors larger or equal 10000 are fatal:
  // the data points specified as input are not sufficient to constrain the
  // unfolding process
  if(unfold.SetInput(histMdetData)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }
  //========================================================================
  // the unfolding is done here
  //
  // scan L curve and find best point
  Int_t nScan=30;
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  // if required, report Info messages (for debugging the L-curve scan)
#ifdef VERBOSE_LCURVE_SCAN
  Int_t oldinfo=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kInfo;
#endif
  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
  iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
  // if required, switch to previous log-level
#ifdef VERBOSE_LCURVE_SCAN
  gErrorIgnoreLevel=oldinfo;
#endif
  //==========================================================================
  // define a correlated systematic error
  // for example, assume there is a 10% correlated error for all reconstructed
  // masses larger than 7
  Double_t SYS_ERROR1_MSTART=6;
  Double_t SYS_ERROR1_SIZE=0.1;
  TH2D *histMdetGenSys1=new TH2D("Mdetgensys1",";mass(det);mass(gen)",
                                 25, 0., 2.*TMath::Pi(), 24, 0., 2.*TMath::Pi() );
  for(Int_t i=0;i<=25+1;i++) {
     if(histMdetData->GetBinCenter(i)>=SYS_ERROR1_MSTART) {
        for(Int_t j=0;j<=24+1;j++) {
           histMdetGenSys1->SetBinContent(i,j,SYS_ERROR1_SIZE);
        }
     }
  }
  // unfold.AddSysError(histMdetGenSysMC,"SYSERROR_MC",TUnfold::kHistMapOutputVert,
  //                    TUnfoldSys::kSysErrModeMatrix);
  unfold.AddSysError(histMdetGenSys1,"SYSERROR1",TUnfold::kHistMapOutputVert,
                     TUnfoldSys::kSysErrModeRelative);
  //==========================================================================
  // print some results
  //
  std::cout<<"tau="<<unfold.GetTau()<<"\n";
  std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
           <<" / "<<unfold.GetNdf()<<"\n";
  std::cout<<"chi**2(sys)="<<unfold.GetChi2Sys()<<"\n";
  //==========================================================================
  // create graphs with one point to visualize the best choice of tau
  //
  Double_t t[1],x[1],y[1];
  logTauX->GetKnot(iBest,t[0],x[0]);
  logTauY->GetKnot(iBest,t[0],y[0]);
  TGraph *bestLcurve=new TGraph(1,x,y);
  TGraph *bestLogTauLogChi2=new TGraph(1,t,x);
  //==========================================================================
  // retreive results into histograms
  // get unfolded distribution
  TH1 *histMunfold=unfold.GetOutput("Unfolded");
  // get unfolding result, folded back
  TH1 *histMdetFold=unfold.GetFoldedOutput("FoldedBack");
  // get error matrix (input distribution [stat] errors only)
  // TH2D *histEmatData=unfold.GetEmatrix("EmatData");
  // get total error matrix:
  //   migration matrix uncorrelated and correlated systematic errors
  //   added in quadrature to the data statistical errors
  TH2 *histEmatTotal=unfold.GetEmatrixTotal("EmatTotal");
  // create data histogram with the total errors
  TH1D *histTotalError=
     new TH1D("TotalError",";mass(gen)",25, 0., 2.*TMath::Pi());
  for(Int_t bin=1;bin<=25;bin++) {
    histTotalError->SetBinContent(bin,histMunfold->GetBinContent(bin));
    histTotalError->SetBinError
       (bin,TMath::Sqrt(histEmatTotal->GetBinContent(bin,bin)));
  }
  // get global correlation coefficients
  // for this calculation one has to specify whether the
  // underflow/overflow bins are included or not
  // default: include all bins
  // here: exclude underflow and overflow bins
  TH1 *histRhoi=unfold.GetRhoItotal("rho_I",
                                    0, // use default title
                                    0, // all distributions
                                    "*[UO]", // discard underflow and overflow bins on all axes
                                    kTRUE, // use original binning
                                    &gHistInvEMatrix // store inverse of error matrix
                                    );
  //=====================================================================
  // plot some histograms
  TCanvas output;
  output.Divide(3,2);
  // Show the matrix which connects input and output
  // There are overflow bins at the bottom, not shown in the plot
  // These contain the background shape.
  // The overflow bins to the left and right contain
  // events which are not reconstructed. These are necessary for proper MC
  // normalisation
  output.cd(1);
  histMdetGenMC->Draw("BOX");
  // draw generator-level distribution:
  //   data (red) [for real data this is not available]
  //   MC input (black) [with completely wrong peak position and shape]
  //   unfolded data (blue)
  output.cd(2);
  // histTotalError->GetYaxis()->SetRangeUser(-20000, 20000);
  histTotalError->GetYaxis()->SetRangeUser(100000, 500000);
  histTotalError->SetLineColor(kBlue);
  histTotalError->Draw("E");
  histMunfold->SetLineColor(kGreen);
  histMunfold->Draw("SAME E1");
  // histDensityGenData->SetLineColor(kRed);
  // histDensityGenData->Draw("SAME");
  // histDensityGenMC->SetLineColor(kRed);
  // histDensityGenMC->Draw("SAME HIST");
  histMgenMC->SetLineColor(kRed);
  histMgenMC->Draw("SAME HIST");
  // show detector level distributions
  //    data (red)
  //    MC (black) [with completely wrong peak position and shape]
  //    unfolded data (blue)
  output.cd(3);
  histMdetFold->SetLineColor(kBlue);
  histMdetFold->Draw();
  histMdetMC->Draw("SAME HIST");
  TH1 *histInput=unfold.GetInput("Minput",";mass(det)");
  histInput->SetLineColor(kRed);
  histInput->Draw("SAME");
  // show correlation coefficients
  output.cd(4);
  histRhoi->Draw();
  // show tau as a function of chi**2
  output.cd(5);
  logTauX->Draw();
  bestLogTauLogChi2->SetMarkerColor(kRed);
  bestLogTauLogChi2->Draw("*");
  // show the L curve
  output.cd(6);
  lCurve->Draw("AL");
  bestLcurve->SetMarkerColor(kRed);
  bestLcurve->Draw("*");
  output.SaveAs("testUnfold1.ps");
  return 0;
}
