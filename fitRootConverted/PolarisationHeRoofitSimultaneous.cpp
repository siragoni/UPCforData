#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TDatime.h"

using namespace RooFit ;


void rf501_simultaneouspdf()
{
  // C r e a t e   m o d e l   f o r   p h y s i c s   s a m p l e
  // -------------------------------------------------------------

  // Create observables
  RooRealVar Angle("Angle","Angle",-7.,7.) ;
  RooRealVar CosTheta("CosTheta","CosTheta",-1.,1.) ;
  RooRealVar Phi("Phi","Phi",-3.1415,3.1415) ;
  RooRealVar TildePhi("TildePhi","TildePhi",0.,2.*3.1415) ;




  // Theta distribution
  RooRealVar LambdaTheta("LambdaTheta","LambdaTheta",0., -1., 1.); //1.67, 1.91);
  RooRealVar LambdaPhi("LambdaPhi","LambdaPhi",0., -1., 1.); //1.67, 1.91);
  RooRealVar LambdaThetaPhi("LambdaThetaPhi","LambdaThetaPhi",0., -1., 1.); //1.67, 1.91);

  // RooGenericPdf  CosThetaPDF("CosThetaPDF","(1. + LambdaTheta*CosTheta*CosTheta)/(3.+LambdaTheta)",                   RooArgSet(CosTheta, LambdaTheta, LambdaPhi));
  // RooGenericPdf  PhiPDF     ("PhiPDF",     "(1. + 2.*LambdaPhi * cos(2.*Phi))/(3.+LambdaTheta)",                      RooArgSet(Phi, LambdaTheta, LambdaPhi));
  // RooGenericPdf  TildePhiPDF("TildePhiPDF","(1. + TMath::Sqrt(2)*LambdaThetaPhi * cos(2.*TildePhi))/(3.+LambdaTheta)",RooArgSet(TildePhi, LambdaTheta, LambdaThetaPhi));
  RooGenericPdf  CosThetaPDF("CosThetaPDF","(1. + LambdaTheta*Angle*Angle)/(3.+LambdaTheta)",                           RooArgSet(Angle, LambdaTheta, LambdaPhi));
  RooGenericPdf  PhiPDF     ("PhiPDF",     "(1. + 2.*LambdaPhi * cos(2.*Angle))/(3.+LambdaTheta)",                      RooArgSet(Angle, LambdaTheta, LambdaPhi));
  RooGenericPdf  TildePhiPDF("TildePhiPDF","(1. + TMath::Sqrt(2)*LambdaThetaPhi * cos(2.*Angle))/(3.+LambdaTheta)",     RooArgSet(Angle, LambdaTheta, LambdaThetaPhi));





  TDatime d;
  // TFile* file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* file1D = 0x0;
  if        ( SignalRangeSelectionMode == 0 ) {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  }
  TH1F* CorrectedCosTheta = (TH1F*) file1D->Get("CorrCosThetaH");
  TH1F* CorrectedPhi      = (TH1F*) file1D->Get("CorrPhiH");
  TH1F* CorrectedTildePhi = (TH1F*) file1D->Get("CorrTildePhiH");

  RooDataHist rCorrectedCosTheta("rCorrectedCosTheta", "rCorrectedCosTheta", CosTheta, Import(*CorrectedCosTheta));
  RooDataHist rCorrectedPhi     ("rCorrectedPhi",      "rCorrectedPhi",      Phi,      Import(*CorrectedPhi)     );
  RooDataHist rCorrectedTildePhi("rCorrectedTildePhi", "rCorrectedTildePhi", TildePhi, Import(*CorrectedTildePhi));






  // C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
  // ---------------------------------------------------------------------------

  // Define category to distinguish physics and control samples events
  RooCategory sample("sample","sample") ;
  sample.defineType("costheta") ;
  sample.defineType("phi") ;
  sample.defineType("tildephi") ;

  // Construct combined dataset in (x,sample)
  RooDataSet combData("combData","combined data",x,Index(sample),Import("physics",*data),Import("control",*data_ctl)) ;



  // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
  // -----------------------------------------------------------------------------------

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

  // Associate model with the physics state and model_ctl with the control state
  simPdf.addPdf(model,"physics") ;
  simPdf.addPdf(model_ctl,"control") ;



  // P e r f o r m   a   s i m u l t a n e o u s   f i t
  // ---------------------------------------------------

  // Perform simultaneous fit of model to data and model_ctl to data_ctl
  simPdf.fitTo(combData) ;



  // P l o t   m o d e l   s l i c e s   o n   d a t a    s l i c e s
  // ----------------------------------------------------------------

  // Make a frame for the physics sample
  RooPlot* frame1 = x.frame(Bins(30),Title("Physics sample")) ;

  // Plot all data tagged as physics sample
  combData.plotOn(frame1,Cut("sample==sample::physics")) ;

  // Plot "physics" slice of simultaneous pdf.
  // NBL You _must_ project the sample index category with data using ProjWData
  // as a RooSimultaneous makes no prediction on the shape in the index category
  // and can thus not be integrated
  simPdf.plotOn(frame1,Slice(sample,"physics"),ProjWData(sample,combData)) ;
  simPdf.plotOn(frame1,Slice(sample,"physics"),Components("px"),ProjWData(sample,combData),LineStyle(kDashed)) ;

  // The same plot for the control sample slice
  RooPlot* frame2 = x.frame(Bins(30),Title("Control sample")) ;
  combData.plotOn(frame2,Cut("sample==sample::control")) ;
  simPdf.plotOn(frame2,Slice(sample,"control"),ProjWData(sample,combData)) ;
  simPdf.plotOn(frame2,Slice(sample,"control"),Components("px_ctl"),ProjWData(sample,combData),LineStyle(kDashed)) ;



  TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",800,400) ;
  c->Divide(2) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;


}
