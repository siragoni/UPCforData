#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TMath.h"
#include "TF1.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "TLatex.h"
using namespace RooFit;
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"

//_____________________________________________________________________________
/* - Function that allows to retrieve the parameters for the fit.
   - The problem is that I have noticed that the templates given by Evegeny
   - Kryshen do not work that well. Especially at the gaussian part of the
   - Crystal Ball, the fit does not reproduce the MC points.
   - Furthermore, the fit does not even reproduce the height of the peaks.
   - I am not sure about what to do. Maybe it is best to just use only the
   - parameters for the background and trust the fit on this one.
   -
 */
std::vector<Double_t> parametersForTheFit(){
  /* - If you open a new TBrowser to check the plots, please mind his coding:

       gl - gamma-gamma-to-mumu continuum starting from 1.8 GeV/c^2
       gh - gamma-gamma-to-mumu continuum starting from 4 GeV/c^2
       1c - coherent J/psi
       2c - coherent psi(2S)
       1i - incoherent J/psi
       2i - incoherent psi(2S)
       fc - coherent feed-down
       fi - incoherent feed-down
       data_noTKLcut - data without SPD veto
       data - data with SPD veto

     - The histograms have the following naming convention:

       set_iy_ipt_rebin
       - iy    = rapidity bin
       - ipt   = pt cut
       - rebin = 1 for MC, 4 for data

     - The rapidity bins and the pt cuts are defined as follows:

       const Int_t nYBins = 10;
       Float_t gYMin[nYBins] = {-4.0, -4.00, -3.75, -3.50, -3.25, -3.00, -2.75, -4.0, -3.5, -3.0};
       Float_t gYMax[nYBins] = {-2.5, -3.75, -3.50, -3.25, -3.00, -2.75, -2.50, -3.5, -3.0, -2.5};
       const Int_t nPtBins = 7;
       const Int_t nPtBinsJpsi = 7;
       Float_t gPtMin[nPtBins] = {   0,   0,   0,     0,    0,      0,     0};
       Float_t gPtMax[nPtBins] = {0.25, 4.0, 0.2, 0.225, 0.25,  0.275,  0.30};

     - Last remark: by default I fit mc from 2 to 8 GeV/c2. Data from 2.2 to 6 GeV/c2.
     -
   */
  TFile* MCPlots = new TFile("MCplotsFromEvgenyKryshen.root");
  /* - This is the background used by Evgeny Kryshen.
     - He has been using ROOT all this time for the fits, so my gut feeling
     - is that he has been using mostly ChiSquare minimization, while my
     - RooFit functions use MINUIT by default.
     -
   */
  // TF1* f2 = new TF1("f2","[0]*[1]*exp(-[1]*x)*(x>=[2] ? 1 : 1+[3]*(x-[2])^2+[4]*(x-[2])^3+[5]*(x-[2])^4)",2,8);
  /* - This is the Crystal Ball used by Evgeny Kryshen.
     -
   */
  // TF1* fCB = new TF1("fCB","[0]*ROOT::Math::crystalball_pdf(x,[3],[4],[2],[1])",2,15);
  /* - The complete fit function used by him is the following:

       TF1 *fSum = new TF1("fSum",fsum,1,15,16);
       double fsum(double *x, double *par){
          Double_t par1s[5];
          Double_t par2s[5];
          for (Int_t i=0;i<5;i++) par1s[i]=par[i+6];
          par2s[0]=par[11]*par1s[0]; // yieldRatio
          par2s[1]=par[12];
          par2s[2]=par[13]*par1s[2]; // widthRatio
          par2s[3]=par[14];
          par2s[4]=par[15];
          return fBG->EvalPar(x,par)+fCB1->EvalPar(x,par1s)+fCB2->EvalPar(x,par2s);
       }
       fSum->SetParameter( 0,intbg*binWidth);
       fSum->SetParameter( 1,fBG->GetParameter(1));
       fSum->FixParameter( 2,fBG->GetParameter(2));
       fSum->FixParameter( 3,fBG->GetParameter(3));
       fSum->FixParameter( 4,fBG->GetParameter(4));
       fSum->FixParameter( 5,fBG->GetParameter(5));
       fSum->SetParameter( 6,intjpsi*binWidth);
       fSum->SetParameter( 7,fCB1->GetParameter(1));
       fSum->SetParameter( 8,fCB1->GetParameter(2));
       fSum->FixParameter( 9,fCB1->GetParameter(3));
       fSum->FixParameter(10,fCB1->GetParameter(4));
       fSum->SetParameter(11,0.02);
       fSum->FixParameter(12,fCB2->GetParameter(1));
       fSum->FixParameter(13,widthRatio);
       fSum->FixParameter(14,fCB2->GetParameter(3));
       fSum->FixParameter(15,fCB2->GetParameter(4));

   */
  TH1D* hM_1c = (TH1D*) MCPlots->Get("1c_0_0_1");
  TH1D* hM_2c = (TH1D*) MCPlots->Get("2c_0_0_1");
  TH1D* hM_gl = (TH1D*) MCPlots->Get("gl_0_0_1");
  TF1* fBG  = (TF1*) hM_gl->GetFunction("f2")->Clone("fBG");
  TF1* fCB1 = (TF1*) hM_1c->GetFunction("fCB")->Clone("fCB1");
  TF1* fCB2 = (TF1*) hM_2c->GetFunction("fCB")->Clone("fCB2");
  Double_t a2      = fBG->GetParameter(3);
  Double_t a3      = fBG->GetParameter(4);
  Double_t a4      = fBG->GetParameter(5);
  Double_t alphaMC = fCB1->GetParameter(3);
  Double_t nMC     = fCB1->GetParameter(4);
  cout << "a2      = " << a2      << endl;
  cout << "a3      = " << a3      << endl;
  cout << "a4      = " << a4      << endl;
  cout << "alphaMC = " << alphaMC << endl;
  cout << "nMC     = " << nMC     << endl;
  std::vector<Double_t> parametersForTheFitVector;
  parametersForTheFitVector.push_back(a2);
  parametersForTheFitVector.push_back(a3);
  parametersForTheFitVector.push_back(a4);
  parametersForTheFitVector.push_back(alphaMC);
  parametersForTheFitVector.push_back(nMC);
  return parametersForTheFitVector;
}


//_____________________________________________________________________________
/* - Fit function
   -
 */
void fitROOFit(const char* AnalysisName, const int selectionFlag){
  /* - There are three cases for the selectionFlag:
     - 1) = 0 ; this implies the traditional pt-integrated plot;
     - 2) = 1 ; this is instead the coherent component;
     - 3) = 2 ; this is the incoherent component;
     - 4) = 3 ; ******************* ;
     -
   */
  TFile* fileList = new TFile(AnalysisName);
  TDirectory* dir = fileList->GetDirectory("MyTask");
  /* - At this level you could check if everything was okay.
   * - We do a dir->ls() to find out! We get:
   *   dir->ls();
   *   TDirectoryFile*		MyTask	MyTask
   *   KEY: TList	MyOutputContainer;1	Doubly linked list
   */
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
  TH1F *fInvariantMassDistributionH = 0x0;
  if      ( selectionFlag == 0 ) fInvariantMassDistributionH = (TH1F*)listings->FindObject("fInvariantMassDistributionH");
  else if ( selectionFlag == 1 ) fInvariantMassDistributionH = (TH1F*)listings->FindObject("fInvariantMassDistributionCoherentH");
  else if ( selectionFlag == 2 ) fInvariantMassDistributionH = (TH1F*)listings->FindObject("fInvariantMassDistributionIncoherentH");
  else                           fInvariantMassDistributionH = (TH1F*)listings->FindObject("fInvariantMassDistributionH");
  fInvariantMassDistributionH->Rebin(5);
  fInvariantMassDistributionH->Draw();







  /* - Given the radiative processes for the muons, the UPC group has always
     - used a Crystal Ball function to model the peak for the J/Psi. This,
     - however, may be a bit too strong. It may be the best option for the
     - decay J/Psi->e+e-, so in the future it may be possible that other pdfs
     - may become object of study.
     -
     - The background is modelled in 2 parts. The most important is the
     - contribution for invariant masses higher than 1 GeV/c, the background
     - due to  gamma+gamma->mu+mu-, that is really exponential. I think for now
     - it is best to give it a simple exponential form and only later alter it,
     - when Monte Carlo become available.
     -
     - The second part is an effect due to the detector. This effect forces to
     - zero the distribution as it reaches null invariant mass. This appears to
     - be caused by the triggers, but it appears that such an effect is not so
     - well understood as of yet.
   */



  /* - What happens now is that we generate three different pdfs to describe
     - respectively the background, the J/Psi Peak and the PsiPrime Peak. I may
     - later on implement a second contribution for the background to account
     - for the detector effects on the lower end of the dimuon invariant mass
     - spectrum. Based on what Kay has told me this contribution is usually
     - modelled as a fourth polinomyal to make it go to zero, but the reason
     - why it really should be this kind of function still eludes me.
     -
   */
  RooRealVar      x           ("x"  ,          "x",          2.0,  6.0         );
  // RooRealVar      lambda      ("lambda",       "lambda",    -0.9, -1, -0.89   );  chi square 2.05
  RooRealVar      lambda      ("lambda",       "lambda",    -1.0, -1.4, -0.6   );
  // RooExponential  Bkg         ("Bkg",          "Bkg",          x, lambda       );
  /* - This will be a slightly different thing from before. I am now changing
     - the values to have them fixed, so, differently from before!
     -
   */
  std::vector<Double_t>  parametersForTheFitVector = parametersForTheFit();
  RooRealVar      bkgFactor   ("bkgFactor",     "bkgFactor", 10.,  0., 1000.   );
  RooRealVar      a2v         ("a2v",           "quadratic height", parametersForTheFitVector[0]*1.);
  RooRealVar      a3v         ("a3v",           "cubic height",     parametersForTheFitVector[1]*1.);
  RooRealVar      a4v         ("a4v",           "quartic height",   parametersForTheFitVector[2]*1.);
  RooGenericPdf   Bkg         ( "Bkg",          "",
                                "bkgFactor*exp( lambda * x ) * ( x > 4 ? 1: 1 + a2v*(x-4)*(x-4) + a3v*(x-4)*(x-4)*(x-4) + a4v*(x-4)*(x-4)*(x-4)*(x-4))",
                                RooArgList( x, bkgFactor, lambda, a2v, a3v, a4v )
                               );



  RooRealVar      mean        ("mean",         "",         3.096, 2.8  ,  3.2  );
  /* - This creates a problem for now without the polynomial on the left side...
     - My guess is that the system tries to compensate the lack of it with a
     - larger gaussian... So for now I will let it be a little freer and we will
     - see!
     -
   */
  RooRealVar      sigma       ("sigma",        "",         0.080, 0.050,  0.100);
  // RooRealVar      sigma       ("sigma",        "",         0.080, 0.050,  1.0  );
  // RooRealVar      alpha       ("alpha",        "",	       1.21,   1.1, 1.3  );
  RooRealVar      alpha       ("alpha",        "",	       parametersForTheFitVector[3]  );
  RooRealVar      n           ("n",            "",	       parametersForTheFitVector[4]/1000000*30.  );
  // RooRealVar      n           ("n",            "",	       2.93,   2,     3  );
  RooCBShape      JPsiPeak    ("JPsiPeak",     "JPsiPeak",     x,  mean,  sigma,  alpha,  n );
  RooRealVar      mean2       ("mean2",        "",         3.686, 3.68,   3.69 );
  RooRealVar      sigma2      ("sigma2",       "",         0.07884             );
  RooRealVar      alpha2      ("alpha2",       "",	       1.099               );
  RooRealVar      n2          ("n2",           "",	      10.0                 );
  RooCBShape      PsiPrimePeak("PsiPrimePeak", "PsiPrimePeak", x,  mean2, sigma2, alpha2, n2);
  /* - Useful reminder: whenever we need a tail on the right side we can use a
     - negative alpha/alpha2!
     -
   */

  /* - The respective values for the global pdf's normalizations
     -
   */
  RooRealVar nJPsiPeak    ("nJPsiPeak",     "",12000.,  10000., 30000.);
  RooRealVar nPsiPrimePeak("nPsiPrimePeak", "",     1,  0.,  1000);
  RooRealVar nBackGround  ("nBackGround",   "",   100,  0., 20000);

  /* - We now create the fit function used to fit the data. Hopefully it should
     - work!
     -
   */
  RooAddPdf* fit = new RooAddPdf( "fit", "",
                                  RooArgList( JPsiPeak,  PsiPrimePeak,  Bkg ),
                                  RooArgList( nJPsiPeak, nPsiPrimePeak, nBackGround)
                                 );
  RooDataHist* fInvariantMassDistributionROOH = new RooDataHist( "fInvariantMassDistributionH",
                                                                 "",
                                                                 x,
                                                                 fInvariantMassDistributionH
                                                                 );
  RooFitResult* fFitResult = fit->fitTo( *fInvariantMassDistributionROOH,
                                          RooFit::Save(),
                                          Range(2.2,6)
                                         );
  new TCanvas( "InvariantMassDimuonFit", "InvariantMassDimuonFit", 900, 800 );
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  /* - This option has been suggested in:
     - https://root-forum.cern.ch/t/tgraph-more-space-for-y-axis-label/5296
     - in order to increase the left margin!
   */
  // gPad->SetLeftMargin(0.30); //15 per cent of pad for left margin, default is 10%



  /* - Beautifying is starting now.
     -
   */
  RooPlot* frame = x.frame();
  frame->SetTitle(  Form(  ";M_{#mu#mu} (GeV/c^{2});Counts / (%.0f MeV/c^{2})",
                           fInvariantMassDistributionH->GetXaxis()->GetBinWidth(1)*1000.  )  );
  fInvariantMassDistributionROOH->plotOn( frame, MarkerStyle(kFullCircle), MarkerSize(1)     );
  fit                           ->plotOn( frame, LineColor(kBlue),         LineStyle(kSolid) );
  fit                           ->plotOn( frame, Components(JPsiPeak),     LineColor(kMagenta), LineStyle(kSolid)  );
  fit                           ->plotOn( frame, Components(PsiPrimePeak), LineColor(kRed),     LineStyle(kSolid)  );
  fit                           ->plotOn( frame, Components(Bkg),          LineColor(kGreen+1), LineStyle(kDashed) );
  frame->GetXaxis()->SetTitleOffset(1.25);
  // frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetYaxis()->SetTitleOffset(1.45);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetXaxis()->SetTitleFont(42);
  frame->GetYaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetLabelFont(42);
  frame->GetYaxis()->SetLabelFont(42);
  frame->GetXaxis()->SetNdivisions(408);
  frame->GetYaxis()->SetRangeUser(0,fInvariantMassDistributionH->GetMaximum()*1.3);
  frame->Draw();
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  latex->DrawLatex(0.55,0.84,"UPC, LHC18q+LHC18r data");
  // latex->DrawLatex(0.55,0.78,"#it{p}_{T} < 0.3 GeV/#it{c}");
  if      ( selectionFlag == 0 ) latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  else if ( selectionFlag == 1 ) latex->DrawLatex(0.55,0.78,"#it{p}_{T} < 0.25 GeV/#it{c}");
  else if ( selectionFlag == 2 ) latex->DrawLatex(0.55,0.78,"#it{p}_{T} > 0.25 GeV/#it{c}");
  else                           latex->DrawLatex(0.55,0.78,"#it{p}_{T}-integrated");
  latex->DrawLatex(0.55,0.72,Form("%.1f < y < %.1f",-4.0,-2.5));


  /* - This is the part where we obtain the actual number of J/Psi, PsiPrime
     - and the background. This is still Kay's original code. I will modify it.
     - Hopefully if everything goes alright, I should have been able to complete
     - it by the time you are reading this.
     -
   */
  x.setRange("signalJPsi",2.75,3.45);
  RooAbsReal* fIntegralJPsiPeak = JPsiPeak.createIntegral(x,NormSet(x),Range("signalJPsi"));
  Double_t pdfIntegralJPsiPeak = fIntegralJPsiPeak->getVal();
  Double_t massEfficiency = pdfIntegralJPsiPeak; //Fraction of dimuons in mass range from CB fit ->used in eff calc
  cout << "pdfIntegralJPsiPeak = " << pdfIntegralJPsiPeak << endl;
  Double_t numberOfTotalJPsi = nJPsiPeak.getVal() * pdfIntegralJPsiPeak;
  cout << "numberOfTotalJPsi = " << numberOfTotalJPsi << endl;

  // Double_t nJpsi = numberOfTotalJPsi;
  // Double_t nJpsiError = nJPsiPeak.getError() * pdfIntegralJPsiPeak;
  // //latex->DrawLatex(0.55,0.66,Form("N_{J/#psi} = %.0f #pm %.0f",n1s.getVal(),n1s.getError()));
  latex->DrawLatex(0.55,0.66,Form("N_{J/#psi} = %.0f #pm %.0f",                numberOfTotalJPsi,      nJPsiPeak.getError() * pdfIntegralJPsiPeak));
  latex->DrawLatex(0.55,0.60,Form("N_{#psi(2S)} = %.0f #pm %.0f",              nPsiPrimePeak.getVal(), nPsiPrimePeak.getError()));
  latex->DrawLatex(0.55,0.54,Form("#sigma_{J/#psi} = %.0f #pm %.0f MeV/c^{2}", sigma.getVal()*1000,    sigma.getError()*1000));
  latex->DrawLatex(0.55,0.48,Form("#sigma_{#psi(2S)} = %.0f MeV/c^{2} fixed",  sigma2.getVal()*1000));

  /* - This part concerns the background of the two signals.
     - Here, we extrapolate the background and compute the significance maybe?
     -
   */
  x.setRange("signalJPsiBis",  2.75, 3.45);
  x.setRange("signalPsiPrime", 3.45, 3.9 );
  RooAbsReal* bkgJPsiPeak                      = Bkg.createIntegral( x, NormSet(x), Range("signalJPsiBis") );
  RooAbsReal* bkgPsiPrimePeak                  = Bkg.createIntegral( x, NormSet(x), Range("signalPsiPrime"));
  Double_t    fValueBkgJPsiPeak                = bkgJPsiPeak->getVal();
  Double_t    fValueBkgPsiPrimePeak            = bkgPsiPrimePeak->getVal();
  Double_t    TotalBackGroundJPsiPeak          = nBackGround.getVal()   * fValueBkgJPsiPeak;
  Double_t    TotalBackGroundPsiPrimePeak      = nBackGround.getVal()   * fValueBkgPsiPrimePeak;
  Double_t    errorTotalBackGroundJPsiPeak     = nBackGround.getError() * fValueBkgJPsiPeak;
  Double_t    errorTotalBackGroundPsiPrimePeak = nBackGround.getError() * fValueBkgPsiPrimePeak;
  latex->DrawLatex(0.55,0.42,Form("N_{BG J/#psi} = %.0f #pm %.0f", TotalBackGroundJPsiPeak , errorTotalBackGroundJPsiPeak));
  latex->DrawLatex(0.55,0.36,Form("N_{BG #psi(2s)} = %.0f #pm %.0f", TotalBackGroundPsiPrimePeak, errorTotalBackGroundPsiPrimePeak));
  //Signal to background ratios
  // Double_t errSB1 = sqrt( pow(( (nJPsiPeak.getError())/(nJPsiPeak.getVal()) ),2) + pow(( (errorTotalBackGroundJPsiPeak)/(TotalBackGroundJPsiPeak) ),2) );
  // Double_t errSB2 = sqrt( pow(( (nPsiPrimePeak.getError())/(nPsiPrimePeak.getVal()) ),2) + pow(( (errorTotalBackGroundPsiPrimePeak)/(TotalBackGroundPsiPrimePeak) ),2) );
  Double_t errSB1 = sqrt( (nJPsiPeak.getError()    *nJPsiPeak.getError()    /( nJPsiPeak.getVal()    *nJPsiPeak.getVal() ))     + (errorTotalBackGroundJPsiPeak    *errorTotalBackGroundJPsiPeak    /( TotalBackGroundJPsiPeak    *TotalBackGroundJPsiPeak)) );
  Double_t errSB2 = sqrt( (nPsiPrimePeak.getError()*nPsiPrimePeak.getError()/( nPsiPrimePeak.getVal()*nPsiPrimePeak.getVal() )) + (errorTotalBackGroundPsiPrimePeak*errorTotalBackGroundPsiPrimePeak/( TotalBackGroundPsiPrimePeak*TotalBackGroundPsiPrimePeak )) );
  latex->DrawLatex(0.55,0.30,Form("S/B for J/#psi = %.2f #pm %.2f", nJPsiPeak.getVal()/(Double_t)TotalBackGroundJPsiPeak, errSB1*nJPsiPeak.getVal()/(Double_t)TotalBackGroundJPsiPeak));
  latex->DrawLatex(0.55,0.24,Form("S/B for #psi(2s) = %.2f #pm %.2f", nPsiPrimePeak.getVal()/(Double_t)TotalBackGroundPsiPrimePeak, errSB2*nPsiPrimePeak.getVal()/(Double_t)TotalBackGroundPsiPrimePeak));
  // Double_t ChiSquareFromTheFit = frame->chiSquare();
  // Double_t ChiSquareFromTheFit = fit->createChi2(*fInvariantMassDistributionROOH)->getVal();
  RooChi2Var chi2_lowstat("chi2_lowstat","chi2",*fit,*fInvariantMassDistributionROOH) ;
  latex->DrawLatex(0.55,0.18,Form("      #tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                  chi2_lowstat.getVal(),
                                  fInvariantMassDistributionROOH->numEntries() - fit->getParameters(*fInvariantMassDistributionROOH)->selectByAttrib("Constant",kFALSE)->getSize(),
                                  chi2_lowstat.getVal()/( fInvariantMassDistributionROOH->numEntries() - fit->getParameters(*fInvariantMassDistributionROOH)->selectByAttrib("Constant",kFALSE)->getSize() )
                                   )
                                 );
  // latex->DrawLatex(0.55,0.18,Form("              #tilde{#chi}^{2} = %.2f ", ChiSquareFromTheFit ));

  if      ( selectionFlag == 0 ) gPad->SaveAs("pngResults/suppressedNjpsi.png",           "RECREATE");
  else if ( selectionFlag == 1 ) gPad->SaveAs("pngResults/suppressedNjpsiCoherent.png",   "RECREATE");
  else if ( selectionFlag == 2 ) gPad->SaveAs("pngResults/suppressedNjpsiIncoherent.png", "RECREATE");
  else                           gPad->SaveAs("pngResults/ops.png",           "RECREATE");

  cout << "nBackGround:           " << nBackGround.getVal() << endl;
  cout << "fValueBkgJPsiPeak:     " << fValueBkgJPsiPeak << endl;
  cout << "fValueBkgPsiPrimePeak: " << fValueBkgPsiPrimePeak << endl;

  // exp = TotalBackGroundJPsiPeak;
  // err = errorTotalBackGroundJPsiPeak;










  //____________________________________________________________________________
  // FIT VALIDATION
  TFile *outputFile = new TFile(Form("test%i.root", selectionFlag), "RECREATE");
  // P l o t   d a t a   a n d   s l i g h t l y   d i s t o r t e d   m o d e l
  // ---------------------------------------------------------------------------

  // Overlay projection of gauss with sigma=3.15 on data with sigma=3.0
  RooPlot* frame1 = x.frame() ;
  frame1->SetTitle(  Form(  ";M_{#mu#mu} (GeV/c^{2});Counts / (%.0f MeV/c^{2})",
                            fInvariantMassDistributionH->GetXaxis()->GetBinWidth(1)*1000.  )  );
  fInvariantMassDistributionROOH->plotOn( frame1, MarkerStyle(kFullCircle), MarkerSize(1),       DataError(RooAbsData::SumW2));
  fit                           ->plotOn( frame1, LineColor(kBlue),         LineStyle(kSolid) );
  fit                           ->plotOn( frame1, Components(JPsiPeak),     LineColor(kMagenta), LineStyle(kSolid)           );
  fit                           ->plotOn( frame1, Components(PsiPrimePeak), LineColor(kRed),     LineStyle(kSolid)           );
  fit                           ->plotOn( frame1, Components(Bkg),          LineColor(kGreen+1), LineStyle(kDashed)          );


  // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
  // -------------------------------------------------------

  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = frame1->residHist() ;

  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = frame1->pullHist() ;

  // Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot* frame2 = x.frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P") ;

  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame3 = x.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P") ;



  TCanvas* c = new TCanvas("Residual and pull distr","Residual and pull distr",900,300) ;
  c->Divide(3) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.6) ; frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.6) ; frame3->Draw() ;
  c->Write();
  hresid->Write();
  hpull->Write();
  // V i s u a l i z e   c o r r e l a t i o n   m a t r i x
  // -------------------------------------------------------

  // Construct 2D color plot of correlation matrix
  gStyle->SetOptStat(0) ;
  TH2* hcorr = fFitResult->correlationHist() ;
  hcorr->Write();

  // Visualize ellipse corresponding to single correlation matrix element
  RooPlot* frame4 = new RooPlot(lambda, mean) ;
  frame->SetTitle("Covariance between lambda and mean") ;
  fFitResult->plotOn(frame4,lambda,mean,"ME12ABHV") ;
  frame4->Draw();
  c->Write();


  // A c c e s s   f i t   r e s u l t   i n f o r m a t i o n
  // ---------------------------------------------------------

  // Access basic information
  cout << "EDM =                " << fFitResult->edm()    << endl ;
  cout << "-log(L) at minimum = " << fFitResult->minNll() << endl ;

  // Access list of final fit parameter values
  cout << "final value of floating parameters" << endl ;
  fFitResult->floatParsFinal().Print("s") ;

  // Access correlation matrix elements
  cout << "correlation between mean and sigma  is " << fFitResult->correlation(mean, sigma) << endl ;
  cout << "correlation between mean and lambda is " << fFitResult->correlation(mean,lambda) << endl ;
  cout << "correlation between nJPsiPeak and nPsiPrimePeak is " << fFitResult->correlation(nJPsiPeak,nPsiPrimePeak) << endl ;
  cout << "correlation between nJPsiPeak and nBackGround is " << fFitResult->correlation(nJPsiPeak,nBackGround) << endl ;

  // Extract covariance and correlation matrix as TMatrixDSym
  const TMatrixDSym& cor = fFitResult->correlationMatrix() ;
  const TMatrixDSym& cov = fFitResult->covarianceMatrix() ;

  // Print correlation, covariance matrix
  cout << "correlation matrix" << endl ;
  cor.Print() ;
  cout << "covariance matrix" << endl ;
  cov.Print() ;
  cor.Write();
  cov.Write();
  outputFile->Close();

}
