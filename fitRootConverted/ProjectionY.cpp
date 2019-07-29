#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TString.h"
using namespace std;
#include <math.h>
#include <vector>


//_____________________________________________________________________________
/* - Projection along Y of the DCA vs Inv Mass.
   -
 */
void ProjectionY(){
  TFile* fileList = new TFile("AnalysisResultsLHC18qr15o29052019noSPD.root");
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
  TH2F *fDCAvsInvMass = (TH2F*)listings->FindObject("fDcaAgainstInvariantMassH");

  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  fDCAvsInvMass->Draw("colZ");


  // Int_t whichBin  = fDCAvsInvMass->GetYaxis()->FindBin(5.9999);   // bump at 4.3
  // Int_t whichBin2 = fDCAvsInvMass->GetYaxis()->FindBin(10.9999);  // bump at 4.3
  // Int_t whichBin  = fDCAvsInvMass->GetYaxis()->FindBin(3.9999);   // best bump at 4.3
  // Int_t whichBin2 = fDCAvsInvMass->GetYaxis()->FindBin(4.3999);   // best bump at 4.3
  Int_t whichBin  = fDCAvsInvMass->GetYaxis()->FindBin(2.9999);   // best bump at 4.3
  Int_t whichBin2 = fDCAvsInvMass->GetYaxis()->FindBin(60.9999);   // best bump at 4.3
  // Int_t whichBin  = fDCAvsInvMass->GetXaxis()->FindBin(4.9999);
  // Int_t whichBin2 = fDCAvsInvMass->GetXaxis()->FindBin(5.9999);

  // TH1D* fProjInvMass = fDCAvsInvMass->ProjectionY("fProjInvMass", whichBin, whichBin2, "");
  TH1D* fProjInvMass = fDCAvsInvMass->ProjectionX("fProjInvMass", whichBin, whichBin2, "");
  TCanvas *c3 = new TCanvas("c3","c3",600,400);
  fProjInvMass->Draw();

  cout << "whichBin  = " << whichBin  << endl;
  cout << "whichBin2 = " << whichBin2 << endl;



  // TGraphErrors *g[2];
  //
  // Double_t x1[8]      = { -0.7, -0.5, -0.3, -0.1, 0.1 , 0.3 , 0.5, 0.7};
  // Double_t y1[8]      = {  5  ,  743, 2289, 3352, 3350, 2091, 682,  12};
  // Double_t x2[8]      = { -0.7, -0.5, -0.3, -0.1, 0.1 , 0.3 , 0.5, 0.7};
  // Double_t y2[8]      = {  15 ,  400, 1241, 1835, 1993, 1366, 483,  12};
  // Double_t y1Error[8] = {  6  ,   41,   71,   87,   89,   71,  41,   7};
  // Double_t y2Error[8] = {  15 ,   95,  169,  210,  222,  210, 134,  36};
  //
  // g[0] = new TGraphErrors(8, x1, y1, 0, y1Error);
  // g[1] = new TGraphErrors(8, x2, y2, 0, y2Error);
  //
  // TMultiGraph *mg = new TMultiGraph();
  // for (int i=0; i<2; i++) {
  //    g[i]->SetMarkerStyle(20);
  //    g[i]->SetMarkerColor(i+2);
  //    g[i]->SetLineColor(i+2);
  //    mg->Add(g[i]);
  // }
  // g[0]->SetTitle("J/Psi");
  // g[1]->SetTitle("GammaGamma");
  //
  // mg->Draw("APL");
  // mg->GetXaxis()->SetTitle("#cos(#theta)");
  // mg->GetYaxis()->SetTitle("Counts");
  // // Change the axis limits
  // gPad->Modified();
  // mg->GetXaxis()->SetLimits(-1., 1.);
  // mg->SetMinimum(0.);
  // mg->SetMaximum(6000.);
  // c2->Print("pngResults/MultiGraph1Dview.png");
  // mg->Draw("a fb l3d");
  // c2->Print("pngResults/MultiGraph2Dview.png");
  // return c2;

}
