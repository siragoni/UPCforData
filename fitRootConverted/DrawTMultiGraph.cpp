#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
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
/* - Drawing function for the invariant mass distributions' retrieved valued,
   - after the fit in terms of CosTheta bins.
   -
 */
void DrawTMultiGraph(){
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  TGraphErrors *g[2];

  Double_t x1[8]      = { -0.7, -0.5, -0.3, -0.1, 0.1 , 0.3 , 0.5, 0.7};
  Double_t y1[8]      = {  5  ,  743, 2289, 3352, 3350, 2091, 682,  12};
  Double_t x2[8]      = { -0.7, -0.5, -0.3, -0.1, 0.1 , 0.3 , 0.5, 0.7};
  Double_t y2[8]      = {  15 ,  400, 1241, 1835, 1993, 1366, 483,  12};
  Double_t y1Error[8] = {  6  ,   41,   71,   87,   89,   71,  41,   7};
  Double_t y2Error[8] = {  15 ,   95,  169,  210,  222,  210, 134,  36};

  g[0] = new TGraphErrors(8, x1, y1, 0, y1Error);
  g[1] = new TGraphErrors(8, x2, y2, 0, y2Error);

  TMultiGraph *mg = new TMultiGraph();
  for (int i=0; i<2; i++) {
     g[i]->SetMarkerStyle(20);
     g[i]->SetMarkerColor(i+2);
     g[i]->SetLineColor(i+2);
     mg->Add(g[i]);
  }
  g[0]->SetTitle("J/Psi");
  g[1]->SetTitle("GammaGamma");

  mg->Draw("APL");
  mg->GetXaxis()->SetTitle("#cos(#theta)");
  mg->GetYaxis()->SetTitle("Counts");
  // Change the axis limits
  gPad->Modified();
  mg->GetXaxis()->SetLimits(-1., 1.);
  mg->SetMinimum(0.);
  mg->SetMaximum(6000.);
  c2->Print("pngResults/MultiGraph1Dview.png");
  mg->Draw("a fb l3d");
  c2->Print("pngResults/MultiGraph2Dview.png");
  return c2;

}
