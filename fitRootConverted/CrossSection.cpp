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
void CrossSection(){
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  TGraphErrors *Coherent;
  TGraphErrors *CoherentAll;
  Double_t DSigmaDy[6] = { 0,0,0,0,0,0};
  Double_t DSigmaDyAll = 0;
  Double_t fI[7]       = { 0.055, 0.060, 0.059, 0.061, 0.052, 0.051, 0.032 };
  Double_t fD          = 0.055;
  Double_t eJPsi[7]    = { 0.120, 0.051, 0.140, 0.204, 0.191, 0.119, 0.029 };
  Double_t LUMI        = 746.26;
  Double_t BR          = 0.05961;
  Double_t NumOfJPsi[7]= { 22063, 0,0,0,0,0,0 };

  for( Int_t iLoop = 0; iLoop < 7; iLoop++ ){
      if( iLoop == 0 ) DSigmaDyAll = NumOfJPsi[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*1.5 );
  }

  Double_t x1[6]      = { -4+1*(4-2.5)/12, -4+3*(4-2.5)/12, -4+5*(4-2.5)/12, -4+7*(4-2.5)/12, -4+9*(4-2.5)/12, -4+11*(4-2.5)/12};
  Double_t y1[6]      = {   2289, 3352, 3350, 2091, 682,  12};
  Double_t x2[1]      = { (-4-2.5)/2 };
  Double_t y2[1]      = { DSigmaDyAll };
  Double_t y1Error[6] = {     71,   87,   89,   71,  41,   7};
  Double_t y2Error[1] = { 0 };
  Double_t x1Error[6] = {  (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12 };
  Double_t x2Error[1] = { (4-2.5)/2 };

  Coherent    = new TGraphErrors(6, x1, y1, x1Error, y1Error);
  CoherentAll = new TGraphErrors(1, x2, y2, x2Error, y2Error);

  TMultiGraph *mg = new TMultiGraph();
  Coherent->SetMarkerStyle(20);
  Coherent->SetMarkerColor(2);
  Coherent->SetLineColor(2);
  mg->Add(Coherent);
  CoherentAll->SetMarkerStyle(20);
  CoherentAll->SetMarkerColor(3);
  CoherentAll->SetLineColor(3);
  mg->Add(CoherentAll);

  Coherent->SetTitle("J/#Psi in rapidity bins");
  CoherentAll->SetTitle("J/#Psi integrated");

  mg->Draw("APL");
  mg->GetXaxis()->SetTitle("y");
  mg->GetYaxis()->SetTitle("d#sigma/dy");
  // Change the axis limits
  gPad->Modified();
  // mg->GetXaxis()->SetLimits(-1., 1.);
  // mg->SetMinimum(0.);
  // mg->SetMaximum(6000.);
  // c2->Print("pngResults/MultiGraph1Dview.png");
  // mg->Draw("a fb l3d");
  // c2->Print("pngResults/MultiGraph2Dview.png");
  // return c2;

}
