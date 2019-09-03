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
void CrossSectionXNXN(){
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  TGraphErrors *Coherent0N0N;
  TGraphErrors *Coherent0NXN;
  TGraphErrors *CoherentXN0N;
  TGraphErrors *CoherentXNXN;
  Double_t DSigmaDy0N0N[3]  = { 0,0,0};
  Double_t DSigmaDy0NXN[3]  = { 0,0,0};
  Double_t DSigmaDyXN0N[3]  = { 0,0,0};
  Double_t DSigmaDyXNXN[3]  = { 0,0,0};
  Double_t fI[3]            = { 0.0595, 0.056,  0.041 };
  Double_t fD               = 0.055;
  Double_t eJPsi[3]         = { (0.051 + 0.140)*0.5, (0.204 + 0.191)*0.5, (0.119 + 0.029)*0.5 };
  Double_t LUMI             = 747.849;
  Double_t BR               = 0.05961;
  Double_t NumOfJPsi0N0N[3] = { 3324, 9766, 4364 };
  Double_t NumOfJPsi0NXN[3] = {  237,  858,  431 };
  Double_t NumOfJPsiXN0N[3] = {  420, 1235,  605 };
  Double_t NumOfJPsiXNXN[3] = {  122,  486,  293 };

  for( Int_t iLoop = 0; iLoop < 3; iLoop++ ){
      DSigmaDy0N0N[iLoop] = NumOfJPsi0N0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*1.5*0.95*1000 );
      DSigmaDy0NXN[iLoop] = NumOfJPsi0NXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*1.5*0.95*1000 );
      DSigmaDyXN0N[iLoop] = NumOfJPsiXN0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*1.5*0.95*1000 );
      DSigmaDyXNXN[iLoop] = NumOfJPsiXNXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*1.5*0.95*1000 );
  }

  Double_t x1[6]      = { -4+1*(4-2.5)/12, -4+3*(4-2.5)/12, -4+5*(4-2.5)/12, -4+7*(4-2.5)/12, -4+9*(4-2.5)/12, -4+11*(4-2.5)/12};
  Double_t y1[6]      = { DSigmaDy[0], DSigmaDy[1], DSigmaDy[2], DSigmaDy[3], DSigmaDy[4], DSigmaDy[5] };
  Double_t x2[1]      = { (-4-2.5)/2 };
  Double_t y2[1]      = { DSigmaDyAll };
  Double_t y1Error[6] = { 45 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                          84 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 ),
                          114/( (1+fI[3]+fD)*eJPsi[3]*BR*LUMI*1.5*0.95*1000 ),
                          97 /( (1+fI[4]+fD)*eJPsi[4]*BR*LUMI*1.5*0.95*1000 ),
                          100/( (1+fI[5]+fD)*eJPsi[5]*BR*LUMI*1.5*0.95*1000 ),
                          55 /( (1+fI[6]+fD)*eJPsi[6]*BR*LUMI*1.5*0.95*1000 )
                          };
  Double_t y2Error[1] = { 223/( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ) };
  Double_t x1Error[6] = {  (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12 };
  Double_t x2Error[1] = {  (4-2.5)/2 };

  Coherent0N0N = new TGraphErrors(3, x1, y1, x1Error, y1Error);
  Coherent0NXN = new TGraphErrors(3, x1, y1, x1Error, y1Error);
  CoherentXN0N = new TGraphErrors(3, x1, y1, x1Error, y1Error);
  CoherentXNXN = new TGraphErrors(3, x1, y1, x1Error, y1Error);

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
  mg->GetYaxis()->SetTitle("d#sigma/dy [mb]");
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
