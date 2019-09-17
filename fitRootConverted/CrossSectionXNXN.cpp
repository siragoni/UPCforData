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
  Double_t fI0N0N[3]        = { 0.0430, 0.027,  0.043 };
  Double_t fI0NXN[3]        = { 0.0840, 0.084,  0.084 };
  Double_t fIXN0N[3]        = { 5.6430, 2.335,  1.511 };
  Double_t fIXNXN[3]        = { 0.5000, 0.485,  0.500 };
  Double_t fD               = 0.055;
  Double_t eJPsi[3]         = { (0.051 + 0.140)*0.5, (0.204 + 0.191)*0.5, (0.119 + 0.029)*0.5 };
  Double_t LUMI             = 747.849;
  Double_t BR               = 0.05961;
  Double_t NumOfJPsi0N0N[3] = { 3404, 9719, 4148 };
  Double_t NumOfJPsi0NXN[3] = {  242,  856,  412 };
  Double_t NumOfJPsiXN0N[3] = {  420, 1232,  592 };
  Double_t NumOfJPsiXNXN[3] = {  123,  485,  290 };

  for( Int_t iLoop = 0; iLoop < 3; iLoop++ ){
      // DSigmaDy0N0N[iLoop] = NumOfJPsi0N0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDy0NXN[iLoop] = NumOfJPsi0NXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDyXN0N[iLoop] = NumOfJPsiXN0N[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      // DSigmaDyXNXN[iLoop] = NumOfJPsiXNXN[iLoop]/( (1+fI[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDy0N0N[iLoop] = NumOfJPsi0N0N[iLoop]/( (1+fI0N0N[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDy0NXN[iLoop] = NumOfJPsi0NXN[iLoop]/( (1+fI0NXN[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDyXN0N[iLoop] = NumOfJPsiXN0N[iLoop]/( (1+fIXN0N[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
      DSigmaDyXNXN[iLoop] = NumOfJPsiXNXN[iLoop]/( (1+fIXNXN[iLoop]+fD)*eJPsi[iLoop]*BR*LUMI*0.5*0.95*1000 );
  }

  Double_t x1[3]          = { -4+1*(4-2.5)/6, -4+3*(4-2.5)/6, -4+5*(4-2.5)/6    };
  Double_t y1Error0N0N[3] = { // 83 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 119/( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 98 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              83 /( (1+fI0N0N[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              119/( (1+fI0N0N[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              98 /( (1+fI0N0N[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t y1Error0NXN[3] = { // 18 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 34 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 32 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              18 /( (1+fI0NXN[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              34 /( (1+fI0NXN[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              32 /( (1+fI0NXN[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t y1ErrorXN0N[3] = { // 29 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 52 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 29 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              29 /( (1+fIXN0N[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              52 /( (1+fIXN0N[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              29 /( (1+fIXN0N[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t y1ErrorXNXN[3] = { // 16 /( (1+fI[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              // 25 /( (1+fI[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              // 24 /( (1+fI[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              16 /( (1+fIXNXN[0]+fD)*eJPsi[0]*BR*LUMI*1.5*0.95*1000 ),
                              25 /( (1+fIXNXN[1]+fD)*eJPsi[1]*BR*LUMI*1.5*0.95*1000 ),
                              24 /( (1+fIXNXN[2]+fD)*eJPsi[2]*BR*LUMI*1.5*0.95*1000 )
                              };
  Double_t x1Error[3]     = {  (4-2.5)/6, (4-2.5)/6, (4-2.5)/6 };

  Coherent0N0N = new TGraphErrors(3, x1, DSigmaDy0N0N, x1Error, y1Error0N0N);
  Coherent0NXN = new TGraphErrors(3, x1, DSigmaDy0NXN, x1Error, y1Error0NXN);
  CoherentXN0N = new TGraphErrors(3, x1, DSigmaDyXN0N, x1Error, y1ErrorXN0N);
  CoherentXNXN = new TGraphErrors(3, x1, DSigmaDyXNXN, x1Error, y1ErrorXNXN);

  TMultiGraph *mg = new TMultiGraph();
  Coherent0N0N->SetMarkerStyle(20);
  Coherent0N0N->SetMarkerColor(2);
  Coherent0N0N->SetLineColor(2);
  mg->Add(Coherent0N0N);
  Coherent0NXN->SetMarkerStyle(20);
  Coherent0NXN->SetMarkerColor(3);
  Coherent0NXN->SetLineColor(3);
  mg->Add(Coherent0NXN);
  CoherentXN0N->SetMarkerStyle(20);
  CoherentXN0N->SetMarkerColor(4);
  CoherentXN0N->SetLineColor(4);
  mg->Add(CoherentXN0N);
  CoherentXNXN->SetMarkerStyle(20);
  CoherentXNXN->SetMarkerColor(6);
  CoherentXNXN->SetLineColor(6);
  mg->Add(CoherentXNXN);

  Coherent0N0N->SetTitle("J/#Psi 0N0N");
  Coherent0NXN->SetTitle("J/#Psi 0NXN");
  CoherentXN0N->SetTitle("J/#Psi XN0N");
  CoherentXNXN->SetTitle("J/#Psi XNXN");

  mg->Draw("APL");
  mg->GetXaxis()->SetTitle("y");
  mg->GetYaxis()->SetTitle("d#sigma/dy [mb]");
  // Change the axis limits
  gPad->BuildLegend();
  gPad->Modified();
  // mg->GetXaxis()->SetLimits(-1., 1.);
  // mg->SetMinimum(0.);
  // mg->SetMaximum(6000.);
  // c2->Print("pngResults/MultiGraph1Dview.png");
  // mg->Draw("a fb l3d");
  // c2->Print("pngResults/MultiGraph2Dview.png");
  // return c2;

  cout << "0N0N : " << endl << DSigmaDy0N0N[0] << endl << DSigmaDy0N0N[1] << endl << DSigmaDy0N0N[2] << endl;
  cout << "0NXN : " << endl << DSigmaDy0NXN[0] << endl << DSigmaDy0NXN[1] << endl << DSigmaDy0NXN[2] << endl;
  cout << "XN0N : " << endl << DSigmaDyXN0N[0] << endl << DSigmaDyXN0N[1] << endl << DSigmaDyXN0N[2] << endl;
  cout << "XNXN : " << endl << DSigmaDyXNXN[0] << endl << DSigmaDyXNXN[1] << endl << DSigmaDyXNXN[2] << endl;

  cout << "0N0N E: " << endl << y1Error0N0N[0] << endl << y1Error0N0N[1] << endl << y1Error0N0N[2] << endl;
  cout << "0NXN E: " << endl << y1Error0NXN[0] << endl << y1Error0NXN[1] << endl << y1Error0NXN[2] << endl;
  cout << "XN0N E: " << endl << y1ErrorXN0N[0] << endl << y1ErrorXN0N[1] << endl << y1ErrorXN0N[2] << endl;
  cout << "XNXN E: " << endl << y1ErrorXNXN[0] << endl << y1ErrorXNXN[1] << endl << y1ErrorXNXN[2] << endl;

}
