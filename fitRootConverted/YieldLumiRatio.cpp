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
void YieldLumiRatio(){
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  TGraphErrors *CoherentMine;
  TGraphErrors *CoherentEvgeny;
  TGraphErrors *CoherentRatio;
  Double_t LUMImine        = 747.849;
  Double_t LUMIevgeny      = 754;
  Double_t NumOfJPsiMine[7]   = { 22141, 968, 3263, 5890, 6459, 4283, 1212 };
  Double_t NumOfJPsiEvgeny[7] = { 21747, 974, 3217, 5769, 6387, 4229, 1190 };
  Double_t YieldLumiRatioMine[7];
  Double_t YieldLumiRatioEvgeny[7];
  Double_t FinalRatio[7];
  for( Int_t iLoop = 0; iLoop < 7; iLoop++ ){
    YieldLumiRatioMine[iLoop]   = NumOfJPsiMine[iLoop]/LUMImine;
    YieldLumiRatioEvgeny[iLoop] = NumOfJPsiEvgeny[iLoop]/LUMIevgeny;
    FinalRatio[iLoop]           = YieldLumiRatioMine[iLoop]/YieldLumiRatioEvgeny[iLoop];
  }

  Double_t x1[7]      = { (-4-2.5)/2, -4+1*(4-2.5)/12, -4+3*(4-2.5)/12, -4+5*(4-2.5)/12, -4+7*(4-2.5)/12, -4+9*(4-2.5)/12, -4+11*(4-2.5)/12};
  // Double_t y1[7]      = { DSigmaDy[0], DSigmaDy[1], DSigmaDy[2], DSigmaDy[3], DSigmaDy[4], DSigmaDy[5] };
  Double_t y1Error[7] = { 223/LUMImine,   45/LUMImine,   84/LUMImine,   114/LUMImine,    97/LUMImine,   100/LUMImine,   55/LUMImine   };
  Double_t y2Error[7] = { 190/LUMIevgeny, 36/LUMIevgeny, 70/LUMIevgeny,  98/LUMIevgeny, 105/LUMIevgeny,  85/LUMIevgeny, 47/LUMIevgeny };
  Double_t y3Error[7] = { 0,0,0,0,0,0,0 };
  Double_t x1Error[7] = {  (4-2.5)/2, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12, (4-2.5)/12 };
  // Double_t x2Error[1] = {  (4-2.5)/2 };

  CoherentMine    = new TGraphErrors(7, x1, YieldLumiRatioMine,   x1Error, y1Error);
  CoherentEvgeny  = new TGraphErrors(7, x1, YieldLumiRatioEvgeny, x1Error, y2Error);
  CoherentRatio   = new TGraphErrors(7, x1, FinalRatio,           x1Error, y3Error);

  TMultiGraph *mg = new TMultiGraph();
  CoherentMine->SetMarkerStyle(20);
  CoherentMine->SetMarkerColor(2);
  CoherentMine->SetLineColor(2);
  mg->Add(CoherentMine);
  CoherentEvgeny->SetMarkerStyle(20);
  CoherentEvgeny->SetMarkerColor(3);
  CoherentEvgeny->SetLineColor(3);
  mg->Add(CoherentEvgeny);
  CoherentRatio->SetMarkerStyle(20);
  CoherentRatio->SetMarkerColor(4);
  CoherentRatio->SetLineColor(4);

  CoherentMine  ->SetTitle("J/#Psi Yield Ratio Mine");
  CoherentEvgeny->SetTitle("J/#Psi Yield Ratio Evgeny");

  mg->Draw("APL");
  mg->GetXaxis()->SetTitle("y");
  mg->GetYaxis()->SetTitle("d#sigma/dy [mb]");
  // Change the axis limits
  gPad->Modified();
  new TCanvas;
  CoherentRatio->Draw();
  // mg->GetXaxis()->SetLimits(-1., 1.);
  // mg->SetMinimum(0.);
  // mg->SetMaximum(6000.);
  // c2->Print("pngResults/MultiGraph1Dview.png");
  // mg->Draw("a fb l3d");
  // c2->Print("pngResults/MultiGraph2Dview.png");
  // return c2;

}
