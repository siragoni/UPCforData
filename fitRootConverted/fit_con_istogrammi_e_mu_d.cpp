#include "TProfile.h"
#include "iostream"
#include "TAxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2F.h"
#include "TFile.h"
#include "TH2.h"
#include "TROOT.h"

TH1 *fPiTemp;
TH1 *fKaTemp;
TH1 *fPrTemp;
TH1 *fMiTemp;

TH1 *feTemp;
TH1 *fMuTemp;
TH1 *fDTemp;


Double_t fTOFResponseT(Double_t* x,Double_t* par)
{

/* TF1* HELP = new TF1
  Double_t val  = par[0]*fPiTemp->Eval(x[0]);
  val += par[1]*fKaTemp->Eval(x[0]);
  val += par[2]*fPrTemp->Eval(x[0]);
//  val + = par[3]*fMiTemp->Eval(x[0]);
*/
 /* double  bin = TMath::Floor(TMath::Abs(x[0])/20);
  if(x[0] < 0.){bin = bin + 1; bin = bin * (-1);}
cout << "BIN:" << bin << endl;*/

 

  //double bin = fPiTemp -> FindBin (x[0]);
  Double_t val  = par[0]*fPiTemp->Interpolate(x[0]);
  val += par[1]*fKaTemp->Interpolate(x[0]);
  val += par[2]*fPrTemp->Interpolate(x[0]);
  val += par[3]*fMiTemp->Interpolate(x[0]);

  
  val += par[4]*feTemp->Interpolate(x[0]);
  val += par[5]*fMuTemp->Interpolate(x[0]);
  val += par[6]*fDTemp->Interpolate(x[0]);

/*

  double bin = fPiTemp -> FindBin (x[0]);
  Double_t val  = par[0]*fPiTemp->GetBinContent(bin);
  val += par[1]*fKaTemp->GetBinContent(bin);
  val += par[2]*fPrTemp->GetBinContent(bin);
*/



//  delete HELP;
  return val;
}


Double_t fTOFResponseT_original(Double_t* x,Double_t* par)
{

/* TF1* HELP = new TF1
  Double_t val  = par[0]*fPiTemp->Eval(x[0]);
  val += par[1]*fKaTemp->Eval(x[0]);
  val += par[2]*fPrTemp->Eval(x[0]);
//  val + = par[3]*fMiTemp->Eval(x[0]);
*/
 /* double  bin = TMath::Floor(TMath::Abs(x[0])/20);
  if(x[0] < 0.){bin = bin + 1; bin = bin * (-1);}
cout << "BIN:" << bin << endl;*/

 

  //double bin = fPiTemp -> FindBin (x[0]);
  Double_t val  = par[0]*fPiTemp->Interpolate(x[0]);
  val += par[1]*fKaTemp->Interpolate(x[0]);
  val += par[2]*fPrTemp->Interpolate(x[0]);
  val += par[3]*fMiTemp->Interpolate(x[0]);

  
//  val += par[4]*feTemp->Interpolate(x[0]);
//  val += par[5]*fMuTemp->Interpolate(x[0]);
  val += par[4]*fDTemp->Interpolate(x[0]);

/*

  double bin = fPiTemp -> FindBin (x[0]);
  Double_t val  = par[0]*fPiTemp->GetBinContent(bin);
  val += par[1]*fKaTemp->GetBinContent(bin);
  val += par[2]*fPrTemp->GetBinContent(bin);
*/



//  delete HELP;
  return val;
}




void combine(){

  //gROOT->Reset();
  //gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  Float_t mean_pi=0, mean_ka=600,mean_pr=2200,sigma=80;

  TFile* f1 = new TFile("all_histo_NEW_TEMPLATE.root");
  TCanvas *c1 = new TCanvas("c1","the fit canvas",500,400);
 // TCanvas *c2 = new TCanvas("c2","the fit canvas",500,400);
 // TCanvas *c3 = new TCanvas("c3","the fit canvas",500,400);
 // TCanvas *c4 = new TCanvas("c4","the fit canvas",500,400);
 // TCanvas *c5 = new TCanvas("c5","the fit canvas",500,400);
 // TCanvas *c6 = new TCanvas("c6","the fit canvas",500,400);


	TF1* f3 = new TF1("myfunc3",fTOFResponseT,-10000,10000,7);
 	f3->SetNpx(5000);

	TF1* f3_original = new TF1("myfunc3_original",fTOFResponseT_original,-10000,10000,5);
 	f3_original->SetNpx(5000);

// N.B.
//  shape_timeEXP_differences [nC] [0]   pioni rispetto ai kaoni
//  shape_timeEXP_differences [nC] [1]   pioni rispetto ai protoni
//  shape_timeEXP_differences [nC] [2]   kaoni rispetto ai pioni
//  shape_timeEXP_differences [nC] [3]   kaoni rispetto ai protoni
//  shape_timeEXP_differences [nC] [4]   protoni rispetto ai pioni
//  shape_timeEXP_differences [nC] [5]   protoni rispetto ai kaoni



	c1->cd();


TH2F* new_template[3];


TH2F* ordinamento[14];
ordinamento[0]=(TH2F*)f1->Get(Form("shape_timeEXP_differencesEXTRA0_%d",0));
ordinamento[1]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",0));
ordinamento[2]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",1));
ordinamento[3]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",2));
ordinamento[4]=(TH2F*)f1->Get(Form("shape_timeEXP_differencesEXTRA0_%d",1));
ordinamento[5]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",3));
ordinamento[6]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",4));
ordinamento[7]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",5));
ordinamento[8]=(TH2F*)f1->Get(Form("shape_timeEXP_differencesEXTRA0_%d",2));
(ordinamento[0])->ProjectionY(Form("bin%d",51),51,51)->SaveAs("HELP.root");


ordinamento[9]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",6));
ordinamento[10]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",7));
ordinamento[11]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",8));

ordinamento[12]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",9));
ordinamento[13]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",10));


TH2F* mismatch[3];
mismatch[0]=(TH2F*)f1->Get(Form("tempMismSpec0_%d",2));
//mismatch[0]->Draw();
//c1->SaveAs("mismatch0.pdf");
mismatch[1]=(TH2F*)f1->Get(Form("tempMismSpec0_%d",3));
//mismatch[1]->Draw();
//c1->SaveAs("mismatch1.pdf");
mismatch[2]=(TH2F*)f1->Get(Form("tempMismSpec0_%d",4));
//mismatch[2]->Draw();
//c1->SaveAs("mismatch2.pdf");

/*
 TCanvas *c2 = new TCanvas("c2","the fit canvas",500,400);
ordinamento[1]->Draw();
  TCanvas *c3 = new TCanvas("c3","the fit canvas",500,400);
ordinamento[2]->Draw();
  TCanvas *c4 = new TCanvas("c4","the fit canvas",500,400);
ordinamento[3]->Draw();
  TCanvas *c5 = new TCanvas("c5","the fit canvas",500,400);
ordinamento[4]->Draw();
  TCanvas *c6 = new TCanvas("c6","the fit canvas",500,400);
ordinamento[5]->Draw();
  TCanvas *c7 = new TCanvas("c7","the fit canvas",500,400);
ordinamento[6]->Draw();
  TCanvas *c8 = new TCanvas("c8","the fit canvas",500,400);
ordinamento[7]->Draw();
  TCanvas *c9 = new TCanvas("c9","the fit canvas",500,400);
ordinamento[8]->Draw();
*/




for(int SPECIE=0; SPECIE<3; ++SPECIE)
{

	TH1F* normalizations[4];
	for (int j=0; j<4; j++) normalizations[j]=new TH1F(Form("normalizations%d",j),"normalizations",100,0, 5);

	TH1F* integral_vs_pt = new TH1F("integral_vs_pt","integral_vs_pt", 100,0,5);


TH1 *h[100] ;//
TH2F* help;

for (int i=0;i<100;i++) 
{ 
  help=(TH2F*)f1->Get(Form("timeEXP_differences0_%d",SPECIE));
  h[i]= help->ProjectionY(Form("bin%d",i+1),i+1,i+1);

  h[i]->SetLineColor(96);
  h[i]->GetXaxis()->SetTitle("time differences");
  h[i]->GetYaxis()->SetTitle("");
  h[i]->SetTitle("FIT shape EXP time differences");


 c1->SaveAs("prova.pdf","recreate");



/*	f2->SetParameter(0, f->GetParameter(0));
	f2->SetParameter(1, f->GetParameter(1));
	f2->SetParameter(2, f->GetParameter(2));
	f2->SetParameter(3, f->GetParameter(9));
	f2->SetLineColor(4);
	f2->Draw("SAME");
	double a = f2->Integral(-4000,4000);
	double b = a / h[i]->GetBinWidth(1);
	integral_vs_pt->SetBinContent(i+1, b);*/





fPiTemp = (ordinamento[3*SPECIE + 0])->ProjectionY(Form("fPiTemp%d",i+1),i+1,i+1);
//fPiTemp->Draw();
fKaTemp = ordinamento[3*SPECIE + 1]->ProjectionY(Form("fKaTemp%d",i+1),i+1,i+1);
fPrTemp = ordinamento[3*SPECIE + 2]->ProjectionY(Form("fPrTemp%d",i+1),i+1,i+1);

fMiTemp = mismatch[SPECIE]->ProjectionY(Form("fMiTemp%d",i+1),i+1,i+1);

feTemp  = ordinamento[9]->ProjectionY(Form("feTemp%d",i+1),i+1,i+1);
fMuTemp = ordinamento[10]->ProjectionY(Form("fMuTemp%d",i+1),i+1,i+1);
fDTemp  = ordinamento[11 + SPECIE]->ProjectionY(Form("fDTemp%d",i+1),i+1,i+1);


/*
double Integral_Pi = 0;
for (int l=0;l<1000;l++) Integral_Pi += fPiTemp -> GetBinContent(l+1);
cout << "Integral_Pi=" << Integral_Pi << endl;

double Integral_Ka = 0;
for (int l=0;l<1000;l++) Integral_Ka += fKaTemp -> GetBinContent(l+1);
cout << "Integral_Ka=" << Integral_Ka << endl;

double Integral_Pr = 0;
for (int l=0;l<1000;l++) Integral_Pr += fPrTemp -> GetBinContent(l+1);
cout << "Integral_Pr=" << Integral_Pr << endl;
*/

double Integral_Pi = fPiTemp->Integral(-10000,10000) ;//* fPiTemp -> GetBinWidth(1);
double Integral_Ka = fKaTemp->Integral(-10000,10000) ;//* fKaTemp -> GetBinWidth(1);
double Integral_Pr = fPrTemp->Integral(-10000,10000) ;//* fPrTemp -> GetBinWidth(1);
double Integral_Mi = fMiTemp->Integral(-10000,10000) ;//* fMiTemp -> GetBinWidth(1);
double Inverse_Pi = 0;
double Inverse_Ka = 0;
double Inverse_Pr = 0;  
double Inverse_Mi = 0;  
if(Integral_Pi != 0){Inverse_Pi  = 1/Integral_Pi;}
else Inverse_Pi=0;
if(Integral_Ka != 0){Inverse_Ka  = 1/Integral_Ka;}
else Inverse_Ka=0;
if(Integral_Pr != 0){Inverse_Pr  = 1/Integral_Pr;}
else Inverse_Pr=0;
if(Integral_Mi != 0){Inverse_Mi  = 1/Integral_Mi;}
else Inverse_Mi=0;
fPiTemp -> Scale(Inverse_Pi);
fKaTemp -> Scale(Inverse_Ka);
fPrTemp -> Scale(Inverse_Pr);
fMiTemp -> Scale(Inverse_Mi);

if(SPECIE == 0){
	if(i < 10){
		h[i]->Fit(f3,"","",-10000,10000);
		   }
	else h[i]->Fit(f3_original,"","",-10000,10000);
}
else   h[i]->Fit(f3_original,"","",-10000,10000);

c1->SaveAs("prova.pdf","recreate");


if(SPECIE == 0){
	if(i < 10){
			 for (int j=0; j<4; j++) normalizations[j]->SetBinContent(i+1,f3->GetParameter(j));
 			 for (int j=0; j<4; j++) normalizations[j]->SetBinError(i+1,f3->GetParError(j));
			}
	else {
			 for (int j=0; j<4; j++) normalizations[j]->SetBinContent(i+1,f3_original->GetParameter(j));
 			 for (int j=0; j<4; j++) normalizations[j]->SetBinError(i+1,f3_original->GetParError(j));
		}
		}
else   {		 for (int j=0; j<4; j++) normalizations[j]->SetBinContent(i+1,f3_original->GetParameter(j));
 			 for (int j=0; j<4; j++) normalizations[j]->SetBinError(i+1,f3_original->GetParError(j));
	}


cout << "CICLO=" << i << endl << "SPECIE=" << SPECIE << endl;

} 

	c1->SaveAs(Form("prova%d.pdf", SPECIE),"recreate");
	TFile * ffiling=new TFile(Form("allfit_NEW_TEMPLATE%d.root",SPECIE),"recreate");
	for (int i=0; i<100; i++) h[i]->Write();
	for (int i=0; i<4;   i++) normalizations[i]->Write();
	ffiling->Close();

}  
}



