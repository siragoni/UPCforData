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
/*
TH1 *fPiTemp;
TH1 *fKaTemp;
TH1 *fPrTemp;
TH1 *fMiTemp;

Double_t fTOFResponseT(Double_t* x,Double_t* par)
{

  TF1* HELP = new TF1
  Double_t val  = par[0]*fPiTemp->Eval(x[0]);
  val += par[1]*fKaTemp->Eval(x[0]);
  val += par[2]*fPrTemp->Eval(x[0]);
//  val + = par[3]*fMiTemp->Eval(x[0]);

  double  bin = TMath::Floor(x[0]/0.05);
  Double_t val  = par[0]*fPiTemp->GetBinContent(bin);
  val += par[1]*fKaTemp->GetBinContent(bin);
  val += par[2]*fPrTemp->GetBinContent(bin);

//  delete HELP;
  return val;
}*/

Double_t fTOFResponseF(Double_t* x,Double_t* par)
{

return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2])* (x[0] < par[1]+par[9]*par[2]) + (x[0] > par[1]+par[9]*par[2])*par[0]*TMath::Exp(-(x[0]-par[1]-par[9]*par[2]*0.5)*par[9]/par[2]) + par[3]*TMath::Exp(-(x[0]-par[4])*(x[0]-par[4])/2/(par[2]*par[2] + par[5]))* (x[0] <par[4]+par[9]*sqrt(par[2]*par[2]+par[5])) + (x[0] >par[4]+par[9]*sqrt(par[2]*par[2]+par[5]))*par[3]*TMath::Exp(-(x[0]-par[4]-par[9]*sqrt(par[2]*par[2]+par[5])*0.5)*par[9]/sqrt(par[2]*par[2]+par[5]))+ par[6]*TMath::Exp(-(x[0]-par[7])*(x[0]-par[7])/2/(par[2]*par[2] + par[8]))* (x[0] <par[7]+par[9]*sqrt(par[2]*par[2]+par[8])) + (x[0] >par[7]+par[9]*sqrt(par[2]*par[2]+par[8]))*par[6]*TMath::Exp(-(x[0]-par[7]-par[9]*sqrt(par[2]*par[2]+par[8])*0.5)*par[9]/sqrt(par[2]*par[2]+par[8])) + par[10] +

par[11]*TMath::Exp(-(x[0]-par[12])*(x[0]-par[12])/2/(par[2]*par[2] + par[13]))* (x[0] <par[12]+par[9]*sqrt(par[2]*par[2]+par[12])) + (x[0] >par[12]+par[9]*sqrt(par[2]*par[2]+par[13]))*par[11]*TMath::Exp(-(x[0]-par[12]-par[9]*sqrt(par[2]*par[2]+par[13])*0.5)*par[9]/sqrt(par[2]*par[2]+par[13]))     +

par[14]*TMath::Exp(-(x[0]-par[15])*(x[0]-par[15])/2/(par[2]*par[2] + par[16]))* (x[0] <par[15]+par[9]*sqrt(par[2]*par[2]+par[16])) + (x[0] >par[15]+par[9]*sqrt(par[2]*par[2]+par[16]))*par[14]*TMath::Exp(-(x[0]-par[15]-par[9]*sqrt(par[2]*par[2]+par[16])*0.5)*par[9]/sqrt(par[2]*par[2]+par[16]))
;
}

Double_t fTOFResponseI(Double_t* x,Double_t* par)
{

return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2])* (x[0] < par[1]+par[3]*par[2]) + (x[0] > par[1]+par[3]*par[2])*par[0]*TMath::Exp(-(x[0]-par[1]-par[3]*par[2]*0.5)*par[3]/par[2]);
}

void combine(){

  //gROOT->Reset();
  //gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  Float_t mean_pi=0, mean_ka=600,mean_pr=2200,sigma=80;

  TFile* f1 = new TFile("all_histo_NEW_TEMPLATE.root");
  TCanvas *c1 = new TCanvas("c1","the fit canvas",500,400);
  TCanvas *c2 = new TCanvas("c2","the fit canvas",500,400);
  TCanvas *c3 = new TCanvas("c3","the fit canvas",500,400);
 // TCanvas *c4 = new TCanvas("c4","the fit canvas",500,400);
 // TCanvas *c5 = new TCanvas("c5","the fit canvas",500,400);
 // TCanvas *c6 = new TCanvas("c6","the fit canvas",500,400);


//	TF1* f = new TF1("myfunc",fTOFResponseF,-10000,10000,10);
	TF1* f = new TF1("myfunc",fTOFResponseF,-10000,10000,17);
 	f->SetNpx(5000);


	TF1* f2 = new TF1("myfunc2",fTOFResponseI,-5000,5000,4);
 	f2->SetNpx(2500);


/*	TF1* f3 = new TF1("myfunc3",fTOFResponseT,-10000,10000,3);
 	f3->SetNpx(5000);*/

// N.B.
//  shape_timeEXP_differences [nC] [0]   pioni rispetto ai kaoni
//  shape_timeEXP_differences [nC] [1]   pioni rispetto ai protoni
//  shape_timeEXP_differences [nC] [2]   kaoni rispetto ai pioni
//  shape_timeEXP_differences [nC] [3]   kaoni rispetto ai protoni
//  shape_timeEXP_differences [nC] [4]   protoni rispetto ai pioni
//  shape_timeEXP_differences [nC] [5]   protoni rispetto ai kaoni


   // pions
   f->SetParameter(0,1); // norm
   f->SetParameter(1,mean_pi); // peak
   f->SetParLimits(1,mean_pi -30,mean_pi+30);
   f->SetParameter(2,sigma); // sigma
   f->SetParLimits(2,50,1000);

   // kaons
   f->SetParameter(3,1); // norm
   f->SetParameter(4,mean_ka); // peak
//   f->SetParLimits(4,mean_ka-100,mean_ka+100);
   f->SetParameter(5,0); // sigma_add
   f->SetParLimits(5,0,500000);

   // protons
   f->SetParameter(6,1); // norm
   f->SetParameter(7,mean_pr); // peak
//   f->SetParLimits(7,mean_pr -100,mean_pr+100);
   f->SetParameter(8,0); // sigma_add
   f->SetParLimits(8,0,500000);

   f->SetParameter(9,1); // tof tail
   f->SetParLimits(9,0.8,5);








   f->SetParameter(11,1); // norm
   f->SetParameter(14,1); // norm
   f->SetParameter(13,0); // sigma_add
   f->SetParameter(16,0); // sigma_add



	c1->cd();


TH2F* new_template[3];


TH2F* ordinamento[9];
//ordinamento[0]=(TH2F*)f1->Get(Form("shape_timeEXP_differencesEXTRA0_%d",0));
ordinamento[1]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",6));
ordinamento[2]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",7));
/*ordinamento[3]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",2));
ordinamento[4]=(TH2F*)f1->Get(Form("shape_timeEXP_differencesEXTRA0_%d",1));
ordinamento[5]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",3));
ordinamento[6]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",4));
ordinamento[7]=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",5));
ordinamento[8]=(TH2F*)f1->Get(Form("shape_timeEXP_differencesEXTRA0_%d",2));
*/



for(int SPECIE=0; SPECIE<1; ++SPECIE)
{
	TH1F* fitted_parameters[17];
	TH1F* normalizations[3];

	for (int j=0; j<17; j++) fitted_parameters[j]=new TH1F(Form("fitted_parameters%d",j),"fitted_parameters",100,0, 5);
	for (int j=0; j<3; j++) normalizations[j]=new TH1F(Form("normalizations%d",j),"normalizations",100,0, 5);

	TH1F* integral_vs_pt = new TH1F("integral_vs_pt","integral_vs_pt", 100,0,5);
	new_template[SPECIE] = new TH2F(Form("new_template%d", SPECIE), Form("Tracks %d", SPECIE), 100, 0, 5, 2000, -1000, 1000);

TH1 *h[100] ;//
TH2F* help;
TH2F* help2;
TH1 *ordinamento_projection[100];
TH1 *ordinamento_projection2[100];
TH1 *shape_timeEXP_differences_projection[100] ;//
TH1 *shape_timeEXP_differences_projection2[100] ;//
for (int i=0;i<100;i++) 
{ 
//  h[i] = ((TH2F*)f1->Get(Form("timeEXP_differences0_%d",SPECIE))->ProjectionY(Form("bin%d",i+1),i+1,i+1);
  help=(TH2F*)f1->Get(Form("timeEXP_differences0_%d",SPECIE));
  h[i]= help->ProjectionY(Form("bin%d",i+1),i+1,i+1);
  //shape_timeEXP_differences_projection[i] = ((TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",2*SPECIE))->ProjectionY(Form("A%d",i+1),i+1,i+1);
help=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",2*SPECIE));
shape_timeEXP_differences_projection[i] = help->ProjectionY(Form("A%d",i+1),i+1,i+1);
help2=(TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",2*SPECIE+1));
  //shape_timeEXP_differences_projection2[i] = ((TH2F*)f1->Get(Form("shape_timeEXP_differences0_%d",2*SPECIE+1))->ProjectionY(Form("B%d",i+1),i+1,i+1);

/*c2 -> cd();
help2 -> Draw();
c2 -> SaveAs("provaB.pdf","recreate");*/
c3 -> cd();
shape_timeEXP_differences_projection2[i] = help2->ProjectionY(Form("B%d",i+1),i+1,i+1);
shape_timeEXP_differences_projection2[i] -> Draw();

ordinamento_projection[i]=ordinamento[1]->ProjectionY(Form("elec%d",i+1),i+1,i+1);
ordinamento_projection2[i]=ordinamento[2]->ProjectionY(Form("muon%d",i+1),i+1,i+1);

c3 -> SaveAs("provaC.pdf","recreate");
c1 -> cd();


   f->SetParameter(4,shape_timeEXP_differences_projection[i]->GetMean()); // peak
   f->SetParLimits(4,shape_timeEXP_differences_projection[i]->GetMean()-100,shape_timeEXP_differences_projection[i]->GetMean()+100);
   f->SetParameter(7,shape_timeEXP_differences_projection2[i]->GetMean()); // peak
   f->SetParLimits(7,shape_timeEXP_differences_projection2[i]->GetMean() -100,shape_timeEXP_differences_projection2[i]->GetMean()+100);
   
   f->SetParameter(12,ordinamento_projection[i]->GetMean()); // peak
   f->SetParLimits(12,ordinamento_projection[i]->GetMean() -100,ordinamento_projection[i]->GetMean()+100);
   f->SetParameter(15,ordinamento_projection2[i]->GetMean()); // peak
   f->SetParLimits(15,ordinamento_projection2[i]->GetMean() -100,ordinamento_projection2[i]->GetMean()+100);


  h[i]->SetLineColor(96);
  h[i]->GetXaxis()->SetTitle("time differences");
  h[i]->GetYaxis()->SetTitle("");
  h[i]->SetTitle("FIT shape EXP time differences");
  h[i]->Fit(f,"","",-6000,6000);


  for (int j=0; j<17; j++) fitted_parameters[j]->SetBinContent(i+1,f->GetParameter(j));
  for (int j=0; j<17; j++) fitted_parameters[j]->SetBinError(i+1,f->GetParError(j));


 c1->SaveAs("prova.pdf","recreate");



	f2->SetParameter(0, f->GetParameter(0));
	f2->SetParameter(1, f->GetParameter(1));
	f2->SetParameter(2, f->GetParameter(2));
	f2->SetParameter(3, f->GetParameter(9));
	f2->SetLineColor(4);
	f2->Draw("SAME");
	double a = f2->Integral(-4000,4000);
	double b = a / h[i]->GetBinWidth(1);
	integral_vs_pt->SetBinContent(i+1, b);



	for(int m=0; m<100000; m++)new_template[SPECIE]->Fill(i*0.05+0.025, f2->GetRandom());


/*

fPiTemp = ordinamento[3*SPECIE + 0]->ProjectionY(Form("bin%d",i+1),i+1,i+1);
fKaTemp = ordinamento[3*SPECIE + 1]->ProjectionY(Form("bin%d",i+1),i+1,i+1);
fPrTemp = ordinamento[3*SPECIE + 2]->ProjectionY(Form("bin%d",i+1),i+1,i+1);
double Integral_Pi = fPiTemp->Integral(-10000,10000);
double Integral_Ka = fKaTemp->Integral(-10000,10000);
double Integral_Pr = fPrTemp->Integral(-10000,10000);
double Inverse_Pi  = 1/Integral_Pi;
double Inverse_Ka  = 1/Integral_Ka;
double Inverse_Pr  = 1/Integral_Pr;
fPiTemp -> Scale(Inverse_Pi);
fKaTemp -> Scale(Inverse_Ka);
fPrTemp -> Scale(Inverse_Pr);

cout << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

h[i]->Fit(f3,"","",-8000,8000);
c1->SaveAs("prova.pdf","recreate");

  for (int j=0; j<3; j++) normalizations[j]->SetBinContent(i+1,f3->GetParameter(j));
  for (int j=0; j<3; j++) normalizations[j]->SetBinError(i+1,f3->GetParError(j));
*/

} 

	c1->SaveAs(Form("prova%d.pdf", SPECIE),"recreate");
	TFile * ffiling=new TFile(Form("allfit%d.root",SPECIE),"recreate");
	for (int i=0; i<100; i++) h[i]->Write();
	for (int i=0; i<100; i++) shape_timeEXP_differences_projection[i]->Write();
	for (int i=0; i<100; i++) shape_timeEXP_differences_projection2[i]->Write();
	for (int i=0; i<100; i++) ordinamento_projection[i]->Write();
	for (int i=0; i<100; i++) ordinamento_projection2[i]->Write();
	for (int j=0; j<17;  j++) fitted_parameters[j]->Write();
	for (int i=0; i<3;   i++) normalizations[i]->Write();



	new_template[SPECIE]->Write();
//	new_template[1]->Write();
//	new_template[2]->Write();

	integral_vs_pt->Write();
	ffiling->Close();

}  
}



