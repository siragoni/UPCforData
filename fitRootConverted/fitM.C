#include "env.h"
#include "TPad.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"

Bool_t ds = 0;

TF1* fCB1;
TF1* fCB2;
TF1* fBG;
Double_t widthRatio;

Double_t m1sPDG = 3.096900;
Double_t m2sPDG = 3.688097;

double bg(double* xx, double* p){
  double x  = xx[0];
  double y = x-p[2];
  return p[0]*p[1]*exp(-p[1]*x)*(y>=0 ? 1 : 1+p[3]*y*y+p[4]*y*y*y+p[5]*y*y*y*y);
}

double cb(double* xx, double* p){
  double x  = xx[0];
  double m  = p[1];
  double s  = p[2];
  double a1 = p[3];
  double n1 = p[4];
  double a2 = p[5];
  double n2 = p[6];
  if (s<=0) return 0;
  if (n1<1) return 0;
  if (n2<1) return 0;
  double z  = (x-m)/s;
  double na1 = n1/a1;
  double na2 = n2/a2;
  double C = n1/a1*1./(n1-1)*exp(-a1*a1/2.);
  double D = sqrt(M_PI/2.)*(1.+ROOT::Math::erf(a1/sqrt(2.)));
  double N = p[0]/s/(C+D);
  if (z<-a1) return N*exp(-0.5*a1*a1)*pow(na1/(na1-a1-z),n1);
//  if (z>+a2) return N*exp(-0.5*a2*a2)*pow(na2/(na2-a2+z),n2);
  return N*exp(-0.5*z*z);
}

double fsum2(double *x, double *par){
  Double_t par1s[7];
  Double_t par2s[7];
  for (Int_t i=0;i<7;i++) par1s[i]=par[i+6];
  par2s[0]=par[13]; // yieldRatio
//  par2s[0]=par[13]*par1s[0]; // yieldRatio
  par2s[1]=par[14]+par1s[1]; // m2s = m1s + m2s_PDG - m1s_PDG
//  par2s[1]=par[14];            // m2s = m2s_PDG
  par2s[2]=par[15]*par1s[2]; // widthRatio
  par2s[3]=par[16];
  par2s[4]=par[17];
  par2s[5]=par[18];
  par2s[6]=par[19];
  return bg(x,par)+cb(x,par1s)+cb(x,par2s);
}

double fsum3(double *x, double *par){
  Double_t par1s[7];
  for (Int_t i=0;i<7;i++) par1s[i]=par[i+6];
  return bg(x,par)+cb(x,par1s);
}

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
//void fitM(TString suffix="gl",Int_t iy=0, Int_t ipt=0, Int_t rebin=1,Bool_t debug=0){
void fitM(TString suffix="data_noTKLcut",Int_t iy=0, Int_t ipt=0, Int_t rebin=4,Bool_t debug=1){


//void fitM(TString suffix="1c_15",Int_t iy=0, Int_t ipt=0, Int_t rebin=1,Bool_t debug=0){
//void fitM(TString suffix="data_noTKLcut",Int_t iy=0, Int_t ipt=0, Int_t rebin=4,Bool_t debug=1){
//void fitM(TString suffix="gl_18",Int_t iy=0, Int_t ipt=0, Int_t rebin=4,Bool_t debug=1){
//void fitM(TString suffix="gl_18wgtData",Int_t iy=0, Int_t ipt=0, Int_t rebin=4,Bool_t debug=1){
//void fitM(TString suffix="1r_18",Int_t iy=0, Int_t ipt=0, Int_t rebin=4,Bool_t debug=1){
//void fitM(TString suffix="gl",Int_t iy=0, Int_t ipt=0, Int_t rebin=1, Bool_t debug=1){


//void fitM(TString suffix="gl",Int_t iy=0, Int_t ipt=0, Int_t rebin=1, Bool_t debug=1){
//void fitM(TString suffix="2c",Int_t iy=4, Int_t ipt=0, Int_t rebin=1, Bool_t debug=1){
//void fitM(TString suffix="1c_15",Int_t iy=2, Int_t ipt=0, Int_t rebin=1, Bool_t debug=1){
//void fitM(TString suffix="data",Int_t iy=0, Int_t ipt=0, Int_t rebin=4,Bool_t debug=1){
//void fitM(TString suffix="data_noTKLcut",Int_t iy=9, Int_t ipt=6, Int_t rebin=4,Bool_t debug=1){
//void fitM(TString suffix="data",Int_t iy=9, Int_t ipt=0, Int_t rebin=4,Bool_t debug=1){
  SetStyle();
  TFile* f = new TFile(Form("templates/%s.root",suffix.Data()));
  TH1D* hM = (TH1D*) f->Get(Form("hM%i%i",iy,ipt));
  TString setName = "";

  Double_t mMin=2.2;
  Double_t mMax=6.;

  hM->GetXaxis()->SetRangeUser(mMin,mMax);
  hM->Sumw2();

  if (suffix.Contains("g")){
    /*
    TFile* fMc   = new TFile("templates/glwgtMc.root");
    TFile* fData = new TFile("templates/glwgtData.root");
    TH1D* hMcM    = (TH1D*) fMc->Get("hM00");
    TH1D* hDataM  = (TH1D*) fData->Get("hM00");
    hDataM->Divide(hDataM,hMcM,1,1,"B");
    hM->Multiply(hDataM);
    */
    hM->GetXaxis()->SetRangeUser(2,8);
    TF1* f1 = new TF1("f1","[0]*[1]*exp(-[1]*x)",suffix.Contains("gh") ? 5 : 3,8);
    f1->SetParameter(0,10000);
    f1->SetParameter(1,1.0);
    hM->Fit(f1,"0NLRQ","");
    TF1* f2 = new TF1("f2","[0]*[1]*exp(-[1]*x)*(x>=[2] ? 1 : 1+[3]*(x-[2])^2+[4]*(x-[2])^3+[5]*(x-[2])^4)",suffix.Contains("gh") ? 5 : 2,8);
    f2->SetParameter(0,f1->GetParameter(0));
    f2->SetParameter(1,f1->GetParameter(1));
    f2->FixParameter(2,4);
    f2->SetNpx(1000);
    hM->Fit(f2,"L","",suffix.Contains("gh") ? 5 : 2,8);
  } else if (suffix.Contains("1r") ||suffix.Contains("1c") || suffix.Contains("2c") || suffix.Contains("1i") || suffix.Contains("2i") ||suffix.Contains("f")){
    Double_t intall  = hM->Integral(hM->GetXaxis()->FindFixBin(2.0),hM->GetXaxis()->FindFixBin(6.0));
    Double_t binWidth = hM->GetBinWidth(1);
    TF1* fCB = new TF1("fCB",cb,2,15,7);
    fCB->SetParameter(0,intall*binWidth);
    fCB->SetParameter(1,(suffix.Contains("2c") || suffix.Contains("2i")) ? 3.686 : 3.096);
    fCB->SetParameter(2,0.072);
    fCB->SetParameter(3,1.0);
    fCB->SetParameter(4,1.6e+7);
    fCB->SetParameter(5,2.5);
    fCB->SetParameter(6,2.0e+7);
    if (!ds){
      fCB->FixParameter(5,10);
      fCB->FixParameter(6,2.0e+7);
    }
    fCB->SetNpx(1000);

    hM->Fit(fCB,"ML","",mMin,mMax);
    printf("%f %f %f\n",fCB->GetParameter(0),fCB->Integral(2.85,3.35,1.e-4)*5,hM->Integral());
    /*
    TF1* fCB = new TF1("fCB","[0]*ROOT::Math::crystalball_pdf(x,[3],[4],[2],[1])",2,15);
    fCB->SetParameter(0,intall*binWidth);
    fCB->SetParameter(1,(suffix.Contains("2c") || suffix.Contains("2i")) ? 3.686 : 3.096);
    fCB->SetParameter(2,0.072);
    fCB->SetParameter(3,1.0);
    fCB->SetParameter(4,1e+7);
    fCB->SetNpx(1000);
//    hM->Rebin(2);
//    hM->GetXaxis()->SetRangeUser(mMin,mMax);
//    if (suffix.Contains("1")) hM->Fit(fCB,"ML","",2.80,3.35);
//    else                      hM->Fit(fCB,"ML","",3.40,3.90);
    if (suffix.Contains("1")) hM->Fit(fCB,"ML","",mMin,mMax);
    else                      hM->Fit(fCB,"ML","",mMin,mMax);
    */
  } else {
    hM->Rebin(rebin);
    Double_t intjpsi = hM->Integral(hM->GetXaxis()->FindFixBin(2.85),hM->GetXaxis()->FindFixBin(3.35));
    Double_t intall  = hM->Integral(hM->GetXaxis()->FindFixBin(2.0),hM->GetXaxis()->FindFixBin(6.0));
    Double_t intbg = intall-intjpsi;
    Double_t binWidth = hM->GetBinWidth(1);

    TFile* fin = new TFile("fitM.root");
    //fin->ls();
    //return;
    TString syear = suffix.Contains("15") ? "_15" : (suffix.Contains("18") ? "_18" : "");
    TH1D* hM_1c = (TH1D*) fin->Get(Form("1c%s_%i_%i_1",syear.Data(),iy,ipt));
    TH1D* hM_2c = (TH1D*) fin->Get(Form("2c%s_%i_%i_1",syear.Data(),iy,ipt));
    TH1D* hM_gl = (TH1D*) fin->Get(Form("gl%s_%i_%i_1",syear.Data(),iy,ipt));
    fBG  = new TF1("fBG" ,bg,2,15,6);
    fCB1 = new TF1("fCB1",cb,2,15,7);
    fCB2 = new TF1("fCB2",cb,2,15,7);
    fBG->SetParameters(hM_gl->GetFunction("f2")->GetParameters());
    fCB1->SetParameters(hM_1c->GetFunction("fCB")->GetParameters());
    fCB2->SetParameters(hM_2c->GetFunction("fCB")->GetParameters());
    fBG->SetNpx(1000);
    fCB1->SetNpx(1000);
    fCB2->SetNpx(1000);
    widthRatio = fCB2->GetParameter(2)/fCB1->GetParameter(2);
    printf("widthRatio = %f\n",widthRatio);

    Double_t m2s = m2sPDG;
    TF1 *fSum = new TF1("fSum2",fsum2,1,15,20);
    fSum->SetParameter( 0,intbg*binWidth);        // Kay: no fix
    fSum->SetParameter( 1,fBG->GetParameter(1));  // Kay: no fix
    fSum->FixParameter( 2,fBG->GetParameter(2));  // Kay: fix
    fSum->FixParameter( 3,fBG->GetParameter(3));  // Kay: fix
    fSum->FixParameter( 4,fBG->GetParameter(4));  // Kay: fix
    fSum->FixParameter( 5,fBG->GetParameter(5));  // Kay: fix
    fSum->SetParameter( 6,intjpsi);             // Kay: no fix
//    fSum->SetParameter( 6,intjpsi*5);             // Kay: no fix
    fSum->SetParameter( 7,fCB1->GetParameter(1)); // Kay: no fix
    fSum->SetParameter( 8,fCB1->GetParameter(2)); // Kay: no fix
    fSum->FixParameter( 9,fCB1->GetParameter(3)); // Kay: fix
    fSum->FixParameter(10,fCB1->GetParameter(4)); // Kay: fix to 10
    fSum->FixParameter(11,fCB1->GetParameter(5)); // 
    fSum->FixParameter(12,fCB1->GetParameter(6)); // 
    fSum->SetParameter(13,0.02*intjpsi);           // Kay: no fix
//    fSum->SetParameter(13,0.02);                  // Kay: no fix
    fSum->FixParameter(14,m2sPDG-m1sPDG);         // Kay: fix 3.686097 +/-0.000025
//    fSum->FixParameter(14,m2sPDG);         // Kay: fix 3.686097 +/-0.000025
//    fSum->FixParameter(15,1.0);            // Kay: fix ratio
    fSum->FixParameter(15,widthRatio*0.9);            // Kay: fix ratio
    fSum->FixParameter(16,fCB2->GetParameter(3)); // Kay: fix
    fSum->FixParameter(17,fCB2->GetParameter(4)); // Kay: fix to 10
    fSum->FixParameter(18,fCB2->GetParameter(5)); // 
    fSum->FixParameter(19,fCB2->GetParameter(6)); // 
//    fSum->SetParLimits(13,0.005,0.060);
    fSum->SetParLimits(13,0.1,10000);
    fSum->SetParLimits(6 ,1,1000000);
    fSum->SetParLimits(7 ,3.0,3.2);
    fSum->SetParLimits(8 ,0.06,0.10);
//    fSum->SetParLimits(9 ,0.1,5.2);
//    fSum->SetParLimits(10,1.1,1e8);

    fSum->SetNpx(1000);
    
    printf("%f %f\n",fBG->GetParameter(1),hM_gl->GetFunction("f2")->GetParError(1));
    
    TFitResultPtr res = hM->Fit(fSum,debug ? "NLM" : "NLQM","",2.2,mMax);
    const Double_t* par = fSum->GetParameters();
    const Double_t* err = fSum->GetParErrors();
    fBG->SetParameters(&(par[0]));
    fCB1->SetParameters(&(par[6]));
    fCB2->SetParameters(&(par[13]));
    fBG->SetParErrors(&(err[0]));
    fCB1->SetParErrors(&(err[6]));
    fCB2->SetParErrors(&(err[6]));
//    fCB2->SetParameter(0,par[6]*par[13]);
    fCB2->SetParameter(0,par[13]);
//    fCB2->SetParameter(1,par[14]);
    fCB2->SetParameter(1,par[7]+par[14]);
    fCB2->SetParameter(2,par[8]*par[15]);
    fCB2->SetParError(0,err[13]);
//    fCB2->SetParError(0,TMath::Sqrt(err[6]*err[6]/par[6]/par[6]+err[13]*err[13]/par[13]/par[13])*par[6]*par[13]);
    fCB2->SetParError(2,err[8]);
    hM->GetListOfFunctions()->Add(fBG);
    hM->GetListOfFunctions()->Add(fCB1);
    hM->GetListOfFunctions()->Add(fCB2);
    hM->GetListOfFunctions()->Add(fSum);

    Double_t chi2ndf = fSum->GetChisquare()/fSum->GetNDF();
    Double_t chi2 = fSum->GetChisquare();
    Int_t dof = fSum->GetNDF();
    printf("jpsi=%.1f chi2ndf=%.2f\n",fCB1->Integral(2.,6.,1e-6)/binWidth,chi2ndf);

    if (debug) {
      Double_t lumi = suffix.Contains("15") ? 216 : (suffix.Contains("18") ? 538 : 754);
      Double_t r       = par[13];
      Double_t r_err   = err[13];
      Double_t n1s     = fCB1->Integral(2.,6.,1e-6)/binWidth;
      Double_t n2s     = fCB2->Integral(2.,6.,1e-6)/binWidth;
      Double_t n1s_err = fCB1->GetParError(0)/fCB1->GetParameter(0)*n1s;
      Double_t n2s_err = fCB2->GetParError(0)/fCB2->GetParameter(0)*n2s;
      Double_t nbg     = fBG->Integral(2.85+0.001,3.35-0.001)/binWidth;
      Double_t nbg_err = fBG->GetParError(0)/fBG->GetParameter(0)*nbg;

//      TCanvas* c = new TCanvas("c","c",285*2,500);
      TCanvas* c = new TCanvas("c","c",285*3,750);
      c->SetRightMargin(0.01);
      c->SetLeftMargin(0.12);
      c->SetTopMargin(0.01);
      c->SetBottomMargin(0.11);
      hM->SetTitle("");
      hM->GetXaxis()->SetTitleOffset(1.1);
      hM->GetYaxis()->SetTitleOffset(1.3);
      hM->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}",binWidth*1000));
      hM->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
      hM->GetXaxis()->SetTitleSize(0.045);
      hM->GetYaxis()->SetTitleSize(0.045);
      hM->SetMarkerStyle(kFullCircle);
      hM->SetMarkerSize(0.7);
      hM->SetMaximum(hM->GetMaximum()*1.2);
      hM->GetXaxis()->SetRangeUser(2.0,6);
      hM->Draw("e0");
      fBG->SetLineColor(kGreen+1);
      fCB1->SetLineColor(kMagenta);
      fCB2->SetLineColor(kRed);
      fSum->SetLineColor(kBlue);
      fBG->SetLineStyle(9);
      fBG->Draw("same");
      fCB1->Draw("same");
      fCB2->Draw("same");
      fSum->Draw("same");
      gPad->RedrawAxis();
      TLatex* latex = new TLatex();
      latex->SetNDC();
      latex->SetTextSize(0.045);
      latex->SetTextAlign(22);
      latex->DrawLatex(0.56,0.94,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
//      latex->DrawLatex(0.56,0.94,"ALICE Performance, Pb-Pb #sqrt{s_{NN}} = 5.02 TeV");
      latex->SetTextSize(0.045);
      latex->SetTextAlign(12);
      latex->DrawLatex(0.50,0.85,Form("UPC, L_{#lower[-0.3]{int}} = %.0f #pm %.0f #mub^{-1}",lumi,lumi*0.05));
      latex->DrawLatex(0.50,0.79,"#it{p}_{#lower[-0.3]{T}} < 0.25 GeV/#it{c}");
      latex->DrawLatex(0.50,0.73,Form("#minus%.2f < #it{y} < #minus%.2f",-gYMin[iy],-gYMax[iy]));
      latex->DrawLatex(0.50,0.67,Form("#it{N}_{J/#psi} = %.0f #pm %.0f",n1s,n1s_err));
      latex->DrawLatex(0.50,0.61,Form("#it{N}_{#psi'} = %.0f #pm %.0f",n2s,n2s_err));
      latex->DrawLatex(0.50,0.55,Form("#chi^{2}/#it{dof} = %.2f (%.1f/%i)",chi2ndf,chi2,dof));
//      latex->DrawLatex(0.50,0.49,Form("N_{#psi(2S)}/N_{J/#psi} = %.4f #pm %.4f\n",r,r_err));
//      latex->DrawLatex(0.50,0.43,Form("N_{bg}(2.85,3.35) = %.0f #pm %.0f\n",nbg,nbg_err));
//      gPad->Print(Form("jpsi_upc_%i.png",iy));
//      gPad->Print(Form("jpsi_upc_%i.eps",iy));
//      gPad->Print(Form("jpsi_upc_%i.pdf",iy));
      gPad->Print(Form("m%i.pdf",iy));
    }
  }
  
  TFile* fout = new TFile("fitM.root","update");
  hM->Write(Form("%s_%i_%i_%i",suffix.Data(),iy,ipt,rebin),TObject::kOverwrite);
  fout->Close();
}

