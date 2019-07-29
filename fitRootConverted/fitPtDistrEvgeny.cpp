#include "env.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
TH1D* hPt;
TH1D* h1c;
TH1D* h1i;
TH1D* hfc;
TH1D* hfi;
TH1D* hgl;
TH1D* hun;
TF1* fgu;

double ffsum(double *x, double *par){
  Double_t pt=x[0];
  Double_t n1c=par[0];
  Double_t n1i=par[1];
  Double_t ngl=par[2];
  Double_t fd =par[3];
  Double_t nun=par[4];
  Double_t sum=0;
  sum+=hgl->GetBinContent(hgl->FindFixBin(pt))*ngl;
  sum+=h1c->GetBinContent(h1c->FindFixBin(pt))*n1c;
  sum+=h1i->GetBinContent(h1i->FindFixBin(pt))*n1i;
  sum+=hfc->GetBinContent(hfc->FindFixBin(pt))*n1c*fd;
  sum+=hfi->GetBinContent(hfi->FindFixBin(pt))*n1i*fd;
  sum+=hun->GetBinContent(hun->FindFixBin(pt))*nun;
  return sum;
}

void fitPt(TString suffix="data_noTKLcut", Int_t iy=6, Int_t im=0, Double_t rSigma2s1s = 0.150, Int_t debug=1){
  SetStyle();
  gStyle->SetLineScalePS(2);
  TFile* feff = new TFile(Form("efficiency%s.root", suffix.Contains("15") ? "15" : (suffix.Contains("18") ? "18" : "")));
  TH1D* hEff1cM = (TH1D*) feff->Get(Form("h1cEffM%i",iy));
  TH1D* hEfffcM = (TH1D*) feff->Get(Form("hfcEffM%i",iy));
  double brfd = 0.614;
  double eff1c_m = hEff1cM->GetBinContent(im+1);
  double efffd_m = hEfffcM->GetBinContent(im+1);
  double fd = rSigma2s1s*efffd_m*brfd/eff1c_m;

  printf("FD in: %.2f<y<%.2f: %.4f\n",gYMin[iy],gYMax[iy],fd);

  TFile* filedata = new TFile(Form("templates/%s.root",suffix.Data()));
  TFile* file1c = new TFile(Form("templates/1c%s.root",suffix.Contains("15") ? "_15" : (suffix.Contains("18") ? "_18" : "")));
  TFile* file1i = new TFile(Form("templates/1i%s.root",suffix.Contains("15") ? "_15" : (suffix.Contains("18") ? "_18" : "")));
  TFile* filefc = new TFile(Form("templates/fc%s.root",suffix.Contains("15") ? "_15" : (suffix.Contains("18") ? "_18" : "")));
  TFile* filefi = new TFile(Form("templates/fi%s.root",suffix.Contains("15") ? "_15" : (suffix.Contains("18") ? "_18" : "")));
  TFile* filegl = new TFile(Form("templates/gl%s.root",suffix.Contains("15") ? "_15" : (suffix.Contains("18") ? "_18" : "")));

  TH1D* hPtSideBandLow  = (TH1D*) filedata->Get(Form("hPt%i%i",iy,nMBins-2));
  TH1D* hPtSideBandHigh = (TH1D*) filedata->Get(Form("hPt%i%i",iy,nMBins-1));
  hPt = (TH1D*) filedata->Get(Form("hPt%i%i",iy,im));
  h1c = (TH1D*) file1c->Get(Form("hPt%i%i",iy,im));
  h1i = (TH1D*) file1i->Get(Form("hPt%i%i",iy,im));
  hfc = (TH1D*) filefc->Get(Form("hPt%i%i",iy,im));
  hfi = (TH1D*) filefi->Get(Form("hPt%i%i",iy,im));
  hgl = (TH1D*) filegl->Get(Form("hPt%i%i",iy,im));

  if (debug==5) hgl = hPtSideBandLow;
  if (debug==6) hgl = hPtSideBandHigh;
  hPt->Rebin(5);
  h1c->Rebin(5);
  h1i->Rebin(5);
  hfc->Rebin(5);
  hfi->Rebin(5);
  hgl->Rebin(5);

  hPt->Sumw2();
//  new TCanvas;
//  hPtSideBandLow->SetLineColor(kBlack);
//  hPtSideBandHigh->SetLineColor(kBlue);
//  hPtSideBandLow->Draw();
//  hPtSideBandHigh->Draw("same");
//
//  hPtSideBandLow->
//  return;

  Int_t nBinsX = hPt->GetNbinsX();

  TF1* fun = new TF1("fun","[0]*x*(1+[1]/[2]*x*x)^(-[2])",0,4);
  fun->SetParameter(0,1);
//  fun->SetParameter(1,debug==4 ? 1.25 : 1.);
//  fun->SetParameter(2,debug==4 ? 6.1 : 1.);
  fun->SetParameter(1,debug==4 ? 1.6 : 1.79);
  fun->SetParameter(2,debug==4 ? 3.58 : 3.58);
  fun->SetNpx(nBinsX);
  hun = (TH1D*) fun->GetHistogram()->Clone("hun");
  for (Int_t ibin=1;ibin<=hun->GetNbinsX();ibin++) hun->SetBinError(ibin,0);

  h1c->Scale(1./h1c->Integral(1,h1c->GetNbinsX()));
  h1i->Scale(1./h1i->Integral(1,h1i->GetNbinsX()));
  hfc->Scale(1./hfc->Integral(1,hfc->GetNbinsX()));
  hfi->Scale(1./hfi->Integral(1,hfi->GetNbinsX()));
  hgl->Scale(1./hgl->Integral(1,hgl->GetNbinsX()));
  hun->Scale(1./hun->Integral(1,hun->GetNbinsX()));

  h1c->SetLineColor(kBlue);
  h1i->SetLineColor(kRed);
  hfc->SetLineColor(kCyan);
  hfi->SetLineColor(kYellow+1);
  hgl->SetLineColor(kGreen);
  hun->SetLineColor(kMagenta);

  h1c->SetMarkerColor(kBlue);
  h1i->SetMarkerColor(kRed);
  hfc->SetMarkerColor(kCyan);
  hfi->SetMarkerColor(kYellow+1);
  hgl->SetMarkerColor(kGreen);
  hun->SetMarkerColor(kMagenta);

  h1c->SetLineWidth(2);
  h1i->SetLineWidth(2);
  hfc->SetLineWidth(2);
  hfi->SetLineWidth(2);
  hgl->SetLineWidth(2);
  hun->SetLineWidth(2);

  // Get expected number of gg events from invariant mass fit (with pt cut)
  TFile* ffitM = new TFile("fitM.root");
  TH1D* hM = (TH1D*) ffitM->Get(Form("%s_%i_0_4",suffix.Data(),iy));
  TF1* fBG = (TF1*) hM->GetListOfFunctions()->FindObject("fBG");
  Double_t n_val = fBG->Integral(gMMin[im]+0.001,gMMax[im]-0.001)/hM->GetBinWidth(1);
  Double_t n_err = fBG->GetParError(0)/fBG->GetParameter(0)*n_val;
  if (debug==2) n_val+=n_err;
  if (debug==3) n_val-=n_err;
  printf("nbg=%f\n",n_val);
  // Estimate fraction of gg events below pt cut
  Double_t ngl = n_val/hgl->Integral(1,hgl->FindFixBin(gPtMax[0]))*hgl->Integral(1,hgl->GetNbinsX());

  TF1* fsum = new TF1("fsum",ffsum,0,4,5);
  fsum->SetNpx(2000);
  fsum->SetParameter(0,4000);
  fsum->SetParameter(1,1000);
  fsum->FixParameter(2,ngl);
  fsum->FixParameter(3,fd);
  fsum->SetParameter(4,100);
  hPt->Fit(fsum,debug?"":"Q");

  Double_t n1c     = fsum->GetParameter(0);
  Double_t n1i     = fsum->GetParameter(1);
  Double_t nun     = fsum->GetParameter(4);
  Double_t n1c_err = fsum->GetParError(0);
  Double_t n1i_err = fsum->GetParError(1);

  Double_t chi2ndf = fsum->GetChisquare()/fsum->GetNDF();
  printf("chi2/ndf=%.2f\n",chi2ndf);
  h1c->Scale(n1c);
  h1i->Scale(n1i);
  hfc->Scale(n1c*fd);
  hfi->Scale(n1i*fd);
  hgl->Scale(ngl);
  hun->Scale(nun);

//  hPt->SetTitle(Form("%s, %.1f < y <%.1f, %.2f < M (GeV/#it{c}^{2}) < %.2f, fd=%.2f\n",suffix.Data(),gYMin[iy],gYMax[iy],gMMin[im],gMMax[im],fd));
  hPt->SetTitle("");
  hPt->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}",1000*hPt->GetBinWidth(1)));
  hPt->GetXaxis()->SetTitle("Dimuon #it{p}_{T} (GeV/#it{c})");

  TFile* fout = new TFile("fitPt.root","update");
  h1c->Write(Form("h1c_%s_%i_%i_%.0f_%s",suffix.Data(),iy,im,rSigma2s1s*1000,debug>1 ? Form("%i",debug) : ""),TObject::kOverwrite);
  h1i->Write(Form("h1i_%s_%i_%i_%.0f_%s",suffix.Data(),iy,im,rSigma2s1s*1000,debug>1 ? Form("%i",debug) : ""),TObject::kOverwrite);
  hfc->Write(Form("hfc_%s_%i_%i_%.0f_%s",suffix.Data(),iy,im,rSigma2s1s*1000,debug>1 ? Form("%i",debug) : ""),TObject::kOverwrite);
  hfi->Write(Form("hfi_%s_%i_%i_%.0f_%s",suffix.Data(),iy,im,rSigma2s1s*1000,debug>1 ? Form("%i",debug) : ""),TObject::kOverwrite);
  hgl->Write(Form("hgl_%s_%i_%i_%.0f_%s",suffix.Data(),iy,im,rSigma2s1s*1000,debug>1 ? Form("%i",debug) : ""),TObject::kOverwrite);
  hun->Write(Form("hun_%s_%i_%i_%.0f_%s",suffix.Data(),iy,im,rSigma2s1s*1000,debug>1 ? Form("%i",debug) : ""),TObject::kOverwrite);
  hPt->Write(Form("hPt_%s_%i_%i_%.0f_%s",suffix.Data(),iy,im,rSigma2s1s*1000,debug>1 ? Form("%i",debug) : ""),TObject::kOverwrite);
//  fsum->Write(Form("fitPt_%s_%i_%i_%i",suffix.Data(),iy,im,ipt),TObject::kOverwrite);
  fout->Close();

  if (!debug) return;

  TCanvas* c = new TCanvas("c","c",285*3,750);
  c->SetRightMargin(0.01);
  c->SetLeftMargin(0.12);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.11);

  hPt->SetLineColor(kBlack);
  hPt->SetMarkerStyle(kFullCircle);
  hPt->SetMarkerSize(0.7);
  fsum->SetLineColor(kBlack);
//
  Bool_t log=1;
  if (log) {
    gPad->SetLogy();
    hPt->GetXaxis()->SetRangeUser(0,2.7);
    hPt->SetMaximum(hPt->GetBinContent(hPt->GetMaximumBin())*3);
    hPt->SetMinimum(2.);
  } else {
    gPad->SetLeftMargin(0.14);
    hPt->GetXaxis()->SetRangeUser(0,1.5);
    hPt->GetYaxis()->SetTitleOffset(1.7);
    hPt->SetMaximum(hPt->GetBinContent(hPt->GetMaximumBin())*1.2);
    hPt->SetMinimum(0.);
  }
  hPt->GetXaxis()->SetTitleOffset(1.17);
  hPt->GetXaxis()->SetTitleSize(0.046);
  hPt->GetYaxis()->SetTitleSize(0.0468);

  hPt->GetListOfFunctions()->Clear();
  hPt->Draw("e0");
  h1c->Scale(1./n1c);
  h1i->Scale(1./n1i);
  hfc->Scale(1./n1c/fd);
  hfi->Scale(1./n1i/fd);
  hgl->Scale(1./ngl);
  hun->Scale(1./nun);
  fsum->DrawClone("same");
  h1c->Scale(n1c);
  h1i->Scale(n1i);
  hfc->Scale(n1c*fd);
  hfi->Scale(n1i*fd);
  hgl->Scale(ngl);
  hun->Scale(nun);

  Double_t int1c = h1c->Integral(1,h1c->FindFixBin(gPtMax[0]-0.0001));
  Double_t int1i = h1i->Integral(1,h1c->FindFixBin(gPtMax[0]-0.0001));
  Double_t intfc = hfc->Integral(1,h1c->FindFixBin(gPtMax[0]-0.0001));
  Double_t intfi = hfi->Integral(1,h1c->FindFixBin(gPtMax[0]-0.0001));
  Double_t intun = hun->Integral(1,h1c->FindFixBin(gPtMax[0]-0.0001));
  Double_t all = int1i+intfc+intfi+intun+int1c;
  Double_t ratio = (int1i+intfc+intfi+intun)/int1c;
  Double_t fI = (int1i+intun)/int1c;
  Double_t fD = (intfc+intfi)/int1c;
  printf("%f\n",int1c);
  printf("%f\n",int1i);
  printf("%f\n",intfc);
  printf("%f\n",intfi);
  printf("%f\n",intun);
  printf("%f\n",all);
  printf("%f\n",fI);
  printf("%f\n",fD);
  printf("%f\n",ratio);

  h1c->Draw("hist same");
  h1i->Draw("hist same");
  hfc->Draw("hist same");
  hfi->Draw("hist same");
  hgl->Draw("hist same");
  hun->Draw("hist same");

  Double_t lumi = suffix.Contains("15") ? 216 : (suffix.Contains("18") ? 533 : 749);

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextAlign(22);
  latex->DrawLatex(0.56,0.95,Form("ALICE, UPC Pb-Pb #sqrt{s_{NN}} = 5.02 TeV, L_{#lower[-0.3]{int}} #approx %.0f #mub^{-1}",lumi));
  latex->SetTextSize(0.035);
  latex->SetTextAlign(12);
  latex->DrawLatex(0.20,0.88,Form("%.2f < y < %.2f",gYMin[iy],gYMax[iy]));
  latex->DrawLatex(0.46,0.88,Form("%.2f < #it{m}_{#mu#mu} < %.2f GeV/#it{c}^{2}",gMMin[im],gMMax[im]));

  if (0) {
    TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
//    TLegend* l = new TLegend(0.55,0.50,0.98,0.84);
    l->SetMargin(0.1);
    l->AddEntry(hPt,"ALICE data");
    l->AddEntry(fsum,Form("Fit: #chi^{2}/NDF=%.1f\n",chi2ndf));
    l->AddEntry(h1c,Form("Coherent J/#psi: %.0f #pm %.0f",n1c,n1c_err));
    l->AddEntry(h1i,Form("Incoherent J/#psi: %.0f #pm %.0f",n1i,n1i_err));
    l->AddEntry(hun,Form("Incoherent dissocitive J/#psi: %.0f",nun));
    l->AddEntry(hfc,"Coherent #psi(2S) feeddown");
    l->AddEntry(hfi,"Incoherent #psi(2S) feeddown");
    l->AddEntry(hgl,Form("Continuum #gamma#gamma #rightarrow #mu#mu: %.0f",ngl));
    l->Draw();
  }  else {
    TLegend* l = new TLegend(0.45,0.55,0.98,0.85);
    l->SetMargin(0.1);
    l->AddEntry(hPt,"ALICE data");
    l->AddEntry(h1c,"Coherent J/#psi");
    l->AddEntry(h1i,"Incoherent J/#psi");
    l->AddEntry(hun,"Incoherent J/#psi with nucleon dissociation");
    l->AddEntry(hfc,"Coherent J/#psi from #psi(2S) decay");
    l->AddEntry(hfi,"Incoherent J/#psi from #psi(2S) decay");
    l->AddEntry(hgl,"Continuum #gamma#gamma #rightarrow #mu#mu");
    l->AddEntry(fsum,"Sum");
    l->Draw();
  }

//  TLatex* latex = new TLatex();
//  latex->SetNDC();
//  latex->DrawLatex(0.1,0.2,);
  gPad->Print(Form("fitPt_%s_%i_%i.png",suffix.Data(),iy,log));
  gPad->Print(Form("pt%i.pdf",iy));
}
