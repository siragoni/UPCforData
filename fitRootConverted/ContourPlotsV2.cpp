#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"


//_____________________________________________________________________________
/* -
 * -
 */
void ContourPlots(){

  TFile* fileHe = new TFile("pngResults/ContourFileHeV2.root");
  TFile* fileCs = new TFile("pngResults/ContourFileCsV2.root");

  TGraph *LThetaVsLPthreeHe  = (TGraph*)fileHe->Get("LThetaVsLPthree");
  TGraph *LThetaVsLPthreeCs  = (TGraph*)fileCs->Get("LThetaVsLPthree");
  TGraph *LThetaVsLPtwoHe    = (TGraph*)fileHe->Get("LThetaVsLPtwo");
  TGraph *LThetaVsLPtwoCs    = (TGraph*)fileCs->Get("LThetaVsLPtwo");
  TGraph *LThetaVsLPoneHe    = (TGraph*)fileHe->Get("LThetaVsLPone");
  TGraph *LThetaVsLPoneCs    = (TGraph*)fileCs->Get("LThetaVsLPone");
  TGraph *LThetaVsLTPthreeHe = (TGraph*)fileHe->Get("LThetaVsLTPthree");
  TGraph *LThetaVsLTPthreeCs = (TGraph*)fileCs->Get("LThetaVsLTPthree");
  TGraph *LThetaVsLTPtwoHe   = (TGraph*)fileHe->Get("LThetaVsLTPtwo");
  TGraph *LThetaVsLTPtwoCs   = (TGraph*)fileCs->Get("LThetaVsLTPtwo");
  TGraph *LThetaVsLTPoneHe   = (TGraph*)fileHe->Get("LThetaVsLTPone");
  TGraph *LThetaVsLTPoneCs   = (TGraph*)fileCs->Get("LThetaVsLTPone");
  TGraph *LPhiVsLTPthreeHe   = (TGraph*)fileHe->Get("LPhiVsLTPthree");
  TGraph *LPhiVsLTPthreeCs   = (TGraph*)fileCs->Get("LPhiVsLTPthree");
  TGraph *LPhiVsLTPtwoHe     = (TGraph*)fileHe->Get("LPhiVsLTPtwo");
  TGraph *LPhiVsLTPtwoCs     = (TGraph*)fileCs->Get("LPhiVsLTPtwo");
  TGraph *LPhiVsLTPoneHe     = (TGraph*)fileHe->Get("LPhiVsLTPone");
  TGraph *LPhiVsLTPoneCs     = (TGraph*)fileCs->Get("LPhiVsLTPone");






  // fThetaMC-> SetLineColor(kRed);
  // fPhiMC  -> SetLineColor(kRed);
  // fTheta  -> SetLineColor(kBlue);
  // fPhi    -> SetLineColor(kBlue);
  // fThetaMC-> SetLineWidth(3);
  // fPhiMC  -> SetLineWidth(3);
  // fTheta  -> SetLineWidth(3);
  // fPhi    -> SetLineWidth(3);











  TCanvas *c3 = new TCanvas("c3","multipads",900,700);
  gStyle->SetOptStat(0);
  c3->Divide(2,3,0,0);
  TH2F h1("h1","test1",10,0,1,20,0,20);
  TH2F h2("h2","test2",10,0,1,20,0,100);
  TH2F h3("h3","test3",10,0,1,20,-1,1);
  TH2F h4("h4","test4",10,0,1,20,0,1000);

  c3->cd(1);
  gPad->SetTickx(2);
  gPad->SetLeftMargin(0.2);
  // gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.12);
  // fThetaMC->GetYaxis()->SetRangeUser(0, 2.*fPhiMC->GetMaximum());
  // fThetaMC->GetYaxis()->SetRangeUser(0, 0.0195);
  // fThetaMC->Draw("L HIST");
  // fTheta  ->Draw("same");

  TLatex* latex3 = new TLatex();
  latex3->SetTextSize(0.045);
  latex3->SetTextFont(42);
  latex3->SetTextAlign(11);
  latex3->SetNDC();
  latex3->DrawLatex(0.25,0.94,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3->SetTextSize(0.036);
  // latex3->DrawLatex(0.55,0.84,"UPC, #it{L} = 235 ub^{-1}");
  // latex3->DrawLatex(0.55,0.84,"UPC, single muons");
  // latex3->DrawLatex(0.55,0.78,"2.9 < M_{#mu#mu} < 3.2 GeV/#it{c^{2}}");
  // latex3->DrawLatex(0.55,0.72,"#it{p}_{T}^{#mu#mu} < 0.2 GeV/#it{c}");
  // latex3->DrawLatex(0.55,0.66,Form("%.1f < y < %.1f",-4.0,-2.5));
  latex3->DrawLatex(0.25,0.84,"UPC, single muons");
  latex3->DrawLatex(0.25,0.78,"2.9 < M_{#mu#mu} < 3.2 GeV/#it{c^{2}}");
  latex3->DrawLatex(0.25,0.72,"#it{p}_{T}^{#mu#mu} < 0.2 GeV/#it{c}");
  latex3->DrawLatex(0.25,0.66,Form("%.1f < y < %.1f",-4.0,-2.5));

  // TLegend* l3 = new TLegend(0.27,0.55,0.5,0.75);
  // // TLegend* l3 = new TLegend(0.45,0.55,0.98,0.85);
  // l3->SetMargin(0.1);
  // l3->SetBorderSize(0);
  // l3->AddEntry(  fThetaMC, "ALICE Simulation");
  // l3->AddEntry(  fTheta  , "ALICE Data"      );
  // l3->Draw();





  c3->cd(2);
  gPad->SetTickx(2);
  // gPad->SetTicky(2);
  // gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  // gPad->SetBottomMargin(0.12);
  // fPhiMC->GetYaxis()->SetRangeUser(0, 0.0195);
  // fPhiMC->GetXaxis()->SetRangeUser(-0.05, 3.14*2.);
  // fPhiMC->GetYaxis()->SetLabelOffset(0.01);
  // // fPhiMC->GetXaxis()->SetTitleSize(42);
  // fPhiMC->SetTitleSize(1.1*fThetaMC->GetXaxis()->GetTitleSize(), "x");
  // fPhiMC->GetXaxis()->SetTitleOffset(0.95);
  // fPhiMC->Draw("L HIST");
  // fPhi  ->Draw("same");



  // TLegend* l4 = new TLegend(0.27,0.55,0.75,0.75);
  TLegend* l4 = new TLegend(0.2,0.15,0.75,0.28);
  // TLegend* l4 = new TLegend(0.15,0.85,0.45,0.75);
  // TLegend* l4 = new TLegend(0.45,0.55,0.98,0.85);
  gStyle->SetLegendFont(42);
  l4->SetMargin(0.1);
  l4->SetBorderSize(0);
  // l4->AddEntry(  fPhiMC, "ALICE Simulation");
  // l4->AddEntry(  fPhi  , "ALICE Data"      );
  l4->Draw();



  //
  //
  //
  // c3->cd(3);
  // h3.Draw();
  //
  // c3->cd(4);
  // gPad->SetTicky(2);
  // h4.Draw();
  //


  // c3->SaveAs("pngResults/ContourPlotsPaper.png",  "RECREATE");













  TH2F *histo=new TH2F("histo","",10, 0, 70, 10, -0.22, 0.22);
  histo->GetXaxis()->SetTitle("Test");
  histo->GetYaxis()->SetTitle("Test");


  // TCanvas *c_super = new TCanvas("c_super","comparing modes",100,100,900,700);
  // c_super->cd();
  // c_super->Divide(1, 2);
  // //c_super->Divide(1,3,1e-17,1e-17);
  // c_super->cd(1);
  // TPad *smallpad = new TPad("smallpad","",0.0,0.0,1,1);
  // smallpad->SetFillStyle(4000);
  // smallpad->Draw();
  // smallpad->Divide(3,1,0,0);
  //
  // for(int m=1; m< 4; m++){
  //
  //     smallpad->cd(m);
  //     if( m == 3)   gPad->SetRightMargin(0.05);
  //     histo->Draw();
  // }
  //
  //
  //
  // c_super->cd(2);
  // TPad *smallpad2 = new TPad("smallpad2","",0.0,0.0,1,1);
  // smallpad2->SetFillStyle(4000);
  // smallpad2->Draw();
  // smallpad2->Divide(3,1,0,0);
  //
  // for(int m=1; m< 4; m++){
  //
  //     smallpad2->cd(m);
  //     if( m == 3)   gPad->SetRightMargin(0.05);
  //     histo->Draw();
  // }
  // TCanvas *c_super = new TCanvas("c_super","comparing modes",100,100,1500,700);
  TCanvas *c_super = new TCanvas("c_super","comparing modes",100,100,1200,700);
  c_super->cd();
  c_super->Divide(3, 1, 0,0 );
  // c_super->Divide(3,1,1e-17,1e-17);
  c_super->cd(1);
  TPad *smallpad = new TPad("smallpad","",0.0,0.0,1,1);
  // smallpad->SetTopMargin(0.3);
  smallpad->SetFillStyle(4000);
  smallpad->Draw();
  smallpad->Divide(1,2,0,0);

  for(int m=1; m< 3; m++){

      smallpad->cd(m);
      gPad->SetRightMargin(0.05);
      histo->Draw();
  }



  c_super->cd(2);
  TPad *smallpad2 = new TPad("smallpad2","",0.0,0.0,1,1);
  // smallpad2->SetTopMargin(0.3);
  smallpad2->SetFillStyle(4000);
  smallpad2->Draw();
  smallpad2->Divide(1,2,0,0);

  for(int m=1; m< 3; m++){

      smallpad2->cd(m);
      gPad->SetRightMargin(0.05);
      histo->Draw();
  }



  c_super->cd(3);
  TPad *smallpad3 = new TPad("smallpad3","",0.0,0.0,1,1);
  // smallpad3->SetTopMargin(0.3);
  smallpad3->SetFillStyle(4000);
  smallpad3->Draw();
  smallpad3->Divide(1,2,0,0);

  for(int m=1; m< 3; m++){

      smallpad3->cd(m);
      gPad->SetRightMargin(0.05);
      histo->Draw();
  }




  c_super ->cd(1);
  smallpad->cd(1);
  // gPad->SetLeftMargin(0.2);
  gPad->SetLeftMargin(0.2);
  // gPad->SetLeftMargin(0.12);
  // gPad->SetRightMargin(0.1);
  gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.1);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  // TLine *line = new TLine(0.4947579,0.002894678,2.484597,0.0003370599);
  // line->SetLineWidth(2);
  // line->SetLineColor(kGray);
  // line->Draw();
  LThetaVsLTPtwoHe->SetTitle("");
  LThetaVsLTPtwoHe->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLTPtwoHe->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LThetaVsLTPtwoHe->GetXaxis()->CenterTitle(true);
  LThetaVsLTPtwoHe->GetYaxis()->CenterTitle(true);
  // LThetaVsLTPtwoHe->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLTPtwoHe->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLTPtwoHe->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLTPtwoHe->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLTPtwoHe->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLTPtwoHe->GetXaxis()->SetTitleSize(0.07);
  TAxis *axis1 = LThetaVsLTPtwoHe->GetXaxis();
  // axis1->SetLimits(0.6,2.0);                 // along X
  // axis1->SetLimits(-1.,2.0);                 // along X
  // axis1->SetLimits(0.5,2.5);                 // along X // OFFICIAL
  axis1->SetLimits(0.5,1.5);                 // along X
  // LThetaVsLTPtwoHe->GetHistogram()->SetMaximum( 0.19);   // along
  // LThetaVsLTPtwoHe->GetHistogram()->SetMinimum(-0.19);  //   Y
  // LThetaVsLTPtwoHe->GetHistogram()->SetMaximum( 0.95);   // along
  // LThetaVsLTPtwoHe->GetHistogram()->SetMinimum(-0.95);  //   Y
  LThetaVsLTPtwoHe->GetHistogram()->SetMaximum( 0.35);   // along
  LThetaVsLTPtwoHe->GetHistogram()->SetMinimum(-0.35);  //   Y
  LThetaVsLTPtwoHe->SetLineWidth(3);
  LThetaVsLTPtwoHe->SetLineColor(4);
  LThetaVsLTPtwoHe->Draw("AL");
  // TLatex* latex5 = new TLatex();
  // latex5->SetTextSize(0.045);
  // latex5->SetTextFont(24);
  // latex5->SetTextAlign(11);
  // latex5->SetNDC();
  // latex5->DrawLatex(0.25,0.94,"ALICE Public");
  // latex5->DrawLatex(0.25,0.86,"PbPb #sqrt{s_{NN}} = 5.02 TeV");
  auto texA = new TLatex(0.7665754,0.1546582,"(a)");
  texA->SetTextSize(0.08);
  texA->SetLineWidth(2);
  texA->Draw();


  // superimpose the second graph by leaving out the axis option "A"
  LThetaVsLTPoneHe->SetTitle("");
  LThetaVsLTPoneHe->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLTPoneHe->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LThetaVsLTPoneHe->GetXaxis()->CenterTitle(true);
  LThetaVsLTPoneHe->GetYaxis()->CenterTitle(true);
  // LThetaVsLTPoneHe->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLTPoneHe->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLTPoneHe->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLTPoneHe->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLTPoneHe->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLTPoneHe->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLTPoneHe->SetLineWidth(3);
  // LThetaVsLTPoneHe->SetMarkerStyle(21);
  LThetaVsLTPoneHe->SetLineColor(2);
  LThetaVsLTPoneHe->Draw("Lsame");

  LThetaVsLTPthreeHe->SetTitle("");
  LThetaVsLTPthreeHe->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLTPthreeHe->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LThetaVsLTPthreeHe->GetXaxis()->CenterTitle(true);
  LThetaVsLTPthreeHe->GetYaxis()->CenterTitle(true);
  // LThetaVsLTPthreeHe->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLTPthreeHe->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLTPthreeHe->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLTPthreeHe->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLTPthreeHe->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLTPthreeHe->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLTPthreeHe->SetLineWidth(3);
  // LThetaVsLTPthreeHe->SetMarkerStyle(21);
  LThetaVsLTPthreeHe->SetLineColor(kBlack);
  LThetaVsLTPthreeHe->Draw("Lsame");
  TPolyMarker *pmThetaVsLTPHe = new TPolyMarker(1);
  pmThetaVsLTPHe->SetPoint(0,0.956,-0.033);
  pmThetaVsLTPHe->SetMarkerStyle(20);
  pmThetaVsLTPHe->SetMarkerColor(kBlack);
  pmThetaVsLTPHe->SetMarkerSize(1.0);
  pmThetaVsLTPHe->Draw("same");
  // TLine *line = new TLine(0.5,0.,2.5,0.); // OFFICIAL
  TLine *line = new TLine(0.5,0.,1.5,0.);
  line->SetLineWidth(2);
  line->SetLineColor(kGray);
  line->Draw();
  pmThetaVsLTPHe->Draw("same");
  LThetaVsLTPthreeHe->Draw("Lsame");
  LThetaVsLTPtwoHe->Draw("Lsame");
  LThetaVsLTPoneHe->Draw("Lsame");
  // TLine *line = new TLine(0.4947579,0.002894678,2.484597,0.0003370599);
  // line->SetLineWidth(2);
  // line->SetLineColor(kGray);
  // line->Draw();





  smallpad->cd(2);
  // gPad->SetLeftMargin(0.2);
  gPad->SetLeftMargin(0.2);
  // gPad->SetLeftMargin(0.12);
  // gPad->SetRightMargin(0.1);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.13);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  LThetaVsLTPtwoCs->SetTitle("");
  LThetaVsLTPtwoCs->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLTPtwoCs->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LThetaVsLTPtwoCs->GetXaxis()->CenterTitle(true);
  LThetaVsLTPtwoCs->GetYaxis()->CenterTitle(true);
  // LThetaVsLTPtwoCs->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLTPtwoCs->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLTPtwoCs->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLTPtwoCs->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLTPtwoCs->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLTPtwoCs->GetXaxis()->SetTitleSize(0.07);
  TAxis *axis2 = LThetaVsLTPtwoCs->GetXaxis();
  // axis2->SetLimits(0.6,2.0);                 // along X
  // axis2->SetLimits(-1.,2.0);                 // along X
  // axis2->SetLimits(0.5,2.5);                 // along X // OFFICIAL
  axis2->SetLimits(0.5,1.5);                 // along X
  // LThetaVsLTPtwoCs->GetHistogram()->SetMaximum( 0.19);   // along
  // LThetaVsLTPtwoCs->GetHistogram()->SetMinimum(-0.19);  //   Y
  // LThetaVsLTPtwoCs->GetHistogram()->SetMaximum( 0.95);   // along
  // LThetaVsLTPtwoCs->GetHistogram()->SetMinimum(-0.95);  //   Y
  LThetaVsLTPtwoCs->GetHistogram()->SetMaximum( 0.35);   // along
  LThetaVsLTPtwoCs->GetHistogram()->SetMinimum(-0.35);  //   Y
  LThetaVsLTPtwoCs->SetLineWidth(3);
  LThetaVsLTPtwoCs->SetLineColor(4);
  LThetaVsLTPtwoCs->Draw("AL");

  // superimpose the second graph by leaving out the axis option "A"
  LThetaVsLTPoneCs->SetTitle("");
  LThetaVsLTPoneCs->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLTPoneCs->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LThetaVsLTPoneCs->GetXaxis()->CenterTitle(true);
  LThetaVsLTPoneCs->GetYaxis()->CenterTitle(true);
  // LThetaVsLTPoneCs->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLTPoneCs->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLTPoneCs->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLTPoneCs->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLTPoneCs->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLTPoneCs->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLTPoneCs->SetLineWidth(3);
  // LThetaVsLTPoneCs->SetMarkerStyle(21);
  LThetaVsLTPoneCs->SetLineColor(2);
  LThetaVsLTPoneCs->Draw("Lsame");


  LThetaVsLTPthreeCs->SetTitle("");
  LThetaVsLTPthreeCs->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLTPthreeCs->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LThetaVsLTPthreeCs->GetXaxis()->CenterTitle(true);
  LThetaVsLTPthreeCs->GetYaxis()->CenterTitle(true);
  LThetaVsLTPthreeCs->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLTPthreeCs->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLTPthreeCs->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLTPthreeCs->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLTPthreeCs->SetLineWidth(3);
  // LThetaVsLTPthreeCs->SetMarkerStyle(21);
  LThetaVsLTPthreeCs->SetLineColor(kBlack);
  LThetaVsLTPthreeCs->Draw("Lsame");
  TPolyMarker *pmThetaVsLTPCs = new TPolyMarker(1);
  pmThetaVsLTPCs->SetPoint(0,0.965,-0.028);
  pmThetaVsLTPCs->SetMarkerStyle(20);
  pmThetaVsLTPCs->SetMarkerColor(kBlack);
  pmThetaVsLTPCs->SetMarkerSize(1.0);
  pmThetaVsLTPCs->Draw("same");
  auto texD = new TLatex(0.7665754,0.1546582,"(d)");
  texD->SetTextSize(0.08);
  texD->SetLineWidth(2);
  texD->Draw();
  // TLine *line2 = new TLine(0.5,0.,2.5,0.); // OFFICIAL
  TLine *line2 = new TLine(0.5,0.,1.5,0.);
  line2->SetLineWidth(2);
  line2->SetLineColor(kGray);
  line2->Draw();
  pmThetaVsLTPCs->Draw("same");
  LThetaVsLTPthreeCs->Draw("Lsame");
  LThetaVsLTPtwoCs->Draw("Lsame");
  LThetaVsLTPoneCs->Draw("Lsame");
























  c_super ->cd(2);
  smallpad2->cd(1);
  // gPad->SetLeftMargin(0.2);
  // gPad->SetLeftMargin(0.15);
  gPad->SetLeftMargin(0.12);
  // gPad->SetRightMargin(0.1);
  gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.13);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  LPhiVsLTPtwoHe->SetTitle("");
  LPhiVsLTPtwoHe->GetXaxis()->SetTitle("#lambda_{#varphi}");
  LPhiVsLTPtwoHe->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LPhiVsLTPtwoHe->GetXaxis()->CenterTitle(true);
  LPhiVsLTPtwoHe->GetYaxis()->CenterTitle(true);
  // LPhiVsLTPtwoHe->GetYaxis()->SetTitleOffset(2.2);
  // LPhiVsLTPtwoHe->GetXaxis()->SetTitleOffset(1.1);
  LPhiVsLTPtwoHe->GetYaxis()->SetTitleOffset(0.8);
  LPhiVsLTPtwoHe->GetXaxis()->SetTitleOffset(0.8);
  LPhiVsLTPtwoHe->GetYaxis()->SetTitleSize(0.07);
  LPhiVsLTPtwoHe->GetXaxis()->SetTitleSize(0.07);
  TAxis *axis3 = LPhiVsLTPtwoHe->GetXaxis();
  // axis3->SetLimits(-0.09,0.19);                    // along X
  // axis3->SetLimits(-1.,2.);                    // along X
  axis3->SetLimits(-0.11,0.19);                    // along X
  // LPhiVsLTPtwoHe->GetHistogram()->SetMaximum( 0.13);   // along
  // LPhiVsLTPtwoHe->GetHistogram()->SetMinimum(-0.17);  //   Y
  // LPhiVsLTPtwoHe->GetHistogram()->SetMaximum( 0.95);   // along
  // LPhiVsLTPtwoHe->GetHistogram()->SetMinimum(-0.95);  //   Y
  LPhiVsLTPtwoHe->GetHistogram()->SetMaximum( 0.35);   // along
  LPhiVsLTPtwoHe->GetHistogram()->SetMinimum(-0.35);  //   Y
  LPhiVsLTPtwoHe->SetLineWidth(3);
  LPhiVsLTPtwoHe->SetLineColor(4);
  LPhiVsLTPtwoHe->Draw("AL");

  // superimpose the second graph by leaving out the axis option "A"
  LPhiVsLTPoneHe->SetTitle("");
  LPhiVsLTPoneHe->GetXaxis()->SetTitle("#lambda_{#varphi}");
  LPhiVsLTPoneHe->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LPhiVsLTPoneHe->GetXaxis()->CenterTitle(true);
  LPhiVsLTPoneHe->GetYaxis()->CenterTitle(true);
  // LPhiVsLTPoneHe->GetYaxis()->SetTitleOffset(2.2);
  // LPhiVsLTPoneHe->GetXaxis()->SetTitleOffset(1.1);
  LPhiVsLTPoneHe->GetYaxis()->SetTitleOffset(0.8);
  LPhiVsLTPoneHe->GetXaxis()->SetTitleOffset(0.8);
  LPhiVsLTPoneHe->GetYaxis()->SetTitleSize(0.07);
  LPhiVsLTPoneHe->GetXaxis()->SetTitleSize(0.07);
  LPhiVsLTPoneHe->SetLineWidth(3);
  // LPhiVsLTPoneHe->SetMarkerStyle(21);
  LPhiVsLTPoneHe->SetLineColor(2);
  LPhiVsLTPoneHe->Draw("Lsame");

  LPhiVsLTPthreeHe->SetTitle("");
  LPhiVsLTPthreeHe->GetXaxis()->SetTitle("#lambda_{#varphi}");
  LPhiVsLTPthreeHe->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LPhiVsLTPthreeHe->GetXaxis()->CenterTitle(true);
  LPhiVsLTPthreeHe->GetYaxis()->CenterTitle(true);
  // LPhiVsLTPthreeHe->GetYaxis()->SetTitleOffset(2.2);
  // LPhiVsLTPthreeHe->GetXaxis()->SetTitleOffset(1.1);
  LPhiVsLTPthreeHe->GetYaxis()->SetTitleOffset(0.8);
  LPhiVsLTPthreeHe->GetXaxis()->SetTitleOffset(0.8);
  LPhiVsLTPthreeHe->GetYaxis()->SetTitleSize(0.07);
  LPhiVsLTPthreeHe->GetXaxis()->SetTitleSize(0.07);
  LPhiVsLTPthreeHe->SetLineWidth(3);
  // LPhiVsLTPthreeHe->SetMarkerStyle(21);
  LPhiVsLTPthreeHe->SetLineColor(kBlack);
  LPhiVsLTPthreeHe->Draw("Lsame");
  TPolyMarker *pmPhiVsLTPthreeHe = new TPolyMarker(1);
  pmPhiVsLTPthreeHe->SetPoint(0,0.040,-0.033);
  pmPhiVsLTPthreeHe->SetMarkerStyle(20);
  pmPhiVsLTPthreeHe->SetMarkerColor(kBlack);
  pmPhiVsLTPthreeHe->SetMarkerSize(1.0);
  pmPhiVsLTPthreeHe->Draw("same");
  auto texB = new TLatex(-0.035,0.1546582,"(b)");
  texB->SetTextSize(0.08);
  texB->SetLineWidth(2);
  texB->Draw();
  TLine *line3 = new TLine(-0.11,0.,0.19,0.);
  line3->SetLineWidth(2);
  line3->SetLineColor(kGray);
  line3->Draw();
  pmPhiVsLTPthreeHe->Draw("same");
  LPhiVsLTPthreeHe->Draw("Lsame");
  LPhiVsLTPtwoHe->Draw("Lsame");
  LPhiVsLTPoneHe->Draw("Lsame");






  smallpad2->cd(2);
  // gPad->SetLeftMargin(0.2);
  // gPad->SetLeftMargin(0.15);
  gPad->SetLeftMargin(0.12);
  // gPad->SetRightMargin(0.1);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.13);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  LPhiVsLTPtwoCs->SetTitle("");
  LPhiVsLTPtwoCs->GetXaxis()->SetTitle("#lambda_{#varphi}");
  LPhiVsLTPtwoCs->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LPhiVsLTPtwoCs->GetXaxis()->CenterTitle(true);
  LPhiVsLTPtwoCs->GetYaxis()->CenterTitle(true);
  // LPhiVsLTPtwoCs->GetYaxis()->SetTitleOffset(2.2);
  // LPhiVsLTPtwoCs->GetXaxis()->SetTitleOffset(1.1);
  LPhiVsLTPtwoCs->GetYaxis()->SetTitleOffset(0.8);
  LPhiVsLTPtwoCs->GetXaxis()->SetTitleOffset(0.8);
  LPhiVsLTPtwoCs->GetYaxis()->SetTitleSize(0.07);
  LPhiVsLTPtwoCs->GetXaxis()->SetTitleSize(0.07);
  TAxis *axis4 = LPhiVsLTPtwoCs->GetXaxis();
  // axis4->SetLimits(-0.09,0.19);                 // along X
  axis4->SetLimits(-0.11,0.19);                 // along X
  // LPhiVsLTPtwoCs->GetHistogram()->SetMaximum( 0.13);   // along
  // LPhiVsLTPtwoCs->GetHistogram()->SetMinimum(-0.17);  //   Y
  // LPhiVsLTPtwoCs->GetHistogram()->SetMaximum( 0.95);   // along
  // LPhiVsLTPtwoCs->GetHistogram()->SetMinimum(-0.95);  //   Y
  LPhiVsLTPtwoCs->GetHistogram()->SetMaximum( 0.35);   // along
  LPhiVsLTPtwoCs->GetHistogram()->SetMinimum(-0.35);  //   Y
  LPhiVsLTPtwoCs->SetLineWidth(3);
  LPhiVsLTPtwoCs->SetLineColor(4);
  LPhiVsLTPtwoCs->Draw("AL");

  // superimpose the second graph by leaving out the axis option "A"
  LPhiVsLTPoneCs->SetTitle("");
  LPhiVsLTPoneCs->GetXaxis()->SetTitle("#lambda_{#varphi}");
  LPhiVsLTPoneCs->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LPhiVsLTPoneCs->GetXaxis()->CenterTitle(true);
  LPhiVsLTPoneCs->GetYaxis()->CenterTitle(true);
  // LPhiVsLTPoneCs->GetYaxis()->SetTitleOffset(2.2);
  // LPhiVsLTPoneCs->GetXaxis()->SetTitleOffset(1.1);
  LPhiVsLTPoneCs->GetYaxis()->SetTitleOffset(0.8);
  LPhiVsLTPoneCs->GetXaxis()->SetTitleOffset(0.8);
  LPhiVsLTPoneCs->GetYaxis()->SetTitleSize(0.07);
  LPhiVsLTPoneCs->GetXaxis()->SetTitleSize(0.07);
  LPhiVsLTPoneCs->SetLineWidth(3);
  // LPhiVsLTPoneCs->SetMarkerStyle(21);
  LPhiVsLTPoneCs->SetLineColor(2);
  LPhiVsLTPoneCs->Draw("Lsame");

  LPhiVsLTPthreeCs->SetTitle("");
  LPhiVsLTPthreeCs->GetXaxis()->SetTitle("#lambda_{#varphi}");
  LPhiVsLTPthreeCs->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  LPhiVsLTPthreeCs->GetXaxis()->CenterTitle(true);
  LPhiVsLTPthreeCs->GetYaxis()->CenterTitle(true);
  // LPhiVsLTPthreeCs->GetYaxis()->SetTitleOffset(2.2);
  // LPhiVsLTPthreeCs->GetXaxis()->SetTitleOffset(1.1);
  LPhiVsLTPthreeCs->GetYaxis()->SetTitleOffset(0.8);
  LPhiVsLTPthreeCs->GetXaxis()->SetTitleOffset(0.8);
  LPhiVsLTPthreeCs->GetYaxis()->SetTitleSize(0.07);
  LPhiVsLTPthreeCs->GetXaxis()->SetTitleSize(0.07);
  LPhiVsLTPthreeCs->SetLineWidth(3);
  // LPhiVsLTPthreeCs->SetMarkerStyle(21);
  LPhiVsLTPthreeCs->SetLineColor(kBlack);
  LPhiVsLTPthreeCs->Draw("Lsame");
  TPolyMarker *pmPhiVsLTPthreeCs = new TPolyMarker(1);
  pmPhiVsLTPthreeCs->SetPoint(0,0.038,-0.028);
  pmPhiVsLTPthreeCs->SetMarkerStyle(20);
  pmPhiVsLTPthreeCs->SetMarkerColor(kBlack);
  pmPhiVsLTPthreeCs->SetMarkerSize(1.0);
  pmPhiVsLTPthreeCs->Draw("same");
  auto texE = new TLatex(-0.035,0.1546582,"(e)");
  texE->SetTextSize(0.08);
  texE->SetLineWidth(2);
  texE->Draw();
  TLine *line4 = new TLine(-0.11,0.,0.19,0.);
  line4->SetLineWidth(2);
  line4->SetLineColor(kGray);
  line4->Draw();
  pmPhiVsLTPthreeCs->Draw("same");
  LPhiVsLTPthreeCs->Draw("Lsame");
  LPhiVsLTPtwoCs->Draw("Lsame");
  LPhiVsLTPoneCs->Draw("Lsame");














  c_super ->cd(3);
  smallpad3->cd(1);
  // gPad->SetLeftMargin(0.2);
  // gPad->SetLeftMargin(0.15);
  gPad->SetLeftMargin(0.12);
  // gPad->SetRightMargin(0.1);
  gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.13);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  LThetaVsLPtwoHe->SetTitle("");
  LThetaVsLPtwoHe->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLPtwoHe->GetYaxis()->SetTitle("#lambda_{#varphi}");
  LThetaVsLPtwoHe->GetXaxis()->CenterTitle(true);
  LThetaVsLPtwoHe->GetYaxis()->CenterTitle(true);
  // LThetaVsLPtwoHe->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLPtwoHe->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLPtwoHe->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLPtwoHe->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLPtwoHe->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLPtwoHe->GetXaxis()->SetTitleSize(0.07);
  TAxis *axis5 = LThetaVsLPtwoHe->GetXaxis();
  // axis5->SetLimits(0.6,2.0);                   // along X
  // axis5->SetLimits(0.5,2.5);                   // along X // OFFICIAL
  axis5->SetLimits(0.5,1.5);                   // along X
  // LThetaVsLPtwoHe->GetHistogram()->SetMaximum( 0.19);   // along
  // LThetaVsLPtwoHe->GetHistogram()->SetMinimum(-0.09);  //   Y
  // LThetaVsLPtwoHe->GetHistogram()->SetMaximum( 0.95);   // along
  // LThetaVsLPtwoHe->GetHistogram()->SetMinimum(-0.95);  //   Y
  LThetaVsLPtwoHe->GetHistogram()->SetMaximum( 0.35);   // along
  LThetaVsLPtwoHe->GetHistogram()->SetMinimum(-0.35);  //   Y
  LThetaVsLPtwoHe->SetLineWidth(3);
  LThetaVsLPtwoHe->SetLineColor(4);
  LThetaVsLPtwoHe->Draw("AL");

  // superimpose the second graph by leaving out the axis option "A"
  LThetaVsLPoneHe->SetTitle("");
  LThetaVsLPoneHe->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLPoneHe->GetYaxis()->SetTitle("#lambda_{#varphi}");
  LThetaVsLPoneHe->GetXaxis()->CenterTitle(true);
  LThetaVsLPoneHe->GetYaxis()->CenterTitle(true);
  // LThetaVsLPoneHe->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLPoneHe->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLPoneHe->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLPoneHe->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLPoneHe->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLPoneHe->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLPoneHe->SetLineWidth(3);
  // LThetaVsLPoneHe->SetMarkerStyle(21);
  LThetaVsLPoneHe->SetLineColor(2);
  LThetaVsLPoneHe->Draw("Lsame");

  LThetaVsLPthreeHe->SetTitle("");
  LThetaVsLPthreeHe->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLPthreeHe->GetYaxis()->SetTitle("#lambda_{#varphi}");
  LThetaVsLPthreeHe->GetXaxis()->CenterTitle(true);
  LThetaVsLPthreeHe->GetYaxis()->CenterTitle(true);
  // LThetaVsLPthreeHe->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLPthreeHe->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLPthreeHe->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLPthreeHe->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLPthreeHe->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLPthreeHe->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLPthreeHe->SetLineWidth(3);
  // LThetaVsLPthreeHe->SetMarkerStyle(21);
  LThetaVsLPthreeHe->SetLineColor(kBlack);
  LThetaVsLPthreeHe->Draw("Lsame");
  TPolyMarker *pmThetaVsLPthreeHe = new TPolyMarker(1);
  pmThetaVsLPthreeHe->SetPoint(0,0.956,0.040);
  pmThetaVsLPthreeHe->SetMarkerStyle(20);
  pmThetaVsLPthreeHe->SetMarkerColor(kBlack);
  pmThetaVsLPthreeHe->SetMarkerSize(1.0);
  pmThetaVsLPthreeHe->Draw("same");
  auto texC = new TLatex(0.7665754,0.1546582,"(c)");
  texC->SetTextSize(0.08);
  texC->SetLineWidth(2);
  texC->Draw();
  // TLine *line5 = new TLine(0.5,0.,2.5,0.); // OFFICIAL
  TLine *line5 = new TLine(0.5,0.,1.5,0.);
  line5->SetLineWidth(2);
  line5->SetLineColor(kGray);
  line5->Draw();
  pmThetaVsLPthreeHe->Draw("same");
  LThetaVsLPthreeHe->Draw("Lsame");
  LThetaVsLPtwoHe->Draw("Lsame");
  LThetaVsLPoneHe->Draw("Lsame");






  smallpad3->cd(2);
  // gPad->SetLeftMargin(0.2);
  // gPad->SetLeftMargin(0.15);
  gPad->SetLeftMargin(0.12);
  // gPad->SetRightMargin(0.1);
  gPad->SetRightMargin(0.01);
  gPad->SetBottomMargin(0.13);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  LThetaVsLPtwoCs->SetTitle("");
  LThetaVsLPtwoCs->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLPtwoCs->GetYaxis()->SetTitle("#lambda_{#varphi}");
  LThetaVsLPtwoCs->GetXaxis()->CenterTitle(true);
  LThetaVsLPtwoCs->GetYaxis()->CenterTitle(true);
  // LThetaVsLPtwoCs->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLPtwoCs->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLPtwoCs->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLPtwoCs->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLPtwoCs->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLPtwoCs->GetXaxis()->SetTitleSize(0.07);
  TAxis *axis6 = LThetaVsLPtwoCs->GetXaxis();
  // axis6->SetLimits(0.6,2.0);                 // along X
  // axis6->SetLimits(-1.,2.0);                 // along X
  // axis6->SetLimits(0.5,2.5);                 // along X // OFFICIAL
  axis6->SetLimits(0.5,1.5);                 // along X
  // LThetaVsLPtwoCs->GetHistogram()->SetMaximum( 0.19);   // along
  // LThetaVsLPtwoCs->GetHistogram()->SetMinimum(-0.09);  //   Y
  // LThetaVsLPtwoCs->GetHistogram()->SetMaximum( 0.95);   // along
  // LThetaVsLPtwoCs->GetHistogram()->SetMinimum(-0.95);  //   Y
  LThetaVsLPtwoCs->GetHistogram()->SetMaximum( 0.35);   // along
  LThetaVsLPtwoCs->GetHistogram()->SetMinimum(-0.35);  //   Y
  LThetaVsLPtwoCs->SetLineWidth(3);
  LThetaVsLPtwoCs->SetLineColor(4);
  LThetaVsLPtwoCs->Draw("AL");

  // superimpose the second graph by leaving out the axis option "A"
  LThetaVsLPoneCs->SetTitle("");
  LThetaVsLPoneCs->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLPoneCs->GetYaxis()->SetTitle("#lambda_{#varphi}");
  LThetaVsLPoneCs->GetXaxis()->CenterTitle(true);
  LThetaVsLPoneCs->GetYaxis()->CenterTitle(true);
  // LThetaVsLPoneCs->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLPoneCs->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLPoneCs->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLPoneCs->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLPoneCs->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLPoneCs->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLPoneCs->SetLineWidth(3);
  // LThetaVsLPoneCs->SetMarkerStyle(21);
  LThetaVsLPoneCs->SetLineColor(2);
  LThetaVsLPoneCs->Draw("Lsame");

  LThetaVsLPthreeCs->SetTitle("");
  LThetaVsLPthreeCs->GetXaxis()->SetTitle("#lambda_{#theta}");
  LThetaVsLPthreeCs->GetYaxis()->SetTitle("#lambda_{#varphi}");
  LThetaVsLPthreeCs->GetXaxis()->CenterTitle(true);
  LThetaVsLPthreeCs->GetYaxis()->CenterTitle(true);
  // LThetaVsLPthreeCs->GetYaxis()->SetTitleOffset(2.2);
  // LThetaVsLPthreeCs->GetXaxis()->SetTitleOffset(1.1);
  LThetaVsLPthreeCs->GetYaxis()->SetTitleOffset(0.8);
  LThetaVsLPthreeCs->GetXaxis()->SetTitleOffset(0.8);
  LThetaVsLPthreeCs->GetYaxis()->SetTitleSize(0.07);
  LThetaVsLPthreeCs->GetXaxis()->SetTitleSize(0.07);
  LThetaVsLPthreeCs->SetLineWidth(3);
  // LThetaVsLPthreeCs->SetMarkerStyle(21);
  LThetaVsLPthreeCs->SetLineColor(kBlack);
  LThetaVsLPthreeCs->Draw("Lsame");
  TPolyMarker *pmThetaVsLPthreeCs = new TPolyMarker(1);
  pmThetaVsLPthreeCs->SetPoint(0,0.965,0.038);
  pmThetaVsLPthreeCs->SetMarkerStyle(20);
  pmThetaVsLPthreeCs->SetMarkerColor(kBlack);
  pmThetaVsLPthreeCs->SetMarkerSize(1.0);
  pmThetaVsLPthreeCs->Draw("same");
  auto texF = new TLatex(0.7665754,0.1546582,"(f)");
  texF->SetTextSize(0.08);
  texF->SetLineWidth(2);
  texF->Draw();
  // TLine *line6 = new TLine(0.5,0.,2.5,0.); // OFFICIAL
  TLine *line6 = new TLine(0.5,0.,1.5,0.); // OFFICIAL
  line6->SetLineWidth(2);
  line6->SetLineColor(kGray);
  line6->Draw();
  pmThetaVsLPthreeCs->Draw("same");
  LThetaVsLPthreeCs->Draw("Lsame");
  LThetaVsLPtwoCs->Draw("Lsame");
  LThetaVsLPoneCs->Draw("Lsame");
  // TLegend *leg = new TLegend(0.5,0.55,0.95,0.79);
  // leg->SetFillStyle(0);
  // leg->SetBorderSize(0);
  // leg->SetTextSize(0.04);
  // // leg->AddEntry("CorrectedTildePhi","Helicity", "ep");
  // // leg->AddEntry("CorrectedTildePhiCS","Collins-Soper #times 1.6", "ep");
  // leg->AddEntry(LThetaVsLPthreeCs,"3 #sigma", "L");
  // leg->AddEntry(LThetaVsLPtwoCs,  "2 #sigma", "L");
  // leg->AddEntry(LThetaVsLPoneCs,  "1 #sigma", "L");
  // leg->Draw();
  // auto leg = new TLatex(1.8,-0.1,"3 #sigma"); // OFFICIAL
  auto leg = new TLatex(1.2,-0.1,"3 #sigma");
  leg->SetTextSize(0.08);
  leg->SetLineWidth(2);
  leg->Draw();
  // auto leg2 = new TLatex(1.8,-0.17,"#color[4]{2 #sigma}"); // OFFICIAL
  auto leg2 = new TLatex(1.2,-0.17,"#color[4]{2 #sigma}");
  leg2->SetTextSize(0.08);
  leg2->SetLineWidth(2);
  leg2->Draw();
  // auto leg3 = new TLatex(1.8,-0.24,"#color[2]{1 #sigma}"); // OFFICIAL
  auto leg3 = new TLatex(1.2,-0.24,"#color[2]{1 #sigma}");
  leg3->SetTextSize(0.08);
  leg3->SetLineWidth(2);
  leg3->Draw();


  c_super->cd();
  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.91,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");




  c_super->SaveAs("pngResults/ContourPlotsPaperV2.png",  "RECREATE");
  c_super->SaveAs("pngResults/ContourPlotsPaperV2.pdf",  "RECREATE");
  // c_super->SaveAs("pngResults/ContourPlotsPaper.C",  "RECREATE");


}
