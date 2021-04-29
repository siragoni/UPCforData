#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TStyle.h"
using namespace std;
#include <vector>


//_____________________________________________________________________________
/* -
 * -
 */
void Plot(){

  TDatime d;
  TFile* file1DHe  = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* file1DCs  = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedCs1Dv2.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );

  TH1F* CorrectedCosThetaHe = (TH1F*) file1DHe->Get("CorrCosThetaH");
  TH1F* CorrectedCosThetaCs = (TH1F*) file1DCs->Get("CorrCosThetaH");
  CorrectedCosThetaHe->Sumw2();
  CorrectedCosThetaCs->Sumw2();


  TH1F* CorrectedPhiHe = (TH1F*) file1DHe->Get("CorrPhiH");
  TH1F* CorrectedPhiCs = (TH1F*) file1DCs->Get("CorrPhiH");
  CorrectedPhiHe->Sumw2();
  CorrectedPhiCs->Sumw2();


  TH1F* CorrectedTildePhiHe = (TH1F*) file1DHe->Get("CorrTildePhiH");
  TH1F* CorrectedTildePhiCs = (TH1F*) file1DCs->Get("CorrTildePhiH");
  CorrectedTildePhiHe->Sumw2();
  CorrectedTildePhiCs->Sumw2();


  TH1F* acceptCosThetaHe = (TH1F*) file1DHe->Get("acceptanceCosTheta");
  TH1F* acceptCosThetaCs = (TH1F*) file1DCs->Get("acceptanceCosTheta");
  acceptCosThetaHe->Sumw2();
  acceptCosThetaCs->Sumw2();


  TH1F* acceptPhiHe = (TH1F*) file1DHe->Get("acceptancePhi");
  TH1F* acceptPhiCs = (TH1F*) file1DCs->Get("acceptancePhi");
  acceptPhiHe->Sumw2();
  acceptPhiCs->Sumw2();


  TH1F* acceptTildePhiHe = (TH1F*) file1DHe->Get("acceptanceTildePhi");
  TH1F* acceptTildePhiCs = (TH1F*) file1DCs->Get("acceptanceTildePhi");
  acceptTildePhiHe->Sumw2();
  acceptTildePhiCs->Sumw2();


  gStyle->SetOptStat(0);


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedCosThetaHe->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedCosThetaHe->GetYaxis()->SetTitleOffset(1.25);
  CorrectedCosThetaHe->GetYaxis()->SetTitleOffset(1.55);
  CorrectedCosThetaHe->GetXaxis()->SetTitleSize(0.045);
  CorrectedCosThetaHe->GetYaxis()->SetTitleSize(0.045);
  CorrectedCosThetaHe->GetXaxis()->SetLabelSize(0.045);
  CorrectedCosThetaHe->GetYaxis()->SetLabelSize(0.045);
  CorrectedCosThetaHe->GetXaxis()->SetTitleFont(42);
  CorrectedCosThetaHe->GetYaxis()->SetTitleFont(42);
  CorrectedCosThetaHe->GetXaxis()->SetLabelFont(42);
  CorrectedCosThetaHe->GetYaxis()->SetLabelFont(42);
  CorrectedCosThetaHe->GetXaxis()->SetNdivisions(408);
  CorrectedCosThetaHe->GetYaxis()->SetRangeUser(0., CorrectedCosThetaHe->GetMaximum()*2);
  // CorrectedCosThetaHe->GetXaxis()->SetRangeUser(2, 6);
  CorrectedCosThetaHe->SetTitle(  Form(  ";cos(#theta); ACCxEFF Corrected Counts / %.3f",
                           CorrectedCosThetaHe->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedCosThetaHe->Draw();
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  latex->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/CosThetaHe.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedCosThetaCs->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedCosThetaCs->GetYaxis()->SetTitleOffset(1.25);
  CorrectedCosThetaCs->GetYaxis()->SetTitleOffset(1.55);
  CorrectedCosThetaCs->GetXaxis()->SetTitleSize(0.045);
  CorrectedCosThetaCs->GetYaxis()->SetTitleSize(0.045);
  CorrectedCosThetaCs->GetXaxis()->SetLabelSize(0.045);
  CorrectedCosThetaCs->GetYaxis()->SetLabelSize(0.045);
  CorrectedCosThetaCs->GetXaxis()->SetTitleFont(42);
  CorrectedCosThetaCs->GetYaxis()->SetTitleFont(42);
  CorrectedCosThetaCs->GetXaxis()->SetLabelFont(42);
  CorrectedCosThetaCs->GetYaxis()->SetLabelFont(42);
  CorrectedCosThetaCs->GetXaxis()->SetNdivisions(408);
  CorrectedCosThetaCs->GetYaxis()->SetRangeUser(0., CorrectedCosThetaCs->GetMaximum()*2);
  // CorrectedCosThetaCs->GetXaxis()->SetRangeUser(2, 6);
  CorrectedCosThetaCs->SetTitle(  Form(  ";cos(#theta); ACCxEFF Corrected Counts / %.3f",
                           CorrectedCosThetaCs->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedCosThetaCs->Draw();
  TLatex* latexa = new TLatex();
  latexa->SetTextSize(0.05);
  latexa->SetTextFont(42);
  latexa->SetTextAlign(11);
  latexa->SetNDC();
  latexa->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latexa->SetTextSize(0.045);
  latexa->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/CosThetaCs.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedPhiHe->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedPhiHe->GetYaxis()->SetTitleOffset(1.25);
  CorrectedPhiHe->GetYaxis()->SetTitleOffset(1.55);
  CorrectedPhiHe->GetXaxis()->SetTitleSize(0.045);
  CorrectedPhiHe->GetYaxis()->SetTitleSize(0.045);
  CorrectedPhiHe->GetXaxis()->SetLabelSize(0.045);
  CorrectedPhiHe->GetYaxis()->SetLabelSize(0.045);
  CorrectedPhiHe->GetXaxis()->SetTitleFont(42);
  CorrectedPhiHe->GetYaxis()->SetTitleFont(42);
  CorrectedPhiHe->GetXaxis()->SetLabelFont(42);
  CorrectedPhiHe->GetYaxis()->SetLabelFont(42);
  CorrectedPhiHe->GetXaxis()->SetNdivisions(408);
  CorrectedPhiHe->GetYaxis()->SetRangeUser(0., CorrectedPhiHe->GetMaximum()*2);
  // CorrectedPhiHe->GetXaxis()->SetRangeUser(2, 6);
  CorrectedPhiHe->SetTitle(  Form(  ";#varphi; ACCxEFF Corrected Counts / %.3f",
                           CorrectedPhiHe->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedPhiHe->Draw();
  TLatex* latex2a = new TLatex();
  latex2a->SetTextSize(0.05);
  latex2a->SetTextFont(42);
  latex2a->SetTextAlign(11);
  latex2a->SetNDC();
  latex2a->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2a->SetTextSize(0.045);
  latex2a->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/PhiHe.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedPhiCs->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedPhiCs->GetYaxis()->SetTitleOffset(1.25);
  CorrectedPhiCs->GetYaxis()->SetTitleOffset(1.55);
  CorrectedPhiCs->GetXaxis()->SetTitleSize(0.045);
  CorrectedPhiCs->GetYaxis()->SetTitleSize(0.045);
  CorrectedPhiCs->GetXaxis()->SetLabelSize(0.045);
  CorrectedPhiCs->GetYaxis()->SetLabelSize(0.045);
  CorrectedPhiCs->GetXaxis()->SetTitleFont(42);
  CorrectedPhiCs->GetYaxis()->SetTitleFont(42);
  CorrectedPhiCs->GetXaxis()->SetLabelFont(42);
  CorrectedPhiCs->GetYaxis()->SetLabelFont(42);
  CorrectedPhiCs->GetXaxis()->SetNdivisions(408);
  CorrectedPhiCs->GetYaxis()->SetRangeUser(0., CorrectedPhiCs->GetMaximum()*2);
  // CorrectedPhiCs->GetXaxis()->SetRangeUser(2, 6);
  CorrectedPhiCs->SetTitle(  Form(  ";#varphi; ACCxEFF Corrected Counts / %.3f",
                           CorrectedPhiCs->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedPhiCs->Draw();
  TLatex* latex2 = new TLatex();
  latex2->SetTextSize(0.05);
  latex2->SetTextFont(42);
  latex2->SetTextAlign(11);
  latex2->SetNDC();
  latex2->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2->SetTextSize(0.045);
  latex2->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/PhiCs.png", "recreate");



  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedTildePhiHe->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedTildePhiHe->GetYaxis()->SetTitleOffset(1.25);
  CorrectedTildePhiHe->GetYaxis()->SetTitleOffset(1.55);
  CorrectedTildePhiHe->GetXaxis()->SetTitleSize(0.045);
  CorrectedTildePhiHe->GetYaxis()->SetTitleSize(0.045);
  CorrectedTildePhiHe->GetXaxis()->SetLabelSize(0.045);
  CorrectedTildePhiHe->GetYaxis()->SetLabelSize(0.045);
  CorrectedTildePhiHe->GetXaxis()->SetTitleFont(42);
  CorrectedTildePhiHe->GetYaxis()->SetTitleFont(42);
  CorrectedTildePhiHe->GetXaxis()->SetLabelFont(42);
  CorrectedTildePhiHe->GetYaxis()->SetLabelFont(42);
  CorrectedTildePhiHe->GetXaxis()->SetNdivisions(408);
  CorrectedTildePhiHe->GetYaxis()->SetRangeUser(0., CorrectedPhiHe->GetMaximum()*2);
  // CorrectedTildePhiHe->GetXaxis()->SetRangeUser(2, 6);
  CorrectedTildePhiHe->SetTitle(  Form(  ";#tilde{#varphi}; ACCxEFF Corrected Counts / %.3f",
                           CorrectedTildePhiHe->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedTildePhiHe->Draw();
  TLatex* latex3a = new TLatex();
  latex3a->SetTextSize(0.05);
  latex3a->SetTextFont(42);
  latex3a->SetTextAlign(11);
  latex3a->SetNDC();
  latex3a->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3a->SetTextSize(0.045);
  latex3a->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/TildePhiHe.png", "recreate");




  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedTildePhiCs->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedTildePhiCs->GetYaxis()->SetTitleOffset(1.25);
  CorrectedTildePhiCs->GetYaxis()->SetTitleOffset(1.55);
  CorrectedTildePhiCs->GetXaxis()->SetTitleSize(0.045);
  CorrectedTildePhiCs->GetYaxis()->SetTitleSize(0.045);
  CorrectedTildePhiCs->GetXaxis()->SetLabelSize(0.045);
  CorrectedTildePhiCs->GetYaxis()->SetLabelSize(0.045);
  CorrectedTildePhiCs->GetXaxis()->SetTitleFont(42);
  CorrectedTildePhiCs->GetYaxis()->SetTitleFont(42);
  CorrectedTildePhiCs->GetXaxis()->SetLabelFont(42);
  CorrectedTildePhiCs->GetYaxis()->SetLabelFont(42);
  CorrectedTildePhiCs->GetXaxis()->SetNdivisions(408);
  CorrectedTildePhiCs->GetYaxis()->SetRangeUser(0., CorrectedPhiHe->GetMaximum()*2);
  // CorrectedTildePhiCs->GetXaxis()->SetRangeUser(2, 6);
  CorrectedTildePhiCs->SetTitle(  Form(  ";#tilde{#varphi}; ACCxEFF Corrected Counts / %.3f",
                           CorrectedTildePhiCs->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedTildePhiCs->Draw();
  TLatex* latex3 = new TLatex();
  latex3->SetTextSize(0.05);
  latex3->SetTextFont(42);
  latex3->SetTextAlign(11);
  latex3->SetNDC();
  latex3->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3->SetTextSize(0.045);
  latex3->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/TildePhiCs.png", "recreate");










  TFile* fileDataRawCosThetaHe =  new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root",       d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* fileDataRawPhiHe      =  new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHEv2/PhiHeFrameV2.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* fileDataRawTildePhiHe =  new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHEv2/TildePhiHeFrameV2.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* fileDataRawCosThetaCs =  new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaCS/CosThetaCsFrame.root",       d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* fileDataRawPhiCs      =  new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiCSv2/PhiCsFrameV2.root",             d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* fileDataRawTildePhiCs =  new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiCSv2/TildePhiCsFrameV2.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TH1F* CosThetaAfterSignalExtractionErrorsRawHe = (TH1F*)fileDataRawCosThetaHe->Get("CosThetaAfterSignalExtractionErrorsH");
  TH1F* PhiAfterSignalExtractionErrorsRawHe      = (TH1F*)fileDataRawPhiHe     ->Get("PhiAfterSignalExtractionErrorsH");
  TH1F* TildePhiAfterSignalExtractionErrorsRawHe = (TH1F*)fileDataRawTildePhiHe->Get("TildePhiAfterSignalExtractionErrorsH");
  TH1F* CosThetaAfterSignalExtractionErrorsRawCs = (TH1F*)fileDataRawCosThetaCs->Get("CosThetaAfterSignalExtractionErrorsH");
  TH1F* PhiAfterSignalExtractionErrorsRawCs      = (TH1F*)fileDataRawPhiCs     ->Get("PhiAfterSignalExtractionErrorsH");
  TH1F* TildePhiAfterSignalExtractionErrorsRawCs = (TH1F*)fileDataRawTildePhiCs->Get("TildePhiAfterSignalExtractionErrorsH");




  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleOffset(1.15);
  // CosThetaAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleOffset(1.25);
  CosThetaAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleOffset(1.55);
  CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->SetLabelSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawHe->GetYaxis()->SetLabelSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleFont(42);
  CosThetaAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleFont(42);
  CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->SetLabelFont(42);
  CosThetaAfterSignalExtractionErrorsRawHe->GetYaxis()->SetLabelFont(42);
  CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->SetNdivisions(408);
  CosThetaAfterSignalExtractionErrorsRawHe->GetYaxis()->SetRangeUser(0., CosThetaAfterSignalExtractionErrorsRawHe->GetMaximum()*2);
  // CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->SetRangeUser(2, 6);
  CosThetaAfterSignalExtractionErrorsRawHe->SetTitle(  Form(  ";cos(#theta); Raw Counts / %.3f",
                           CosThetaAfterSignalExtractionErrorsRawHe->GetXaxis()->GetBinWidth(1)  )  );
  CosThetaAfterSignalExtractionErrorsRawHe->Draw();
  TLatex* latexb = new TLatex();
  latexb->SetTextSize(0.05);
  latexb->SetTextFont(42);
  latexb->SetTextAlign(11);
  latexb->SetNDC();
  latexb->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latexb->SetTextSize(0.045);
  latexb->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/RawCosThetaHe.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleOffset(1.15);
  // CosThetaAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleOffset(1.25);
  CosThetaAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleOffset(1.55);
  CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->SetLabelSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawCs->GetYaxis()->SetLabelSize(0.045);
  CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleFont(42);
  CosThetaAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleFont(42);
  CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->SetLabelFont(42);
  CosThetaAfterSignalExtractionErrorsRawCs->GetYaxis()->SetLabelFont(42);
  CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->SetNdivisions(408);
  CosThetaAfterSignalExtractionErrorsRawCs->GetYaxis()->SetRangeUser(0., CosThetaAfterSignalExtractionErrorsRawCs->GetMaximum()*2);
  // CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->SetRangeUser(2, 6);
  CosThetaAfterSignalExtractionErrorsRawCs->SetTitle(  Form(  ";cos(#theta); Raw Counts / %.3f",
                           CosThetaAfterSignalExtractionErrorsRawCs->GetXaxis()->GetBinWidth(1)  )  );
  CosThetaAfterSignalExtractionErrorsRawCs->Draw();
  TLatex* latexab = new TLatex();
  latexab->SetTextSize(0.05);
  latexab->SetTextFont(42);
  latexab->SetTextAlign(11);
  latexab->SetNDC();
  latexab->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latexab->SetTextSize(0.045);
  latexab->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/RawCosThetaCs.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleOffset(1.15);
  // PhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleOffset(1.25);
  PhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleOffset(1.55);
  PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleSize(0.045);
  PhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleSize(0.045);
  PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetLabelSize(0.045);
  PhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetLabelSize(0.045);
  PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleFont(42);
  PhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleFont(42);
  PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetLabelFont(42);
  PhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetLabelFont(42);
  PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetNdivisions(408);
  PhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetRangeUser(0., PhiAfterSignalExtractionErrorsRawHe->GetMaximum()*2);
  // PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetRangeUser(2, 6);
  PhiAfterSignalExtractionErrorsRawHe->SetTitle(  Form(  ";#varphi; Raw Counts / %.3f",
                           PhiAfterSignalExtractionErrorsRawHe->GetXaxis()->GetBinWidth(1)  )  );
  PhiAfterSignalExtractionErrorsRawHe->Draw();
  TLatex* latex2ab = new TLatex();
  latex2ab->SetTextSize(0.05);
  latex2ab->SetTextFont(42);
  latex2ab->SetTextAlign(11);
  latex2ab->SetNDC();
  latex2ab->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2ab->SetTextSize(0.045);
  latex2ab->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/RawPhiHe.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  PhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleOffset(1.15);
  // PhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleOffset(1.25);
  PhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleOffset(1.55);
  PhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleSize(0.045);
  PhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleSize(0.045);
  PhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetLabelSize(0.045);
  PhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetLabelSize(0.045);
  PhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleFont(42);
  PhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleFont(42);
  PhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetLabelFont(42);
  PhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetLabelFont(42);
  PhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetNdivisions(408);
  PhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetRangeUser(0., PhiAfterSignalExtractionErrorsRawCs->GetMaximum()*2);
  // CorrectedPhiCs->GetXaxis()->SetRangeUser(2, 6);
  PhiAfterSignalExtractionErrorsRawCs->SetTitle(  Form(  ";#varphi; Raw Counts / %.3f",
                           PhiAfterSignalExtractionErrorsRawCs->GetXaxis()->GetBinWidth(1)  )  );
  PhiAfterSignalExtractionErrorsRawCs->Draw();
  TLatex* latex2b = new TLatex();
  latex2b->SetTextSize(0.05);
  latex2b->SetTextFont(42);
  latex2b->SetTextAlign(11);
  latex2b->SetNDC();
  latex2b->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2b->SetTextSize(0.045);
  latex2b->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/RawPhiCs.png", "recreate");



  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleOffset(1.15);
  // TildePhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleOffset(1.25);
  TildePhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleOffset(1.55);
  TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetLabelSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetLabelSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetTitleFont(42);
  TildePhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetTitleFont(42);
  TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetLabelFont(42);
  TildePhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetLabelFont(42);
  TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetNdivisions(408);
  TildePhiAfterSignalExtractionErrorsRawHe->GetYaxis()->SetRangeUser(0., TildePhiAfterSignalExtractionErrorsRawHe->GetMaximum()*2);
  // TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->SetRangeUser(2, 6);
  TildePhiAfterSignalExtractionErrorsRawHe->SetTitle(  Form(  ";#tilde{#varphi}; Raw Counts / %.3f",
                           TildePhiAfterSignalExtractionErrorsRawHe->GetXaxis()->GetBinWidth(1)  )  );
  TildePhiAfterSignalExtractionErrorsRawHe->Draw();
  TLatex* latex3ab = new TLatex();
  latex3ab->SetTextSize(0.05);
  latex3ab->SetTextFont(42);
  latex3ab->SetTextAlign(11);
  latex3ab->SetNDC();
  latex3ab->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3ab->SetTextSize(0.045);
  latex3ab->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/RawTildePhiHe.png", "recreate");




  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleOffset(1.15);
  // TildePhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleOffset(1.25);
  TildePhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleOffset(1.55);
  TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetLabelSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetLabelSize(0.045);
  TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetTitleFont(42);
  TildePhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetTitleFont(42);
  TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetLabelFont(42);
  TildePhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetLabelFont(42);
  TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetNdivisions(408);
  TildePhiAfterSignalExtractionErrorsRawCs->GetYaxis()->SetRangeUser(0., TildePhiAfterSignalExtractionErrorsRawCs->GetMaximum()*2);
  // TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->SetRangeUser(2, 6);
  TildePhiAfterSignalExtractionErrorsRawCs->SetTitle(  Form(  ";#tilde{#varphi}; Raw Counts / %.3f",
                           TildePhiAfterSignalExtractionErrorsRawCs->GetXaxis()->GetBinWidth(1)  )  );
  TildePhiAfterSignalExtractionErrorsRawCs->Draw();
  TLatex* latex3b = new TLatex();
  latex3b->SetTextSize(0.05);
  latex3b->SetTextFont(42);
  latex3b->SetTextAlign(11);
  latex3b->SetNDC();
  latex3b->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3b->SetTextSize(0.045);
  latex3b->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/RawTildePhiCs.png", "recreate");





















  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  acceptCosThetaHe->GetXaxis()->SetTitleOffset(1.15);
  // acceptCosThetaHe->GetYaxis()->SetTitleOffset(1.25);
  acceptCosThetaHe->GetYaxis()->SetTitleOffset(1.55);
  acceptCosThetaHe->GetXaxis()->SetTitleSize(0.045);
  acceptCosThetaHe->GetYaxis()->SetTitleSize(0.045);
  acceptCosThetaHe->GetXaxis()->SetLabelSize(0.045);
  acceptCosThetaHe->GetYaxis()->SetLabelSize(0.045);
  acceptCosThetaHe->GetXaxis()->SetTitleFont(42);
  acceptCosThetaHe->GetYaxis()->SetTitleFont(42);
  acceptCosThetaHe->GetXaxis()->SetLabelFont(42);
  acceptCosThetaHe->GetYaxis()->SetLabelFont(42);
  acceptCosThetaHe->GetXaxis()->SetNdivisions(408);
  acceptCosThetaHe->GetYaxis()->SetRangeUser(0., acceptCosThetaHe->GetMaximum()*2);
  // acceptCosThetaHe->GetXaxis()->SetRangeUser(2, 6);
  acceptCosThetaHe->SetTitle(  ";cos(#theta); ACCxEFF [a.u.]");
  acceptCosThetaHe->Draw();
  TLatex* latexc = new TLatex();
  latexc->SetTextSize(0.05);
  latexc->SetTextFont(42);
  latexc->SetTextAlign(11);
  latexc->SetNDC();
  latexc->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latexc->SetTextSize(0.045);
  latexc->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/AccCosThetaHe.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  acceptCosThetaCs->GetXaxis()->SetTitleOffset(1.15);
  // acceptCosThetaCs->GetYaxis()->SetTitleOffset(1.25);
  acceptCosThetaCs->GetYaxis()->SetTitleOffset(1.55);
  acceptCosThetaCs->GetXaxis()->SetTitleSize(0.045);
  acceptCosThetaCs->GetYaxis()->SetTitleSize(0.045);
  acceptCosThetaCs->GetXaxis()->SetLabelSize(0.045);
  acceptCosThetaCs->GetYaxis()->SetLabelSize(0.045);
  acceptCosThetaCs->GetXaxis()->SetTitleFont(42);
  acceptCosThetaCs->GetYaxis()->SetTitleFont(42);
  acceptCosThetaCs->GetXaxis()->SetLabelFont(42);
  acceptCosThetaCs->GetYaxis()->SetLabelFont(42);
  acceptCosThetaCs->GetXaxis()->SetNdivisions(408);
  acceptCosThetaCs->GetYaxis()->SetRangeUser(0., acceptCosThetaCs->GetMaximum()*2);
  // acceptCosThetaCs->GetXaxis()->SetRangeUser(2, 6);
  acceptCosThetaCs->SetTitle( ";cos(#theta); ACCxEFF [a.u.]" );
  acceptCosThetaCs->Draw();
  TLatex* latexan = new TLatex();
  latexan->SetTextSize(0.05);
  latexan->SetTextFont(42);
  latexan->SetTextAlign(11);
  latexan->SetNDC();
  latexan->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latexan->SetTextSize(0.045);
  latexan->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/AccCosThetaCs.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  acceptPhiHe->GetXaxis()->SetTitleOffset(1.15);
  // acceptPhiHe->GetYaxis()->SetTitleOffset(1.25);
  acceptPhiHe->GetYaxis()->SetTitleOffset(1.55);
  acceptPhiHe->GetXaxis()->SetTitleSize(0.045);
  acceptPhiHe->GetYaxis()->SetTitleSize(0.045);
  acceptPhiHe->GetXaxis()->SetLabelSize(0.045);
  acceptPhiHe->GetYaxis()->SetLabelSize(0.045);
  acceptPhiHe->GetXaxis()->SetTitleFont(42);
  acceptPhiHe->GetYaxis()->SetTitleFont(42);
  acceptPhiHe->GetXaxis()->SetLabelFont(42);
  acceptPhiHe->GetYaxis()->SetLabelFont(42);
  acceptPhiHe->GetXaxis()->SetNdivisions(408);
  acceptPhiHe->GetYaxis()->SetRangeUser(0., acceptPhiHe->GetMaximum()*2);
  // acceptPhiHe->GetXaxis()->SetRangeUser(2, 6);
  acceptPhiHe->SetTitle( ";#varphi; ACCxEFF [a.u.]" );
  acceptPhiHe->Draw();
  TLatex* latex2an = new TLatex();
  latex2an->SetTextSize(0.05);
  latex2an->SetTextFont(42);
  latex2an->SetTextAlign(11);
  latex2an->SetNDC();
  latex2an->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2an->SetTextSize(0.045);
  latex2an->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/AccPhiHe.png", "recreate");


  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  acceptPhiCs->GetXaxis()->SetTitleOffset(1.15);
  // acceptPhiCs->GetYaxis()->SetTitleOffset(1.25);
  acceptPhiCs->GetYaxis()->SetTitleOffset(1.55);
  acceptPhiCs->GetXaxis()->SetTitleSize(0.045);
  acceptPhiCs->GetYaxis()->SetTitleSize(0.045);
  acceptPhiCs->GetXaxis()->SetLabelSize(0.045);
  acceptPhiCs->GetYaxis()->SetLabelSize(0.045);
  acceptPhiCs->GetXaxis()->SetTitleFont(42);
  acceptPhiCs->GetYaxis()->SetTitleFont(42);
  acceptPhiCs->GetXaxis()->SetLabelFont(42);
  acceptPhiCs->GetYaxis()->SetLabelFont(42);
  acceptPhiCs->GetXaxis()->SetNdivisions(408);
  acceptPhiCs->GetYaxis()->SetRangeUser(0., acceptPhiCs->GetMaximum()*2);
  // acceptPhiCs->GetXaxis()->SetRangeUser(2, 6);
  acceptPhiCs->SetTitle( ";#varphi; ACCxEFF [a.u.]" );
  acceptPhiCs->Draw();
  TLatex* latex2n = new TLatex();
  latex2n->SetTextSize(0.05);
  latex2n->SetTextFont(42);
  latex2n->SetTextAlign(11);
  latex2n->SetNDC();
  latex2n->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2n->SetTextSize(0.045);
  latex2n->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/AccPhiCs.png", "recreate");



  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  acceptTildePhiHe->GetXaxis()->SetTitleOffset(1.15);
  // acceptTildePhiHe->GetYaxis()->SetTitleOffset(1.25);
  acceptTildePhiHe->GetYaxis()->SetTitleOffset(1.55);
  acceptTildePhiHe->GetXaxis()->SetTitleSize(0.045);
  acceptTildePhiHe->GetYaxis()->SetTitleSize(0.045);
  acceptTildePhiHe->GetXaxis()->SetLabelSize(0.045);
  acceptTildePhiHe->GetYaxis()->SetLabelSize(0.045);
  acceptTildePhiHe->GetXaxis()->SetTitleFont(42);
  acceptTildePhiHe->GetYaxis()->SetTitleFont(42);
  acceptTildePhiHe->GetXaxis()->SetLabelFont(42);
  acceptTildePhiHe->GetYaxis()->SetLabelFont(42);
  acceptTildePhiHe->GetXaxis()->SetNdivisions(408);
  acceptTildePhiHe->GetYaxis()->SetRangeUser(0., acceptTildePhiHe->GetMaximum()*2);
  // acceptTildePhiHe->GetXaxis()->SetRangeUser(2, 6);
  acceptTildePhiHe->SetTitle( ";#tilde{#varphi}; ACCxEFF [a.u.]"  );
  acceptTildePhiHe->Draw();
  TLatex* latex3an = new TLatex();
  latex3an->SetTextSize(0.05);
  latex3an->SetTextFont(42);
  latex3an->SetTextAlign(11);
  latex3an->SetNDC();
  latex3an->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3an->SetTextSize(0.045);
  latex3an->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  gPad->SaveAs("pngResults/AccTildePhiHe.png", "recreate");




  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  acceptTildePhiCs->GetXaxis()->SetTitleOffset(1.15);
  // acceptTildePhiCs->GetYaxis()->SetTitleOffset(1.25);
  acceptTildePhiCs->GetYaxis()->SetTitleOffset(1.55);
  acceptTildePhiCs->GetXaxis()->SetTitleSize(0.045);
  acceptTildePhiCs->GetYaxis()->SetTitleSize(0.045);
  acceptTildePhiCs->GetXaxis()->SetLabelSize(0.045);
  acceptTildePhiCs->GetYaxis()->SetLabelSize(0.045);
  acceptTildePhiCs->GetXaxis()->SetTitleFont(42);
  acceptTildePhiCs->GetYaxis()->SetTitleFont(42);
  acceptTildePhiCs->GetXaxis()->SetLabelFont(42);
  acceptTildePhiCs->GetYaxis()->SetLabelFont(42);
  acceptTildePhiCs->GetXaxis()->SetNdivisions(408);
  acceptTildePhiCs->GetYaxis()->SetRangeUser(0., acceptTildePhiCs->GetMaximum()*2);
  // acceptTildePhiCs->GetXaxis()->SetRangeUser(2, 6);
  acceptTildePhiCs->SetTitle(  ";#tilde{#varphi}; ACCxEFF [a.u.]" );
  acceptTildePhiCs->Draw();
  TLatex* latex3n = new TLatex();
  latex3n->SetTextSize(0.05);
  latex3n->SetTextFont(42);
  latex3n->SetTextAlign(11);
  latex3n->SetNDC();
  latex3n->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3n->SetTextSize(0.045);
  latex3n->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, CS");
  gPad->SaveAs("pngResults/AccTildePhiCs.png", "recreate");

}
