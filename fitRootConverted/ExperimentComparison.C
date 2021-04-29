#include <TROOT.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TBox.h>
#include <TWbox.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TMarker.h>
#include <stdio.h>
#include <iostream>

void drawAverage(double mean, double nErr, double pErr) {

  TBox* band = new TBox(mean-nErr,0.0,mean+pErr,1.0);
  band->SetFillColor(188);
  band->Draw();
  TLine* center = new TLine(mean, 0.0, mean, 1.0);
  center->SetLineStyle(2); // dash
  center->SetLineWidth(4);
  center->SetLineColor(9);
  center->Draw();
  TLine* left = new TLine(mean-nErr, 0.0, mean-nErr, 1.0);
  left->SetLineStyle(1);
  left->SetLineWidth(2);
  left->SetLineColor(9);
  left->Draw();
  TLine* right = new TLine(mean+pErr, 0.0, mean+pErr, 1.0);
  right->SetLineStyle(1);
  right->SetLineWidth(2);
  right->SetLineColor(9);
  right->Draw();

  return;
}

void drawBkgBox(int i, int aux[6], double vstep, TH2F* histo, int debug=0) {

  double lowY = (i+1)*vstep; // margins are 1 vstep each
  double uppY = (i+2)*vstep;

  // double lowX = histo->GetBinLowEdge(1);
  double lowX = histo->GetXaxis()->GetBinLowEdge(1);
  // double uppX = histo->GetBinLowEdge(histo->GetNbinsX()) +
  //               histo->GetBinWidth(histo->GetNbinsX());
  double uppX = histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()) +
                histo->GetXaxis()->GetBinWidth(histo->GetNbinsX());

  TBox* band = new TBox(lowX,lowY,uppX,uppY);
  band->SetFillColor(aux[2]);
  band->Draw();

  if (debug!=0) {

    TLine* l1 = new TLine(lowX,lowY,uppX,lowY); l1->Draw();
    TLine* l2 = new TLine(lowX,uppY,uppX,uppY); l2->Draw();
    TLine* l3 = new TLine(lowX,lowY,lowX,uppY); l3->Draw();
    TLine* l4 = new TLine(uppX,lowY,uppX,uppY); l4->Draw();
  }

  return;
}

void drawMeasurement(int i, double m[5], char label[2][100], int aux[6],
		     double vstep, TH2F* histo, TCanvas* canvas) {

  // Text size [SetTextSize(Float_t)] is expressed in percentage of
  // the current pad/canvas size (isn't that smart?).  So, we have to
  // implement a correction in order to have the same text size for
  // different canvae size ... but it's more screwed up than this:
  // 'canvas size' is its height for horizontal canvae and is its
  // width for vertical canvae!!!
  double coef = 50.0/TMath::Min(canvas->GetWindowWidth(),
				canvas->GetWindowHeight());

  double lowY = (i+1)*vstep;
  double uppY = (i+2)*vstep;

  double lowX = histo->GetXaxis()->GetBinLowEdge(1);
  double uppX = histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()) +
                histo->GetXaxis()->GetBinWidth(histo->GetNbinsX());
  double widthX = uppX - lowX;

  // REMEMBER: Y-size of the histogram does not scale.  It is defined
  // in absolute terms: 50points*(NUM+2) [2 for upper and lower margins]
  // ... code accordingly
  // there are exactly (NUM+2 vsteps). Take advanatge of the fact that
  // y-range of the histogram is [0...1]

  double startX = lowX + 0.01*widthX;

  // if (i==0 || i == 3) { // HFAG Label
  if (i==10) { // HFAG Label

    TPaveText* textA = new TPaveText(startX, lowY+0.1*vstep,
				     startX+0.085*widthX, lowY+0.35*vstep,
				     "BR");
    textA->SetTextAlign(22);
    textA->SetFillColor(10);
    textA->SetTextColor(1);
    textA->SetLineColor(1);
    textA->SetBorderSize(1);
    TText* tA = textA->AddText("2     0     2     1");
    tA->SetTextSize(0.17*coef);
    tA->SetTextFont(42);
    textA->Draw();

    TPaveText* textB = new TPaveText(startX, lowY+0.35*vstep,
				     startX+0.085*widthX, lowY+0.90*vstep,
				     "BR");
    textB->SetTextAlign(22);
    textB->SetFillColor(1);
    textB->SetTextColor(10);
    textB->SetLineColor(1);
    textB->SetBorderSize(1);
    TText* tB = textB->AddText("ALICE");
    tB->SetTextSize(0.45*coef);
    tB->SetTextFont(62);
    textB->Draw();

    TPaveText* textC = new TPaveText(startX+0.09*widthX, lowY,
				     startX+0.09*widthX, uppY, "BR");
    textC->SetTextAlign(12);
    // textC->SetFillColor(aux[2]);
    textC->SetTextColor(1); // or aux[1], if you like colors
    textC->SetLineColor(1);
    textC->SetBorderSize(0);
    // TText* tC = textC->AddText("UPC");
    TText* tC = 0x0;
    if      (i==0) { tC = textC->AddText("(a)"); }
    else if (i==3) { tC = textC->AddText("(b)"); }
    else           { tC = textC->AddText("UPC"); }
    tC->SetTextSize(0.5*coef);
    tC->SetTextFont(aux[3]);
    textC->Draw();
  }

  else { // normal user set-able labels

    TPaveText* text = new TPaveText(startX, lowY,
				    startX, uppY, "BR");
    text->SetTextAlign(12);
    // text->SetFillColor(aux[2]);
    text->SetFillColor(0);
    text->SetTextColor(aux[1]);
    text->SetLineColor(1);
    text->SetBorderSize(0);
    TText* t0 = text->AddText("  ");
    t0->SetTextSize(0.01*coef);
    t0->SetTextFont(aux[3]);
    TText* t1 = text->AddText(label[0]);
    t1->SetTextSize(0.5*coef);
    t1->SetTextFont(aux[3]);
    TText* t2 = text->AddText(label[1]);
    t2->SetTextSize(0.35*coef);
    t2->SetTextFont(aux[3]);
    text->Draw();
  }

  double ypos = 0.5*(lowY+uppY);
  double mean = m[0];
  double nErr1 = m[1];
  double pErr1 = m[2];
  double nErr2 = sqrt(m[1]*m[1]+m[3]*m[3]);
  double pErr2 = sqrt(m[2]*m[2]+m[4]*m[4]);

  // draw TGraphAsymmErrors 1 (stat only) |---*---|
  TMarker* measurement = new TMarker(mean, ypos, aux[5]);
  measurement->SetMarkerColor(aux[1]);
  measurement->SetMarkerStyle(aux[5]);
  // not clear what are the units for the marker
  //measurement->SetMarkerSize(30.0*coef);
  measurement->SetMarkerSize(vstep/0.04);
  //measurement->SetMarkerSize(1.75);
  measurement->Draw();

  double vsizeErr1 = 0.09*vstep;
  TLine* l1 = new TLine(mean, ypos, mean-nErr1, ypos);
  l1->SetLineWidth(aux[4]);
  l1->SetLineColor(aux[1]);
  l1->Draw();
  TLine* l2 = new TLine(mean, ypos, mean+pErr1, ypos);
  l2->SetLineWidth(aux[4]);
  l2->SetLineColor(aux[1]);
  l2->Draw();
  TLine* l3 = new TLine(mean-nErr1, ypos-vsizeErr1,
			mean-nErr1, ypos+vsizeErr1);
  l3->SetLineWidth(aux[4]);
  l3->SetLineColor(aux[1]);
  l3->Draw();
  TLine* l4 = new TLine(mean+pErr1, ypos-vsizeErr1,
			mean+pErr1, ypos+vsizeErr1);
  l4->SetLineWidth(aux[4]);
  l4->SetLineColor(aux[1]);
  l4->Draw();

  // overlay TGraphAsymmErrors 2 (stat+syst) |----*-----|
  double vsizeErr2 = 0.12*vstep;
  TLine* l5 = new TLine(mean, ypos, mean-nErr2, ypos);
  l5->SetLineWidth(aux[4]);
  l5->SetLineColor(aux[1]);
  l5->Draw();
  TLine* l6 = new TLine(mean, ypos, mean+pErr2, ypos);
  l6->SetLineWidth(aux[4]);
  l6->SetLineColor(aux[1]);
  l6->Draw();
  TLine* l7 = new TLine(mean-nErr2, ypos-vsizeErr2,
			mean-nErr2, ypos+vsizeErr2);
  l7->SetLineWidth(aux[4]);
  l7->SetLineColor(aux[1]);
  l7->Draw();
  TLine* l8 = new TLine(mean+pErr2, ypos-vsizeErr2,
			mean+pErr2, ypos+vsizeErr2);
  l8->SetLineWidth(aux[4]);
  l8->SetLineColor(aux[1]);
  l8->Draw();

  // draw measurement label "XXX+/-YY+/-ZZ"
  TPaveText* num = new TPaveText(uppX-0.32*widthX,
				 lowY, uppX-0.02*widthX, uppY, "BR");
  num->SetTextAlign(12);
  // num->SetFillColor(aux[2]);
  num->SetFillColor(0);
  num->SetTextColor(aux[1]);
  num->SetLineColor(aux[1]);
  num->SetBorderSize(0);
  TString str;
  char s[100];
  if (aux[0]==2) { // 2 decimal digits
    // sprintf(s, "%4.2f#color[%d]{X}", m[0], aux[2]); str +=s;
    sprintf(s, "%4.2f", m[0]); str +=s;
    if (m[1]==m[2]) { // sym. stat. errors
      sprintf(s, "#pm%4.2f",  m[1]); str +=s;
    } else {
      sprintf(s, "^{+%4.2f}",  m[2]); str +=s;
      sprintf(s, "_{-#color[%d]{|}%4.2f}", aux[2], m[1]); str +=s;
    }
    if (m[3]!=0.0 || m[4]!=0.0) {
      if (m[3]==m[4]) { // sym. syst. errors
  // sprintf(s, "#color[%d]{X}#pm%4.2f", aux[2], m[3]); str +=s;
	sprintf(s, "#pm%4.2f", m[3]); str +=s;
      } else {
  // sprintf(s, "#color[%d]{X}", aux[2]); str +=s;
	sprintf(s, "^{+%4.2f}",  m[4]); str +=s;
	sprintf(s, "_{-#color[%d]{|}%4.2f}", aux[2], m[3]); str +=s;
      }
    }
  }
  if (aux[0]==3) { // 3 decimal digits
    // sprintf(s, "%5.3f#color[%d]{X}", m[0], aux[2]); str +=s;
    sprintf(s, "%5.3f", m[0]); str +=s;
    if (m[1]==m[2]) { // sym. stat. errors
      sprintf(s, "#pm%5.3f",  m[1]); str +=s;
    } else {
      sprintf(s, "^{+%5.3f}",  m[2]); str +=s;
      sprintf(s, "_{-#color[%d]{|}%5.3f}", aux[2], m[1]); str +=s;
    }
    if (m[3]!=0.0 || m[4]!=0.0) {
      if (m[3]==m[4]) { // sym. syst. errors
  // sprintf(s, "#color[%d]{X}#pm%5.3f", aux[2], m[3]); str +=s;
	sprintf(s, "#pm%5.3f", m[3]); str +=s;
      } else {
  // sprintf(s, "#color[%d]{X}", aux[2]); str +=s;
	sprintf(s, "^{+%5.3f}",  m[4]); str +=s;
	sprintf(s, "_{-#color[%d]{|}%5.3f}", aux[2], m[3]); str +=s;
      }
    }
  }
  TText* n0 = num->AddText(str);
  n0->SetTextSize(0.5*coef);
  n0->SetTextFont(62);
  num->Draw();

  return;
}

// =====================================================================
// =====================================================================
// === this is where the magic happens =================================
// =====================================================================

void bs() {

  gStyle->SetFrameLineWidth(1); // RENE look, try setting this to 4
                                // to see more vividly

  // const int NUM = 11;
  // // mean, -stat, +stat, -syst, +syst
  // double m[NUM][5] = {
  //   1.461, 0.057, 0.057, 0.0, 0.0,
  //   1.363, 0.100, 0.100, 0.010, 0.007,
  //   1.42, 0.13, 0.14, 0.03, 0.03,
  //   1.53, 0.15, 0.16, 0.07, 0.07,
  //   1.36, 0.09, 0.09, 0.05, 0.06,
  //   1.34, 0.19, 0.23, 0.05, 0.05,
  //   1.72, 0.19, 0.20, 0.17, 0.18,
  //   1.50, 0.15, 0.16, 0.04, 0.04,
  //   1.47, 0.14, 0.14, 0.08, 0.08,
  //   1.60, 0.26, 0.26, 0.15, 0.13,
  //   1.54, 0.13, 0.14, 0.04, 0.04
  // };
  const int NUM = 6;


  // OLD
  // Double_t LambdaThetaHeALICE = 1.208;
  // Double_t LambdaPhiHeALICE   = 0.049;
  // Double_t StatThetaHeALICE   = 0.155;
  // Double_t StatPhiHeALICE     = 0.026;
  // Double_t SysThetaHeALICE    = 0.251;
  // Double_t SysPhiHeALICE      = 0.020;

  Double_t LambdaThetaHeALICE = 0.956;
  Double_t LambdaPhiHeALICE   = 0.040;
  Double_t StatThetaHeALICE   = 0.050;
  Double_t StatPhiHeALICE     = 0.024;
  Double_t SysThetaHeALICE    = 0.068;
  Double_t SysPhiHeALICE      = 0.012;



  Double_t r00_04_HeALICE      = (1.-LambdaThetaHeALICE)/(3.+LambdaThetaHeALICE);
  Double_t r1minus1_04_HeALICE = 0.5*LambdaPhiHeALICE*(1.+r00_04_HeALICE);


  // Double_t Stat_r00_HeALICE    = r00_04_HeALICE*( StatThetaHeALICE*StatThetaHeALICE/( (1.-LambdaThetaHeALICE)*(1.-LambdaThetaHeALICE) )   +   )
  Double_t Stat_r00_HeALICE    = 4.*StatThetaHeALICE/( (3.+LambdaThetaHeALICE)*(3.+LambdaThetaHeALICE) );
  Double_t Sys_r00_HeALICE     = 4.*SysThetaHeALICE /( (3.+LambdaThetaHeALICE)*(3.+LambdaThetaHeALICE) );


  Double_t Stat_r11_HeALICE    = r1minus1_04_HeALICE*(  StatPhiHeALICE*StatPhiHeALICE/(LambdaPhiHeALICE*LambdaPhiHeALICE) +  Stat_r00_HeALICE*Stat_r00_HeALICE/( (1.+r00_04_HeALICE)*(1.+r00_04_HeALICE) )  );
  Double_t Sys_r11_HeALICE     = r1minus1_04_HeALICE*(  SysPhiHeALICE*SysPhiHeALICE  /(LambdaPhiHeALICE*LambdaPhiHeALICE) +  Sys_r00_HeALICE *Sys_r00_HeALICE /( (1.+r00_04_HeALICE)*(1.+r00_04_HeALICE) )  );


  // mean, -stat, +stat, -syst, +syst
  double m[NUM][5] = {
    r00_04_HeALICE, Stat_r00_HeALICE, Stat_r00_HeALICE, Sys_r00_HeALICE, Sys_r00_HeALICE,
    0.003, 0.039, 0.039, 0.028, 0.028,
    0.120, 0.080, 0.080, 0.015, 0.013,
    r1minus1_04_HeALICE, Stat_r11_HeALICE, Stat_r11_HeALICE, Sys_r11_HeALICE, Sys_r11_HeALICE,
   -0.011, 0.036, 0.036, 0.030, 0.030,
    0.34,  0.09,  0.09,  0.06,  0.03
  };


  // label text, sub-label text
  // Note: TString does not work, b/c oce cannot pass an array of TStrings
  // as an argument to a function
  // char label[NUM][2][100] = {
  //   "2003 World Average", "",
  //   "CDF Preliminary", "2004 (not included in the average)",
  //   "DELPHI", "2000",
  //   "DELPHI", "2000",
  //   "CDF"   , "1999",
  //   "CDF"   , "1998",
  //   "OPAL"  , "1998",
  //   "OPAL"  , "1998",
  //   "ALEPH" , "1998",
  //   "DELPHI", "1996",
  //   "ALEPH" , "1996"
  // };
  char label[NUM][2][100] = {
    "ALICE" , "(a)",
    "H1" ,    "(a)",
    "ZEUS" ,  "(a)",
    "ALICE" , "(b)",
    "H1" ,    "(b)",
    "ZEUS" ,  "(b)"
    // "ALICE" , ""
  };


  // format:
  // # dec. digits, fgColor, bgColor, fontStyle, lineWidth, markerStyle
  int aux[NUM][6] = {
    3, 4, 179, 62, 4, 20,
    3, 2,  10, 62, 2, 23,
    3, 1,  18, 62, 2, 22,
    3, 4, 179, 62, 4, 20,
    3, 2,  10, 62, 2, 23,
    3, 1,  18, 62, 2, 22


    // 3, 4, 179, 62, 4, 20,
    // 3, 2,  10, 62, 2, 23,
    // 2, 1,  18, 42, 2, 22,
    // 3, 4, 179, 62, 4, 20,
    // 3, 2,  10, 62, 2, 23,
    // 2, 1,  18, 42, 2, 22



    // 2, 1,  10, 42, 2, 22,
    // 2, 1, 196, 42, 2, 22,
    // 2, 1,  10, 42, 2, 22,
    // 2, 1,  10, 42, 2, 22,
    // 2, 1,  10, 42, 2, 22,
    // 2, 1,  10, 42, 2, 22,
    // 2, 1,  10, 42, 2, 22,
    // 2, 1,  10, 42, 2, 22
  };

  // determning the x size of the plot
  double lowX = 1.0e+32;
  double uppX = -1.0e+32;
  // determining the x range spanned by the measurements w/ error bars
  for (int i=0; i!=NUM; ++i) {

    double l = m[i][0] - sqrt(m[i][1]*m[i][1]+m[i][3]*m[i][3]);
    double u = m[i][0] + sqrt(m[i][2]*m[i][2]+m[i][4]*m[i][4]);

    if (l<lowX) lowX = l;
    if (u>uppX) uppX = u;
  }
  lowX = lowX - 0.52*(uppX-lowX); // for measurement text label
  uppX = uppX + 0.48*(uppX-lowX); // for measurement number+/-uncert.

  TH2F* lft = new TH2F("lft", "", 50, lowX, uppX, 1, 0.0, 1.0);

  // height = NUM*unitHeight + 2*spacers + 1*bottomMargin + 0.5*topMargin
  const double unitHeight = 50.0; // make it an even number
  const double height = (double(NUM)+2.0)*unitHeight+100.0+30.0;
  // how much to step each time to cover the vertical range of the histo in
  // exactly NUM+2 steps
  const double vstep  = unitHeight/(height-100.0-30.0);
  const double width = 800.0;

  printf("Canvas: width=%d, height=%d\n",
	 TMath::Nint(width), TMath::Nint(height));
  printf("Y-step = %6.4f\n", vstep);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 200, 0,
				TMath::Nint(width), TMath::Nint(height));
  canvas->SetFillColor(10);
  canvas->SetRightMargin(10.0/width);
  canvas->SetLeftMargin(10.0/width);
  canvas->SetBottomMargin(100.0/height);
  // canvas->SetTopMargin(30.0/height);
  canvas->SetTopMargin(10.0/height);
  canvas->Draw();
  canvas->cd();

  //printf("TopMargin   : %6.4f\n", canvas->GetTopMargin());
  //printf("BottomMargin: %6.4f\n", canvas->GetBottomMargin());

  // TickLength is a mess.  Don't touch it.
  // Pray that it does not exceed the margins
  //gStyle->SetTickLength(0.2*vstep);

  TAxis* xaxis = lft->GetXaxis();
  TAxis* yaxis = lft->GetYaxis();
  //xaxis->CenterTitle(kTRUE);  // x axis label is put
  //xaxis->SetTitleSize(0.06);  // "by hand" because
  //xaxis->SetTitleFont(42);    // we use different
  //xaxis->SetTitleOffset(1.0); // fonts in it
  xaxis->SetLabelOffset(0.01);
  xaxis->SetLabelSize(0.03);
  xaxis->SetLabelFont(42);
  yaxis->SetLabelSize(0.0);
  yaxis->SetNdivisions(-1);

  lft->SetXTitle("");
  lft->SetYTitle("");
  lft->SetStats(kFALSE);
  lft->SetTitle("");
  lft->Draw();

  // x axis label (combining different font types/sizes)
  double averageInNDC = (m[0][0]-lowX)/(uppX-lowX) *
    (1.0 - canvas->GetLeftMargin() - canvas->GetRightMargin()) +
    canvas->GetLeftMargin();
  double coef = 70.0/TMath::Min(canvas->GetWindowWidth(),
				canvas->GetWindowHeight());
  // TPaveText* title1 = new TPaveText(averageInNDC-0.05, 0.0, averageInNDC,
	// 			    0.7*canvas->GetBottomMargin(), "NDCBR");
  // title1->SetTextAlign(32);
  // title1->SetFillColor(10);
  // title1->SetTextColor(1);
  // title1->SetLineColor(1);
  // title1->SetBorderSize(0);
  // // title1->SetTextSize(0.55*coef);
  // title1->SetTextSize(0.055);
  // // title1->AddText("#tau_{B_{s}}, ");
  // title1->AddText("(a) = r_{00}^{04}, (b) = r_{1, -1}^{04} ");
  // title1->Draw();

  // TPaveText* title2 = new TPaveText(averageInNDC, 0.01,
	// 			    averageInNDC + 0.05,
	// 			    0.7*canvas->GetBottomMargin()+0.01, "NDCBR");
  // title2->SetTextAlign(12);
  // title2->SetFillColor(10);
  // title2->SetTextColor(1);
  // title2->SetLineColor(1);
  // title2->SetBorderSize(0);
  // title2->SetTextSize(0.45*coef);
  // //title2->AddText("#times10^{-12}s");
  // title2->AddText("picoseconds");
  // title2->Draw();

  // mirror the x-axis to the top of the plot
  TGaxis* upax = new TGaxis(lowX, 1.0, uppX, 1.0, lowX, uppX,
			    xaxis->GetNdivisions(), "-");
  upax->SetLabelOffset(-0.01);
  upax->SetLabelSize(0.03);
  upax->SetLabelFont(42);
  // upax->Draw();



  // for (int i=0; i!=NUM; ++i) {
  //
  //   drawBkgBox(i, aux[i], vstep, lft); // RENE look, function is defined at the
  //                                      // beginning of the file
  // }
  // drawAverage(m[0][0], m[0][1], m[0][2]);
  canvas->RedrawAxis(); // RENE, look
  gPad->RedrawAxis();  // do both, niether seems to help
  // upax->Draw();

  for (int i=0; i!=NUM; ++i) {

    drawMeasurement(i, m[i], label[i], aux[i], vstep, lft, canvas);
  }


  TLatex* latex = new TLatex();
  latex->SetTextSize(0.055);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  // latex->DrawLatex(0.37,0.1,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->DrawLatex(0.2,0.92,"ALICE Public, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  // latex->DrawLatex(0.05,0.1,"(a) = r_{00}^{04}, (b) = r_{1, -1}^{04}");
  latex->DrawLatex(0.37,0.1,"(a) = r_{00}^{04}, (b) = r_{1, -1}^{04}");


  // canvas->Print("eps/bs-hfag2004.eps");
  canvas->SaveAs("pngResults/H1ZEUSv2.pdf");
}
