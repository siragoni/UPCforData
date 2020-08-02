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



void ReadHistoCounter(){

  // TFile* fileList = new TFile("AnalysisResultsLHC18l7Michal.root");
  TFile* fileList = new TFile("AnalysisResultsMichal.root");
  // TFile* fileList = new TFile("MCtrainResults/2019-09-17/kCohJpsiToMu/AnalysisResults.root");
  TDirectory* dir = fileList->GetDirectory("MyTask");
  TList* listings;
  dir->GetObject("MyOutputContainer", listings);
  TH1F* fCounterH   = (TH1F*)listings->FindObject("fCosThetaCsFrameTwentyfiveBinsH");
  Double_t counters[25];
  for ( Int_t i = 0; i < fCounterH->GetNbinsX(); i++ ) {
    counters[i] = fCounterH->GetBinContent(i+1);
    cout << i << " & " << counters[i] << " \\ " << endl;
  }
}
