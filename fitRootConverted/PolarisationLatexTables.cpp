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

#include "TDatime.h"


#include "TH2.h"

//_____________________________________________________________________________
/* - Fit function for the ZNC plots.
 * -
 */
void PolarisationLatexTables( const Int_t selectionFlag ){

  TDatime d;
  TFile* fileDataRawCosTheta = 0x0;
  TFile* fileDataRawPhi      = 0x0;
  TFile* fileDataRawTildePhi = 0x0;
  if (        selectionFlag == 0 ){
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaCS/CosThetaCSFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiCS/PhiCsFrame.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiCS/TildePhiCsFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( selectionFlag == 1 ){
    fileDataRawCosTheta = new TFile( Form("pngResults/%d-%2.2d-%2.2d/CosThetaHE/CosThetaHeFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawPhi      = new TFile( Form("pngResults/%d-%2.2d-%2.2d/PhiHE/PhiHeFrame.root",           d.GetYear(), d.GetMonth(), d.GetDay() ) );
    fileDataRawTildePhi = new TFile( Form("pngResults/%d-%2.2d-%2.2d/TildePhiHE/TildePhiHeFrame.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  }

  TH1F* CosThetaAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  CosThetaAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawCosThetaH = (TH1F*) CosThetaAfterSignalExtractionErrorsRawH->Clone("RawCosThetaH");

  TH1F* PhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  PhiAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawPhiH = (TH1F*) PhiAfterSignalExtractionErrorsRawH->Clone("RawPhiH");

  TH1F* TildePhiAfterSignalExtractionErrorsRawH = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
  TildePhiAfterSignalExtractionErrorsRawH->Sumw2();
  TH1F* RawTildePhiH = (TH1F*) TildePhiAfterSignalExtractionErrorsRawH->Clone("RawTildePhiH");

  Double_t dNdCosTheta[25];
  Double_t dNdCosThetaErr[25];
  Double_t dNdPhi[25];
  Double_t dNdPhiErr[25];
  Double_t dNdTildePhi[25];
  Double_t dNdTildePhiErr[25];
  for (size_t i = 1; i < RawCosThetaH->GetNbinsX(); i++) {
    dNdCosTheta[i]    = RawCosThetaH->GetBinContent(i);
    dNdPhi[i]         = RawPhiH     ->GetBinContent(i);
    dNdTildePhi[i]    = RawTildePhiH->GetBinContent(i);
    dNdCosThetaErr[i] = RawCosThetaH->GetBinError(i);
    dNdPhiErr[i]      = RawPhiH     ->GetBinError(i);
    dNdTildePhiErr[i] = RawTildePhiH->GetBinError(i);
  }

  cout << endl;
  cout << endl;
  cout << "\\begin{table}[h!]"          << endl;
  cout << "  \\begin{center}"           << endl;
  cout << "     \\caption{$\\frac{\\text{d}N}{\\text{d}\\cos\\theta}$ distribution.}" << endl;
  cout << "     \\label{tab:table1}" << endl;
  cout << "     \\begin{tabular}{l|c|r} " << endl;
  cout << "        \\textbf{Bin number} & \\textbf{$\\frac{\\text{d}N}{\\text{d}\\cos\\theta}$} & \\textbf{Error (stat.)}\\\\ " << endl;
  cout << "        \\hline" << endl;
  for( size_t i = 1; i < RawCosThetaH->GetNbinsX(); i++ ){
    cout << i << " & " << dNdCosTheta[i] << " & " << dNdCosThetaErr[i] << " \\\\ " << endl;
  }
  cout << "     \\end{tabular}" << endl;
  cout << "   \\end{center}"    << endl;
  cout << "\\end{table}"        << endl;


  cout << endl;
  cout << endl;
  cout << "\\begin{table}[h!]"          << endl;
  cout << "  \\begin{center}"           << endl;
  cout << "     \\caption{$\\frac{\\text{d}N}{\\text{d}\\phi}$ distribution.}" << endl;
  cout << "     \\label{tab:table1}" << endl;
  cout << "     \\begin{tabular}{l|c|r} " << endl;
  cout << "        \\textbf{Bin number} & \\textbf{$\\frac{\\text{d}N}{\\text{d}\\phi}$} & \\textbf{Error (stat.)}\\\\ " << endl;
  cout << "        \\hline" << endl;
  for( size_t i = 1; i < RawCosThetaH->GetNbinsX(); i++ ){
    cout << i << " & " << dNdPhi[i] << " & " << dNdPhiErr[i] << " \\\\ " << endl;
  }
  cout << "     \\end{tabular}" << endl;
  cout << "   \\end{center}"    << endl;
  cout << "\\end{table}"        << endl;


  cout << endl;
  cout << endl;
  cout << "\\begin{table}[h!]"          << endl;
  cout << "  \\begin{center}"           << endl;
  cout << "     \\caption{$\\frac{\\text{d}N}{\\text{d}\\widetilde{\\phi}}$ distribution.}" << endl;
  cout << "     \\label{tab:table1}" << endl;
  cout << "     \\begin{tabular}{l|c|r} " << endl;
  cout << "        \\textbf{Bin number} & \\textbf{$\\frac{\\text{d}N}{\\text{d}\\widetilde{\\phi}}$} & \\textbf{Error (stat.)}\\\\ " << endl;
  cout << "        \\hline" << endl;
  for( size_t i = 1; i < RawCosThetaH->GetNbinsX(); i++ ){
    cout << i << " & " << dNdTildePhi[i] << " & " << dNdTildePhiErr[i] << " \\\\ " << endl;
  }
  cout << "     \\end{tabular}" << endl;
  cout << "   \\end{center}"    << endl;
  cout << "\\end{table}"        << endl;


  // TFile f("pngResults/PolarisationCorrectedCs1D.root", "recreate");
  // acceptanceCosTheta->Write();
  // CorrCosThetaH     ->Write();
  // acceptancePhi     ->Write();
  // CorrPhiH          ->Write();
  // acceptanceTildePhi->Write();
  // CorrTildePhiH     ->Write();
  // // AccErrors         ->Write();
  // // EntErrors         ->Write();
  // // ReconTheta        ->Write();
  // f.Close();
}
