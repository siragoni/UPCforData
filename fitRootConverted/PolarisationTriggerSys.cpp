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
#include "TGraphErrors.h"
#include "TGraph.h"

//_____________________________________________________________________________
/* -
 * -
 */
void PolarisationTriggerSys( const char* AnalysisName, Int_t TriggerMode ){

  /* - Open the file containing all the TLists.
   * -
   * - Recreate new files with single outpust lists.
   * - To have the bash script go through them.
   * -
   */

  TFile*      fileList;
  TDirectory* dir;
  TList*      listings;

  fileList = new TFile(AnalysisName);
  dir      = fileList->GetDirectory("MyTask");
  dir->GetObject( Form("MyOutputContainer%d", TriggerMode) , listings );


  TFile* fileOutput     = new TFile( Form("Polarisation_%d.root", TriggerMode ), "recreate" );
  TDirectory* dirOutput = fileOutput->mkdir("MyTask");
  dirOutput     ->cd();
  TList* listingsOutput = (TList*) listings->Clone("MyOutputContainer");
  listingsOutput->Write("MyOutputContainer", 1);
  fileList      ->Close();
  fileOutput    ->Close();






}
