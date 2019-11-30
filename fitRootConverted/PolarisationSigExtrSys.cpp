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
/* - Compute the Signal extraction systematics.
 * - Compute the Fit Range Variation systematics.
 * -
 */
void PolarisationSigExtrSystematics(){

  /* - Open all files.
   * -
   * - 6 Signal Extraction bins.
   * - 3 Range Variation bins.
   * -
   */
  TFile* FitResultFile[6][3];
  for ( Int_t SigExBin = 0; SigExBin < 6; SigExBin++ ) {
    FitResultFile[6][3]
  }

}
