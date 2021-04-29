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


void Conversion_QT( Int_t mode = 0 ) {


  // Double_t LambdaThetaHe      = 1.208;
  // Double_t LambdaPhiHe        = 0.049;
  // Double_t LambdaThetaPhiHe   =-0.032;
  // Double_t StatThetaHe        = 0.155;
  // Double_t StatPhiHe          = 0.026;
  // Double_t StatThetaPhiHe     = 0.037;
  // Double_t SysThetaHe         = 0.251;
  // Double_t SysPhiHe           = 0.023;
  // Double_t SysThetaPhiHe      = 0.006;
  //
  //
  // Double_t LambdaThetaCs      = 1.326;
  // Double_t LambdaPhiCs        = 0.052;
  // Double_t LambdaThetaPhiCs   =-0.038;
  // Double_t StatThetaCs        = 0.159;
  // Double_t StatPhiCs          = 0.026;
  // Double_t StatThetaPhiCs     = 0.037;
  // Double_t SysThetaCs         = 0.226;
  // Double_t SysPhiCs           = 0.025;
  // Double_t SysThetaPhiCs      = 0.010;

  Double_t LambdaThetaHe      = 0.956;
  Double_t LambdaPhiHe        = 0.048;
  Double_t LambdaThetaPhiHe   =-0.030;
  Double_t StatThetaHe        = 0.050;
  Double_t StatPhiHe          = 0.024;
  Double_t StatThetaPhiHe     = 0.034;
  Double_t SysThetaHe         = 0.068;
  Double_t SysPhiHe           = 0.012;
  Double_t SysThetaPhiHe      = 0.006;


  Double_t LambdaThetaCs      = 0.965;
  Double_t LambdaPhiCs        = 0.038;
  Double_t LambdaThetaPhiCs   =-0.028;
  Double_t StatThetaCs        = 0.050;
  Double_t StatPhiCs          = 0.023;
  Double_t StatThetaPhiCs     = 0.035;
  Double_t SysThetaCs         = 0.086;
  Double_t SysPhiCs           = 0.013;
  Double_t SysThetaPhiCs      = 0.009;



  if        ( mode == 1 ) {
    LambdaThetaHe += TMath::Sqrt(StatThetaHe*StatThetaHe + SysThetaHe*SysThetaHe);
    LambdaThetaCs += TMath::Sqrt(StatThetaCs*StatThetaCs + SysThetaCs*SysThetaCs);
  } else if ( mode == 2 ) {
    LambdaThetaHe -= TMath::Sqrt(StatThetaHe*StatThetaHe + SysThetaHe*SysThetaHe);
    LambdaThetaCs -= TMath::Sqrt(StatThetaCs*StatThetaCs + SysThetaCs*SysThetaCs);
  } else if ( mode == 3 ) {
    LambdaPhiHe += TMath::Sqrt(StatPhiHe*StatPhiHe + SysPhiHe*SysPhiHe);
    LambdaPhiCs += TMath::Sqrt(StatPhiCs*StatPhiCs + SysPhiCs*SysPhiCs);
  } else if ( mode == 4 ) {
    LambdaPhiHe -= TMath::Sqrt(StatPhiHe*StatPhiHe + SysPhiHe*SysPhiHe);
    LambdaPhiCs -= TMath::Sqrt(StatPhiCs*StatPhiCs + SysPhiCs*SysPhiCs);
  } else if ( mode == 5 ) {
    LambdaThetaPhiHe += TMath::Sqrt(StatThetaPhiHe*StatThetaPhiHe + SysThetaPhiHe*SysThetaPhiHe);
    LambdaThetaPhiCs += TMath::Sqrt(StatThetaPhiCs*StatThetaPhiCs + SysThetaPhiCs*SysThetaPhiCs);
  } else if ( mode == 6 ) {
    LambdaThetaPhiHe -= TMath::Sqrt(StatThetaPhiHe*StatThetaPhiHe + SysThetaPhiHe*SysThetaPhiHe);
    LambdaThetaPhiCs -= TMath::Sqrt(StatThetaPhiCs*StatThetaPhiCs + SysThetaPhiCs*SysThetaPhiCs);
  }



  Double_t ModifiedLambdaThetaHe = (-1.)*LambdaThetaHe/(3.*TMath::Sqrt(3) + TMath::Sqrt(3)*LambdaThetaHe);
  Double_t ModifiedLambdaThetaCs = (-1.)*LambdaThetaCs/(3.*TMath::Sqrt(3) + TMath::Sqrt(3)*LambdaThetaCs);

  Double_t ModifiedStatLambdaThetaHe = TMath::Sqrt(3)*StatThetaHe/((3.+LambdaThetaHe)*(3.+LambdaThetaHe));
  Double_t ModifiedStatLambdaThetaCs = TMath::Sqrt(3)*StatThetaCs/((3.+LambdaThetaCs)*(3.+LambdaThetaCs));
  Double_t ModifiedSysLambdaThetaHe  = TMath::Sqrt(3)*SysThetaHe/((3.+LambdaThetaHe)*(3.+LambdaThetaHe));
  Double_t ModifiedSysLambdaThetaCs  = TMath::Sqrt(3)*SysThetaCs/((3.+LambdaThetaCs)*(3.+LambdaThetaCs));


  Double_t ModifiedLambdaThetaPhiHe = LambdaThetaPhiHe/(3. + LambdaThetaHe);
  Double_t ModifiedLambdaThetaPhiCs = LambdaThetaPhiCs/(3. + LambdaThetaCs);

  Double_t ModifiedStatLambdaThetaPhiHe = TMath::Sqrt(3)*StatThetaHe/((3.+LambdaThetaHe)*(3.+LambdaThetaHe));
  Double_t ModifiedStatLambdaThetaPhiCs = TMath::Sqrt(3)*StatThetaCs/((3.+LambdaThetaCs)*(3.+LambdaThetaCs));




  Double_t F_He = (6.*LambdaPhiHe+2.*LambdaThetaHe)/(3.*LambdaThetaHe+9.);
  Double_t F_Cs = (6.*LambdaPhiCs+2.*LambdaThetaCs)/(3.*LambdaThetaCs+9.);

  Double_t F_He_stat = 2.*F_He*TMath::Sqrt( StatPhiHe*StatPhiHe + (StatThetaHe*(1.-LambdaPhiHe)/(LambdaThetaHe+3.))*(StatThetaHe*(1.-LambdaPhiHe)/(LambdaThetaHe+3.)) );
  Double_t F_He_sys  = 2.*F_He*TMath::Sqrt( SysPhiHe *SysPhiHe  + (SysThetaHe *(1.-LambdaPhiHe)/(LambdaThetaHe+3.))*(SysThetaHe *(1.-LambdaPhiHe)/(LambdaThetaHe+3.)) );
  Double_t F_Cs_stat = 2.*F_Cs*TMath::Sqrt( StatPhiHe*StatPhiCs + (StatThetaCs*(1.-LambdaPhiCs)/(LambdaThetaCs+3.))*(StatThetaCs*(1.-LambdaPhiCs)/(LambdaThetaCs+3.)) );
  Double_t F_Cs_sys  = 2.*F_Cs*TMath::Sqrt( SysPhiHe *SysPhiCs  + (SysThetaCs *(1.-LambdaPhiCs)/(LambdaThetaCs+3.))*(SysThetaCs *(1.-LambdaPhiCs)/(LambdaThetaCs+3.)) );



  std::cout << "F_He = " << F_He << " +/- " << F_He_stat << " (stat.) +/- " << F_He_sys << " (sys.)" << '\n';
  std::cout << "F_Cs = " << F_Cs << " +/- " << F_Cs_stat << " (stat.) +/- " << F_Cs_sys << " (sys.)" << '\n';













  /* -
   * - QT matrix.
   */
  TMatrixD QuantumMatrix_He(3, 3);
  TMatrixD QuantumMatrix_Cs(3, 3);
  TArrayD  NullValues(9);
  for (Int_t i = 0; i < 9; i++) {
     NullValues[i] = 0;
  }
  QuantumMatrix_He.SetMatrixArray(NullValues.GetArray());
  QuantumMatrix_Cs.Print();

  TArrayD  qt_ij_He(9);
  TArrayD  qt_ij_Cs(9);

  qt_ij_He[0]  = (-6.*LambdaPhiHe+2.*LambdaThetaHe)/(3.*LambdaThetaHe+9.);
  qt_ij_He[1]  = 0.;
  qt_ij_He[2]  = (-2.)*LambdaThetaPhiHe/(LambdaThetaHe+3.);
  qt_ij_He[3]  = 0.;
  qt_ij_He[4]  = F_He;
  qt_ij_He[5]  = 0.;
  qt_ij_He[6]  = qt_ij_He[2];
  qt_ij_He[7]  = 0.;
  qt_ij_He[8]  = (-4.)*LambdaThetaHe/(3.*LambdaThetaHe+9.);


  qt_ij_Cs[0]  = (-6.*LambdaPhiCs+2.*LambdaThetaCs)/(3.*LambdaThetaCs+9.);
  qt_ij_Cs[1]  = 0.;
  qt_ij_Cs[2]  = (-2.)*LambdaThetaPhiCs/(LambdaThetaCs+3.);
  qt_ij_Cs[3]  = 0.;
  qt_ij_Cs[4]  = F_Cs;
  qt_ij_Cs[5]  = 0.;
  qt_ij_Cs[6]  = qt_ij_Cs[2];
  qt_ij_Cs[7]  = 0.;
  qt_ij_Cs[8]  = (-4.)*LambdaThetaCs/(3.*LambdaThetaCs+9.);

  QuantumMatrix_He.SetMatrixArray( qt_ij_He.GetArray() );
  QuantumMatrix_Cs.SetMatrixArray( qt_ij_Cs.GetArray() );
  QuantumMatrix_He.Print();
  QuantumMatrix_Cs.Print();





  TMatrixD RhoX_He(3, 3);
  TMatrixD RhoX_Cs(3, 3);
  TArrayD  Identity(9);
  for (Int_t i = 0; i < 9; i++) {
    if ( i == 0 || i == 4 || i == 8 ) {
      Identity[i] = 1./3.;
    } else {
      Identity[i] = 0.;
    }
  }
  RhoX_He.SetMatrixArray(Identity.GetArray());
  RhoX_Cs.SetMatrixArray(Identity.GetArray());
  RhoX_He.Print();
  RhoX_Cs.Print();



  RhoX_He += QuantumMatrix_He;
  RhoX_Cs += QuantumMatrix_Cs;
  RhoX_He.Print();
  RhoX_Cs.Print();

}
