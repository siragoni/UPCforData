//_____________________________________________________________________________
/* - The following are code snippets adapted from the AliAODDimuon class.
   - The problem is that that class was adapted specifically for the
   - inclusive people's analysis, hence it is not fit for the UPC...
   -
 */
Double_t AliAnalysisTaskUPCforward::CosThetaCollinsSoper( TLorentzVector muonPositive,
                                                          TLorentzVector muonNegative,
                                                          TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  /* - Determine the CS angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaCS = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaCS;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforward::CosPhiCollinsSoper( TLorentzVector muonPositive,
                                                        TLorentzVector muonNegative,
                                                        TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper PHI for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  TVector3 yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
  TVector3 xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();

  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
  return   phi;
}
//_____________________________________________________________________________
/* -
   - Quantum Tomography for CS.
   -
 */
Double_t AliAnalysisTaskUPCforwardMC::CosThetaQuantumTomCS( TLorentzVector muonPositive,
                                                            TLorentzVector muonNegative,
                                                            TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  TLorentzVector VectorMeson = possibleJPsi;
  // cout << "Vector meson Px: " << VectorMeson.Px() << "    X: " << VectorMeson.X() << endl;
  // cout << "Vector meson Py: " << VectorMeson.Py() << "    Y: " << VectorMeson.Y() << endl;
  // cout << "Vector meson Pz: " << VectorMeson.Pz() << "    Z: " << VectorMeson.Z() << endl;
  // cout << "pProjCM Px: " << pProjCM.Px() << "    X: " << pProjCM.X() << endl;
  // cout << "pProjCM Py: " << pProjCM.Py() << "    Y: " << pProjCM.Y() << endl;
  // cout << "pProjCM Pz: " << pProjCM.Pz() << "    Z: " << pProjCM.Z() << endl;
  // cout << "pTargCM Px: " << pTargCM.Px() << "    X: " << pTargCM.X() << endl;
  // cout << "pTargCM Py: " << pTargCM.Py() << "    Y: " << pTargCM.Y() << endl;
  // cout << "pTargCM Pz: " << pTargCM.Pz() << "    Z: " << pTargCM.Z() << endl;

  // TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) ).Unit();
  TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) );
  // QuantumZ *= ( 1./QuantumZ.Mag2() );

  Double_t       xQuantumY =  pProjCM.Pz()*pTargCM.E()*VectorMeson.Py()-pProjCM.E()*pTargCM.Pz()*VectorMeson.Py();
  Double_t       yQuantumY = -pProjCM.Pz()*pTargCM.E()*VectorMeson.Px()+pProjCM.E()*pTargCM.Pz()*VectorMeson.Px();
  TLorentzVector QuantumYnotnormalised(xQuantumY, yQuantumY, 0., 0.);
  // TLorentzVector QuantumY = QuantumYnotnormalised.Unit();
  TLorentzVector QuantumY = QuantumYnotnormalised;
  // QuantumY *= ( 1./QuantumY.Mag2() );
  // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) ).Unit();
  // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) );
  // QuantumX *= ( 1./QuantumX.Mag2() );
  TLorentzVector QuantumX = VectorMeson - pProjCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pProjCM)) )) - pTargCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pTargCM)) ));
  // QuantumX *= ( 1./QuantumX.Mag2() );



  TLorentzVector DifferenceMuons   = muonPositive - muonNegative;
  Double_t       xFinalVector      = QuantumX.Dot(DifferenceMuons);
  Double_t       yFinalVector      = QuantumY.Dot(DifferenceMuons);
  Double_t       zFinalVector      = QuantumZ.Dot(DifferenceMuons);
  TVector3       FinalVector( xFinalVector, yFinalVector, zFinalVector );

  Double_t       CosThetaQuantumCS = ( FinalVector.Unit() ).Z();
  Double_t       SinThetaQuantumCS = TMath::Sqrt( 1 - CosThetaQuantumCS*CosThetaQuantumCS );
  Double_t       PhiQuantumCS      = TMath::ASin( ( (FinalVector.Unit()).Y() )/ SinThetaQuantumCS );


  return CosThetaQuantumCS;

}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUPCforwardMC::PhiQuantumTomogrCS( TLorentzVector muonPositive,
                                                          TLorentzVector muonNegative,
                                                          TLorentzVector possibleJPsi )
{
  /*    - Quantum Tomography for CS.
        -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  TLorentzVector VectorMeson = possibleJPsi;

  // TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) ).Unit();
  TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) );
  // QuantumZ *= ( 1./QuantumZ.Mag2() );

  Double_t       xQuantumY = -pProjCM.Pz()*pTargCM.E()*VectorMeson.Py()+pProjCM.E()*pTargCM.Pz()*VectorMeson.Py();
  Double_t       yQuantumY =  pProjCM.Pz()*pTargCM.E()*VectorMeson.Px()+pProjCM.E()*pTargCM.Pz()*VectorMeson.Px();
  TLorentzVector QuantumYnotnormalised(xQuantumY, yQuantumY, 0., 0.);
  // TLorentzVector QuantumY = QuantumYnotnormalised.Unit();
  TLorentzVector QuantumY = QuantumYnotnormalised;
  // QuantumY *= ( 1./QuantumY.Mag2() );
  // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) ).Unit();
  // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) );
  // QuantumX *= ( 1./QuantumX.Mag2() );
  TLorentzVector QuantumX = VectorMeson - pProjCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pProjCM)) )) - pTargCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pTargCM)) ));
  // QuantumX *= ( 1./QuantumX.Mag2() );



  TLorentzVector DifferenceMuons   = muonPositive - muonNegative;
  Double_t       xFinalVector      = QuantumX.Dot(DifferenceMuons);
  Double_t       yFinalVector      = QuantumY.Dot(DifferenceMuons);
  Double_t       zFinalVector      = QuantumZ.Dot(DifferenceMuons);
  TVector3       FinalVector( xFinalVector, yFinalVector, zFinalVector );

  Double_t       CosThetaQuantumCS = ( FinalVector.Unit() ).Z();
  Double_t       SinThetaQuantumCS = TMath::Sqrt( 1 - CosThetaQuantumCS*CosThetaQuantumCS );
  Double_t       PhiQuantumCS      = TMath::ASin( ( (FinalVector.Unit()).Y() )/ SinThetaQuantumCS );


  return PhiQuantumCS;

}
//_____________________________________________________________________________
