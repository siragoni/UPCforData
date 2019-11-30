#!/bin/bash

ROOTfile=$1

echo "OK0"
echo $ROOTfile
#aliroot -b -l -q fitZNC.cpp($ROOTfile)


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCosThetaHE.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK1"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCosThetaCS.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK2"


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationPhiHE.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK3"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationPhiCS.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK4"


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationTildePhiHE.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK5"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationTildePhiCS.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK6"




if [ -d "pngResults/$(date +%F)" ]; then rm -rf pngResults/$(date +%F); fi
mkdir -p pngResults/$(date +%F)
mkdir -p pngResults/$(date +%F)/CosThetaHE
mkdir -p pngResults/$(date +%F)/CosThetaCS
mkdir -p pngResults/$(date +%F)/PhiHE
mkdir -p pngResults/$(date +%F)/PhiCS
mkdir -p pngResults/$(date +%F)/TildePhiHE
mkdir -p pngResults/$(date +%F)/TildePhiCS
mkdir -p pngResults/$(date +%F)/1Dresults
mkdir -p pngResults/$(date +%F)/2DHE
mkdir -p pngResults/$(date +%F)/2DCS
mkdir -p pngResults/$(date +%F)/2Dresults


mv pngResults/CosThetaHe_*                         pngResults/$(date +%F)/CosThetaHE
mv pngResults/CosThetaHeFrame*.root                pngResults/$(date +%F)/CosThetaHE
mv pngResults/CosThetaCs_*                         pngResults/$(date +%F)/CosThetaCS
mv pngResults/CosThetaCsFrame*.root                pngResults/$(date +%F)/CosThetaCS
mv pngResults/PhiHe_*                              pngResults/$(date +%F)/PhiHE
mv pngResults/PhiHeFrame*.root                     pngResults/$(date +%F)/PhiHE
mv pngResults/PhiCs_*                              pngResults/$(date +%F)/PhiCS
mv pngResults/PhiCsFrame*.root                     pngResults/$(date +%F)/PhiCS
mv pngResults/TildePhiHe_*                         pngResults/$(date +%F)/TildePhiHE
mv pngResults/TildePhiHeFrame*.root                pngResults/$(date +%F)/TildePhiHE
mv pngResults/TildePhiCs_*                         pngResults/$(date +%F)/TildePhiCS
mv pngResults/TildePhiCsFrame*.root                pngResults/$(date +%F)/TildePhiCS


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribHe1D.cpp
PolarisationCorrectedDistribHe1D($value);
EOF
done
echo "OK5"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribCs1D.cpp
PolarisationCorrectedDistribCs1D($value);
EOF
done
echo "OK6"


mv pngResults/PolarisationCorrectedHe1D*.root      pngResults/$(date +%F)/1Dresults
mv pngResults/PolarisationCorrectedCs1D*.root      pngResults/$(date +%F)/1Dresults


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/Polarisation2DHE.cpp
CreateCosThetaTh2("$ROOTfile", $value);
EOF
done
echo "OK7"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/Polarisation2DCS.cpp
CreateCosThetaTh2("$ROOTfile", $value);
EOF
done
echo "OK8"


mv pngResults/2DHe_*                              pngResults/$(date +%F)/2DHE
mv pngResults/Polarisation2DHE*.root               pngResults/$(date +%F)/2DHE
mv pngResults/2DCs_*                              pngResults/$(date +%F)/2DCS
mv pngResults/Polarisation2DCs*.root               pngResults/$(date +%F)/2DCS


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribHe2D.cpp
PolarisationCorrectedDistribHe2D($value);
EOF
done
echo "OK7"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribCs2D.cpp
PolarisationCorrectedDistribCs2D($value);
EOF
done
echo "OK8"


mv pngResults/PolarisationCorrectedHe2D*.root      pngResults/$(date +%F)/2DHE
mv pngResults/PolarisationCorrectedCs2D*.root      pngResults/$(date +%F)/2DCS


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeFitWithRoot2D.cpp
PolarisationHeFitWithRoot2D();
EOF

echo "OK9"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCsFitWithRoot2D.cpp
PolarisationCsFitWithRoot2D();
EOF

echo "OK9"

mv pngResults/PolFitWithRootHe2D.png              pngResults/$(date +%F)/2Dresults
mv pngResults/PolFitWithRootCs2D.png              pngResults/$(date +%F)/2Dresults


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
for value2 in {0..2}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeMinuit1D.cpp
PolarisationHeMinuit1D($value, $value2);
EOF
done
done
echo "OK10"

for value in {0..5}
do
for value2 in {0..2}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCsMinuit1D.cpp
PolarisationCsMinuit1D($value, $value2);
EOF
done
done
echo "OK11"

mv pngResults/*Minuit.png                         pngResults/$(date +%F)/1Dresults
mv pngResults/Parameters_SigEx_*                  pngResults/$(date +%F)/1Dresults


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeMinuit2D.cpp
PolarisationHeMinuit2D();
EOF

echo "OK10"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCsMinuit2D.cpp
PolarisationCsMinuit2D();
EOF

echo "OK11"

mv pngResults/*Minuit.png                         pngResults/$(date +%F)/2Dresults






#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCosThetaHE_0N0N.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK1"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCosThetaCS_0N0N.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK2"


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationPhiHE_0N0N.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK3"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationPhiCS_0N0N.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK4"


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationTildePhiHE_0N0N.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK5"

for value in {0..5}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationTildePhiCS_0N0N.cpp
CreateCosThetaTh1("$ROOTfile", $value);
EOF
done
echo "OK6"




mkdir -p pngResults/$(date +%F)/CosThetaHE_0N0N
mkdir -p pngResults/$(date +%F)/CosThetaCS_0N0N
mkdir -p pngResults/$(date +%F)/PhiHE_0N0N
mkdir -p pngResults/$(date +%F)/PhiCS_0N0N
mkdir -p pngResults/$(date +%F)/TildePhiHE_0N0N
mkdir -p pngResults/$(date +%F)/TildePhiCS_0N0N
mkdir -p pngResults/$(date +%F)/1Dresults_0N0N
mkdir -p pngResults/$(date +%F)/2DHE_0N0N
mkdir -p pngResults/$(date +%F)/2DCS_0N0N
mkdir -p pngResults/$(date +%F)/2Dresults_0N0N


mv pngResults/CosThetaHe_*                         pngResults/$(date +%F)/CosThetaHE_0N0N
mv pngResults/CosThetaHeFrame_0N0N*.root           pngResults/$(date +%F)/CosThetaHE_0N0N
mv pngResults/CosThetaCs_*                         pngResults/$(date +%F)/CosThetaCS_0N0N
mv pngResults/CosThetaCsFrame_0N0N*.root           pngResults/$(date +%F)/CosThetaCS_0N0N
mv pngResults/PhiHe_*                              pngResults/$(date +%F)/PhiHE_0N0N
mv pngResults/PhiHeFrame_0N0N*.root                pngResults/$(date +%F)/PhiHE_0N0N
mv pngResults/PhiCs_*                              pngResults/$(date +%F)/PhiCS_0N0N
mv pngResults/PhiCsFrame_0N0N*.root                pngResults/$(date +%F)/PhiCS_0N0N
mv pngResults/TildePhiHe_*                         pngResults/$(date +%F)/TildePhiHE_0N0N
mv pngResults/TildePhiHeFrame_0N0N*.root           pngResults/$(date +%F)/TildePhiHE_0N0N
mv pngResults/TildePhiCs_*                         pngResults/$(date +%F)/TildePhiCS_0N0N
mv pngResults/TildePhiCsFrame_0N0N*.root           pngResults/$(date +%F)/TildePhiCS_0N0N


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribHe1D_0N0N.cpp
PolarisationCorrectedDistribHe1D_0N0N();
EOF

echo "OK5"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribCs1D_0N0N.cpp
PolarisationCorrectedDistribCs1D_0N0N();
EOF

echo "OK6"


mv pngResults/PolarisationCorrectedHe1D_0N0N.root      pngResults/$(date +%F)/1Dresults_0N0N
mv pngResults/PolarisationCorrectedCs1D_0N0N.root      pngResults/$(date +%F)/1Dresults_0N0N


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/Polarisation2DHE_0N0N.cpp
CreateCosThetaTh2("$ROOTfile");
EOF

echo "OK7"

aliroot -b -l <<EOF
.L fitRootConverted/Polarisation2DCS_0N0N.cpp
CreateCosThetaTh2("$ROOTfile");
EOF

echo "OK8"


mv pngResults/2DHe_*                              pngResults/$(date +%F)/2DHE_0N0N
mv pngResults/Polarisation2DHE_0N0N.root          pngResults/$(date +%F)/2DHE_0N0N
mv pngResults/2DCs_*                              pngResults/$(date +%F)/2DCS_0N0N
mv pngResults/Polarisation2DCs_0N0N.root          pngResults/$(date +%F)/2DCS_0N0N


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribHe2D_0N0N.cpp
PolarisationCorrectedDistribHe2D_0N0N();
EOF

echo "OK7"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribCs2D_0N0N.cpp
PolarisationCorrectedDistribCs2D_0N0N();
EOF

echo "OK8"


mv pngResults/PolarisationCorrectedHe2D_0N0N.root      pngResults/$(date +%F)/2DHE_0N0N
mv pngResults/PolarisationCorrectedCs2D_0N0N.root      pngResults/$(date +%F)/2DCS_0N0N


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeFitWithRoot2D_0N0N.cpp
PolarisationHeFitWithRoot2D_0N0N();
EOF

echo "OK9"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCsFitWithRoot2D_0N0N.cpp
PolarisationCsFitWithRoot2D_0N0N();
EOF

echo "OK9"

mv pngResults/PolFitWithRootHe2D_0N0N.png         pngResults/$(date +%F)/2Dresults_0N0N
mv pngResults/PolFitWithRootCs2D_0N0N.png         pngResults/$(date +%F)/2Dresults_0N0N


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeMinuit1D_0N0N.cpp
PolarisationHeMinuit1D_0N0N();
EOF

echo "OK10"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCsMinuit1D_0N0N.cpp
PolarisationCsMinuit1D_0N0N();
EOF

echo "OK11"

mv pngResults/*Minuit_0N0N.png                    pngResults/$(date +%F)/1Dresults_0N0N


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeMinuit2D_0N0N.cpp
PolarisationHeMinuit2D_0N0N();
EOF

echo "OK10"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCsMinuit2D_0N0N.cpp
PolarisationCsMinuit2D_0N0N();
EOF

echo "OK11"

mv pngResults/*Minuit_0N0N.png                    pngResults/$(date +%F)/2Dresults_0N0N
