#!/bin/bash

ROOTfile=$1

echo "OK0"
echo $ROOTfile
#aliroot -b -l -q fitZNC.cpp($ROOTfile)


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCosThetaHE.cpp
CreateCosThetaTh1("$ROOTfile");
EOF

echo "OK1"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCosThetaCS.cpp
CreateCosThetaTh1("$ROOTfile");
EOF

echo "OK2"


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationPhiHE.cpp
CreateCosThetaTh1("$ROOTfile");
EOF

echo "OK3"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationPhiCS.cpp
CreateCosThetaTh1("$ROOTfile");
EOF

echo "OK4"


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationTildePhiHE.cpp
CreateCosThetaTh1("$ROOTfile");
EOF

echo "OK5"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationTildePhiCS.cpp
CreateCosThetaTh1("$ROOTfile");
EOF

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

mv pngResults/CosThetaHe_*                        pngResults/$(date +%F)/CosThetaHE
mv pngResults/CosThetaHeFrame.root                pngResults/$(date +%F)/CosThetaHE
mv pngResults/CosThetaCs_*                        pngResults/$(date +%F)/CosThetaCS
mv pngResults/CosThetaCsFrame.root                pngResults/$(date +%F)/CosThetaCS
mv pngResults/PhiHe_*                             pngResults/$(date +%F)/PhiHE
mv pngResults/PhiHeFrame.root                     pngResults/$(date +%F)/PhiHE
mv pngResults/PhiCs_*                             pngResults/$(date +%F)/PhiCS
mv pngResults/PhiCsFrame.root                     pngResults/$(date +%F)/PhiCS
mv pngResults/TildePhiHe_*                        pngResults/$(date +%F)/TildePhiHE
mv pngResults/TildePhiHeFrame.root                pngResults/$(date +%F)/TildePhiHE
mv pngResults/TildePhiCs_*                        pngResults/$(date +%F)/TildePhiCS
mv pngResults/TildePhiCsFrame.root                pngResults/$(date +%F)/TildePhiCS


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedHe1D.cpp
PolarisationCorrectedHe1D();
EOF

echo "OK5"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedCs1D.cpp
PolarisationCorrectedCs1D();
EOF

echo "OK6"


mv pngResults/PolarisationCorrectedHe1D.root      pngResults/$(date +%F)/1Dresults
mv pngResults/PolarisationCorrectedCs1D.root      pngResults/$(date +%F)/1Dresults


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/Polarisation2DHE.cpp
CreateCosThetaTh2("$ROOTfile");
EOF

echo "OK7"

aliroot -b -l <<EOF
.L fitRootConverted/Polarisation2DCS.cpp
CreateCosThetaTh2("$ROOTfile");
EOF

echo "OK8"


mv pngResults/2DHe_*                              pngResults/$(date +%F)/2DHE
mv pngResults/Polarisation2DHE.root               pngResults/$(date +%F)/2DHE
mv pngResults/2DCs_*                              pngResults/$(date +%F)/2DCS
mv pngResults/Polarisation2DCs.root               pngResults/$(date +%F)/2DCS


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribHe2D.cpp
PolarisationCorrectedDistribHe2D();
EOF

echo "OK7"

aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCorrectedDistribCs2D.cpp
PolarisationCorrectedDistribHe2D();
EOF

echo "OK8"


mv pngResults/PolarisationCorrectedHe2D.root      pngResults/$(date +%F)/2DHE
mv pngResults/PolarisationCorrectedCs2D.root      pngResults/$(date +%F)/2DCS
