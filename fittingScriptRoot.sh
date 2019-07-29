#!/bin/bash

ROOTfile=$1
MCfile=$2

echo "OK0"
echo $ROOTfile
echo $MCfile
#aliroot -b -l -q fitZNC.cpp($ROOTfile)


#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/fitZNCkaedenRoofitPlot.cpp
fitZNC("$ROOTfile");
EOF

echo "OK1"

aliroot -b -l <<EOF
.L fitRootConverted/fitZNAkaedenRoofitPlot.cpp
fitZNA("$ROOTfile");
EOF

echo "OK2"

# aliroot -b -l <<EOF
# .L fitRootConverted/fitPtDistr.cpp
# fitPtDistr("$ROOTfile");
# EOF
#
# echo "OK3"

for value in {0..4}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitPtDistr0.cpp
fitPtDistr("$ROOTfile",$value);
EOF
done

echo "OK3bis"

aliroot -b -l <<EOF
.L fitRootConverted/fitBinMigrationWithRootSimple.cpp
fitBinMigration("$MCfile");
EOF

echo "OK4"

var1=$(aliroot -b -l <<EOF
.L fitRootConverted/fitSuppressedNjpsi.cpp
fitROOFit("$ROOTfile", 0);
EOF)

echo "OK5"

var2=$(aliroot -b -l <<EOF
.L fitRootConverted/fitSuppressedNjpsi.cpp
fitROOFit("$ROOTfile", 1);
EOF)

echo "OK6"

var3=$(aliroot -b -l <<EOF
.L fitRootConverted/fitSuppressedNjpsi.cpp
fitROOFit("$ROOTfile", 2);
EOF)

echo "$var1" > outputSuppressedSelFlag0.txt
echo "$var2" > outputSuppressedSelFlag1.txt
echo "$var3" > outputSuppressedSelFlag2.txt

echo "OK7"

var=$(aliroot -b -l <<EOF
.L fitRootConverted/readHisto.cpp
readHisto("$ROOTfile");
EOF)

echo "$var" > output.txt

echo "OK8"

aliroot -b -l <<EOF
.L fitRootConverted/fitHelicityShapeWithRootSimple.cpp
fitHelicityShape("$ROOTfile","$MCfile");
EOF

echo "OK9"

aliroot -b -l <<EOF
.L fitRootConverted/fitCosThetaHeOVB.cpp
fitHelicityShape("$ROOTfile","$MCfile");
EOF

echo "OK10"

aliroot -b -l <<EOF
.L fitRootConverted/fitPhiHeOVB.cpp
fitHelicityShape("$ROOTfile","$MCfile");
EOF

echo "OK11"

for value in {0..7}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitCosThetaHeOVBrapidity.cpp
fitHelicityShape("$ROOTfile","$MCfile",$value);
EOF
done

echo "OK12"

for value in {0..7}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitPhiHeOVBrapidity.cpp
fitHelicityShape("$ROOTfile","$MCfile",$value);
EOF
done

echo "OK13"

var4=$(aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiTemplate.cpp
fitJPsiTemplate("$ROOTfile", 0);
EOF)

echo "OK14"

var5=$(aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiTemplate.cpp
fitJPsiTemplate("$ROOTfile", 1);
EOF)

echo "OK15"

var6=$(aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiTemplate.cpp
fitJPsiTemplate("$ROOTfile", 2);
EOF)

echo "$var4" > outputJPsiTemplate0.txt
echo "$var5" > outputJPsiTemplate1.txt
echo "$var6" > outputJPsiTemplate2.txt

var4=$(aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiCrystalBallV2.cpp
fitJPsiTemplate("$ROOTfile", 0);
EOF)

echo "OK14"

var5=$(aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiCrystalBallV2.cpp
fitJPsiTemplate("$ROOTfile", 1);
EOF)

echo "OK15"

var6=$(aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiCrystalBallV2.cpp
fitJPsiTemplate("$ROOTfile", 2);
EOF)

echo "$var4" > outputJPsiTemplate0.txt
echo "$var5" > outputJPsiTemplate1.txt
echo "$var6" > outputJPsiTemplate2.txt


for value in {0..7}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitXNXNTemplate.cpp
fitJPsiTemplate("$ROOTfile", $value);
EOF
done

for value in {0..9}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiCosThetaBinsTemplate.cpp
fitJPsiTemplate("$ROOTfile", $value);
EOF
done

for value in {0..9}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitCosThetaHeOVBrapidity10bins.cpp
fitHelicityShape("$ROOTfile","$MCfile",$value);
EOF
done

echo "OK16"

for value in {0..9}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitPhiHeOVBrapidity10bins.cpp
fitHelicityShape("$ROOTfile","$MCfile",$value);
EOF
done

echo "OK17"

aliroot -b -l <<EOF
.L fitRootConverted/DrawTMultiGraph.cpp
DrawTMultiGraph();
EOF

echo "OK18"

aliroot -b -l <<EOF
.L fitRootConverted/fitHelicityHeFrame2D.cpp
fitHelicityShape("$ROOTfile", "AnalysisResultsMC18042019.root");
EOF

echo "OK19"

aliroot -b -l <<EOF
.L fitRootConverted/fitHelicityHeFrame2DinclusiveBinning.cpp
fitHelicityShape("$ROOTfile", "AnalysisResultsMC18042019.root");
EOF

echo "OK19"

aliroot -b -l <<EOF
.L fitRootConverted/fitHelicityHeFrame2DonlyGeneratedLevel.cpp
fitHelicityShape("$ROOTfile", "AnalysisResultsMC18042019.root");
EOF

echo "OK20"

aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiTemplateSignalExtractionV2.cpp
CreateTH2("$ROOTfile");
EOF

echo "OK21"



echo "These are the results   of: $ROOTfile" > ReadMeAnalysisResults.md
# echo "with the MC simulations of: $MCfile  " > ReadMeAnalysisResults.md
sort -n output.txt > outputSorted.txt

rm output.txt

if [ -d "pngResults/$(date +%F)" ]; then rm -rf pngResults/$(date +%F); fi
mkdir -p pngResults/$(date +%F)
mkdir -p pngResults/$(date +%F)/HelicityFrame
mkdir -p pngResults/$(date +%F)/CollinsSoperFrame
mkdir -p pngResults/$(date +%F)/XNXNandPTdistr
mkdir -p pngResults/$(date +%F)/FitInvMass
mkdir -p pngResults/$(date +%F)/JPsiCosThetaBins
mkdir -p pngResults/$(date +%F)/HelicityFrame10bins
mkdir -p pngResults/$(date +%F)/Fit2Dpolarisation
mkdir -p pngResults/$(date +%F)/SignalExtraction

mv ReadMeAnalysisResults.md       pngResults/$(date +%F)
mv outputSorted.txt               pngResults/$(date +%F)
mv outputSuppressedSelFlag0.txt   pngResults/$(date +%F)
mv outputSuppressedSelFlag1.txt   pngResults/$(date +%F)
mv outputSuppressedSelFlag2.txt   pngResults/$(date +%F)
mv outputJPsiTemplate0.txt        pngResults/$(date +%F)
mv outputJPsiTemplate1.txt        pngResults/$(date +%F)
mv outputJPsiTemplate2.txt        pngResults/$(date +%F)
mv test*.root                     pngResults/$(date +%F)
mv pngResults/CosThetaHeO*.png    pngResults/$(date +%F)/HelicityFrame
mv pngResults/PhiHeO*.png         pngResults/$(date +%F)/HelicityFrame
mv pngResults/CosThetaHe10Bi*.png pngResults/$(date +%F)/HelicityFrame10bins
mv pngResults/PhiHe10Bi*.png      pngResults/$(date +%F)/HelicityFrame10bins
# mv pngResults/CosThetaHeO*.png  pngResults/$(date +%F)/CollinsSoperFrame
# mv pngResults/PhiHeO*.png       pngResults/$(date +%F)/CollinsSoperFrame
mv pngResults/Coh0N*.png          pngResults/$(date +%F)/XNXNandPTdistr
mv pngResults/CohXN*.png          pngResults/$(date +%F)/XNXNandPTdistr
mv pngResults/Incoh0N*.png        pngResults/$(date +%F)/XNXNandPTdistr
mv pngResults/IncohXN*.png        pngResults/$(date +%F)/XNXNandPTdistr
mv pngResults/fitPtDistrALL.png   pngResults/$(date +%F)/XNXNandPTdistr
mv pngResults/fitPtDistr0*.png    pngResults/$(date +%F)/XNXNandPTdistr
mv pngResults/fitPtDistrX*.png    pngResults/$(date +%F)/XNXNandPTdistr
mv pngResults/Coherent*.png       pngResults/$(date +%F)/FitInvMass
mv pngResults/PtInt*.png          pngResults/$(date +%F)/FitInvMass
mv pngResults/Incoherent*.png     pngResults/$(date +%F)/FitInvMass
mv pngResults/suppressed*.png     pngResults/$(date +%F)/FitInvMass
mv pngResults/JPsiCosThetaB*.png  pngResults/$(date +%F)/JPsiCosThetaBins
mv pngResults/fit2D*.png          pngResults/$(date +%F)/Fit2Dpolarisation
mv pngResults/fitSignalExtra*.png pngResults/$(date +%F)/SignalExtraction
mv pngResults/TH2*.root           pngResults/$(date +%F)/SignalExtraction
mv pngResults/*.png               pngResults/$(date +%F)
