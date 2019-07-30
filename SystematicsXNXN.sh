#!/bin/bash

ROOTfile=$1

echo "OK0"
echo $ROOTfile
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

for value in {1..10}
do
for value2 in {0..4}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiCrystalBallSystematics.cpp
fitJPsiTemplate("$ROOTfile",$value,$value2);
EOF
rm pngResults/Systematics.root
mv pngResults/Systematics2.root pngResults/Systematics.root
done
done

for value in {1..10}
do
for value2 in {0..4}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiCrystalBallSystematicsHalfBin.cpp
fitJPsiTemplate("$ROOTfile",$value,$value2);
EOF
rm pngResults/SystematicsHalfBin.root
mv pngResults/SystematicsHalfBin2.root pngResults/SystematicsHalfBin.root
done
done

for value in {1..10}
do
for value2 in {0..4}
do
aliroot -b -l <<EOF
.L fitRootConverted/fitJPsiCrystalBallSystematicsHalfHalfBin.cpp
fitJPsiTemplate("$ROOTfile",$value,$value2);
EOF
rm pngResults/SystematicsHalfHalfBin.root
mv pngResults/SystematicsHalfHalfBin2.root pngResults/SystematicsHalfHalfBin.root
done
done

# for value in {0..4}
# do
# aliroot -b -l <<EOF
# .L fitRootConverted/fitPtDistr0.cpp
# fitPtDistr("$ROOTfile",$value);
# EOF
# done


# for value in {0..7}
# do
# aliroot -b -l <<EOF
# .L fitRootConverted/fitXNXNTemplate.cpp
# fitJPsiTemplate("$ROOTfile", $value);
# EOF
# done


if [ -d "pngResults/$(date +%F)" ]; then rm -rf pngResults/$(date +%F); fi
mkdir -p pngResults/$(date +%F)
mkdir -p pngResults/$(date +%F)/XNXN
mkdir -p pngResults/$(date +%F)/XNXNhalfBin
mkdir -p pngResults/$(date +%F)/XNXNhalfHalfBin
mkdir -p pngResults/$(date +%F)/PTdistr
mkdir -p pngResults/$(date +%F)/FitInvMass
mkdir -p pngResults/$(date +%F)/FitInvMassHalfBin
mkdir -p pngResults/$(date +%F)/FitInvMassHalfHalfBin

mv pngResults/InvMassSystematics_1_*              pngResults/$(date +%F)/FitInvMass
mv pngResults/InvMassSystematics_2_*              pngResults/$(date +%F)/FitInvMass
mv pngResults/InvMassSystematics_3_*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematics_4_*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematics_5_*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematics_6_*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematics_7_*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematics_8_*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematics_9_*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematics_10*              pngResults/$(date +%F)/XNXN
mv pngResults/InvMassSystematicsHalfBin_1_*       pngResults/$(date +%F)/FitInvMassHalfBin
mv pngResults/InvMassSystematicsHalfBin_2_*       pngResults/$(date +%F)/FitInvMassHalfBin
mv pngResults/InvMassSystematicsHalfBin_3_*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfBin_4_*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfBin_5_*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfBin_6_*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfBin_7_*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfBin_8_*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfBin_9_*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfBin_10*       pngResults/$(date +%F)/XNXNhalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_1_*   pngResults/$(date +%F)/FitInvMassHalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_2_*   pngResults/$(date +%F)/FitInvMassHalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_3_*   pngResults/$(date +%F)/XNXNhalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_4_*   pngResults/$(date +%F)/XNXNhalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_5_*   pngResults/$(date +%F)/XNXNhalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_6_*   pngResults/$(date +%F)/XNXNhalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_7_*   pngResults/$(date +%F)/XNXNhalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_8_*   pngResults/$(date +%F)/XNXNhalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_9_*   pngResults/$(date +%F)/XNXNhalfHalfBin
mv pngResults/InvMassSystematicsHalfHalfBin_10*   pngResults/$(date +%F)/XNXNhalfHalfBin
