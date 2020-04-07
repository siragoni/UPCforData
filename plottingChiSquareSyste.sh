#!/bin/bash

ROOTfile=$1

echo "OK0"
echo $ROOTfile
#aliroot -b -l -q fitZNC.cpp($ROOTfile)




#before I had always used aliroot -b -l -q , but it doesn't seem to work...
for value in {0..5}
do
for value2 in {0..2}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeMinuit1D_chi2.cpp
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
.L fitRootConverted/PolarisationCsMinuit1D_chi2.cpp
PolarisationCsMinuit1D($value, $value2);
EOF
done
done
echo "OK11"

mv pngResults/*Minuit.png                         pngResults/$(date +%F)/1Dresults
mv pngResults/*.png                               pngResults/$(date +%F)/1Dresults
mv pngResults/Parameters_SigEx_*                  pngResults/$(date +%F)/1Dresults
