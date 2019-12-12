#!/bin/bash

ROOTfile=$1

echo "OK0"
echo $ROOTfile


for value in {1..7}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationTriggerSys.cpp
PolarisationTriggerSys("$ROOTfile", $value);
EOF
done
echo "OK1"

. polarisation_small.sh Polarisation_1.root 1
. polarisation_small.sh Polarisation_2.root 2
. polarisation_small.sh Polarisation_3.root 3
. polarisation_small.sh Polarisation_4.root 4
. polarisation_small.sh Polarisation_5.root 5
. polarisation_small.sh Polarisation_6.root 6
. polarisation_small.sh Polarisation_7.root 7
