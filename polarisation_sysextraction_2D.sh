#!/bin/bash


echo "OK0"
#aliroot -b -l -q fitZNC.cpp($ROOTfile)


for value in {0..5}
do
for value2 in {0..2}
do
#before I had always used aliroot -b -l -q , but it doesn't seem to work...
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationHeMinuit2D.cpp
PolarisationHeMinuit2D($value, $value2);
EOF
done
done

echo "OK10"

for value in {0..5}
do
for value2 in {0..2}
do
aliroot -b -l <<EOF
.L fitRootConverted/PolarisationCsMinuit2D.cpp
PolarisationCsMinuit2D($value, $value2);
EOF
done
done


echo "OK11"
