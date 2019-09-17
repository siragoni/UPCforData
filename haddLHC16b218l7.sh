#!/usr/bin/env bash




if [ -d "MCtrainResults/$(date +%F)" ]; then rm -rf MCtrainResults/$(date +%F); fi
mkdir -p MCtrainResults/$(date +%F)

LHC16b2=$(date +%F)
LHC16b2+=_LHC16b2
LHC18l7=$(date +%F)
LHC18l7+=_LHC18l7

echo $kCohJpsiToMu
mkdir -p MCtrainResults/$(date +%F)/kCohJpsiToMu
hadd     MCtrainResults/$(date +%F)/kCohJpsiToMu/AnalysisResults.root        MCtrainResults/$LHC16b2/kCohJpsiToMu/AnalysisResults.root        MCtrainResults/$LHC18l7/kCohJpsiToMu/AnalysisResults.root

echo $kIncohJpsiToMu
mkdir -p MCtrainResults/$(date +%F)/kIncohJpsiToMu
hadd     MCtrainResults/$(date +%F)/kIncohJpsiToMu/AnalysisResults.root      MCtrainResults/$LHC16b2/kIncohJpsiToMu/AnalysisResults.root      MCtrainResults/$LHC18l7/kIncohJpsiToMu/AnalysisResults.root

echo $kTwoGammaToMuMedium
mkdir -p MCtrainResults/$(date +%F)/kTwoGammaToMuMedium
hadd     MCtrainResults/$(date +%F)/kTwoGammaToMuMedium/AnalysisResults.root MCtrainResults/$LHC16b2/kTwoGammaToMuMedium/AnalysisResults.root MCtrainResults/$LHC18l7/kTwoGammaToMuMedium/AnalysisResults.root

echo $kTwoGammaToMuHigh
mkdir -p MCtrainResults/$(date +%F)/kTwoGammaToMuHigh
hadd     MCtrainResults/$(date +%F)/kTwoGammaToMuHigh/AnalysisResults.root   MCtrainResults/$LHC16b2/kTwoGammaToMuHigh/AnalysisResults.root   MCtrainResults/$LHC18l7/kTwoGammaToMuHigh/AnalysisResults.root

echo $kCohPsi2sToMu
mkdir -p MCtrainResults/$(date +%F)/kCohPsi2sToMu
hadd     MCtrainResults/$(date +%F)/kCohPsi2sToMu/AnalysisResults.root       MCtrainResults/$LHC16b2/kCohPsi2sToMu/AnalysisResults.root       MCtrainResults/$LHC18l7/kCohPsi2sToMu/AnalysisResults.root

echo $kCohPsi2sToMuPi
mkdir -p MCtrainResults/$(date +%F)/kCohPsi2sToMuPi
hadd     MCtrainResults/$(date +%F)/kCohPsi2sToMuPi/AnalysisResults.root     MCtrainResults/$LHC16b2/kCohPsi2sToMuPi/AnalysisResults.root     MCtrainResults/$LHC18l7/kCohPsi2sToMuPi/AnalysisResults.root

echo $kIncohPsi2sToMu
mkdir -p MCtrainResults/$(date +%F)/kIncohPsi2sToMu
hadd     MCtrainResults/$(date +%F)/kIncohPsi2sToMu/AnalysisResults.root     MCtrainResults/$LHC16b2/kIncohPsi2sToMu/AnalysisResults.root     MCtrainResults/$LHC18l7/kIncohPsi2sToMu/AnalysisResults.root

echo $kIncohPsi2sToMuPi
mkdir -p MCtrainResults/$(date +%F)/kIncohPsi2sToMuPi
hadd     MCtrainResults/$(date +%F)/kIncohPsi2sToMuPi/AnalysisResults.root   MCtrainResults/$LHC16b2/kIncohPsi2sToMuPi/AnalysisResults.root   MCtrainResults/$LHC18l7/kIncohPsi2sToMuPi/AnalysisResults.root

echo "CIAO!"
