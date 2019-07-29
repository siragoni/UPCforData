#!/usr/bin/env bash

# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/462_20190320-1856/merge/AnalysisResults.root  AnalysisResultsTrying.root

# - NB: the train label should be something along the lines of:
# --- Train number
# --- +
# --- Date the train ran
# --- +
# --- Child number
# - that being said an example is:
# --- 457_20190317-1224_child_1


TRAINlabel=$1
echo $TRAINlabel

if [ -d "MCtrainResults/$(date +%F)" ]; then rm -rf MCtrainResults/$(date +%F); fi
mkdir -p MCtrainResults/$(date +%F)

BASEpath=alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/
LASTpartPath=merge/AnalysisResults.root
kCohJpsiToMu=$BASEpath
kCohJpsiToMu+=$TRAINlabel
kCohJpsiToMu+=_child_1/
kCohJpsiToMu+=$LASTpartPath
echo $kCohJpsiToMu
mkdir -p MCtrainResults/$(date +%F)/kCohJpsiToMu
alien_cp $kCohJpsiToMu  MCtrainResults/$(date +%F)/kCohJpsiToMu/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_1/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kCohJpsiToMu/AnalysisResults.root

kIncohJpsiToMu=$BASEpath
kIncohJpsiToMu+=$TRAINlabel
kIncohJpsiToMu+=_child_2/
kIncohJpsiToMu+=$LASTpartPath
echo $kIncohJpsiToMu
mkdir -p MCtrainResults/$(date +%F)/kIncohJpsiToMu
alien_cp $kIncohJpsiToMu  MCtrainResults/$(date +%F)/kIncohJpsiToMu/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_2/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kIncohJpsiToMu/AnalysisResults.root

kTwoGammaToMuMedium=$BASEpath
kTwoGammaToMuMedium+=$TRAINlabel
kTwoGammaToMuMedium+=_child_3/
kTwoGammaToMuMedium+=$LASTpartPath
echo $kTwoGammaToMuMedium
mkdir -p MCtrainResults/$(date +%F)/kTwoGammaToMuMedium
alien_cp $kTwoGammaToMuMedium  MCtrainResults/$(date +%F)/kTwoGammaToMuMedium/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_3/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kTwoGammaToMuMedium/AnalysisResults.root

kTwoGammaToMuHigh=$BASEpath
kTwoGammaToMuHigh+=$TRAINlabel
kTwoGammaToMuHigh+=_child_4/
kTwoGammaToMuHigh+=$LASTpartPath
echo $kTwoGammaToMuHigh
mkdir -p MCtrainResults/$(date +%F)/kTwoGammaToMuHigh
alien_cp $kTwoGammaToMuHigh  MCtrainResults/$(date +%F)/kTwoGammaToMuHigh/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_4/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kTwoGammaToMuHigh/AnalysisResults.root

kCohPsi2sToMu=$BASEpath
kCohPsi2sToMu+=$TRAINlabel
kCohPsi2sToMu+=_child_5/
kCohPsi2sToMu+=$LASTpartPath
echo $kCohPsi2sToMu
mkdir -p MCtrainResults/$(date +%F)/kCohPsi2sToMu
alien_cp $kCohPsi2sToMu  MCtrainResults/$(date +%F)/kCohPsi2sToMu/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_5/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kCohPsi2sToMu/AnalysisResults.root

kCohPsi2sToMuPi=$BASEpath
kCohPsi2sToMuPi+=$TRAINlabel
kCohPsi2sToMuPi+=_child_6/
kCohPsi2sToMuPi+=$LASTpartPath
echo $kCohPsi2sToMuPi
mkdir -p MCtrainResults/$(date +%F)/kCohPsi2sToMuPi
alien_cp $kCohPsi2sToMuPi  MCtrainResults/$(date +%F)/kCohPsi2sToMuPi/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_6/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kCohPsi2sToMuPi/AnalysisResults.root

kIncohPsi2sToMu=$BASEpath
kIncohPsi2sToMu+=$TRAINlabel
kIncohPsi2sToMu+=_child_7/
kIncohPsi2sToMu+=$LASTpartPath
echo $kIncohPsi2sToMu
mkdir -p MCtrainResults/$(date +%F)/kIncohPsi2sToMu
alien_cp $kIncohPsi2sToMu  MCtrainResults/$(date +%F)/kIncohPsi2sToMu/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_7/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kIncohPsi2sToMu/AnalysisResults.root

kIncohPsi2sToMuPi=$BASEpath
kIncohPsi2sToMuPi+=$TRAINlabel
kIncohPsi2sToMuPi+=_child_8/
kIncohPsi2sToMuPi+=$LASTpartPath
echo $kIncohPsi2sToMuPi
mkdir -p MCtrainResults/$(date +%F)/kIncohPsi2sToMuPi
alien_cp $kIncohPsi2sToMuPi  MCtrainResults/$(date +%F)/kIncohPsi2sToMuPi/AnalysisResults.root
# alien_cp alien:///alice/cern.ch/user/a/alitrain/PWGUD/UD_PbPb_AOD/$(TRAINlabel)_child_8/merge/AnalysisResults.root  MCtrainResults/$(date +%F)/kIncohPsi2sToMuPi/AnalysisResults.root

echo "CIAO!"
