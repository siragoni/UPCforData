void runUnfoldSimple()
{
//Load In Monte Carlo Truth and detector level
// TFile* fResponse = new TFile("AnalysisResultsLHC18l7_coherent_30112021.root");
  TFile* fResponse = new TFile("AnalysisResults.root");
  TDirectory* dir = fResponse->GetDirectory("MyTask");
  TList* lResponse;
  dir->GetObject("MyOutputContainer", lResponse);

  TH1F* JCMeasured = (TH1F*) lResponse->FindObject("fPhiHelicityFrameTwentyfiveBinsHv2_restrict");
  TH1F* JCTruth = (TH1F*) lResponse->FindObject("fMCPhiHelicityFrameTwentyfiveBinsHv2_restrict");

//Normalise Monte Carlo Truth and detector level
  //JCMeasured->Scale(1.0/(JCMeasured->Integral()));
  //JCTruth->Scale(1.0/(JCTruth->Integral()));

//Load in Response matrix
  TH2F* JCResponse = (TH2F*) lResponse->FindObject("fPhiRecVsGenHelicityH");
JCResponse->Rebin2D(24,24);
//Load in Data
  TFile* fData = new TFile("pngResults/2021-09-21/PhiHEv2/PhiHeFrameV2.root");
  TH1F* JCData = (TH1F*)fData->Get("PhiAfterSignalExtractionErrorsH");


  // JCData->Scale(1.0/JCData->Integral());


  RooUnfoldResponse response(JCMeasured ,JCTruth,JCResponse);

  // RooUnfoldBayes unfold(&response, JCData,2);
  RooUnfoldBayes unfold(&response, JCMeasured,2);


  TH1D* hist_reco= (TH1D*) unfold.Hreco();

  // hist_reco->Scale(1.0/hist_reco->Integral());
  // JCTruth->Scale(1.0/(JCTruth->Integral()));


  //Draw the particle level histogram
  JCTruth->SetLineColor(kRed);
  hist_reco->Draw("");
  JCTruth->Draw("SAME");


  JCTruth->Chi2Test(hist_reco, "p");

}
