TString fTreeDir = "/mu2e/data/projects/mu2eii_snowmass/";
int fEvaluateMVAs = 0; //SU2020 TrkQual MVA currently gives worse resolution, possibly due to DAR training applied to PAR tracks
int fMomCorrMode  = 0; //Nominal: 0; Mu2e era tracker straws: 1;

using namespace mu2eii;

int process_file(TString fname, TString dname, int nentries = -1) {
  TFile* file = TFile::Open(fname.Data(), "READ");
  if(!file) return 1;
  TTree* tree = (TTree*) file->Get("TrkAnaNeg/trkana");
  if(!tree) {
    cout << "TrkAna Tree for dataset " << dname.Data() << " not found\n";
    return 2;
  }
  cout << "Processing dataset " << dname.Data() << endl;
  TrkAnaHist analyzer(dname.Data(), dname.Data());
  analyzer.fEvaluateMVAs = fEvaluateMVAs;
  analyzer.fMomCorrMode  = fMomCorrMode;
  TString outname = dname; outname.ToLower();
  outname = Form("mu2eiisnowmass.%s.trkanahist.1011.hist", outname.Data());
  analyzer.Analyze(tree, outname.Data(), nentries);
  delete tree;
  file->Close();
  return 0;
}

int process_trkana(TString tag = "", int nentries = -1) {
  fMomCorrMode = 0;

  vector<TString> datasets = {"CE", "CE_UT", "CE_Mix", "CE_Mix_UT", "CE_Mix_NT", "DIO"};
  for(TString dataset : datasets) {
    if(tag != "" && !dataset.Contains(tag)) continue;
    TString fname;
    if(dataset == "CE_Mix") {
      fname = "/mu2e/data/users/mmackenz/mu2eii/snowmass/tree/nts.mu2e.CeEndpointMixTriggered.Mu2eIIa2b.root";
    } else if(dataset == "CE_Mix_UT") {
      fname = "/mu2e/data/users/mmackenz/mu2eii/snowmass/tree/nts.mu2e.CeEndpointMixUntriggered.TrkAnaReco.root";
    } else if(dataset == "CE_Mix_NT") {
      fname = "/mu2e/data/users/mmackenz/mu2eii/snowmass/tree/nts.mu2e.CeEndpointMix.TrkAnaReco.root";
    } else if(dataset == "CE_UT") {
      fname = fTreeDir + "CE_TrkAna_Untrig_Mu2eII_NewGeom.root";
    } else {
      fname = fTreeDir + dataset + "_TrkAna_Trig_Mu2eII_NewGeom.root";
    }
    process_file(fname, dataset, nentries);
  }
  return 0;
}

//Process the no-mix ntuple with the Mu2e era tracker straw parameters
int process_old_tracker(int nentries = -1) {
  fMomCorrMode = 1; //different momentum offset
  TString fname = "/mu2e/data/users/mmackenz/mu2eii/snowmass/tree/nts.mu2e.CeEndpointDigiTriggered_TrkAnaReco.Mu2eIIbOldTrk.root";
  return process_file(fname, "CE_OldTrk", nentries);
}

//Process the no-mix ntuple without the IPA
int process_no_ipa(int nentries = -1) {
  fMomCorrMode = 0;
  TString fname = "/mu2e/data/users/mmackenz/mu2eii/snowmass/tree/mcs.mu2e.CeEndpointDigiTriggered.Mu2eIIbNoIPA.root";
  return process_file(fname, "CE_NoIPA", nentries);
}

int process_all() {
  int status(0);
  status += process_trkana();
  status += process_old_tracker();
  status += process_no_ipa();
  return status;
}
