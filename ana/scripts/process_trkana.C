TString fTreeDir = "/mu2e/data/projects/mu2eii_snowmass/";
int fEvaluateMVAs = 0; //SU2020 TrkQual MVA currently gives worse resolution, possibly due to DAR training applied to PAR tracks

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
  TString outname = dname; outname.ToLower();
  outname = Form("mu2eiisnowmass.%s.trkanahist.1011.hist", outname.Data());
  analyzer.Analyze(tree, outname.Data(), nentries);
  delete tree;
  file->Close();
  return 0;
}

int process_trkana(int nentries = -1) {
  vector<TString> datasets = {"CE", "DIO"};

  for(TString dataset : datasets) {
    TString fname = fTreeDir + dataset + "_TrkAna_Trig_Mu2eII_NewGeom.root";
    process_file(fname, dataset, nentries);
  }
  return 0;
}

int process_old_tracker(int nentries = -1) {
  TString fname = "/mu2e/data/users/mmackenz/mu2eii/snowmass/tree/nts.mu2e.CeEndpointDigiTriggered_TrkAnaReco.Mu2eIIbOldTrk.root";
  return process_file(fname, "CE_OldTrk", nentries);
}
