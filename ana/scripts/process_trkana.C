TString fTreeDir = "/mu2e/data/projects/mu2eii_snowmass/";

using namespace mu2eii;

int process_trkana(int nentries = -1) {
  vector<TString> datasets = {"CE", "DIO"};

  for(TString dataset : datasets) {
    TString fname = fTreeDir + dataset + "_TrkAna_Trig_Mu2eII_NewGeom.root";
    TFile* file = TFile::Open(fname.Data(), "READ");
    if(!file) continue;
    TTree* tree = (TTree*) file->Get("TrkAnaNeg/trkana");
    if(!tree) {
      cout << "TrkAna Tree for dataset " << dataset.Data() << " not found\n";
      continue;
    }
    cout << "Processing dataset " << dataset.Data() << endl;
    // tree->Print();
    TrkAnaHist analyzer(dataset.Data(), dataset.Data());
    TString outname = dataset; outname.ToLower();
    outname = Form("mu2eiisnowmass.%s.trkanahist.1011.hist", outname.Data());
    analyzer.Analyze(tree, outname.Data(), nentries);
    delete tree;
    file->Close();
  }
  return 0;
}
