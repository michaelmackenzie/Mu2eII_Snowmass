#ifndef TRKANA_HIST
#define TRKANA_HIST

//local includes
#include "mu2eii/ana/ana/TrackHist_t.hh"
#include "mu2eii/ana/ana/TrackPar_t.hh"
#include "mu2eii/ana/ana/Utilities.hh"
#include "mu2eii/ana/ana/mva_data.hh"

//Offline/TrkAna includes
#include "TrkAna/inc/TrkInfo.hh"
#include "TrkAna/inc/TrkCaloHitInfo.hh"
#include "TrkAna/inc/RecoQualInfo.hh"
// #include "TrkAna/inc/TrkQualInfo.hh"
#include "TrkAna/inc/SimInfo.hh"
#include "TrkAna/inc/EventInfo.hh"

//ROOT includes
#include "TNamed.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"

namespace mu2eii {
  // struct TrkQualInfo { //Currently using mu2e::TrkQualInfo is broken due to compiling issues, so use a local definition for now
  //   Float_t _trkqualvars[10]; //mu2e::TrkQualDetail::n_vars];
  //   Float_t _mvaout;
  //   Int_t _mvastat;
  // };

  class TrkAnaHist : public TNamed {
  public:
    TrkAnaHist(const char* Name, const char* Title);

    void InitBranches(TTree* tree);
    void InitTrackPar(TrackPar_t& tp);
    void InitHistograms(TrackHist_t& Hist, int ihist);
    void FillHistograms(TrackHist_t& Hist, TrackPar_t& tp, Double_t weight = 1.);
    void BeginJob();
    void EndJob();

    int Analyze(TTree* tree, const char* fileout = "trkana.hist", int nevents = -1);

    enum {kNHists = 1000};

    TrackPar_t           fTp;
    mu2e::TrkInfo        fTrkInfo;
    mu2e::TrkFitInfo     fTrkFitInfo;
    mu2e::TrkInfoMC      fTrkInfoMC;
    mu2e::TrkInfoMCStep  fTrkInfoMCStep;
    mu2e::TrkCaloHitInfo fTrkCaloHitInfo;
    // mu2e::TrkQualInfo    fTrkQualInfo;
    // TrkQualInfo          fTrkQualInfo;
    mu2e::RecoQualInfo   fRecoQualInfo;
    mu2e::SimInfo        fSimInfo;
    mu2e::EventInfo      fEventInfo;
    TrackHist_t          fHists[kNHists];
    TString*             fHistTitles[kNHists];
    TFile*               fOutFile;
    TDirectory*          fTopDir;
    TDirectory*          fDirs[kNHists];

    Int_t                fEvaluateMVAs;
    mva_data*            fTrkQualMVA; //for on-the-fly MVA evaluations
    mva_data*            fPIDMVA;

};
}
#endif
