#ifndef TRACKHIST_T
#define TRACKHIST_T
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

struct TrackHist_t {
  TH1F* hP[3];
  TH1F* hT0;
  TH2F* hPvsT0;
  TH3F* hPvsT0vsMVA;
  TH1F* hPFront;
  TH1F* hDpf;
  TH1F* hDpGen;
  TH2F* hPvsPFront;
  TH1F* hPExit;
  TH1F* hPMCExit;
  TH1F* hDpExit;
  TH1F* hPErr;
  TH1F* hT0Err;
  TH1F* hD0;
  TH1F* hRMax;
  TH1F* hTDip;
  TH1F* hNHits;
  TH1F* hNActiveHits;
  TH1F* hNActiveFraction;
  TH1F* hNDoublets;
  TH1F* hNActiveDoublets;
  TH1F* hNActDoubletFrac;
  TH1F* hNHitsAmbZero;
  TH1F* hNHitsAmbZeroFrac;
  TH1F* hNMat;
  TH1F* hNActiveMat;
  TH1F* hNActiveMatFrac;
  TH1F* hChiSq;
  TH1F* hChiSqR;
  TH1F* fLogFitCon;
  TH1F* hECluster;
  TH1F* hEoverP;
  TH1F* hNCrystals;
  TH1F* hSeedFrac;
  TH1F* hTchDt;
  TH1F* hTchDz;
  TH1F* hTchDr;
  TH1F* hPath;

  TH1F* hTrkQual;
  TH1F* hPIDScore;

  TH1F* hNPOT;

  TH1F* hGenEnergy;
  TH1F* hGenCode;
};

#endif
