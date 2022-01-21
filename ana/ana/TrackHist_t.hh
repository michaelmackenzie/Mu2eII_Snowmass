#ifndef TRACKHIST_T
#define TRACKHIST_T
#include "TH1.h"
#include "TH2.h"

struct TrackHist_t {
  TH1F* hP[3];
  TH1F* hT0;
  TH2F* hPvsT0;
  TH1F* hPFront;
  TH1F* hDpf;
  TH2F* hPvsPFront;
  TH1F* hPErr;
  TH1F* hT0Err;
  TH1F* hD0;
  TH1F* hRMax;
  TH1F* hTDip;
  TH1F* hNHits;
  TH1F* hChiSq;
  TH1F* hChiSqR;
  TH1F* hEoverP;
  TH1F* hECluster;
  TH1F* hTrkQual;
  TH1F* hGenEnergy;
  TH1F* hGenCode;
};

#endif
