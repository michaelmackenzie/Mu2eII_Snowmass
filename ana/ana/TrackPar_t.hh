#ifndef TRACKPAR_T
#define TRACKPAR_T

struct TrackPar_t {
  double fP;
  double fPFront;
  double fT0;
  double fPErr;
  double fT0Err;
  double fD0;
  double fRMax;
  double fTanDip;
  int    fNHits;
  double fChiSq;
  int    fNDoF;
  double fECluster;
  double fTrkQual;

  double fGenEnergy;
  double fGenCosTh;
  int    fGenCode;
  double fEventWeight;
  double fDIOWeight;

  int    fTrackID[20];
};

#endif
