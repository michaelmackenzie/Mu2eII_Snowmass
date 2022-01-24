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
  int    fNActiveHits;
  double fNActiveFraction;
  int    fNDoublets;
  int    fNActiveDoublets;
  int    fNHitsAmbZero;
  int    fNMat;
  int    fNActiveMat;
  double fFitCons;
  double fChiSq;
  int    fNDoF;
  double fECluster;
  double fEOverP;
  double fTrkQual;
  double fPIDScore;
  double fDefaultTrkQual;

  int    fNPOT;

  //PID variables
  int    fNCrystals;
  double fSeedFr;
  double fTchDt;
  double fTchDz;
  double fTchDr;
  double fPath;

  double fGenEnergy;
  double fGenP;
  double fGenCosTh;
  int    fGenCode;
  double fEventWeight;
  double fDIOWeight;

  int    fTrackID[20];
};

#endif
