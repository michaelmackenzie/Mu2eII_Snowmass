//make some plots of relevant parameters
#include "mu2eii/ana/ana/TrkAnaHist.hh"

namespace mu2eii {
  TrkAnaHist::TrkAnaHist(const char* Name, const char* Title) : TNamed(Name, Title) {
    fEvaluateMVAs = 0; //whether or not to re-evaluate MVAs on the fly

    //SU2020 MVAs
    fTrkQualMVA = new mva_data("fele2s51b1",1070);
    fPIDMVA     = new mva_data("ele00s61b0",1000);
  }

  void TrkAnaHist::BeginJob() {

    for(int ihist = 0; ihist < kNHists; ++ihist) fHistTitles[ihist] = 0;

    fHistTitles[  0] = new TString("all events");
    fHistTitles[  1] = new TString("all tracks");
    fHistTitles[  2] = new TString("all tracks with weights");
    fHistTitles[100] = new TString("all tracks passing track ID");
    fHistTitles[101] = new TString("all tracks passing track ID with weights");
    fHistTitles[102] = new TString("all tracks passing track ID except failing TrkQual");
    fHistTitles[103] = new TString("all tracks passing track ID without TrkQual cut");
    fHistTitles[110] = new TString("all tracks passing track ID, T0 > 700 ns");
    fHistTitles[111] = new TString("all tracks passing track ID with weights, T0 > 700 ns");
    fHistTitles[120] = new TString("all tracks passing track ID, T0 > 640 ns");
    fHistTitles[121] = new TString("all tracks passing track ID with weights, T0 > 640 ns");
    fHistTitles[130] = new TString("all tracks passing track ID, T0 > 680 ns");
    fHistTitles[131] = new TString("all tracks passing track ID with weights, T0 > 680 ns");

    fHistTitles[200] = new TString("all tracks passing track ID with Offline TrkQual cut");
    fHistTitles[210] = new TString("all tracks passing track ID with Offline TrkQual cut, T0 > 700 ns");
    fHistTitles[220] = new TString("all tracks passing track ID with Offline TrkQual cut, T0 > 640 ns");
    fHistTitles[230] = new TString("all tracks passing track ID with Offline TrkQual cut, T0 > 680 ns");

    fTopDir = fOutFile->mkdir(GetName());
    for(int ihist = 0; ihist < kNHists; ++ihist) {
      if(fHistTitles[ihist]) InitHistograms(fHists[ihist], ihist);
      else fDirs[ihist] = 0;
    }
  }

  void TrkAnaHist::EndJob() {
    printf("TrkAnaHist::%s: Writing out histogram file to %s\n", __func__, fOutFile->GetName());
    fTopDir->cd();
    // fTopDir->Write();
    // for(int ihist = 0; ihist < kNHists; ++ihist) {
    //   if(fDirs[ihist]) fDirs[ihist]->Write();
    // }
    fOutFile->Write();
    fOutFile->Close();
  }

  void TrkAnaHist::InitHistograms(TrackHist_t& Hist, int ihist) {
    TDirectory* dir = fTopDir->mkdir(Form("trk_%i", ihist));
    dir->SetTitle(fHistTitles[ihist]->Data());
    dir->cd();
    fDirs[ihist] = dir;
    Hist.hP[0]      = new TH1F("p_0", "Track P", 400, 0, 200);
    Hist.hP[1]      = new TH1F("p_1", "Track P", 400, 80, 120);
    Hist.hP[2]      = new TH1F("p_2", "Track P", 400, 0, 200);
    Hist.hT0        = new TH1F("t0", "Track T0", 200, 0, 2000);
    Hist.hPvsT0     = new TH2F("p_vs_t0", "Track P vs T0", 800, 80, 120, 200, 0, 2000);
    Hist.hPFront    = new TH1F("pfront", "Track P(Front MC)", 1000, 0, 200);
    Hist.hDpf       = new TH1F("dpf", "Track P - Track P(Front MC)", 200, -10, 10);
    Hist.hPvsPFront = new TH2F("p_vs_front", "Track P vs P(Front MC)", 200, 50, 150, 200, 50, 150);
    Hist.hPErr      = new TH1F("p_err", "Track P uncertainty", 100, 0, 5);
    Hist.hT0Err     = new TH1F("t0_err", "Track T0 uncertainty", 100, 0, 10);
    Hist.hD0        = new TH1F("d0", "Track D0", 100, -250, 250);
    Hist.hRMax      = new TH1F("rmax", "Track R(max)", 100, 0, 1000);
    Hist.hTDip      = new TH1F("tdip", "Track tan(dip)", 100, 0, 3);
    Hist.hNHits             = new TH1F("nhits", "Track N(hits)", 200, 0, 200);
    Hist.hNActiveHits       = new TH1F("nactvhits", "Track N(active hits)", 200, 0, 200);
    Hist.hNActiveFraction   = new TH1F("nactvfrac", "Track N(active hits) / N(hits)", 110, 0., 1.1);
    Hist.hNDoublets         = new TH1F("ndoublets", "Track N(doublets)", 200, 0, 200);
    Hist.hNActiveDoublets   = new TH1F("nactvdoublets", "Track N(active doublets)", 200, 0, 200);
    Hist.hNActDoubletFrac   = new TH1F("nactvdfrac", "Track N(active doublets) / N(active hits)", 110, 0., 1.1);
    Hist.hNHitsAmbZero      = new TH1F("nnullhits", "Track N(null ambiguity hits)", 200, 0, 200);
    Hist.hNHitsAmbZeroFrac  = new TH1F("nnullfrac", "Track N(null hits) / N(active hits)", 110, 0., 1.1);
    Hist.hNMat              = new TH1F("nmat", "Track N(mat)", 200, 0, 200);
    Hist.hNActiveMat        = new TH1F("nactivemat", "Track N(active mat)", 200, 0, 200);
    Hist.hNActiveMatFrac    = new TH1F("nactvmatfrac", "Track N(active mat) / N(mat)", 110, 0., 1.1);
    Hist.hChiSq     = new TH1F("chisq", "Track #chi^2", 100, 0, 100);
    Hist.hChiSqR    = new TH1F("chisq_r", "Track reduced #chi^2", 100, 0, 10);
    Hist.hECluster  = new TH1F("e_cluster", "Cluster energy", 300, 0, 150);
    Hist.hEoverP    = new TH1F("e_over_p", "Cluster E over Track P", 100, 0, 2);
    Hist.hNCrystals = new TH1F("n_crystals", "N(Crystals)", 10, 0, 10);
    Hist.hSeedFrac  = new TH1F("seed_frac", "Seed Fraction", 110, 0., 1.1);
    Hist.hTchDt     = new TH1F("tch_dt", "TrkCaloHit #DeltaT0", 100, -50, 50.);
    Hist.hTchDz     = new TH1F("tch_dz", "TrkCaloHit #Deltaz", 200, -250, 250.);
    Hist.hTchDr     = new TH1F("tch_dr", "TrkCaloHit #DeltaR", 200, -250, 250.);
    Hist.hPath      = new TH1F("tch_path", "TrkCaloHit Path", 200, 00, 600.);

    Hist.hTrkQual   = new TH1F("trkqual", "Track ANN quality score", 110, 0., 1.1);
    Hist.hPIDScore  = new TH1F("pid", "Track PID score", 260, -1.5, 1.1);

    Hist.hNPOT      = new TH1F("npot", "N(POT)", 150, 0, 1.5e8);
    Hist.hGenEnergy = new TH1F("gen_e", "Gen energy", 200, 0, 200);
    Hist.hGenCode   = new TH1F("gen_id", "Gen ID", 200, 0, 200);
  }

  void TrkAnaHist::FillHistograms(TrackHist_t& Hist, TrackPar_t& tp, Double_t weight) {
    Hist.hP[0]      ->Fill(tp.fP , weight);
    Hist.hP[1]      ->Fill(tp.fP , weight);
    Hist.hP[2]      ->Fill(tp.fP , weight);
    Hist.hT0        ->Fill(tp.fT0 , weight);
    Hist.hPvsT0     ->Fill(tp.fP, tp.fT0 , weight);
    Hist.hPFront    ->Fill(tp.fPFront , weight);
    Hist.hDpf       ->Fill(tp.fP - tp.fPFront , weight);
    Hist.hPvsPFront ->Fill(tp.fP, tp.fPFront , weight);
    Hist.hPErr      ->Fill(tp.fPErr , weight);
    Hist.hT0Err     ->Fill(tp.fT0Err , weight);
    Hist.hD0        ->Fill(tp.fD0 , weight);
    Hist.hRMax      ->Fill(tp.fRMax , weight);
    Hist.hTDip      ->Fill(tp.fTanDip , weight);
    Hist.hNHits     ->Fill(tp.fNHits , weight);
    Hist.hNActiveHits       ->Fill(tp.fNActiveHits       , weight);
    Hist.hNActiveFraction   ->Fill(tp.fNActiveFraction, weight);
    Hist.hNDoublets         ->Fill(tp.fNDoublets         , weight);
    Hist.hNActiveDoublets   ->Fill(tp.fNActiveDoublets   , weight);
    Hist.hNActDoubletFrac   ->Fill(tp.fNActiveDoublets / (1.*tp.fNActiveHits), weight);
    Hist.hNHitsAmbZero      ->Fill(tp.fNHitsAmbZero      , weight);
    Hist.hNHitsAmbZeroFrac  ->Fill(tp.fNHitsAmbZero / (1.*tp.fNActiveHits), weight);
    Hist.hNMat              ->Fill(tp.fNMat              , weight);
    Hist.hNActiveMat        ->Fill(tp.fNActiveMat        , weight);
    Hist.hNActiveMatFrac    ->Fill(tp.fNActiveMat / (1.*tp.fNMat), weight);


    Hist.hChiSq     ->Fill(tp.fChiSq , weight);
    Hist.hChiSqR    ->Fill((tp.fNDoF > 0) ? tp.fChiSq/tp.fNDoF : -1., weight);

    Hist.hECluster  ->Fill(tp.fECluster , weight);
    Hist.hEoverP    ->Fill(tp.fEOverP   , weight);
    Hist.hNCrystals ->Fill(tp.fNCrystals, weight);
    Hist.hSeedFrac  ->Fill(tp.fSeedFr   , weight);
    Hist.hTchDt     ->Fill(tp.fTchDt    , weight);
    Hist.hTchDz     ->Fill(tp.fTchDz    , weight);
    Hist.hTchDr     ->Fill(tp.fTchDr    , weight);
    Hist.hPath      ->Fill(tp.fPath     , weight);

    Hist.hTrkQual   ->Fill(tp.fTrkQual  , weight);
    Hist.hPIDScore  ->Fill(tp.fPIDScore , weight);

    Hist.hNPOT      ->Fill(tp.fNPOT     , weight);
    Hist.hGenEnergy ->Fill(tp.fGenEnergy, weight);
    Hist.hGenCode   ->Fill(tp.fGenCode  , weight);
  }

  void TrkAnaHist::InitBranches(TTree* tree) {
    TString type = "de";
    tree->SetBranchAddress(Form("%s"       , type.Data()),  &fTrkInfo   );
    tree->SetBranchAddress(Form("%sent"    , type.Data()),  &fTrkFitInfo);
    tree->SetBranchAddress(Form("%smcpri"  , type.Data()),  &fSimInfo   );
    tree->SetBranchAddress(Form("%smcent"  , type.Data()),  &fTrkInfoMCStep);
    tree->SetBranchAddress(Form("%stch"    , type.Data()),  &fTrkCaloHitInfo);
    // tree->SetBranchAddress(Form("%strkqual", type.Data()),  &fTrkQualInfo);
    tree->SetBranchAddress(Form("%squal"   , type.Data()),  &fRecoQualInfo);
  }

  void TrkAnaHist::InitTrackPar(TrackPar_t& tp) {
    tp.fT0    = fTrkInfo._t0;
    tp.fT0Err = fTrkInfo._t0err;
    tp.fNHits = fTrkInfo._nhits;
    tp.fNActiveHits = fTrkInfo._nactive;
    tp.fNDoublets = fTrkInfo._ndouble;
    tp.fNActiveDoublets = fTrkInfo._ndactive;
    tp.fNHitsAmbZero = fTrkInfo._nnullambig;
    tp.fNMat = fTrkInfo._nmat;
    tp.fNActiveMat = fTrkInfo._nmatactive;
    tp.fNDoF  = fTrkInfo._ndof;
    tp.fChiSq = fTrkInfo._chisq;
    tp.fFitCons = fTrkInfo._fitcon;
    tp.fP     = fTrkFitInfo._fitmom;
    tp.fPErr  = fTrkFitInfo._fitmomerr;
    tp.fPFront = fTrkInfoMCStep._mom;
    tp.fTanDip = fTrkFitInfo._fitpar._td;
    tp.fD0     = fTrkFitInfo._fitpar._d0;
    tp.fRMax   = fabs(tp.fD0 + 2./fTrkFitInfo._fitpar._om);
    tp.fGenEnergy = fSimInfo._mom;
    tp.fGenCode   = fSimInfo._gen;
    tp.fECluster = (fTrkCaloHitInfo._did < 0) ? 0. : fTrkCaloHitInfo._edep;
    tp.fDefaultTrkQual = fRecoQualInfo._qualsAndCalibs[2];

    tp.fNPOT = fEventInfo._nprotons;

    //PID variables  FIXME: Initialize PID variables
    tp.fEOverP = tp.fECluster / tp.fP;
    tp.fTchDt = tp.fT0 - fTrkCaloHitInfo._t0; //FIXME: Check that this is the Dt definition, likely missing a dT offset
    tp.fTchDr = 0.;
    tp.fTchDz = 0.;
    tp.fSeedFr = 1.;

    //Track Quality variables
    tp.fNActiveFraction = ((double) tp.fNActiveHits) / (tp.fNHits);

    if(fEvaluateMVAs) {

      //PID MVA
      fPIDMVA->fVar[0] = tp.fEOverP;
      fPIDMVA->fVar[1] = tp.fNCrystals;
      fPIDMVA->fVar[2] = tp.fSeedFr;
      fPIDMVA->fVar[3] = tp.fTchDt;
      fPIDMVA->fVar[4] = tp.fTchDz;
      fPIDMVA->fVar[5] = tp.fTchDr; //DOCA
      fPIDMVA->fVar[6] = tp.fPath; //FIXME: Path (ds) is not apparent in the TrkCaloHitInfo branch

      if ((fabs(tp.fTchDt) <  20.) && (fabs(tp.fTchDr) < 100. ) && //Cuts to define the PID selection
          (tp.fTchDz       > -50.) && (tp.fTchDz       < 250. ) &&
          (tp.fEOverP      >   0.) && (tp.fEOverP      < 1.05)    ) {
        tp.fPIDScore = 1.; //FIXME: For now assume all pass PID, which is 99.2% efficient for CE in SU2020, after applying cuts that make the PID selection well defined
      } else {
        tp.fPIDScore = 0.;
      }

      //Track Quality MVA
      double na(tp.fNActiveHits), nm(tp.fNMat);
      fTrkQualMVA->fVar[0] = tp.fNActiveHits;
      fTrkQualMVA->fVar[1] = tp.fNActiveFraction;
      fTrkQualMVA->fVar[2] = std::log10(tp.fFitCons);
      fTrkQualMVA->fVar[3] = tp.fPErr;
      fTrkQualMVA->fVar[4] = tp.fT0Err;
      fTrkQualMVA->fVar[5] = tp.fNActiveDoublets/na;
      fTrkQualMVA->fVar[6] = tp.fNHitsAmbZero/na;
      fTrkQualMVA->fVar[7] = tp.fNActiveMat/nm;


      tp.fTrkQual   = fTrkQualMVA->Eval();

    } else { //read in MVA scores from the TrkAna Tree
      // tp.fTrkQual  = (fTrkQualInfo._mvastat != 2) ? -1. : fTrkQualInfo._mvaout;
      //In RecoQualInfo, entries are: TrkPID, TrkPIDCalib, TrkQual, TrkQualCalib
      tp.fTrkQual  = fRecoQualInfo._qualsAndCalibs[2];
      tp.fPIDScore = fRecoQualInfo._qualsAndCalibs[0];
    }

    tp.fDIOWeight = Util::DioWeightAl_LL(tp.fGenEnergy);
    tp.fEventWeight = (tp.fGenCode == 166) ? tp.fDIOWeight : 1.;

    for(int i = 0; i < 20; ++i) tp.fTrackID[i] = 0;
    tp.fTrackID[0] += (1 << 0)*!(fabs(tp.fD0) < 100.); //D0 consistent with the stopping target
    tp.fTrackID[0] += (1 << 1)*!(tp.fECluster > 10. && tp.fEOverP > 0. && tp.fEOverP < 1.05); //has a cluster with reasonable E/P
    tp.fTrackID[0] += (1 << 2)*!(0.5 < tp.fTanDip && tp.fTanDip < 1.); //reject cosmic backgrounds/beam electrons
    tp.fTrackID[0] += (1 << 3)*!(tp.fNHits >= 20); //enough hits to reconstruct the helix well
    tp.fTrackID[0] += (1 << 4)*!(tp.fRMax < 680.); //avoid interactions with tracker edge
    tp.fTrackID[0] += (1 << 5)*!(tp.fT0Err < 0.9); //ignore poor fits, especially inconsistent with cluster timing
    if(fEvaluateMVAs) {
      tp.fTrackID[0] += (1 << 6)*!(tp.fTrkQual > fTrkQualMVA->CutValue()); //apply the ANN track quality
      tp.fTrackID[0] += (1 << 9)*!(tp.fPIDScore > 0.5); //passes PID selection
    } else {
      tp.fTrackID[0] += (1 << 6)*!(tp.fTrkQual > 0.8); //apply the ANN track quality
      // FIXME: Investigate low PID efficiency with default score
      // tp.fTrackID[0] += (1 << 9)*!(tp.fPIDScore > 0.); //passes PID selection
    }
    tp.fTrackID[0] += (1 << 7)*!(tp.fT0 < 1650.); //apply the upper limit time cut
    tp.fTrackID[0] += (1 << 8)*!(tp.fP > 100.); //must be above the DIO generation minimum energy
  }

  int TrkAnaHist::Analyze(TTree* tree, const char* fileout, int nevents) {
    int status(0);
    if(!tree) {
      printf("TrkAnaHist::%s: Input tree is not defined\n", __func__);
      return 1;
    }
    fOutFile = new TFile(fileout, "RECREATE");
    if(!fOutFile) {
      printf("TrkAnaHist::%s: Output file cannot be created\n", __func__);
      return 2;
    }
    BeginJob();
    InitBranches(tree);

    //Main processing loop
    if(nevents < 0) nevents = tree->GetEntriesFast();
    for(int ievent = 0; ievent < nevents; ++ievent) {
      if(ievent % 10000 == 0) printf("TrkAnaHist::%s: Processing event %10i (%6.2f%%)...\n", __func__, ievent, (100.*ievent)/nevents);
      tree->GetEntry(ievent);
      InitTrackPar(fTp);
      //All events
      FillHistograms(fHists[0], fTp);
      //All events with tracks
      if(fTp.fNHits > 0) {
        FillHistograms(fHists[1], fTp);
        //With event weights (if applicable)
        FillHistograms(fHists[2], fTp, fTp.fEventWeight);
        //Track selection cut
        if(fTp.fTrackID[0] == 0) {
          FillHistograms(fHists[100], fTp);
          FillHistograms(fHists[101], fTp, fTp.fEventWeight);
          if(fTp.fT0 > 700.) { //CD3 timing cut
            FillHistograms(fHists[110] ,fTp);
            FillHistograms(fHists[111], fTp, fTp.fEventWeight);
          }
          if(fTp.fT0 > 640.) { //SU2020 timing cut
            FillHistograms(fHists[120], fTp);
            FillHistograms(fHists[121], fTp, fTp.fEventWeight);
          }
          if(fTp.fT0 > 680.) { //Mu2e-II timing cut
            FillHistograms(fHists[130], fTp);
            FillHistograms(fHists[131], fTp, fTp.fEventWeight);
          }
        } else if((fTp.fTrackID[0] & (~(1 << 6))) == 0) { //fails only TrkQual
          FillHistograms(fHists[102], fTp);
        }
        if((fTp.fTrackID[0] & (~(1 << 6))) == 0) { //Track ID without TrkQual cut
          FillHistograms(fHists[103], fTp);
          if(fTp.fDefaultTrkQual > 0.8) { //Offline TrkQual
            FillHistograms(fHists[200], fTp);
            if(fTp.fT0 > 700.) { //CD3 timing cut
              FillHistograms(fHists[210], fTp);
            }
            if(fTp.fT0 > 640.) { //SU2020 timing cut
              FillHistograms(fHists[220], fTp);
            }
            if(fTp.fT0 > 680.) { //Mu2e-II timing cut
              FillHistograms(fHists[230], fTp);
            }
          }
        }
      }
    }
    EndJob();
    return status;
  }
}
