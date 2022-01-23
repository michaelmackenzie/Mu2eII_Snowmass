///////////////////////////////////////////////////////////////////////////////
// use numbers from
// to be compiled

#include "mu2eii/stat/stat/channel.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "TEnv.h"
#include "TF2.h"

namespace mu2eii {

  //-----------------------------------------------------------------------------
  // free function -- functional parameterization of the 2D cosmic background
  //-----------------------------------------------------------------------------
  double f_bgr_cosmics(double* X, double* P) {
    //-----------------------------------------------------------------------------
    // to begin with, assume cosmics uniform
    // missing scale factor to CORSIKA - set to 1, add a 5% uncertainty
    // assume P = X[0] ; T = X[1]
    //-----------------------------------------------------------------------------
    const float p        = X[0];
    const float t        = X[1];
    float rho_per_mev    = (14.57+0.3199*(p-105.))/5./110.73 ;     // 2025lo background per MeV, integrated over time 700-1695 = 995 ns
    rho_per_mev         += 2./50./3.70;                            // 2 events from 2025hi (positrons+electrons = 1 event)/2
    // account for T>1650... - correction on the tail "eats up" 2.4%,
    // renormalize the integral
    if (t > 1650) rho_per_mev = rho_per_mev*(1.-(t-1650.)/40.)*1.024;
    if (t > 1690) rho_per_mev = 0;

    const float time_window    = 995.; //width of the time window used to make the initial estimate
    const float rho_per_mev_ns = rho_per_mev/time_window; //density per MeV/c per ns

    return rho_per_mev_ns;
  }

  //-----------------------------------------------------------------------------
  // free function -- plot the 2D cosmic background using the functional parameterization
  //-----------------------------------------------------------------------------
  void plot_bgr_cosmics(float pmin = 100., float pmax = 110., float tmin = 500., float tmax = 1700.) {

    TF2* bgr_cosmics_2D = new TF2("bgr_cosmics_2D",f_bgr_cosmics,pmin, pmax, tmin, tmax);

    bgr_cosmics_2D->Draw();
  }

  //-----------------------------------------------------------------------------
  // constructor that doesn't initialize the histograms
  //-----------------------------------------------------------------------------
  channel::channel(const char* ChannelName, int Mode, int Verbose) : TNamed(ChannelName,ChannelName),
                                                                     fVerbose(Verbose) {
    const std::string name = ChannelName;
    const int isMu2eIIDataset = (name == "DIO" || name == "CE") * (Mode % 1000) / 100;
    _constants = constants(Mode, isMu2eIIDataset);
    fHist      = nullptr;

    kLumiSF_1B = _constants.npot_1b();
    kLumiSF_2B = _constants.npot_2b();

    const char* hist_dir = gEnv->GetValue("mu2e.HistDir","/projects/mu2e/hist");
    const char* project  = gEnv->GetValue("mu2e.HistProject","su2020");
    fHistDir = Form("%s/%s",hist_dir, project);
  };


  //-----------------------------------------------------------------------------
  channel::channel(const char* ChannelName, float ExtraSF, int Mode, int Verbose) : channel(ChannelName, Mode, Verbose) {
    // HistName: generic: "time", "mom", "time_vs_mom"
    // allow additional scaling of the histograms, by default, ExtraSF = 1

    const float PMin =  103.60;
    const float PMax =  104.90;
    const float TMin =  640.00;
    const float TMax = 1650.00;

    if(fVerbose > 0)
      printf("ChannelName, SF1B, SF2B, ExtraSF: %-10s  %9.3e %9.3e %9.3e\n",ChannelName,kLumiSF_1B,kLumiSF_2B,ExtraSF);

    const double NPOT_1B = _constants.npot_1b();
    const double NPOT_2B = _constants.npot_2b();

    TString channel_name = ChannelName;

    if (channel_name == "DIO" && (Mode % 1000) / 100 == 0) { //SU2020 histogram
      const char* dsid    = "su2020.fele2s51b1";
      const char* ana_job = "su2020_track_ana.1011";
      double      ngen (1.e7), erange(10.);            // properties of the dataset

      double sf1b = NPOT_1B*_constants.muon_stop_rate()*(1.-_constants.muon_capture())*erange/ngen*ExtraSF;
      double sf2b = NPOT_2B*_constants.muon_stop_rate()*(1.-_constants.muon_capture())*erange/ngen*ExtraSF;

      if(fVerbose > 0)
        printf("  DIO : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2006/p_vs_time")->Clone("DIO_t_vs_p");
      double sum1 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2007/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);

      if(fVerbose > 0)
        printf("  DIO      : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);

      delete h;
    }
    else if (channel_name == "DIO" && (Mode % 1000) / 100 == 1) { //Mu2e-II histogram
      const char* dsid    = "mu2eiisnowmass.dio";
      const char* ana_job = "trkanahist.1011";
      double      ngen (1.e7), erange(104.97 - 100.), frac_tail(3.4456778e-13); // properties of the dataset
      //frac_tail is the integral of the DIO LL spectrum on Al from 100 to 104.97 MeV

      double sf1b = NPOT_1B*_constants.muon_stop_rate()*(1.-_constants.muon_capture())*frac_tail*erange/ngen*ExtraSF;
      double sf2b = NPOT_2B*_constants.muon_stop_rate()*(1.-_constants.muon_capture())*frac_tail*erange/ngen*ExtraSF;

      if(fVerbose > 0)
        printf("  DIO : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

      // fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2006/p_vs_time")->Clone("DIO_t_vs_p");
      TFile* f = TFile::Open(Form("%s/%s.%s.hist", GetHistDir(), dsid, ana_job));
      if(!f) throw 20;
      fTimeVsMom = (TH2F*) f->Get("DIO/trk_100/p_vs_t0")->Clone("DIO_t_vs_p");
      double sum1 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      TH2F* h = (TH2F*) f->Get("DIO/trk_100/p_vs_t0")->Clone("tmp"); //FIXME: Add two-batch DIO histogram
      // TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2007/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);

      if(fVerbose > 0)
        printf("  DIO      : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);

      delete h;
    }
    else if (channel_name == "CE" && (Mode % 1000) / 100 == 0) { //SU2020 histogram
      const char* dsid    = "su2020.cele0s61b1";
      const char* ana_job = "su2020_track_ana.1011";
      double      ngen (1.e6);

      double sf1b = NPOT_1B*_constants.muon_stop_rate()*_constants.muon_capture()/ngen*ExtraSF;
      double sf2b = NPOT_2B*_constants.muon_stop_rate()*_constants.muon_capture()/ngen*ExtraSF;

      if(fVerbose > 0)
        printf("  CE : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2004/p_vs_time")->Clone("CE_t_vs_p");
      double sum1 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2005/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);

      if(fVerbose > 0)
        printf("  CE       : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);
      delete h;
    }
    else if (channel_name == "CE" && (Mode % 1000) / 100 == 1) { //Mu2e-II histogram
      const char* dsid    = "mu2eiisnowmass.ce";
      const char* ana_job = "trkanahist.1011";
      double      ngen (1.e6);

      double sf1b = NPOT_1B*_constants.muon_stop_rate()*_constants.muon_capture()/ngen*ExtraSF;
      double sf2b = NPOT_2B*_constants.muon_stop_rate()*_constants.muon_capture()/ngen*ExtraSF;

      if(fVerbose > 0)
        printf("  CE : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

      // fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2004/p_vs_time")->Clone("CE_t_vs_p");
      TFile* f = TFile::Open(Form("%s/%s.%s.hist", GetHistDir(), dsid, ana_job));
      if(!f) throw 20;
      fTimeVsMom = (TH2F*) f->Get("CE/trk_100/p_vs_t0")->Clone("CE_t_vs_p");
      double sum1 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      TH2F* h = (TH2F*) f->Get("CE/trk_100/p_vs_t0")->Clone("tmp"); //FIXME: Add two-batch CE histogram
      // TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2005/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);

      if(fVerbose > 0)
        printf("  CE       : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);
      delete h;
    }
    else if (channel_name == "PbarAnni") {
      //-----------------------------------------------------------------------------
      // instead of the luminosity re-weighing, introduce a CE-based efficiency scale
      // factor of 0.98
      //-----------------------------------------------------------------------------
      const char* dsid    = "su2020.pbar0s61b0";
      const char* ana_job = "su2020_pbar_ana.1010";

      double ngen       (2.e7);         // number of pbars generated at ST
      double nsa_per_pot(4.84e-18);     // number of pbars stopped in ST per POT

      double sf1b = NPOT_1B*nsa_per_pot/ngen*ExtraSF;
      double sf2b = NPOT_2B*nsa_per_pot/ngen*ExtraSF*0.98;

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("PbarAnni_t_vs_p");
      double sum1 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);

      if(fVerbose > 0) {
        printf("  PbarAnni : NPOT_1B, NPOT_2B, sf1b, sf2b = %12.5e %12.5e %12.5e %12.5e\n",NPOT_1B,NPOT_2B,sf1b,sf2b);
        printf("  PbarAnni : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);
      }
      delete h;
    }
    else if (channel_name == "PbarRPCe") {
      const char* dsid    = "su2020.rpce3s41b0";
      const char* ana_job = "su2020_pbar_ana.1010";

      double ngen(5.e7);                 // number of generated RPCe events
      double nsp_per_pot(1.87e-16);      // number of stopped pi^- produced by pbars, per POT
      double br_rpc     (2.15e-2);       // BR(RPC)

      double sf1b = NPOT_1B*nsp_per_pot/ngen*br_rpc*ExtraSF;
      double sf2b = NPOT_2B*nsp_per_pot/ngen*br_rpc*ExtraSF*0.98;

      if(fVerbose > 0)
        printf("  PbarRPCe : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("PbarRPCe_t_vs_p");
      double sum1 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      if(fVerbose > 0)
        printf("  PbarRPCe : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);
      delete h;
    }
    else if (channel_name == "PbarRPCi") {
      const char* dsid    = "su2020.rpci3s41b0";
      const char* ana_job = "su2020_pbar_ana.1010";

      double ngen(5.e5);                 // number of generated RPCe events
      double nsp_per_pot(1.87e-16);      // number of stopped pi^- produced by pbars, per POT
      double br_rpc     (2.15e-2);       // BR(RPC)
      double rho        (6.94e-3);       // internal conversion fraction, RPC

      double sf1b = NPOT_1B*nsp_per_pot/ngen*br_rpc*rho*ExtraSF;
      double sf2b = NPOT_2B*nsp_per_pot/ngen*br_rpc*rho*ExtraSF*0.98;

      if(fVerbose > 0)
        printf("  PbarRPCi : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("PbarRPCi_t_vs_p");
      double sum1 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      if(fVerbose > 0)
        printf("  PbarRPCi : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);

      delete h;
    }
    else if (channel_name == "PbarRPC") {
      //-----------------------------------------------------------------------------
      // sum of internal and external , just assume gPbarRPCe and gPbarRPCi are already initialized
      //-----------------------------------------------------------------------------
      if(fVerbose > 0)
        printf("\n");
      channel pbar_rpci("PbarRPCi",ExtraSF,Mode,fVerbose);
      channel pbar_rpce("PbarRPCe",ExtraSF,Mode,fVerbose);
      fTimeVsMom = (TH2F*) pbar_rpci.fTimeVsMom->Clone("PbarRPC_t_vs_p");
      fTimeVsMom->Add(pbar_rpce.fTimeVsMom);
    }
    else if (channel_name == "PbarTOT") {
      //-----------------------------------------------------------------------------
      // sum of internal and external , just assume gPbarRPCe and gPbarRPCi are already initialized
      //-----------------------------------------------------------------------------
      if(fVerbose > 0)
        printf("\n");
      channel pbar_anni("PbarAnni",ExtraSF,Mode,fVerbose);
      channel pbar_rpc ("PbarRPC" ,ExtraSF,Mode,fVerbose);
      fTimeVsMom = (TH2F*) pbar_anni.fTimeVsMom->Clone("PbarTOT_t_vs_p");
      fTimeVsMom->Add(pbar_rpc.fTimeVsMom);
    }
    else if (channel_name == "Cosmics") {
      //-----------------------------------------------------------------------------
      // kludge - initialize using a parameterization
      //-----------------------------------------------------------------------------
      double x[2], p[2];

      const char* dsid    = "su2020.cry33s51b0";
      const char* ana_job = "su2020_cosmic_ana.1010";

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_CosmicAna","trk_2000/p_vs_time")->Clone("Cosmics_t_vs_p");
      fTimeVsMom->Reset();

      double xmin = fTimeVsMom->GetXaxis()->GetXmin();
      double binx = fTimeVsMom->GetXaxis()->GetBinWidth(1);

      double ymin = fTimeVsMom->GetYaxis()->GetXmin();
      double biny = fTimeVsMom->GetYaxis()->GetBinWidth(1);

      if(fVerbose > 0)
        printf("  Cosmics: binx = %10.4f biny = %10.4f\n",binx,biny);

      int nbx = fTimeVsMom->GetNbinsX();
      int nby = fTimeVsMom->GetNbinsY();
      for(int ix=1; ix<=nbx; ix++) {
        x[0] = xmin+(ix-0.5)*binx;
        for(int iy=1; iy<=nby; iy++) {
          x[1] = ymin+(iy-0.5)*biny;
          //-----------------------------------------------------------------------------
          // f_bgr_cosmics returns normalized background
          //-----------------------------------------------------------------------------
          double bgr = f_bgr_cosmics(x,p)*binx*biny*ExtraSF;
          fTimeVsMom->SetBinContent(ix,iy,bgr);
          fTimeVsMom->SetBinError  (ix,iy,bgr*0.1);
        }
      }
      //-----------------------------------------------------------------------------
      // for cosmics, determine only one number - normalized to the time...
      //----------------------------------------------------------------------------------------
      double sum3 = GetIntegral(PMin+1.e-3, PMax-1.e-3, TMin+1.e-3, TMax-1.e-3);
      if(fVerbose > 0)
        printf("  Cosmics : BGR(total): %12.5e\n",sum3);
    }
  }

  //-----------------------------------------------------------------------------
  // create momentum histogram with given  timing cutoff
  //-----------------------------------------------------------------------------
  TH1D* channel::CreateMomHist(double TMin, double TMax) {

    if(!fTimeVsMom || !fTimeVsMom->GetYaxis()) {
      std::cout << __func__ << ": Undefined TimeVsMom histogram/axis!\n";
      return nullptr;
    }
    double tmin = fTimeVsMom->GetYaxis()->GetXmin();
    double tmax = fTimeVsMom->GetYaxis()->GetXmax();
    double bin  = fTimeVsMom->GetYaxis()->GetBinWidth(1);

    int iy1(0), iy2(-1); // defaults

    if (TMin >  tmin) iy1 = (TMin+1.e-6*bin - tmin)/bin + 1;
    if (TMax <= tmax) iy2 = (TMax-1.e-6*bin - tmin)/bin + 1;
    if(fVerbose > 0)
      printf("iy1, iy2 = %5i %5i \n",iy1,iy2);

    TH1D* htemp   = fTimeVsMom->ProjectionX("htemp",iy1,iy2);

    TString name  = Form("h_%s_mom_hist",GetName());
    TH1D* hist    = (TH1D*) htemp->Clone(name.Data());
    delete htemp;

    return hist;
  }

  //-----------------------------------------------------------------------------
  // create momentum histogram with given  timing cutoff
  //-----------------------------------------------------------------------------
  TH1D* channel::CreateTimeHist(double PMin, double PMax) {

    if(!fTimeVsMom || !fTimeVsMom->GetXaxis()) {
      std::cout << __func__ << ": Undefined TimeVsMom histogram/axis!\n";
      return nullptr;
    }
    double pmin = fTimeVsMom->GetXaxis()->GetXmin();
    double pmax = fTimeVsMom->GetXaxis()->GetXmax();
    double bin  = fTimeVsMom->GetXaxis()->GetBinWidth(1);

    int ix1(0), ix2(-1); // defaults

    if(fVerbose > 0)
      printf ("PMin, PMax, pmin, pmax, bin = %10.3f  %10.3f  %10.3f  %10.3f %10.4f\n",
              PMin, PMax, pmin, pmax,bin);

    if (PMin >  pmin ) ix1 = (PMin+1.e-6*bin - pmin)/bin + 1;
    if (PMax <= pmax ) ix2 = (PMax-1.e-6*bin - pmin)/bin + 1;

    if(fVerbose > 0)
      printf("ix1, ix2 = %5i %5i \n",ix1,ix2);

    TH1D* htemp  = fTimeVsMom->ProjectionY("htemp",ix1,ix2);
    TString name = Form("h_%s_mom_hist",GetName());
    TH1D* hist   = (TH1D*) htemp->Clone(name);
    delete htemp;

    return hist;
  }

  //-----------------------------------------------------------------------------
  // create running integral (T_i, TMax) histogram with given timing cutoff for
  // the time distribution within given (PMin,PMax) momentum band
  //-----------------------------------------------------------------------------
  TH1D* channel::CreateTimeIntegralHist(double PMin, double PMax, double TMin, double TMax) {

    if(!fTimeVsMom || !fTimeVsMom->GetXaxis()) {
      std::cout << __func__ << ": Undefined TimeVsMom histogram/axis!\n";
      return nullptr;
    }
    double pmin = fTimeVsMom->GetXaxis()->GetXmin();
    double pmax = fTimeVsMom->GetXaxis()->GetXmax();
    double binp = fTimeVsMom->GetXaxis()->GetBinWidth(1);

    double tmin = fTimeVsMom->GetYaxis()->GetXmin();
    double tmax = fTimeVsMom->GetYaxis()->GetXmax();
    double bint = fTimeVsMom->GetYaxis()->GetBinWidth(1);
    int    nbt  = fTimeVsMom->GetYaxis()->GetNbins();

    int ix1(0), ix2(-1), iy1(0), iy2(nbt); // defaults

    if(fVerbose > 0)
      printf ("PMin, PMax, pmin, pmax, binp = %10.3f  %10.3f  %10.3f  %10.3f %10.4f\n",
              PMin, PMax, pmin, pmax,binp);

    if (PMin > pmin) ix1 = (PMin+1.e-6*binp - pmin)/binp + 1;
    if (PMax < pmax) ix2 = (PMax-1.e-6*binp - pmin)/binp + 1;

    if (TMin > tmin) iy1 = (TMin+1.e-6*bint - tmin)/bint + 1;
    if (TMax < tmax) iy2 = (TMax-1.e-6*bint - tmin)/bint + 1;

    if(fVerbose > 0) {
      printf("ix1, ix2 = %5i %5i \n",ix1,ix2);
      printf("iy1, iy2 = %5i %5i \n",iy1,iy2);
    }

    TH1D* htemp  = fTimeVsMom->ProjectionY("htemp",ix1,ix2);
    TString name = Form("h_%s_time_integral_hist",GetName());
    TH1D* hist   = (TH1D*) htemp->Clone(name);

    hist->Reset();

    for (int iy = iy1; iy<=iy2; iy++) {
      double sy  = 0;
      double sw = 0;
      for (int iyy = iy; iyy<=iy2; iyy++) {
        sy += htemp->GetBinContent(iyy);
        sw += htemp->GetBinError(iyy)*fTimeVsMom->GetBinError(iyy);
        //        printf("iyy, y, sy, err, sw : %5i %12.5e  %12.5e  %12.5e  %12.5e\n",
        //     iyy, fTimeVsMom->GetBinContent(iyy),sy,fTimeVsMom->GetBinError(iyy), sw);
      }
      hist->SetBinContent(iy,sy);
      hist->SetBinError  (iy,sqrt(sw));
      // printf("--------- iy, y, sw : %5i %12.5e  %12.5e\n",iy,sy,sw);
    }

    // delete htemp;
    return hist;
  }

  //-----------------------------------------------------------------------------
  // create running integral (T_i, TMax) histogram with given timing cutoff for
  // the time distribution within given (PMin,PMax) momentum band
  //-----------------------------------------------------------------------------
  TH1D* channel::CreateMomIntegralHist(double PMin, double PMax, double TMin, double TMax) {

    if(!fTimeVsMom || !fTimeVsMom->GetXaxis()) {
      std::cout << __func__ << ": Undefined TimeVsMom histogram/axis!\n";
      return nullptr;
    }
    double pmin = fTimeVsMom->GetXaxis()->GetXmin();
    double pmax = fTimeVsMom->GetXaxis()->GetXmax();
    double binp = fTimeVsMom->GetXaxis()->GetBinWidth(1);
    //  int    nbp  = fTimeVsMom->GetXaxis()->GetNbins();

    double tmin = fTimeVsMom->GetYaxis()->GetXmin();
    double tmax = fTimeVsMom->GetYaxis()->GetXmax();
    double bint = fTimeVsMom->GetYaxis()->GetBinWidth(1);
    int    nbt  = fTimeVsMom->GetYaxis()->GetNbins();

    int ix1(0), ix2(-1), iy1(0), iy2(nbt); // defaults

    if(fVerbose > 0)
      printf ("PMin, PMax, pmin, pmax, binp = %10.3f  %10.3f  %10.3f  %10.3f %10.4f\n",
              PMin, PMax, pmin, pmax,binp);

    if (PMin > pmin) ix1 = (PMin+1.e-6*binp - pmin)/binp + 1;
    if (PMax < pmax) ix2 = (PMax-1.e-6*binp - pmin)/binp + 1;

    if (TMin > tmin) iy1 = (TMin+1.e-6*bint - tmin)/bint + 1;
    if (TMax < tmax) iy2 = (TMax-1.e-6*bint - tmin)/bint + 1;

    if(fVerbose > 0) {
      printf("ix1, ix2 = %5i %5i \n",ix1,ix2);
      printf("iy1, iy2 = %5i %5i \n",iy1,iy2);
    }

    TH1D*   htemp = fTimeVsMom->ProjectionX("htemp",iy1,iy2);
    TString name  = Form("h_%s_mom_integral_hist",GetName());
    TH1D*   hist  = (TH1D*) htemp->Clone(name);
    hist->Reset();

    for (int ix = ix1; ix<=ix2; ix++) {
      double sx  = 0;
      double sw = 0;
      for (int ixx = ix; ixx<=ix2; ixx++) {
        sx += htemp->GetBinContent(ixx);
        sw += htemp->GetBinError(ixx)*fTimeVsMom->GetBinError(ixx);
        //        printf("iyy, y, sy, err, sw : %5i %12.5e  %12.5e  %12.5e  %12.5e\n",
        //     iyy, fTimeVsMom->GetBinContent(iyy),sy,fTimeVsMom->GetBinError(iyy), sw);
      }
      hist->SetBinContent(ix,sx);
      hist->SetBinError  (ix,sqrt(sw));
      // printf("--------- iy, y, sw : %5i %12.5e  %12.5e\n",iy,sy,sw);
    }

    delete htemp;
    return hist;
  }

  //-----------------------------------------------------------------------------
  double channel::GetIntegral(float PMin, float PMax, float TMin, float TMax) {
    double error;
    return GetIntegral(PMin, PMax, TMin, TMax, error);
  }

  //-----------------------------------------------------------------------------
  double channel::GetIntegral(float PMin, float PMax, float TMin, float TMax, double& error) {

    if(!fTimeVsMom || !fTimeVsMom->GetXaxis()) {
      std::cout << __func__ << ": Undefined TimeVsMom histogram/axis!\n";
      return 0.;
    }
    TAxis* ax = fTimeVsMom->GetXaxis();
    TAxis* ay = fTimeVsMom->GetYaxis();
    //push the bounds to ensure they're within the correct bin
    PMin += 1.e-3;
    PMax -= 1.e-3;
    TMin += 1.e-3;
    TMax -= 1.e-3;
    return fTimeVsMom->IntegralAndError(ax->FindBin(PMin), ax->FindBin(PMax), ay->FindBin(TMin), ay->FindBin(TMax), error);

    // int nbp   = ax->GetNbins();
    // int nbt   = ay->GetNbins();

    // float sw  = 0;

    // for (int ix=0; ix<nbp; ix++) {
    //   float p = ax->GetBinCenter(ix);
    //   if ((p < PMin) or (p > PMax))   continue;
    //   for (int iy=0; iy<nbt; iy++) {
    //     float t = ay->GetBinCenter(iy);
    //     if ((t < TMin) or (t > TMax))   continue;
    //     float w = fTimeVsMom->GetBinContent(ix,iy);
    //     sw += w;
    //   }
    // }

    // return sw;
  }
};
