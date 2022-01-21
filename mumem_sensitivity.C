///////////////////////////////////////////////////////////////////////////////
// use example:
// ------------
// mumem* m3 = new mumem(3)
// m3->scan_pmin(104.5,105.0,700.,1650)
//
// DatasetConfigCode : signal_code+cosmics_code+dio_code+pbar_code+rpc_code (10000)
// background dataset names are arranged alphabetically
// Mode = 1 : default benchmarking
// Mode = 2 : cosmics = 0.21
// Mode = 3 : default, the deadtime and the pulse intensity cut-offs accounted for
// Mode = 4 : cosmics = 0.21, the deadtime and the pulse intensity cut-offs accounted for
//
// channel initialization expects mu2e.HistDir to be defined in .rootrc and be pointing
// to the histogram area.
// On mu2egpvm* the official su2020 histogram area is located at
//                                  /mu2e/data/projects/su2020/hist/su2020
///////////////////////////////////////////////////////////////////////////////
// const char* FiguresDir         = "/projects/mu2e/notes/mu2e-37812.su2020_sens_mumem/figures";
// const char* su2020_HistDir     = "/projects/mu2e/hist/su2020" ;
#ifndef __mu2eii_stat_mumem_sensitivity_C__
#define __mu2eii_stat_mumem_sensitivity_C__

#include "Stntuple/scripts/stn_catalog.hh"

#include "Stntuple/scripts/plot_data.hh"
#include "Stntuple/scripts/plot_hist_1D.C"
#include "Stntuple/scripts/plot_hist_2D.C"
#include "Stntuple/alg/TFeldmanCousinsA.hh"

#include "mu2eii/stat/stat/channel_rpc.hh"
#include "mu2eii/stat/stat/analysis.hh"
#include "mu2eii/stat/stat/constants.hh"

class mumem : public mu2eii::analysis {
public:
  double  fSFSignal;            // try to rescale to account for tan(dip) cut
  double  fSFCosmics;           //
  double  fPMin;                // signal window for plotting histograms
  double  fPMax;
  double  fTMin;
  double  fTMax;
  //-----------------------------------------------------------------------------
  // constructors and destructor
  //-----------------------------------------------------------------------------
  mumem(int Mode = 1, int DatasetConfigCode = 0, const char* Name = "mu-e-", const char* Title = "mu- --> e- conversion");
  //-----------------------------------------------------------------------------
  // other functions
  //-----------------------------------------------------------------------------
  int  init_channels(int Mode, int DatasetConfigCode);
  void plot          (int Figure, int Plot = -1);
  int  scan_pmin     (double P0=103.85, double PMax=105.0, double TMin=700, double TMax = 1700, int nsteps = 20);
  int  scan_tmin     (double P0=103.85, double PMax=105.0, double TMin=700, double TMax = 1700, int nsteps = 20);
  int  scan_tmax     (double P0=103.85, double PMax=105.0, double TMin=700, double TMax = 1700, int nsteps = 20);
};

// constructor
mumem::mumem(int Mode, int DatasetConfigCode, const char* Name, const char* Title) : mu2eii::analysis(Name,Title) {
  fSFSignal = 1.; fSFCosmics = 1.;
  init_channels(Mode,DatasetConfigCode) ;
  fPMin = 103.85;
  fPMax = 105.10;
  fTMin =  700;
  fTMax = 1700;

  // plot_data_t::fgFiguresDir = ".";
}

//-----------------------------------------------------------------------------
// channel names - predefined
// 1. "CE",
// 2. "Cosmics",
// 3. "DIO",
// 4. "PbarTOT" = "PbarAnni" + "PbarRPCe" + "PbarPRCi",
// 4. "PbarRPC" = "PbarRPCe" + "PbarPRCi",
// 5. "RPC"     = "RPCe" + "RPCi"
//
// any change in the names needs to be propagated to the channel initialization code
// for now, the DatasetConfigCode is used only for RPC
//-----------------------------------------------------------------------------
int mumem::init_channels(int Mode, int DatasetConfigCode) {
  mu2eii::constants cnsts(Mode);
  const double npot_1b = cnsts.npot_1b();
  const double npot_2b = cnsts.npot_2b();
  const bool isMu2eII = (Mode % 100) / 10 > 0;
  fSFCosmics = cnsts.cosmics_scale();

  printf("%s: Effective N(POT) = %.3e (%.3e + %.3e)\n", __func__, npot_1b+npot_2b, npot_1b, npot_2b);

  //-----------------------------------------------------------------------------
  // 1. CE signal. For CE, specify signal in units of acceptance ?
  //-----------------------------------------------------------------------------
  // Rmue comes independently
  float Rmue = (isMu2eII) ? 1.e-16 : 1.e-15;
  mu2eii::channel* ce = new mu2eii::channel("CE",fSFSignal,Mode);
  SetSignal(ce,Rmue);
  //-----------------------------------------------------------------------------
  // 2. cosmics: fit , NPOT numbers are not used
  //-----------------------------------------------------------------------------
  mu2eii::channel* cosmics = new mu2eii::channel("Cosmics",fSFCosmics,Mode);
  AddBgrChannel(cosmics);
  //-----------------------------------------------------------------------------
  // 3. DIO
  //-----------------------------------------------------------------------------
  mu2eii::channel* dio = new mu2eii::channel("DIO",fSFSignal,Mode);
  AddBgrChannel(dio);
  //-----------------------------------------------------------------------------
  // 4. pbars
  //-----------------------------------------------------------------------------
  mu2eii::channel* pbar = new mu2eii::channel("PbarTOT",fSFSignal*cnsts.pbar_scale(),Mode);
  AddBgrChannel(pbar);

  //-----------------------------------------------------------------------------
  // 5. RPC - for fixed hist ID's, initialization of the RPC channel requires four datasets
  //-----------------------------------------------------------------------------
  int rpc_config_code = DatasetConfigCode % 10 ;   // the least significant digit

  printf("rpc_config_code = %i\n",rpc_config_code);

  mu2eii::channel* rpc(nullptr);
  if (rpc_config_code == 0)  rpc = new mu2eii::channel_rpc("RPC",fSFSignal*cnsts.rpc_scale(),Mode,rpc_config_code);
  if (rpc_config_code == 1)  rpc = new mu2eii::channel_rpc("RPC",fSFSignal*cnsts.rpc_scale(),Mode,rpc_config_code);

  AddBgrChannel(rpc);

  //-----------------------------------------------------------------------------
  // 6. RPC out-of-time: create channels, but don't add
  //-----------------------------------------------------------------------------
  printf("---------------------------------------------------------------\n");
  mu2eii::channel* rpc_oot(nullptr);
  rpc_oot = new mu2eii::channel_rpc("RPC_OOT",fSFSignal,Mode,rpc_config_code);
  return 0;
}

//-----------------------------------------------------------------------------
int mumem::scan_pmin(double PMin, double PMax, double TMin, double TMax, int nsteps) {
  //  init_channels();

  TFeldmanCousinsA fc("mumem::Print",-1);
  fc.SetNExp(100000); // enough for defining 5

  double smin(1.), smax(20.);  // range of tested signal values
  int npt(50);
  double s[npt], msign[npt];

  double qbgr[100];

  int nbgr = GetNBgrChannels();

  double pmax = PMax;

  printf("      pmin       pmax       tmin       tmax    %-10s     %-10s %-10s      %-10s Total        Signal         s5          SES         Rmue\n",
         GetBgrChannel(0)->GetName(),
         GetBgrChannel(1)->GetName(),
         GetBgrChannel(2)->GetName(),
         GetBgrChannel(3)->GetName()
         );

  printf("--------------------------------------------------------------------------------");
  printf("-------------------------------------------------------------------------------\n");

  for (int i=0; i<nsteps; i++) {
    double pmin = PMin-i*0.050;
    double sig = GetSigIntegral(pmin,pmax,TMin,TMax);
    double ses = GetSES        (pmin,pmax,TMin,TMax);

    double bgr_tot = 0;
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      mu2eii::channel* bgr = GetBgrChannel(ibgr);
      qbgr[ibgr] = bgr->GetIntegral(pmin,pmax,TMin,TMax);
      bgr_tot   += qbgr[ibgr] ;
    }
    //-----------------------------------------------------------------------------
    // for given bgr_tot, calculate the expected mean 5-sigma sensitivity
    // assume bgr_tot < 0.5, [1,11] covers the range in all cases
    //-----------------------------------------------------------------------------
    fc.DiscoveryProbMean(bgr_tot, smin, smax, 20, s, msign);

    double sigd(5.), s5;                    // signal coresponding to 1-sided "5-sigma"
    fc.SolveFor(sigd,s,msign,npt,&s5);
    printf("%10.3f %10.3f %10.3f %10.3f",pmin,PMax,TMin,TMax);
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      printf(" %12.5e", qbgr[ibgr]);
    }
    printf(" %12.5e %12.5e %12.5e %12.5e %12.5e\n",bgr_tot,sig,s5,ses,s5*ses);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int mumem::scan_tmin(double PMin, double PMax, double TMin,double TMax, int nsteps) {
  //  init_channels();

  TFeldmanCousinsA fc("mumem::Print",-1);
  fc.SetNExp(100000); // enough for defining 5

  double smin(1.), smax(20.);  // range of tested signal values
  int npt(50);
  double s[npt], msign[npt];

  double qbgr[100];

  int nbgr = GetNBgrChannels();

  printf("      pmin       pmax       tmin       tmax    %-10s     %-10s %-10s      %-10s Total        Signal         s5          SES         Rmue \n",
         GetBgrChannel(0)->GetName(),
         GetBgrChannel(1)->GetName(),
         GetBgrChannel(2)->GetName(),
         GetBgrChannel(3)->GetName()
         );

  printf("---------------------------------------------------------------------------------");
  printf("-------------------------------------------------------------------------------\n");

  for (int i=0; i<nsteps; i++) {
    double tmin = TMin-i*10;            // the bin width is 10 ns
    double sig = GetSigIntegral(PMin,PMax,tmin,TMax);
    double ses = GetSES        (PMin,PMax,tmin,TMax);

    double bgr_tot = 0;
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      mu2eii::channel* bgr = GetBgrChannel(ibgr);
      qbgr[ibgr] = bgr->GetIntegral(PMin,PMax,tmin,TMax);
      bgr_tot   += qbgr[ibgr] ;
    }
    //-----------------------------------------------------------------------------
    // for given bgr_tot, calculate the expected mean 5-sigma sensitivity
    // assume bgr_tot < 0.5, [1,11] covers the range in all cases
    //-----------------------------------------------------------------------------
    fc.DiscoveryProbMean(bgr_tot, smin, smax, npt, s, msign);

    double s5;                         // signal coresponding to 1-sided "5-sigma"
    fc.SolveFor(5,s,msign,npt,&s5);
    printf("%10.3f %10.3f %10.3f %10.3f",PMin,PMax,tmin,TMax);
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      printf(" %12.5e", qbgr[ibgr]);
    }
    printf(" %12.5e %12.5e %12.5e %12.5e %12.5e\n",bgr_tot,sig,s5,ses,s5*ses);
  }
  return 0;
}


//-----------------------------------------------------------------------------
int mumem::scan_tmax(double PMin, double PMax, double TMin,double TMax, int nsteps) {
  //  init_channels();

  TFeldmanCousinsA fc("mumem::Print",-1);
  fc.SetNExp(100000); // enough for defining 5

  double smin(1.), smax(20.);  // range of tested signal values
  int npt(50);
  double s[npt], msign[npt];

  double qbgr[100];

  int nbgr = GetNBgrChannels();

  printf("      pmin       pmax       tmin       tmax    %-10s     %-10s %-10s      %-10s Total        Signal         s5          SES         Rmue \n",
         GetBgrChannel(0)->GetName(),
         GetBgrChannel(1)->GetName(),
         GetBgrChannel(2)->GetName(),
         GetBgrChannel(3)->GetName()
         );

  printf("---------------------------------------------------------------------------------");
  printf("-------------------------------------------------------------------------------\n");

  for (int i=0; i<nsteps; i++) {
    double tmax = TMax-i*10;            // the bin width is 10 ns
    double sig = GetSigIntegral(PMin,PMax,TMin,tmax);
    double ses = GetSES        (PMin,PMax,TMin,tmax);

    double bgr_tot = 0;
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      mu2eii::channel* bgr = GetBgrChannel(ibgr);
      qbgr[ibgr] = bgr->GetIntegral(PMin,PMax,TMin,tmax);
      bgr_tot   += qbgr[ibgr] ;
    }
    //-----------------------------------------------------------------------------
    // for given bgr_tot, calculate the expected mean 5-sigma sensitivity
    // assume bgr_tot < 0.5, [1,11] covers the range in all cases
    //-----------------------------------------------------------------------------
    fc.DiscoveryProbMean(bgr_tot, smin, smax, npt, s, msign);

    double s5;                         // signal coresponding to 1-sided "5-sigma"
    fc.SolveFor(5,s,msign,npt,&s5);
    printf("%10.3f %10.3f %10.3f %10.3f",PMin,PMax,TMin,tmax);
    for (int ibgr=0; ibgr<nbgr; ibgr++) {
      printf(" %12.5e", qbgr[ibgr]);
    }
    printf(" %12.5e %12.5e %12.5e %12.5e %12.5e\n",bgr_tot,sig,s5,ses,s5*ses);
  }
  return 0;
}


//-----------------------------------------------------------------------------
void mumem::plot(int Figure, int Print = 0) {

  if (Figure == 1) {
    //-----------------------------------------------------------------------------
    // momentum distribution :
    //-----------------------------------------------------------------------------
    int nbgr = GetNBgrChannels();

    plot_data_t  p(nbgr+1);

    mu2eii::channel* ce  = GetSignal();

    TH1D* h1             = ce->CreateMomHist(fTMin,fTMax);
    h1->SetName("CE");

    h1->Scale(fRmue);
    p.hd[0]              = hist_data_t(h1);
    p.hd[0].fLabel       = "CE";
    p.hd[0].fLineColor   = kRed+1;
    p.hd[0].fMarkerColor = kRed+1;
    p.hd[0].fMarkerStyle = 20;
    p.hd[0].fMarkerSize  = 0.8;
    p.hd[0].fStats       = 0;

    int col[4] = { EColor::kBlue-1, EColor::kGreen+2, EColor::kMagenta, EColor::kRed+2 };

    for (int i=0; i<nbgr; i++) {
      int ih = i+1;                     // signal first

      mu2eii::channel* bgr   = GetBgrChannel(i);
      TH1D* h2               = bgr->CreateMomHist(fTMin,fTMax);
      h2->SetName(bgr->GetName());

      p.hd[ih]                = hist_data_t(h2);

      TString s               = bgr->GetName();
      if (s == "PbarTOT") s = "Pbar";
      p.hd[ih].fLabel         = s;

      p.hd[ih].fLineColor     = col[i];
      p.hd[ih].fMarkerColor   = col[i];
      p.hd[ih].fMarkerStyle   = 20;
      p.hd[ih].fMarkerSize    = 0.8;
      p.hd[ih].fStats         = 0;
    }
    //-----------------------------------------------------------------------------
    // common for the plot opt stat: 1110001
    // by default, assume that all histograms have the same binning,
    // so the rebinning factor can be specified just once per plot
    //-----------------------------------------------------------------------------
    // Y-scale and rebinning
    p.fYLogScale     = 0;
    p.fRebin         = 1;

    p.fXMin          = 102.5;
    p.fXMax          = 105.199;
    p.fXAxisTitle    = "e^{-} momentum, MeV/c";

    p.fYMin          = 1e-5;
    p.fYMax          = 0.30;

    p.fStatBoxYMin   = 0.83;         // up to 0.90, dY=0.07
    p.fOptStat       = 1000001;

    p.fLabel         = Form("SU2020: R = %12.5e, T0 in [%5.0f,%5.0f] ns",fRmue,fTMin,fTMax);
    p.fLabelFontSize = 0.04;

    p.fLegendXMin    = 0.25; p.fLegendXMax  = 0.40; p.fLegendYMin  = 0.50; p.fLegendYMax  = 0.70;

    p.fCanvasName    = Form("Figure_%04i",Figure);
    p.fName          = Form("figure_%05i_mumem_total_p",Figure);

    plot_hist_1d(&p,-1);

    TLegend* leg2 = new TLegend(0.35,0.50,0.50,0.70);
    leg2->SetBorderSize(0);
    leg2->AddEntry("",Form(": %5.3f",p.hd[0].fHist->Integral(451,504))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[1].fHist->Integral(451,504))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[2].fHist->Integral(451,504))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[3].fHist->Integral(451,504))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[4].fHist->Integral(451,504))," ");
    leg2->Draw();

    p.fCanvas->Modified();
    p.fCanvas->Update();

    if (Print == 1) p.print();
  }

  if (Figure == 2) {
    //-----------------------------------------------------------------------------
    // timing distribution : 103.5-105.1 MeV
    //-----------------------------------------------------------------------------
    int nbgr = GetNBgrChannels();

    plot_data_t  p(nbgr+1);

    mu2eii::channel* ce  = GetSignal();

    TH1D* h1             = ce->CreateTimeHist(fPMin,fPMax);
    h1->SetName("CE");

    h1->Scale(fRmue);
    p.hd[0]              = hist_data_t(h1);

    float sum1 = h1->Integral(51,165);  // 500-1650 ns

    p.hd[0].fLabel       = Form("CE");
    p.hd[0].fLineColor   = kRed+1;
    p.hd[0].fMarkerColor = kRed+1;
    p.hd[0].fMarkerStyle = 20;
    p.hd[0].fMarkerSize  = 0.8;
    p.hd[0].fStats       = 0;

    int col[4] = { EColor::kBlue-1, EColor::kGreen+2, EColor::kMagenta, EColor::kRed+2 };

    for (int i=0; i<nbgr; i++) {
      int ih = i+1;                     // signal first

      mu2eii::channel* bgr   = GetBgrChannel(i);
      TH1D* h2               = bgr->CreateTimeHist(fPMin,fPMax);

      float sum = h2->Integral(51,165);  // 500-1650 ns

      p.hd[ih]                = hist_data_t(h2);

      TString s               = bgr->GetName();
      if (s == "PbarTOT") s = "Pbar";
      p.hd[ih].fLabel         = s;

      p.hd[ih].fNewName       = s; //    h2->SetName(s.Data());

      p.hd[ih].fLineColor     = col[i];
      p.hd[ih].fMarkerColor   = col[i];
      p.hd[ih].fMarkerStyle   = 20;
      p.hd[ih].fMarkerSize    = 0.8;
      p.hd[ih].fStats         = 0;
    }
    //-----------------------------------------------------------------------------
    // common for the plot opt stat: 1110001
    // by default, assume that all histograms have the same binning,
    // so the rebinning factor can be specified just once per plot
    //-----------------------------------------------------------------------------
    // Y-scale and rebinning
    p.fYLogScale     = 1;
    p.fRebin         = 1;

    p.fXMin          = 500;
    p.fXMax          = 1650;
    p.fXAxisTitle    = "reconstucted T_{0}, ns";

    p.fYMin          = 1e-5;
    p.fYMax          = 1e3;

    p.fStatBoxYMin   = 0.83;         // up to 0.90, dY=0.07
    p.fOptStat       = 1000001;

    p.fLabel         = Form("SU2020: R = %12.5e, P in [%7.2f,%7.2f] ns",fRmue,fPMin,fPMax);
    p.fLabelFontSize = 0.04;

    p.fLegendXMin    = 0.25; p.fLegendXMax  = 0.40; p.fLegendYMin  = 0.50; p.fLegendYMax  = 0.70;

    p.fCanvasName    = Form("Figure_%04i",Figure);
    p.fName          = Form("figure_%05i_mumem_total_time",Figure);

    plot_hist_1d(&p,-1);

    TLegend* leg2 = new TLegend(0.35,0.50,0.50,0.70);
    leg2->SetBorderSize(0);
    leg2->AddEntry("",Form(": %5.3f",p.hd[0].fHist->Integral(51,165))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[1].fHist->Integral(51,165))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[2].fHist->Integral(51,165))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[3].fHist->Integral(51,165))," ");
    leg2->AddEntry("",Form(": %5.3f",p.hd[4].fHist->Integral(51,165))," ");
    leg2->Draw();

    p.fCanvas->Modified();
    p.fCanvas->Update();
    if (Print == 1) p.print();
  }
}
#endif
