///////////////////////////////////////////////////////////////////////////////
// use numbers from
// to be compiled
// dsid[2] and dsid[3] - out of time RPC
///////////////////////////////////////////////////////////////////////////////
#include "mu2eii/stat/stat/channel_rpc.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "TEnv.h"
#include "TF2.h"

namespace mu2eii {

  channel_rpc::channel_rpc(const char* ChannelName, float ExtraSF, int Mode, int DatasetConfigMode, int Verbose) :
    channel(ChannelName,Mode,Verbose)
  {
    const char* dsid_0 [] = { "su2020.rpce0s51b1.su2020_rpc_ana.1011.hist",
                              "su2020.rpci0s51b1.su2020_rpc_ana.1011.hist",
                              "su2020.rpce1s61b0.su2020_rpc_ana.1010.hist",
                              "su2020.rpci1s61b0.su2020_rpc_ana.1010.hist" };

    const char* dsid_1 [] = { "su2020.rpce0s51b1.su2020_rpc_ana.1011.wo_3_cuts.hist",
                              "su2020.rpci0s51b1.su2020_rpc_ana.1011.wo_3_cuts.hist",
                              "su2020.rpce1s61b0.su2020_rpc_ana.1010.wo_3_cuts.hist",
                              "su2020.rpci1s61b0.su2020_rpc_ana.1010.wo_3_cuts.hist" };

    // allow additional scaling of the histograms, by default, ExtraSF = 1

    TString channel_name = ChannelName;

    const char** dsid(nullptr);

    if      (DatasetConfigMode == 0) dsid = dsid_0;
    else if (DatasetConfigMode == 1) dsid = dsid_1;
    else {
      printf("dont know what to do: DatasetConfigMode = %i\n",DatasetConfigMode);
    }

    const double NPOT_1B = _constants.npot_1b();
    const double NPOT_2B = _constants.npot_2b();

    if (channel_name == "RPCe") {
      //-----------------------------------------------------------------------------
      // external in-time RPC (HisfFN[0]), generation parameters
      // ($npot/$npot_gen)*($nstops/$eff_ngen)*$br_rpc*$fr
      //-----------------------------------------------------------------------------
      double npot_gen(1.e8);          // N(POT) for pi^- beam transport
      double eff_450 (2.423e-2);      // efficiency of the T>450 ns cut-off

      double nstops  (210551);        // N(stopped pi^-)

      double ngen_ph(1.e8);           // N(generated photons), from ST
      double br_rpc (2.15e-2);        // BR(RPC)

      double eff_ngen = ngen_ph/eff_450; // effective number of generated photons

      double sf1b = (NPOT_1B/npot_gen)*(nstops/eff_ngen)*br_rpc*ExtraSF;
      double sf2b = (NPOT_2B/npot_gen)*(nstops/eff_ngen)*br_rpc*ExtraSF*0.98; // additional efficiency scale factor
      if(fVerbose > 0)
        printf("RPCe : sf1b, sf2b, dsid = %12.5e %12.5e %s\n",sf1b,sf2b,dsid[0]);

      fTimeVsMom  = (TH2F*) gh2(Form("%s/%s",GetHistDir(),dsid[0]),"su2020_RPCAna","trk_2004/p_vs_time")->Clone("RPCe_t_vs_p");
      double sum1 = GetIntegral(103.85,105.1,700,1700);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(103.85,105.1,700,1700);
      //-----------------------------------------------------------------------------
      // for RPC in-time datasets, because of the limited statistics, the luminosity
      // reweighing produces a biased (underestimated) answer
      // use the efficiency scale factor (0.98), derived from the 1B and 2B conversion
      // electron (CE) datasets instead trk_2005 --> trk_2004 * 0.98
      //-----------------------------------------------------------------------------
      //    TH2F* h     = (TH2F*) gh2(Form("%s/%s",GetHistDir(),dsid[0]),"su2020_RPCAna","trk_2005/p_vs_time")->Clone("tmp");
      TH2F* h     = (TH2F*) gh2(Form("%s/%s",GetHistDir(),dsid[0]),"su2020_RPCAna","trk_2004/p_vs_time")->Clone("tmp");

      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(103.85,105.1,700,1700);
      if(fVerbose > 0)
        printf("sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);

      delete h;
    }
    else if (channel_name == "RPCi") {
      //-----------------------------------------------------------------------------
      // internal in-time RPC: HistFn[1], generation parameters
      //-----------------------------------------------------------------------------
      double npot_gen(1.e8);          // N(POT) for pi^- beam transport
      double eff_450 (2.423e-2);      // efficiency of the T>450 ns cut-off
      double nstops  (210551);        // N(stopped pi^-)

      double ngen_ph (1.e6);          // N(generated photons), from ST
      double rho     (6.9e-3);        // internal conversion fraction
      double br_rpc  (2.15e-2);       // BR(RPC)

      double eff_ngen = ngen_ph/eff_450; // effective number of generated photons

      double sf1b = (NPOT_1B/npot_gen)*(nstops/eff_ngen)*br_rpc*rho*ExtraSF;
      double sf2b = (NPOT_2B/npot_gen)*(nstops/eff_ngen)*br_rpc*rho*ExtraSF;

      if(fVerbose > 0)
        printf("RPCi : sf1b, sf2b,dsid = %12.5e %12.5e %s\n",sf1b,sf2b,dsid[1]);

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s",GetHistDir(),dsid[1]),"su2020_RPCAna","trk_2004/p_vs_time")->Clone("RPCi_t_vs_p");
      double sum1 = GetIntegral(103.85,105.1,700,1700);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(103.85,105.1,700,1700);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s",GetHistDir(),dsid[1]),"su2020_RPCAna","trk_2005/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(103.85,105.1,700,1700);
      if(fVerbose > 0)
        printf("sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);

      delete h;
    }
    else if (channel_name == "RPC") {
      //-----------------------------------------------------------------------------
      // total: sum of internal and external
      //-----------------------------------------------------------------------------
      if(fVerbose > 0)
        printf("\n");
      channel_rpc rpce("RPCe",ExtraSF,Mode,DatasetConfigMode);
      channel_rpc rpci("RPCi",ExtraSF,Mode,DatasetConfigMode);

      fTimeVsMom = (TH2F*) rpci.fTimeVsMom->Clone("RPC_t_vs_p");
      fTimeVsMom->Add(rpce.fTimeVsMom);
    }
    else if (channel_name == "RPCe_OOT") {
      //-----------------------------------------------------------------------------
      // factor of 0.98
      //-----------------------------------------------------------------------------
      const char* dsid    = "su2020.rpce1s61b0";
      const char* ana_job = "su2020_rpc_ana.1010";

      double npot_gen(1.e8);          // N(POT) for pi^- beam transport
      double nstops  (210551);        // N(stopped pi^-)

      double ngen_ph (1.e8);         // number of pbars generated at ST
      double br_rpc (2.15e-2);        // BR(RPC)

      double ext_fact(_constants.extinction());        // proton extinction factor

      double sf1b = NPOT_1B*(nstops/npot_gen)*(1./ngen_ph)*ext_fact*br_rpc*ExtraSF;
      double sf2b = NPOT_2B*(nstops/npot_gen)*(1./ngen_ph)*ext_fact*br_rpc*ExtraSF*0.98; // additional efficiency scale factor

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2004/p_vs_time")->Clone("RPCe_OOT_t_vs_p");
      double sum1 = GetIntegral(103.85,105.1,700,1700);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(103.85,105.1,700,1700);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2004/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(103.85,105.1,700,1700);

      if(fVerbose > 0) {
        printf("RPCe_OOT : NPOT_1B, NPOT_2B, sf1b, sf2b = %12.5e %12.5e %12.5e %12.5e\n",NPOT_1B,NPOT_2B,sf1b,sf2b);
        printf("RPCe_OOT : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);
      }
      delete h;
    }
    else if (channel_name == "RPCi_OOT") {
      //-----------------------------------------------------------------------------
      // factor of 0.98
      //-----------------------------------------------------------------------------
      const char* dsid    = "su2020.rpci1s61b0";
      const char* ana_job = "su2020_rpc_ana.1010";

      double npot_gen(1.e8);          // N(POT) for pi^- beam transport
      double nstops  (210551);        // N(stopped pi^-)

      double ngen_ph (1.e6);         // number of pbars generated at ST
      double rho     (6.9e-3);        // internal conversion fraction
      double br_rpc (2.15e-2);        // BR(RPC)

      double ext_fact(_constants.extinction());        // proton extinction factor

      double sf1b = NPOT_1B*(nstops/npot_gen)*(1./ngen_ph)*ext_fact*br_rpc*rho*ExtraSF;
      double sf2b = NPOT_2B*(nstops/npot_gen)*(1./ngen_ph)*ext_fact*br_rpc*rho*ExtraSF*0.98; // additional efficiency scale factor

      fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2004/p_vs_time")->Clone("RPCi_OOT_t_vs_p");
      double sum1 = GetIntegral(103.85,105.1,700,1700);
      fTimeVsMom->Scale(sf1b);
      double sum2 = GetIntegral(103.85,105.1,700,1700);
      TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2004/p_vs_time")->Clone("tmp");
      fTimeVsMom->Add(h,sf2b);
      double sum3 = GetIntegral(103.85,105.1,700,1700);

      if(fVerbose > 0) {
        printf("RPCi_OOT : NPOT_1B, NPOT_2B, sf1b, sf2b = %12.5e %12.5e %12.5e %12.5e\n",NPOT_1B,NPOT_2B,sf1b,sf2b);
        printf("RPCi_OOT : sw(signal), BGR(one-batch), BGR(2-batch), BGR(total): %12.5e %12.5e %12.5e %12.5e\n",sum1,sum2,sum3-sum2,sum3);
      }
      delete h;
    }
    else if (channel_name == "RPC_OOT") {
      //-----------------------------------------------------------------------------
      // total: sum of internal and external
      //-----------------------------------------------------------------------------
      if(fVerbose > 0)
        printf("\n");
      channel_rpc rpce_oot("RPCe_OOT",ExtraSF,Mode,DatasetConfigMode);
      channel_rpc rpci_oot("RPCi_OOT",ExtraSF,Mode,DatasetConfigMode);

      fTimeVsMom = (TH2F*) rpci_oot.fTimeVsMom->Clone("RPC_t_vs_p");
      fTimeVsMom->Add(rpce_oot.fTimeVsMom);
    }
  }

};
