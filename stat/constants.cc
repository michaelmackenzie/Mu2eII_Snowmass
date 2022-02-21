//
#include "mu2eii/stat/stat/constants.hh"
#include <cstdio>

namespace mu2eii {

  constants::constants() : _muon_capture(0.609), _muon_stop_rate(1.59e-3),
                           _npot_1b(2.86e19), _npot_2b(9.03e18), _cosmics_scale(1.),
                           _pbar_scale(1.), _rpc_scale(1.), _rpc_oot_scale(1.),_extinction(1.e-10),
                           _mixed_scale(1.), _unmixed_scale(1.) {
  }

  /**
     Mode definiton: 10,000 * use fit DIO PDF + 1000 * RPC OOT mode + 100 * use Mu2e-II datasets + 10 * Mu2e-II mode + 1 * Overall mode
     Mu2e-II mode:
     0 : Mu2e era values
     1 : use 4.4e22 POT and 2-batch mode SU2020 histograms
     2 : use 4.4e22 POT and 1-batch mode SU2020 histograms
     3 : Mu2e-II mode 1 + suppress RPC by a factor of 10
     4 : Mu2e-II mode 1 + suppress RPC to 0
     5 : Mu2e-II mode 1 with 3.6e22 POT (3-year run)
     6 : Mu2e-II mode 1 with 5/3*3.6e22 POT (5-year run)
     7 : Mu2e-II mode 1 with RPC suppressed by a factor of 34
     Overall Mode:
     0 : No additional scales
     3 : Apply deadtime, upstream rejection, luminosity cut-off, and trigger efficiency
     5 : 3 but with 3.6e20 POT, with Run 1 1-batch mode POT and the remaining in 2-batch mode
   **/
  constants::constants(int Mode, int isMu2eIIDataset) : constants() {
    float eff_ur_1b     (0.994);
    float eff_ur_2b     (0.986);
    float deadtime_1b   (1.67927e-02);
    float deadtime_2b   (7.62028e-02);
    float lumi_cutoff_1b(1-0.99342);
    float lumi_cutoff_2b(1-0.88785);
    float trig_eff      (0.98);

    const int mu2eIIMode = (Mode %  100) /  10;
    Mode %= 10;

    if(mu2eIIMode) { //Mu2e-II parameters
      _mixed_scale = (0.95/0.70)*(0.95/0.90); //remove the high inefficiency from the mixed samples
      _unmixed_scale = (0.95)*(0.95); //apply the trigger + selection pileup inefficiency
      _muon_stop_rate = 9.1e-5;
      _rpc_scale = (9.6e-5)/(2.1e-3); //change in the pion stopping rate, ignoring arrival time changes
      _rpc_oot_scale = (9.6e-5)/(2.1e-3); //change in the pion stopping rate, should be independent of in time vs out-of-time
      _cosmics_scale = 1.; //Mu2e-II era scaling is from switching 2D density functions for cosmics
      _pbar_scale = 0.;
      _extinction = 1.e-11;
      deadtime_1b = 0.05; //assume 5% deadtime from beam induced clusters in the CRV
      deadtime_2b = 0.05;
      lumi_cutoff_1b = 0.; //ignore the high intensity cut-off loss
      lumi_cutoff_2b = 0.;
      if(mu2eIIMode == 2) { //use one-batch mode parameters, originally assumed POT
        _npot_1b = 4.4e22;
        _npot_2b = 0.;
      } else if(mu2eIIMode == 5) { //3 year running nominal, from Frank
        _npot_2b = 3.6e22;
      } else if(mu2eIIMode == 6) { //running for 5 years instead of 3
        _npot_2b = 3.6e22*5./3.;
        _cosmics_scale *= 5./3.;
      } else { //originally assumed POT (Mode 1/3/4)
        _npot_1b = 0.;
        _npot_2b = 4.4e22;
      }
      if(mu2eIIMode == 3) _rpc_scale *= 0.1; //suppress RPC by an additional factor
      if(mu2eIIMode == 4) _rpc_scale = 0.; //remove in-time RPC
      if(mu2eIIMode == 7) _rpc_scale *= 1./34.; //suppress by the ratio of the rate in a 65 ns shifted window
    // else {
    //   printf("constants::%s: Unknown Mu2e-II mode %i\n",__func__,mu2eIIMode);
    //   _npot_1b = 0.;
    //   _npot_2b = 4.4e22;
    // }
    }

    if(isMu2eIIDataset) {
      trig_eff = 1.; //Already split into triggered/untriggered data streams
    }

    if(mu2eIIMode == 0 && Mode == 5) { // assume full Mu2e running
      _npot_2b = (3.6e20 - _npot_1b);
      _cosmics_scale *= 5.86; //multiply (Run 1 livetime + (2-batch time / N(2-batch POT) * (3.6e20 - Run 1 POT))) / (Run 1 livetime) = (1.06e7 + 5.01e7) / 1.e06e7
    }

    if(Mode == 3 || Mode == 5) { //apply a series of efficiencies
      _npot_1b *= (1-deadtime_1b)*eff_ur_1b*(1-lumi_cutoff_1b/2.)*trig_eff;
      _npot_2b *= (1-deadtime_2b)*eff_ur_2b*(1-lumi_cutoff_2b/2.)*trig_eff;
      _cosmics_scale *= trig_eff;
    }
  }
};
