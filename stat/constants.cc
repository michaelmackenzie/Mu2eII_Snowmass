//
#include "mu2eii/stat/stat/constants.hh"
#include <cstdio>

namespace mu2eii {

  constants::constants() : _muon_capture(0.609), _muon_stop_rate(1.59e-3),
                           _npot_1b(2.86e19), _npot_2b(9.03e18), _cosmics_scale(1.),
                           _pbar_scale(1.), _rpc_scale(1.), _extinction(1.e-10) {
  }
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
      _muon_stop_rate = 8.9e-5;
      _rpc_scale = (9.6e-5)/(2.1e-3); //change in the pion stopping rate, ignoring arrival time changes
      if(mu2eIIMode == 1 || mu2eIIMode == 3) { //use two-batch mode parameters
        _npot_1b = 0.;
        _npot_2b = 4.4e22;
        if(mu2eIIMode == 3)_rpc_scale *= 0.1; //suppress RPC by an additional factor
      } else if(mu2eIIMode == 2) { //use one-batch mode parameters
        _npot_1b = 4.4e22;
        _npot_2b = 0.;
      } else {
        printf("constants::%s: Unknown Mu2e-II mode %i\n",__func__,mu2eIIMode);
        _npot_1b = 0.;
        _npot_2b = 4.4e22;
      }
      _cosmics_scale = 3.;
      _pbar_scale = 0.;
      _extinction = 1.e-11;
    }

    if(isMu2eIIDataset) {
      trig_eff = 1.; //Already split into triggered/untriggered data streams
    }

    if(Mode == 3) { //apply a series of efficiencies
      _npot_1b *= (1-deadtime_1b)*eff_ur_1b*(1-lumi_cutoff_1b/2.)*trig_eff;
      _npot_2b *= (1-deadtime_2b)*eff_ur_2b*(1-lumi_cutoff_2b/2.)*trig_eff;
      _cosmics_scale *= trig_eff;
    }
  }
};
