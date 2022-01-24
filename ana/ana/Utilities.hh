#ifndef MU2EII_UTILITIES
#define MU2EII_UTILITIES

#include <cmath>
#include "TH1.h"

namespace mu2eii {
  class Util {
  public:
    //-----------------------------------------------------------------------------
    // Taken from STNTUPLE:
    // parameterization of the DIO spectrum on Al with LL radiative corrections 
    // from mu2e-6309 (by R.Szafron)
    // 'emu' - energy of the muon bound in Al nucleus, not the muon mass
    //-----------------------------------------------------------------------------
    static double DioWeightAl_LL(double E) {

      const double a5(8.9), a6(1.17169), a7(-1.06599e-2), a8(8.14251e-3);
      const double emu(105.194), mAl(25133.), mmu(105.658), me(0.511);
      const double alpha(1./137.036) ;  // alpha EM

      double de, de5, w;

      de  = emu-E-E*E/(2*mAl);

      if (de < 0.) w = 0;

      de5 = de*de*de*de*de;

      double f = 2.*std::log((mmu/me)*(1-de/mmu))-2+2.*std::log(2);

      w   = 1.e-17*std::pow(de/mmu,(alpha/M_PI)*f)*de5*(a5 + de*(a6+de*(a7+a8*de)));


      return w;
    }

    //Get FWHM of a TH1
    static double th1_fwhm(TH1* h) {
      if(!h) return -1.;
      const double max_val = h->GetMaximum();
      const int bin1 = h->FindFirstBinAbove(max_val/2.);
      const int bin2 = h->FindLastBinAbove(max_val/2.);
      const double width = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
      return width;
    }
  };
}
#endif
