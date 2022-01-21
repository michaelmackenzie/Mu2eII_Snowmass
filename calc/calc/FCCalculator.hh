#ifndef FCSYS_FCCALCULATOR
#define FCSYS_FCCALCULATOR

//local includes
#include "mu2eii/calc/calc/Var.hh"
#include "mu2eii/calc/calc/Poisson.hh"

//ROOT includes
#include "TH1.h"
#include "TRandom3.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"

//c++ includes
#include <set>

namespace FCSys {

  //////////////////////////////////////////////////////////////
  // class to perform Feldman-Cousins limits
  //////////////////////////////////////////////////////////////

  class FCCalculator {
  public:
    FCCalculator(Poisson_t& model, var_t& poi, TRandom3& rnd, double cl = 0.9, int verbose = 0) :
      model_(model), poi_(poi), rnd_(rnd), cl_(cl), verbose_(verbose), res_(1.e-3) {
      double oldval = poi.val_;
      poi.val_ = 0.;
      hNull_ = model.GeneratePDF(rnd);
      hNull_->SetName("FC_NULL");
      null_mu_ = hNull_->GetMean();
      if(verbose > 0) printf("%s: Null mean is %.3e events\n", __func__, null_mu_);
      poi.val_ = oldval;
      CalculateIndividualInterval(hNull_, null_min_, null_max_);
      if(null_min_ > 0) {
        printf("!!! %s: N(obs) = 0 is not contained within the NULL hypothesis!\n", __func__);
      }
      if(verbose_ > 0) printf("%s: Null interval at %.3f CL: %i - %i\n", __func__, cl, null_min_, null_max_);
    }

    //Number of events needed to be seen to be >x sigma on the right tail
    int NSigmaThreshold(TH1D* hPDF, double nsigma);

    //Get the median of the distribution
    int GetMedian(TH1D* hPDF);

    //Find the minimum value of the POI that has a median of n
    double FindForMedianN(int n);

    //For a given PDF, construct the FC interval in N(observed)
    void CalculateIndividualInterval(TH1D* hPDF, int& nmin, int& nmax);

    //Get the upper or lower limit for the model given an observation
    double FindLimit(int nobs, bool upperLimit);

    //Calculate the interval for a given observation
    void CalculateInterval(int nobs, double& mu_min, double& mu_max);

    Poisson_t& model_;
    var_t& poi_;
    TRandom3& rnd_;
    double cl_;
    int verbose_;
    TH1D* hNull_;
    double null_mu_;
    int null_min_;
    int null_max_;
    double res_;
  };
}
#endif
