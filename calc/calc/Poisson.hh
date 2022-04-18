#ifndef FCSYS_POISSON
#define FCSYS_POISSON

//local includes
#include "mu2eii/calc/calc/Var.hh"

//ROOT includes
#include "TH1.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"

//c++ includes
#include <vector>
#include <map>
#include <iostream>

namespace FCSys {

  //////////////////////////////////////////////////////////////
  //class for a Poission PDF with systematics
  //////////////////////////////////////////////////////////////
  class Poisson_t {
  public:

    Poisson_t(TString name, var_t& obs, std::vector<var_t*> mu,
              std::vector<var_t*> sys = {}, std::vector<var_t*> spectators = {}) : name_(name), obs_(obs), mean_(new var_t("total_mean", 0., 0., obs.max_)),
                                                                                   mu_(mu), sys_(sys), spec_(spectators), spectate_(false),
                                                                                   specBins_(100), verbose_(0),
                                                                                   ngen_(1e5), nmax_((int) obs.max_) {
    }

    void SetVerbose(int verbose) {verbose_ = verbose;}

    void InitSpectators();

    void FillSpectators();

    double GetMean();

    double GetNominalMean();

    double Eval(int n);

    void RandomSys(TRandom3& rnd);

    int RandomSample(TRandom3& rnd);

    TH1D* GeneratePDF(TRandom3& rnd);

    void Print();
    void PrintSpectators(const char* path);

    TString name_;
    var_t& obs_;
    var_t* mean_; //total mean
    std::vector<var_t*> mu_;
    std::vector<var_t*> sys_;
    std::vector<var_t*> spec_;
    std::map<var_t*, TH1*> specMap_;
    bool spectate_;
    int specBins_;
    int verbose_;
    int ngen_;
    int nmax_;
  };
}
#endif
