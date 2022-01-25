#include "mu2eii/calc/calc/FCCalculator.hh"

namespace FCSys {

  //Number of events needed to be seen to be >x sigma on the right tail
  int FCCalculator::NSigmaThreshold(TH1D* hPDF, double nsigma) {
    if(nsigma < 0.) return 0; //ignore this region
    //p-value for the N sigma deviation
    const double psigma = ROOT::Math::gaussian_cdf(-nsigma);
    double p = 1.;
    const int nbins = hPDF->GetNbinsX();
    int bin = 0;
    if(verbose_ > 1) printf("%s: Printing threshold calculation:\n", __func__);
    //loop through bins until the right-side integral is <= the desired p-value
    do { //ensure it enters the loop at least once
      ++bin;
      p = hPDF->Integral(bin, nbins);
      if(verbose_ > 1) {
        const double binc = hPDF->GetBinCenter(bin);
        const int nevents = binc + 0.1; //ensure it rounds correctly
        printf(" n = %2i P(n' >= n) = %.3e\n", nevents, p);
      }
    } while(p > psigma && bin <= nbins);

    //get the N(events) value for the bin
    const double binc = hPDF->GetBinCenter(bin);
    //Two reasonable scenarios:
    //Bin low edge = n --> bin center = n + 0.5 --> rounds down to n correctly
    //Bin center = n --> ben center either rounds to n or n-1, depending on accuracy
    //Adding 0.1 won't change this in any of the cases
    const int nevents = binc + 0.1;
    return nevents;
  }

  //Get the median of the distribution
  int FCCalculator::GetMedian(TH1D* hPDF) {
    double p = 0.;
    int bin = 0;
    while(p < 0.5) {
      ++bin;
      p += hPDF->GetBinContent(bin);
    }
    const double binc = hPDF->GetBinCenter(bin);
    const int nevents = binc + 0.1;
    return nevents;
  }

  //Find the minimum value of the POI that has a median of n
  double FCCalculator::FindForMedianN(int n) {
    double mu_max = poi_.max_;
    double mu_min = poi_.min_;
    while(fabs(mu_max - mu_min)/(mu_max+mu_min) > res_) {
      const double mu = (mu_max + mu_min) / 2.;
      poi_.val_ = mu;
      TH1D* h = model_.GeneratePDF(rnd_);
      int median = GetMedian(h);
      if(median >= n) mu_max = mu; //mu satisfies criteria, so set as maximum
      else            mu_min = mu; //mu fails, so must be larger
      delete h;
    }
    return (mu_min + mu_max) / 2.;
  }

  //For a given PDF, construct the FC interval in N(observed)
  void FCCalculator::CalculateIndividualInterval(TH1D* hPDF, int& nmin, int& nmax) {
    nmin = hPDF->GetBinCenter(hPDF->GetNbinsX()) + 0.1;
    nmax = hPDF->GetBinCenter(1) + 0.1;
    double p = 0.;
    //mu_bkg for the denominator of the likelihood ratio ordering parameter
    const double mu = null_mu_;//hNull_->GetMean() - 0.5; //bins are centered at n + 0.5
    if(hPDF->Integral() < cl_) {
      printf("!!! %s: PDF doesn't have the range to calculate the confidence level precision!\n", __func__);
      nmax = nmin;
      nmin = 0;
      return;
    }
    std::set<int> ns;
    if(verbose_ > 2) printf("%s: Printing interval construction:\n", __func__);
    while(p < cl_) {
      int nbest = -1; //N(events) for highest R
      int bbest =  0; //Bin for highest R
      double rbest = -1.; //Highest R
      for(int bin = 1; bin <= hPDF->GetNbinsX(); ++bin) { //check each histogram bin
        const double binc = hPDF->GetBinCenter(bin);
        const int nevents = binc + 0.1;
        if(ns.count(nevents)) continue; //skip if already added to the interval
        //FIXME: could generate an array of PDFs for mu + mu_s = n for n = (int) mu - nmax to use for the denominator
        const double r = hPDF->GetBinContent(bin) / ROOT::Math::poisson_pdf(nevents, std::max(mu, (double) nevents));
        if(r > rbest) {
          rbest = r;
          nbest = nevents;
          bbest = bin;
        }
      }
      if(nbest < 0) {
        printf("!!! Error: FCCalculator::%s: No next interval addition for the PDF (likely due to integral of PDF < 1)! p(current) = %.3e, p(target) = %.3e, nmin = %i, nmax = %i",
               __func__, p, cl_, nmin, nmax);
      }
      //Add the highest R N(events) to the interval
      p += hPDF->GetBinContent(bbest);
      nmin = std::min(nbest, nmin);
      nmax = std::max(nbest, nmax);
      ns.insert(nbest);
      if(verbose_ > 9) printf(" n = %i, r = %.3e, p = %.3f\n", nbest, rbest, p);
    }
    if(verbose_ > 2) printf(" Final interval: %i - %i\n", nmin, nmax);
  }

  //Get the upper or lower limit for the model given an observation
  double FCCalculator::FindLimit(int nobs, bool upperLimit) {
    int attempts = 0;
    const int maxAttempts = 100;
    double mu_min = 0.;
    double mu_max = poi_.max_;
    double mu_range = poi_.max_ - poi_.min_;
    while(fabs((mu_max - mu_min) / mu_range) > res_ && attempts < maxAttempts) {
      ++attempts;
      double mu = (mu_max + mu_min) / 2.;
      poi_.val_ = mu;
      if(verbose_ > 2) model_.Print();
      TH1D* h = model_.GeneratePDF(rnd_);
      int nmin, nmax;
      CalculateIndividualInterval(h, nmin, nmax);
      if(upperLimit) {
        if(nmin > nobs) mu_max = mu; //gone past allowed
        else mu_min = mu; //still allowed
      } else {
        if(nmax < nobs) mu_min = mu; //gone past allowed
        else mu_max = mu; //still allowed
      }
      delete h;
      if(verbose_ > 1) printf("%s: Attempt %i has bounds: %10.4f - %10.4f (%.3e vs res_ requirement of %.3e)\n",
                              __func__, attempts, mu_min, mu_max, fabs((mu_max - mu_min) / mu_range), res_);
    }
    if(attempts == maxAttempts) printf("!!! %s: Hit maximum limit finding attempts for N(obs) = %i, mu in %.3f - %.3f (upper = %o)\n",
                                       __func__, nobs, mu_min, mu_max, upperLimit);
    return (mu_max + mu_min) / 2.;
  }

  //Calculate the interval for a given observation
  void FCCalculator::CalculateInterval(int nobs, double& mu_min, double& mu_max) {
    //first test if mu = 0 is included in the interval
    if(nobs > null_max_) {
      mu_min = FindLimit(nobs, false);
    } else {
      mu_min = 0.;
    }
    mu_max = FindLimit(nobs, true);
  }

} //end FCSys namespace
