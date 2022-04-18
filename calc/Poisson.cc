#include "mu2eii/calc/calc/Poisson.hh"

namespace FCSys {

  void Poisson_t::InitSpectators() {
    for(auto spectator : spec_) {
      if(specMap_.find(spectator) != specMap_.end() && specMap_[spectator]) {
        delete specMap_[spectator];
        specMap_[spectator] = 0;
      }
      const double xmin = spectator->min_;
      const double xmax = (spectator->nom_ > 1.e-5) ? std::min(20.*spectator->nom_, spectator->max_) : spectator->max_; //betas have nominal of 0, rest use 20*nominal
      specMap_[spectator] = new TH1D(Form("h%s", spectator->name_.Data()), Form("%s samples", spectator->name_.Data()),
                                     specBins_, xmin, xmax);
    }
    //FIXME: Set max to model related parameter
    specMap_[mean_] = new TH1D("htotal_mean", "Total mean samples", specBins_, mean_->min_, 20.); //mean_->max_);
  }

  void Poisson_t::FillSpectators() {
    for(auto spectator : spec_) {
      specMap_[spectator]->Fill(spectator->get_val());
    }
    specMap_[mean_]->Fill(GetMean());
  }

  double Poisson_t::GetMean() {
    double mu = 0.;
    for(var_t* var : mu_) {
      mu += var->get_val();
      if(verbose_ > 9) std::cout << __func__ << ": mu(" << var->name_.Data() << ") = " << var->get_val() << std::endl;
    }
    if(verbose_ > 9) std::cout << __func__ << ": mu = " << mu << std::endl;
    return mu;
  }

  double Poisson_t::GetNominalMean() {
    double mu = 0.;
    for(var_t* var : sys_) {
      var->val_ = var->nom_;
    }
    for(var_t* var : mu_) {
      var->val_ = var->nom_;
      mu += var->get_val();
    }
    if(verbose_ > 9) std::cout << __func__ << ": mu = " << mu << std::endl;
    return mu;
  }

  double Poisson_t::Eval(int n) {
    const double mu = GetMean();
    const double p = ROOT::Math::poisson_pdf(n, mu);
    if(verbose_ > 9) std::cout << __func__ << ": p = " << p << std::endl;
    return p;
  }

  void Poisson_t::RandomSys(TRandom3& rnd) {
    for(var_t* var : sys_) {
      var->set_rnd_val(rnd);
      if(verbose_ > 9) {std::cout << __func__ << ":\n"; var->print();}
    }
  }

  int Poisson_t::RandomSample(TRandom3& rnd) {
    //first sample the systematics
    RandomSys(rnd);
    //next sample the poisson distribution
    const double mu = GetMean();
    const int n = rnd.Poisson(mu);
    if(verbose_ > 9) std::cout << __func__ << ": n = " << n << std::endl;
    return n;
  }

  TH1D* Poisson_t::GeneratePDF(TRandom3& rnd) {
    //sample the nuisance parameters to define a mean, then add a poisson PDF for this
    int nbins = nmax_;
    TH1D* hpdf = new TH1D("hpdf", "PDF", nbins, -0.5, nbins - 0.5);
    const int nattempts = ngen_;
    if(spectate_) InitSpectators();
    for(int attempt = 0; attempt < nattempts; ++attempt) {
      RandomSys(rnd);
      const double mu = GetMean();
      for(int n = 0; n < nbins; ++n) {
        hpdf->Fill(n, ROOT::Math::poisson_pdf(n, mu));
        if(nattempts == 1) hpdf->SetBinError(n+1, 0.);
      }
      if(spectate_) FillSpectators();
    }
    hpdf->Scale(1. / nattempts);
    return hpdf;
  }

  void Poisson_t::Print() {
    printf(" %s: obs = %s, mu = {", name_.Data(), obs_.name_.Data());
    for(var_t* var : mu_) {
      printf("%s", var->name_.Data());
      if(var != mu_.back()) printf(", ");
    }
    printf("} sys = {");
    for(var_t* var : sys_) {
      printf("%s", var->name_.Data());
      if(var != sys_.back()) printf(", ");
    }
    printf("} spectators = {");
    for(var_t* var : spec_) {
      printf("%s", var->name_.Data());
      if(var != spec_.back()) printf(", ");
    }
    printf("}\n");
  }

  void Poisson_t::PrintSpectators(const char* path) {
    gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", path, path));
    for(auto spectator : spec_) {
      TCanvas* c = new TCanvas();
      specMap_[spectator]->Draw();
      c->Print(Form("%s/c_%s.png", path, spectator->name_.Data()));
      delete c;
    }
    TCanvas* c = new TCanvas();
    specMap_[mean_]->Draw();
    c->Print(Form("%s/c_total_mean.png", path));
    delete c;
  }
}
