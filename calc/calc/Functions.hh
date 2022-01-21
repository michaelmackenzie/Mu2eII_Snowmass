#ifndef FCSYS_FUNCTIONS
#define FCSYS_FUNCTIONS

//local includes
#include "mu2eii/calc/calc/Var.hh"

//ROOT includes
#include "TString.h"

//c++ includes
#include <vector>
#include <cmath>

namespace FCSys {

  //Exponential function wrapper
  class exp_t : public var_t {
  public:
    exp_t(TString name = "Default") : var_t(name) {}
    exp_t(TString name, double nom) : var_t(name, nom) {}
    ~exp_t() {};

    double get_val() override {
      double val = var_t::get_val();
      return std::exp(val);
    }
    void set_dependents(std::vector<var_t*> add, std::vector<var_t*> mul = {}, std::vector<var_t*> pow = {}) {
      var_t::set_dependents(add, mul, pow);
    }
    void print() {
      printf(" %s: %.3e = e^{((%.3e", name_.Data(), get_val(), val_);
      for(unsigned index = 0; index < add_.size(); ++index) {
        printf(" + %s", add_[index]->name_.Data());
      }
      printf(")");
      for(unsigned index = 0; index < mul_.size(); ++index) {
        printf(" * %s", mul_[index]->name_.Data());
      }
      printf(")");
      if(pow_.size() > 0) {
        printf("^{");
        for(unsigned index = 0; index < pow_.size(); ++index) {
          printf("%s", pow_[index]->name_.Data());
          if(index < pow_.size() - 1) printf(" * ");
        }
        printf("}");
      }
      printf("}\n");
    }
  };

  //Polynomial function
  class pow_t : public var_t {
  public:
    pow_t(TString name = "Default") : var_t(name) {}
    pow_t(TString name, double nom) : var_t(name, nom) {}
    ~pow_t() {};

    double get_val() override {
      double val(0.);
      for(unsigned index = 0; index < powers_.size(); ++index) {
        val += coeff_[index]->get_val()*std::pow(vars_[index]->get_val(), powers_[index]);
      }
      return val;
    }
    void set_dependents(std::vector<var_t*> vars, std::vector<var_t*> coeff, std::vector<int> powers) {
      if(vars.size() != coeff.size() || vars.size() != powers.size()) {
        std::cout << "!!! Error: pow::" << __func__ << " polynomial parameter lists not the same length!\n";
        throw 10;
      }
      vars_ = vars;
      coeff_ = coeff;
      powers_ = powers;
    }
    void print() {
      printf(" %s: %.3e = ", name_.Data(), get_val());
      for(unsigned index = 0; index < vars_.size(); ++index) {
        printf(" (%s) * (%s)^{%i}", coeff_[index]->name_.Data(), vars_[index]->name_.Data(), powers_[index]);
        if(index < vars_.size() - 1) printf(" + ");
      }
      printf("\n");
    }
    std::vector<var_t*> vars_;
    std::vector<var_t*> coeff_;
    std::vector<int> powers_;
  };
}
#endif
