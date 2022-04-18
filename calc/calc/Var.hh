#ifndef FCSYS_VAR
#define FCSYS_VAR

//ROOT includes
#include "TString.h"
#include "TRandom3.h"

//c++ includes
#include <vector>
#include <iostream>

namespace FCSys {

  //////////////////////////////////////////////////////////////
  //a variable class:
  //knows about its:
  // bounds
  // nominal value
  // additive (applied first) variables
  // multiplicative (second) variables
  // power (third) variables
  // PDF describing its uncertainty
  //////////////////////////////////////////////////////////////
  class var_t {
  public:
    var_t(TString name = "Default") : name_(name), min_(-1.), max_(1.), nom_(0.), val_(0.), constant_(false), verbose_(0), retry_(false) {}
    var_t(TString name, double nom) : name_(name), min_(nom-1.), max_(nom+1.), nom_(nom), val_(nom), constant_(false), verbose_(0), retry_(false) {}
    var_t(TString name, double nom, double min, double max, TString pdf = "Gauss") :
      name_(name), min_(min), max_(max), nom_(nom), val_(nom), pdf_(pdf), constant_(false), verbose_(0), retry_(false) {}
    virtual ~var_t();

    //To ignore calls to set_rnd_val, remains constant
    void set_constant(bool constant = true) {
      constant_ = constant;
    }

    //initialize the variable's dependence on other variables
    void set_dependents(std::vector<var_t*> add, std::vector<var_t*> mul, std::vector<var_t*> pow = {}) {
      add_ = add;
      mul_ = mul;
      pow_ = pow;
    }

    //set the base value of the parameter, before the evaluation of other variables
    void set_val(double val) {
      val_ = std::min(max_, std::max(min_, val));
    }

    //set the base value to a random value following its PDF
    void set_rnd_val(TRandom3& rnd) {
      if(constant_) return;
      double val = nom_;
      do {
        if(pdf_ == "Gauss") {
          val = rnd.Gaus(nom_);
        }
        else if(pdf_ == "Flat") {
          val = min_ + (max_-min_)*rnd.Uniform();
        }
      } while(retry_ && (val < min_ || val > max_));
      set_val(val);
    }

    //Evaluate the value by taking the base value and applying the dependent variable values
    virtual double get_val() {
      double val = val_;
      if(verbose_ > 2) printf("Variable %s has starting value %.3e\n", name_.Data(), val);
      for(var_t* var : add_) {
        val += var->get_val();
        if(verbose_ > 2) printf("Variable %s add %s = %.3e --> %.3e\n", name_.Data(), var->name_.Data(), var->get_val(), val);
      }
      for(var_t* var : mul_) {
        val *= var->get_val();
        if(verbose_ > 2) printf("Variable %s mul %s = %.3e --> %.3e\n", name_.Data(), var->name_.Data(), var->get_val(), val);
      }
      for(var_t* var : pow_) {
        val = std::pow(val, var->get_val());
        if(verbose_ > 2) printf("Variable %s pow %s = %.3e --> %.3e\n", name_.Data(), var->name_.Data(), var->get_val(), val);
      }
      if(val > max_ || val < min_) {
        if(verbose_ > 0) std::cout << "Variable " << name_.Data() << " value outside of bounds! Returning bound...\n";
        return std::min(max_, std::max(min_, val));
      }
      return val;
    }

    //print information about the variable
    void print() {
      printf(" %s: %.3e (%.3e) [%.3e - %.3e]", name_.Data(), get_val(), nom_, min_, max_);
      if(add_.size() > 0) {
        printf(" add = {");
        for(var_t* var : add_) {
          printf("%s", var->name_.Data());
          if(var != add_[add_.size()-1]) printf(", ");
        }
        printf("}");
      }
      if(mul_.size() > 0) {
        printf(" mul = {");
        for(var_t* var : mul_) {
          printf("%s", var->name_.Data());
          if(var != mul_[mul_.size()-1]) printf(", ");
        }
        printf("}");
      }
      if(pow_.size() > 0) {
        printf(" pow = {");
        for(var_t* var : pow_) {
          printf("%s", var->name_.Data());
          if(var != pow_[pow_.size()-1]) printf(", ");
        }
        printf("}");
      }
      printf("\n");
    }

    //Fields for the variable, all public for convenience
    TString name_;
    double min_;
    double max_;
    double nom_;
    double val_;
    TString pdf_;
    bool constant_;
    int verbose_;
    bool retry_; //re-evaluate random drawing if out of bounds
    std::vector<var_t*> add_;
    std::vector<var_t*> mul_;
    std::vector<var_t*> pow_;
  };
}
#endif
