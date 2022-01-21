#ifndef FCSYS_MODELBUILDER
#define FCSYS_MODELBUILDER

//local includes
#include "mu2eii/calc/calc/Var.hh"
#include "mu2eii/calc/calc/Functions.hh"
#include "mu2eii/calc/calc/Poisson.hh"

//ROOT includes
#include "TNamed.h"

//c++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

namespace FCSys {
  //-----------------------------------------------------------------------------
  class ModelBuilder: public TNamed {
  public:
    ModelBuilder(const char* Name, const char* Title);
    ~ModelBuilder();

    int LoadModel(const char* file);
    Poisson_t* GetModel() {return _model;}
    var_t& GetPOI() { return _r;}

  private:
    int ParseLine(std::string line);
    int DefineSources(std::string line);
    int AddClasses(std::string line);
    int AddCategories(std::string line);
    int AddRateSources(std::string line);
    int AddRates(std::string line);
    int AddSystematic(std::string line);
    int AddModifier(std::string line);

    int InitializeVariables();
    int ApplyModifications();

  public:
    std::string TrimLine(std::string line) {
      std::string newline = line;
      //Trim whitespace
      boost::trim_left(newline);
      boost::trim_right(newline);
      return newline;
    }


    void SetVerbose(int verbose) {_verbose = verbose;}

  private:
    Poisson_t* _model;
    int _verbose;
    var_t _obs;
    var_t _r; //signal strength modifier
    const double _NOEFFECT = -1.e50; //ignore for modifiers and systematics

    std::vector<var_t>  _sources;
    std::vector<std::string>  _sourceNames;
    std::vector<int>    _classes;
    std::vector<double> _rates;
    // std::vector<std::string> _categories;
    // std::map<std::string, std::vector<std::string>> _catSources; //map from the category to the source names
    // std::map<std::string, std::vector<std::string>> _catRates; //map from the category to the source rates
    std::vector<var_t> _systematics;
    std::vector<var_t> _systematicBetas;
    std::vector<std::string> _systematicTypes;
    std::vector<std::vector<double>>  _systematicValues; //vector of systematics values
    std::vector<std::vector<var_t>>  _systematicVars; //vector of systematics application variables
    std::vector<std::vector<double>>  _modifiers; //map of modifiers and their values
    ClassDef(ModelBuilder,0);
  };
}
#endif
