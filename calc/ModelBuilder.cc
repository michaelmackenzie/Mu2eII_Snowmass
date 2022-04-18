#include "mu2eii/calc/calc/ModelBuilder.hh"

namespace FCSys {

  ModelBuilder::ModelBuilder(const char* Name, const char* Title) : TNamed(Name, Title), _model(nullptr), _verbose(0) {}
  ModelBuilder::~ModelBuilder() {}

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::LoadModel(const char* file) {
    int status(0);
    ///////////////////////////////
    // Open the data card
    ///////////////////////////////
    std::fstream infile;
    infile.open(file, std::ios::in);
    if(!infile.is_open()) {
      printf("ModelBuilder::%s Error! File %s is not found", __func__, file);
      return 10;
    }
    std::string line;
    std::vector<std::string> comments = {"#", "//", "!"};
    if(_verbose > 2) printf("ModelBuilder::%s: Printing datacard %s\n", __func__, file);

    ///////////////////////////////
    // Read the data card
    ///////////////////////////////
    while(getline(infile, line)) {
      if(_verbose > 3) std::cout << line << std::endl;
      //Remove comments
      for(auto comment : comments) {
        line = line.substr(0, line.find(comment));
      }
      if(line.size() == 0) continue;

      line = TrimLine(line);
      if(line.size() == 0) continue;

      if(_verbose == 3) std::cout << line << std::endl;

      status = ParseLine(line);
      if(status) {
        printf("ModelBuilder::%s: Error reading line %s, status = %i\n", __func__, line.c_str(), status);
        return status;
      };
    }
    //Apply the effects given in the card, such as systematic uncertainties
    status = ApplyModifications();
    if(status) return status;

    ///////////////////////////////
    // Create the model
    ///////////////////////////////
    var_t& signal_mu = GetPOI();
    if(_verbose > 2) {
      printf("Signal strength modifier information:\n");
      signal_mu.print();
    }
    //initialize an observable
    _obs.name_ = "Number observed";
    _obs.nom_ = 0.;
    _obs.val_ = 0.;
    _obs.min_ = 0.;
    _obs.max_ = 1e6;

    if(_verbose > 4) printf("ModelBuilder::%s: Printing variable information before model creation:\n", __func__);
    if(_verbose > 6) {
      for(std::vector<var_t>& vars : _systematicVars) {
        std::cout << " Printing next var list\n";
        for(var_t& var : vars) {
          std::cout << "var  " << var.name_.Data() << " address = " << &(var) << std::endl;
          var.print();
        }
      }
    }


    //initialize a list of source contributions
    std::vector<var_t*> vars;
    std::vector<var_t*> spectators;
    for(var_t& var : _sources) {
      vars.push_back(&var);
      spectators.push_back(&var);
      if(_verbose > 6) {
        var.verbose_ = 3;
        printf(" Printing information for variable %s:\n", var.name_.Data());
      }
      if(_verbose > 4) var.print();
      if(_verbose > 6) var.verbose_ = 0;
      for(var_t* dep : var.mul_) {
        if(!std::count(spectators.begin(), spectators.end(), dep)) spectators.push_back(dep);
      }
    }

    //initialize a list of systematic betas
    std::vector<var_t*> betas;
    for(var_t& var : _systematicBetas) {
      betas.push_back(&var);
      spectators.push_back(&var);
      if(_verbose > 6) {
        var.verbose_ = 3;
        printf(" Printing information for variable %s:\n", var.name_.Data());
      }
      if(_verbose > 4) var.print();
      if(_verbose > 6) var.verbose_ = 0;
    }

    //initialize the model
    _model = new Poisson_t("Counting model", _obs, vars, betas, spectators);
    if(_verbose > 2) _model->Print();


    return status;
  }

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::ParseLine(std::string line) {
    int status(0);
    std::string label = line.substr(0, line.find(" "));
    if(_verbose > 9) printf("ModelBuilder::%s: Line label %s\n", __func__, label.c_str());
    if(line == "") {
      printf("ModelBuilder::%s: Error! Empty line\n", __func__);
      return 1;
    }
    if(label == "source") {
      if(_sources.size() == 0)      return DefineSources(line);
      else                          return AddRateSources(line);
    } else if(label == "class"   ){ return AddClasses(line);
    } else if(label == "category"){ return AddCategories(line);
    } else if(label == "rate"    ){ return AddRates(line);
    } else if(label == "sys"     ){ return AddSystematic(line);
    } else if(label == "mod"     ){ return AddModifier(line);
    } else {
      printf("ModelBuilder::%s: Error! Unrecognized line label %s for line = %s\n", __func__, label.c_str(), line.c_str());
      return 2;
    }
    return status;
  }

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::DefineSources(std::string line){
    int status(0);
    std::string label = line.substr(0, line.find(" "));
    if(label != "source") {
      printf("ModelBuilder::%s: Error! Defining sources but line label = %s for line %s\n", __func__, label.c_str(), line.c_str());
      return 1;
    }
    if(_sources.size() != 0) {
      printf("ModelBuilder::%s: Error! Defining sources but sources are already defined for line %s\n", __func__, line.c_str());
      return 2;
    }
    std::string remainder = line;
    remainder.erase(0, remainder.find(" ")); //remove the label
    remainder = TrimLine(remainder);
    std::size_t pos(0);
    while(remainder.size() > 0) {
      pos = ((pos = remainder.find(" ")) != std::string::npos) ? pos : remainder.size();
      std::string next = remainder.substr(0, pos);
      if(std::count(_sourceNames.begin(), _sourceNames.end(), next)) {
        printf("ModelBuilder::%s: Error! Duplicate source %s in line %s\n", __func__, next.c_str(), line.c_str());
      }
      _sources.push_back(var_t(next.c_str()));
      _sourceNames.push_back(next);
      remainder.erase(0, pos); //remove the current token
      remainder = TrimLine(remainder);
      if(_verbose > 4) printf("ModelBuilder::%s: Adding source %s\n", __func__, next.c_str());
      if(_verbose > 5) printf("ModelBuilder::%s: Remaining line to process %s\n", __func__, next.c_str());
    }
    if(_sources.size() < 2) {
      printf("ModelBuilder::%s: Error! There must be at least 2 sources, %li found for line %s\n", __func__, _sources.size(), line.c_str());
      return 2;
    }
    return status;
  }

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::AddClasses(std::string line){
    int status(0);
    std::string label = line.substr(0, line.find(" "));
    if(label != "class") {
      printf("ModelBuilder::%s: Error! Defining classes but line label = %s for line %s\n", __func__, label.c_str(), line.c_str());
      return 1;
    }
    if(_classes.size() != 0) {
      printf("ModelBuilder::%s: Error! Defining classes but classes are already defined for line %s\n", __func__, line.c_str());
      return 2;
    }
    if(_sources.size() == 0) {
      printf("ModelBuilder::%s: Error! Defining classes but sources are not yet defined for line %s\n", __func__, line.c_str());
      return 3;
    }
    std::string remainder = line;
    remainder.erase(0, remainder.find(" ")); //remove the label
    remainder = TrimLine(remainder);
    std::size_t pos(0);
    int signal(0), background(0);
    while(remainder.size() > 0) {
      pos = ((pos = remainder.find(" ")) != std::string::npos) ? pos : remainder.size();
      std::string next = remainder.substr(0, pos);
      try {
        int val = std::stoi(next);
        if(val == 0) ++signal;
        else ++background;
        _classes.push_back(val);
      } catch(const std::invalid_argument& e) {
        printf("ModelBuilder::%s: Value conversion error on token %s in line %s\n", __func__, next.c_str(), line.c_str());
        throw e;
      }
      remainder.erase(0, pos); //remove the current token
      remainder = TrimLine(remainder);
    }
    if(_classes.size() != _sources.size()) {
      printf("ModelBuilder::%s: Error! Number of class definitions (%li) doesn't match the sources (%li) for line %s\n",
             __func__, _classes.size(), _sources.size(), line.c_str());
      return 4;
    }
    if(!background || !signal) {
      printf("ModelBuilder::%s: Error! A background and a signal was not found for line %s\n",
             __func__, line.c_str());
      return 5;
    }
    if(signal != 1) {
      printf("ModelBuilder::%s: Error! Currently only 1 signal is supported, but %i found for line %s\n",
             __func__, signal, line.c_str());
      return 6;
    }
    return status;
  }

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::AddCategories(std::string line){return 0;}

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::AddRateSources(std::string line){return 0;}

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::AddRates(std::string line){
    int status(0);
    std::string label = line.substr(0, line.find(" "));
    if(label != "rate") {
      printf("ModelBuilder::%s: Error! Defining rates but line label = %s for line %s\n", __func__, label.c_str(), line.c_str());
      return 1;
    }
    if(_rates.size() != 0) {
      printf("ModelBuilder::%s: Error! Defining rates but rates are already defined for line %s\n", __func__, line.c_str());
      return 2;
    }
    if(_classes.size() == 0 || _sources.size() == 0) {
      printf("ModelBuilder::%s: Error! Sources (N = %li) and classes (N = %li) must be defined before rates for line %s\n",
             __func__, _sources.size(), _classes.size(), line.c_str());
      return 3;
    }
    std::string remainder = line;
    remainder.erase(0, remainder.find(" ")); //remove the label
    remainder = TrimLine(remainder);
    std::size_t pos(0);
    while(remainder.size() > 0) {
      pos = ((pos = remainder.find(" ")) != std::string::npos) ? pos : remainder.size();
      std::string next = remainder.substr(0, pos);
      try {
        double val = std::stod(next);
        _rates.push_back(val);
      } catch(const std::invalid_argument& e) {
        printf("ModelBuilder::%s: Value conversion error on token %s in line %s\n", __func__, next.c_str(), line.c_str());
        throw e;
      }
      remainder.erase(0, pos); //remove the current token
      remainder = TrimLine(remainder);
    }
    if(_rates.size() != _sources.size()) {
      printf("ModelBuilder::%s: Error! Number of rate definitions (%li) doesn't match the sources (%li) for line %s\n",
             __func__, _rates.size(), _sources.size(), line.c_str());
      return 4;
    }

    //Now that the source, class, and rates are defined, can initialize the variables
    InitializeVariables();

    return status;
  }

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::AddSystematic(std::string line){
    //Remove the "sys" label from the line defining a systematic
    std::string label = line.substr(0, line.find(" "));
    if(label != "sys") {
      printf("ModelBuilder::%s: Error! Defining systematic but line label = %s for line %s\n", __func__, label.c_str(), line.c_str());
      return 1;
    }
    std::string remainder = line;
    remainder.erase(0, remainder.find(" ")); //remove the label
    remainder = TrimLine(remainder);

    if(remainder.size() == 0) {
      printf("ModelBuilder::%s: Error! Defining systematic but no definition information for line %s\n", __func__, line.c_str());
      return 2;
    }

    //Get the name of the systematic
    std::string name = remainder.substr(0, remainder.find(" "));
    remainder.erase(0, remainder.find(" ")); //remove the name
    remainder = TrimLine(remainder);
    _systematics.push_back(var_t(name));
    // var_t& sys = _systematics.back();

    //Get the systematic type, e.g. log-normal (lnN)
    std::string type = remainder.substr(0, remainder.find(" "));
    remainder.erase(0, remainder.find(" ")); //remove the type
    remainder = TrimLine(remainder);
    _systematicTypes.push_back(type);

    //Retrieve the size of the systematic for each process
    std::vector<double> values;
    std::size_t pos(0);
    while(remainder.size() > 0) {
      pos = ((pos = remainder.find(" ")) != std::string::npos) ? pos : remainder.size();
      std::string next = remainder.substr(0, pos);
      if(next == "-") { //Does not apply to this process
        values.push_back(_NOEFFECT);
      } else {
        try {
          const double val = std::stod(next);
          values.push_back(val);
        } catch(const std::invalid_argument& e) {
          printf("ModelBuilder::%s: Value conversion error on token %s in line %s\n", __func__, next.c_str(), line.c_str());
          throw e;
        }
      }
      remainder.erase(0, pos); //remove the current token
      remainder = TrimLine(remainder);
    }

    //Check that the values size agrees with the number of processes
    if(values.size() != _sources.size()) {
      printf("ModelBuilder::%s: Error! Number of systematic values (%li) doesn't match the sources (%li) for line %s\n",
             __func__, values.size(), _sources.size(), line.c_str());
      return 4;
    }

    //store the list of values associated with this systematic
    _systematicValues.push_back(values);
    return 0;
  }

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::AddModifier(std::string line){return 0;}

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::InitializeVariables(){
    if(_verbose > 3) printf("ModelBuilder::%s: Priniting initial variable information:\n", __func__);
    _r.nom_ = 1.;
    _r.val_ = 1.;
    _r.name_ = "Signal POI (r)";
    _r.min_  = 0.;
    _r.max_  = 20.;
    for(unsigned index = 0; index < _sources.size(); ++index) {
      var_t& var = _sources[index];
      var.nom_ = _rates[index];
      var.val_ = _rates[index];
      var.min_ = 0.;
      if(_classes[index]) {
        var.max_ = std::max(100.*var.nom_, 1.);
      } else { //signal
        var.max_ = std::max(10.*var.nom_, 50.);
        var.set_dependents({}, {&_r}, {});
      }
      if(_verbose > 3) var.print();
    }
    return 0;
  }

  //-----------------------------------------------------------------------------------------
  int ModelBuilder::ApplyModifications(){
    ///////////////////////////////////////////////////
    //First create the systematics
    for(unsigned index = 0; index < _systematics.size(); ++index) {
      var_t& sys = _systematics[index];
      _systematicBetas.push_back(var_t(Form("%s_beta",sys.name_.Data())));
      _systematicVars.push_back(std::vector<var_t>()); //create an empty list of variables
      std::vector<var_t>& vars = _systematicVars.back();
      //add a systematic variable for each process
      for(unsigned isource = 0; isource < _sources.size(); ++isource) {
        vars.push_back(var_t(Form("%s_%s_kappa", _sourceNames[isource].c_str(), sys.name_.Data())));
      }
    }

    ///////////////////////////////////////////////////
    //Now configure the systematics
    for(unsigned index = 0; index < _systematics.size(); ++index) {
      var_t& sys = _systematics[index];
      std::string type = _systematicTypes[index];
      std::vector<double> values = _systematicValues[index];

      //Create a normally distributed beta to control the uncertainty fluctuations
      var_t& beta = _systematicBetas[index];
      beta.nom_ = 0.; beta.val_ = 0.; beta.min_ = -10.; beta.max_ = 10.; beta.pdf_ = "Gauss"; beta.constant_ = false;
      if(_verbose > 5) {
        printf("ModelBuilder::%s: Retrieving systematic beta for %s:\n", __func__, sys.name_.Data());
        beta.print();
      }
      std::vector<var_t>& vars = _systematicVars[index];

      //apply the effect to each source with the source specific size
      for(unsigned isource = 0; isource < _sources.size(); ++isource) {
        var_t& source       = _sources[isource]; //process this systematic effects
        const double value  = values[isource]; //nominal rate
        var_t& var          = vars[isource]; //the systematic variable for this process
        if(_verbose > 6) std::cout << "Variable " << var.name_.Data() << " address = " << &var << std::endl;
        if(value == _NOEFFECT) continue; //does not effect this source
        if(type == "lnN") { //Log-normal
          var.nom_ = 1. + value; var.val_ = 1. + value; var.min_ = 0.; var.max_ = 1+10.*value;
          var.set_dependents({}, {}, {&beta});
          source.mul_.push_back(&var);
        } else if(type == "Gaus") { //Gaussian
          //source = nominal + sys; sys = Gaussian(mean = 0, width = uncertainty) = uncertainty*Gaussian(0, 1)
          var.nom_ = 0.; var.val_ = 0.; var.min_ = -10.*value; var.max_ = 10.*value;
          var.set_dependents({}, {&beta}, {});
          source.add_.push_back(&var);
        } else {
          printf("ModelBuilder::%s: Error! Unknown uncertainty type %s for systematic %s", __func__, type.c_str(), sys.name_.Data());
          return 1;
        }
        if(_verbose > 5) {
          printf("ModelBuilder::%s: Creating systematic kappa for %s and source %s:\n", __func__, sys.name_.Data(), source.name_.Data());
          var.print();
          source.print();
        }
      }
    }
    if(_verbose > 6) {
      printf("ModelBuilder::%s: Printing the systematic kappas again\n", __func__);
      for(std::vector<var_t>& vars : _systematicVars) {
        std::cout << " Printing next var list\n";
        for(var_t& var : vars) {
          std::cout << "var  " << var.name_.Data() << " address = " << &(var) << std::endl;
          var.print();
        }
      }
    }

    return 0;
  }

}
