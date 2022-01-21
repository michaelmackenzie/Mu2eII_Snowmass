//
#include "mu2eii/stat/stat/analysis.hh"

namespace mu2eii {

analysis::analysis(const char* Name, const char* Title) : TNamed(Name,Title) {
  fSignal = nullptr;
  fListOfBgrChannels = new TObjArray();
}

//-----------------------------------------------------------------------------
analysis::~analysis() {
  if (fSignal) delete fSignal;
  
  fListOfBgrChannels->Delete();
  delete fListOfBgrChannels;
}

  
double analysis::GetSigIntegral(float PMin, float PMax, float TMin, float TMax) {

  float sw  = fSignal->GetIntegral(PMin,PMax,TMin,TMax); 
  sw        = sw*fRmue;

  return sw;
}


double analysis::GetBgrIntegral(float PMin, float PMax, float TMin, float TMax) {

  float sw  = 0;

  int nbgr = fListOfBgrChannels->GetEntriesFast();
  
  for (int ibgr=0; ibgr<nbgr; ibgr++) {
    channel* bgr = GetBgrChannel(ibgr);
    sw += bgr->GetIntegral(PMin,PMax,TMin,TMax);
  }
					//no additional scaling needed
  return sw;
}
};
