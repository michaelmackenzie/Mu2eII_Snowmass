//
#ifndef mu2eii_sens_analysis_hh
#define mu2eii_sens_analysis_hh

#include "TNamed.h"
#include "TObjArray.h"

#include "mu2eii/stat/stat/channel.hh"


namespace mu2eii {
//-----------------------------------------------------------------------------
class analysis: public TNamed {
public:
//-----------------------------------------------------------------------------
// signal normalized to a certain proton flux
//-----------------------------------------------------------------------------
  channel*   fSignal;
  float      fRmue;
  
  TObjArray* fListOfBgrChannels;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  analysis(const char* Name, const char* Title);
  
  ~analysis();
  
  void       AddBgrChannel(channel* Channel) {
    fListOfBgrChannels->Add(Channel);
  };

  void       SetSignal(channel* Signal, float Rmue) { fSignal = Signal; fRmue = Rmue; }

  // accessors

  int      GetNBgrChannels() { return fListOfBgrChannels->GetEntriesFast(); }

  channel* GetBgrChannel(const char* Name) { return (channel*) fListOfBgrChannels->FindObject(Name); }

  channel* GetBgrChannel(int         I   ) { return (channel*) fListOfBgrChannels->UncheckedAt(I)  ; }

  channel* GetSignal() { return fSignal; }

  double   GetSES   (float PMin, float PMax, float TMin, float TMax) {
    double ses = 1./fSignal->GetIntegral(PMin,PMax,TMin,TMax);
    return ses;
  }

  double   GetBgrIntegral(float PMin, float PMax, float TMin, float TMax);
  double   GetSigIntegral(float PMin, float PMax, float TMin, float TMax);

  ClassDef(analysis,0);
};
}
#endif
