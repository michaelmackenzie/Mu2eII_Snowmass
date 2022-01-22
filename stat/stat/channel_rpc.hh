///////////////////////////////////////////////////////////////////////////////
// the only difference from the base class - initialization
//-----------------------------------------------------------------------------
#ifndef __mu2eii_stat_channel_h__
#define __mu2eii_stat_channel_h__

#include "TH1.h"
#include "TH2.h"

#include "mu2eii/stat/stat/channel.hh"

namespace mu2eii {
class  channel_rpc : public channel {
public:
// -----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  channel_rpc(const char* ChannelName, float ExtraSF = 1., int Mode = 0, int DatasetConfigCode = 0, int Verbose = 0);
};
}
#endif
