//
#ifndef mu2eii_constants_hh
#define mu2eii_constants_hh


namespace mu2eii {
  //-----------------------------------------------------------------------------
  class constants {
  public:
    constants();
    constants(int Mode, int isMu2eIIDataset = 0);
    ~constants() {}

    double muon_capture()   { return _muon_capture  ;}
    double muon_stop_rate() { return _muon_stop_rate;}
    double npot_1b()        { return _npot_1b       ;}
    double npot_2b()        { return _npot_2b       ;}
    double cosmics_scale()  { return _cosmics_scale ;}
    double pbar_scale()     { return _pbar_scale    ;}
    double rpc_scale()      { return _rpc_scale     ;}
    double rpc_oot_scale()  { return _rpc_oot_scale ;}
    double extinction()     { return _extinction    ;}
    double mixed_scale()    { return _mixed_scale   ;}
    double unmixed_scale()  { return _unmixed_scale ;}

    double _muon_capture  ;
    double _muon_stop_rate;
    double _npot_1b       ;
    double _npot_2b       ;
    double _cosmics_scale ;
    double _pbar_scale    ;
    double _rpc_scale     ;
    double _rpc_oot_scale ;
    double _extinction    ;
    double _mixed_scale   ;
    double _unmixed_scale ;
  };
}
#endif
