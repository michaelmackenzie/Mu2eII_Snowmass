using namespace FCSys;

TRandom3* rnd_ = new TRandom3(90);
bool comparePDF_ = false; //make a plot comparing full random sample to semi-random sampling
bool doConstraints_ = false; //use systematic uncertainties in limit calculation
bool useLogNormal_ = true; //use log-normal systematic uncertainty PDFs
double scaleLuminosity_ = 1.; //scale luminosity for testing purposes
bool useCRLightYieldUnc_ = false; //consider our uncertainty on future light-yield in cosmics estimate
bool useLYSafetyFactor_ = false; //use the cosmics background estimate with the lower light yield
double scaleUncertainty_ = 1.; //artificially change the systematic uncertainties
double scalePrecision_ = 10.; //increase or decrease numerical precision
int verbose_ = 0;
bool printSpectators_ = false;
bool reduceRPC_ = true; //add an additional RPC suppression factor

//results, stored for external use
double upperLimit_;
double discovery_;
double nominalBkg_;
double nullBkg_;

void Mu2eII_model() {

  ///////////////////////////////////////////////////////
  // Initialize experimental parameters
  ///////////////////////////////////////////////////////

  //optimization output
  //    pmin      pmax      tmin      tmax   Cosmics       DIO       PbarTOT        RPC       Total       Signal        s5         SES        Rmue
  // ----------------------------------------------------------------------------------------------------------------------------------------------------
  //   104.000   104.900   680.0     1640.0  9.39716e-02  2.20949e-01 0.00000e+00  4.15179e-02  3.56438e-01  1.76602e+01  6.53214e+00  5.66245e-18  3.69879e-17

  //lower light yield
  //  103.500    104.850   620.0    1680.0  2.43785e-1   6.55531e-2  1.15306e-2  3.85590e-2  3.59428e-1  4.74246e-1  6.54537e0  2.10861e-16 1.38016e-15

  //All in terms of nominal and fractional uncertainty on the nominal
  const double dio_bkg          = (reduceRPC_) ? 2.202949e-1 : 2.20949e-1;
  const double dio_frac_unc     = scaleUncertainty_*0.70; //using systematic upper estimate with statistical error
  const double rpc_bkg          = 4.15179e-2;
  const double rpc_frac_unc     = scaleUncertainty_*0.29; //using systematic upper estimate with statistical error
  const double extinction       = 0.; //in units of 10^-10, FIXME: Get expected value
  const double rpc_oot_bkg      = 14.5e-4*extinction; //see docdb-36503 section 11.3, scaled 13.9 by (tmax - tmin) / (1695 - 700)
  const double rpc_oot_frac_unc = scaleUncertainty_*0.12; //using systematic upper estimate with statistical error
  const double pbar_bkg         = 0.;
  const double pbar_frac_unc    = scaleUncertainty_*1.; //100% uncertainty quoted - careful with Gaussian mode!

  const double cr_lo_bkg        = 9.39716e-2;
  const double cr_lo_frac_unc   = scaleUncertainty_*0.20;
  const double cr_hi_bkg        = 0.; //0.02; //see docdb-38052 slide 20: COMBINED INTO CRV_LO FOR NOW
  const double cr_hi_frac_unc   = scaleUncertainty_*0.50;

  //consider our current uncertainty on the light-yield at data-taking time as a systematic on the measurement
  const double cr_lo_ly = cr_lo_bkg;
  const double cr_hi_ly = 0.5; //FIXME: increased to not be negative with high value

  const double one_batch_pot = 0.; //one-batch mode POT, for separate batch uncertainties
  const double two_batch_pot = 4.4e22; //two-batch mode POT, for separate batch uncertainties
  const double tot_pot = one_batch_pot + two_batch_pot;

  const double muon_stop_rate = 8.9e-5;
  double ses = 5.66245e-18; //for signal strength -> R_mue, SES if acceptance was 100%
  const double signal_acceptance = 1/(tot_pot*muon_stop_rate*0.609*ses);
  ses *= signal_acceptance; //SES if 100% acceptance
  const double sig_frac_unc = scaleUncertainty_*0.04; //taken as the maximum from the momentum scale, but should be two-sided and correlated with DIO

  const double lumi_frac_unc = scaleUncertainty_*0.1; //muon stop normalization uncertainty
  const double deadtime_nom = 0.; //Deadtime already in background estimates
  const double deadtime_frac_unc = scaleUncertainty_*0.5; //assume 10 +- 5% deadtime
  const double trig_eff = 1.00; //trigger efficiency included in the background estimates already
  const double two_batch_ineff = 0.; //inefficiency due to intensity cut-off already in background estimates
  const double two_batch_unc = scaleUncertainty_*1.00; //assume two-batch mode loss is 6% +- 6%

  ///////////////////////////////////////////////////////
  // Initialize model parameters
  ///////////////////////////////////////////////////////

  /////////////////////////////////////
  // Deadtime variables
  var_t deadtime_unc("Deadtime uncertainty", deadtime_frac_unc);
  var_t deadtime_beta("Deadtime beta", 0., -10., 10.);
  var_t deadtime_var("Deadtime variation", 0., -1., 10.);
  var_t deadtime("Deadtime", deadtime_nom, 0., 1.);
  if(useLogNormal_) {
    //Log-normal via value = nominal * (1 + uncertainty = kappa) ^ (Normal distribution = beta)
    deadtime_var.nom_ = 1.; deadtime_var.val_ = 1.;
    deadtime_var.set_dependents({&deadtime_unc}, {}, {&deadtime_beta});
    deadtime.set_dependents({}, {&deadtime_var}, {});
  } else {
    deadtime_var.set_dependents({&deadtime_unc}, {&deadtime_beta}); //add uncertainty * unit width gaussian to deadtime
    deadtime.set_dependents({&deadtime_var}, {});
  }
  pow_t livetime_frac("Livetime fraction");
  var_t one("one", 1.); //place-holder
  var_t neg_one("neg_one", -1.); //place-holder
  livetime_frac.set_dependents({&one, &deadtime}, {&one, &neg_one}, {0, 1}); //set livetime = 1*1^0 + (-1)*deadtime^1 = 1 - deadtime

  /////////////////////////////////////
  // Intensity cut-off variables
  // Two-batch mode generation had a cutoff of the intensity distribution, where 12% of the
  // log-normal distribution was not included. To account for this, assume the 2-batch
  // integrated intensity has a -6% loss with +- 6% uncertainty on this (with a max of 0)

  var_t intensity_cutoff_2b_unc("Intensity cut-off two-batch uncertainty", two_batch_unc);
  var_t intensity_cutoff_2b_beta("Intensity cut-off two-batch beta", 0., -10., 10.);
  var_t intensity_cutoff_2b_var("Intensity cut-off two-batch variation", 0., -1., 10.);
  var_t intensity_cutoff_2b("Intensity cut-off 2b effect", 1. - two_batch_ineff, 0., 1.);
  var_t intensity_cutoff_1b("Intensity cut-off 1b effect", 1.); //no real effect for 1-batch mode
  if(useLogNormal_) {
    intensity_cutoff_2b_var.nom_ = 1.; intensity_cutoff_2b_var.val_ = 1.;
    intensity_cutoff_2b_var.set_dependents({&intensity_cutoff_2b_unc}, {}, {&intensity_cutoff_2b_beta});
    intensity_cutoff_2b.set_dependents({}, {&intensity_cutoff_2b_var}, {});
  } else {
    intensity_cutoff_2b_var.set_dependents({&intensity_cutoff_2b_unc}, {&intensity_cutoff_2b_beta});
    intensity_cutoff_2b.set_dependents({&intensity_cutoff_2b_var}, {});
  }
  var_t intensity_frac_1b("Intensity in one-batch mode", one_batch_pot/tot_pot);
  var_t intensity_frac_2b("Intensity in two-batch mode", two_batch_pot/tot_pot);
  pow_t intensity_frac("Intensity fraction due to cutoff effect");
  // intensity = 1-batch intensity * cutoff effect + 2-batch intensity * cutoff effect
  intensity_frac.set_dependents({&intensity_cutoff_1b, &intensity_cutoff_2b}, {&intensity_frac_1b, &intensity_frac_2b}, {1,1});

  /////////////////////////////////////
  // Trigger efficiency variables
  var_t trigger_eff("Trigger efficiency", trig_eff);

  /////////////////////////////////////
  // Luminosity variables
  var_t lumi_unc("Luminosity uncertainty", lumi_frac_unc);
  var_t lumi_beta("Luminosity beta", 0., -10., 10.);
  var_t lumi_var("Luminosity variation", 0., -1., 10.);
  var_t lumi("Luminosity", (scaleLuminosity_ > 0.) ? scaleLuminosity_ : 1., 0., 3.);
  if(useLogNormal_) {
    lumi_var.nom_ = 1.; lumi_var.val_ = 1.;
    lumi_var.set_dependents({&lumi_unc}, {}, {&lumi_beta}); //multiply (1 + uncertainty) ^ unit width gaussian to lumi prediction
    lumi.set_dependents({}, {&lumi_var, &livetime_frac, &intensity_frac, &trigger_eff}, {}); // lumi = nominal * kappa ^ beta; kappa = 1 + uncertainty, beta = normal, width 1
    //{add}, {multiply}, {powers} --> value = ((sum of adds) * (product of mul)) ^ {product of powers}
  } else {
    lumi_var.set_dependents({&lumi_unc}, {&lumi_beta}); //add uncertainty * unit width gaussian to lumi prediction
    lumi.set_dependents({&lumi_var}, {&livetime_frac, &intensity_frac, &trigger_eff});
  }

  var_t dio_unc("DIO uncertainty", dio_bkg*dio_frac_unc);
  var_t dio_beta("DIO beta", 0., -10., 10.);
  var_t dio_var("DIO variation", 0., -1.*dio_bkg, 5.);
  var_t dio("DIO expectation", dio_bkg, 0., std::max(1., dio_bkg + 10.*sqrt(dio_bkg)));
  if(useLogNormal_) {
    dio_var.nom_ = 1.; dio_var.val_ = 1.;
    dio_var.set_dependents({&dio_unc}, {}, {&dio_beta});
    dio.set_dependents({}, {&dio_var, &lumi});
  } else {
    dio_var.set_dependents({&dio_unc}, {&dio_beta});
    dio.set_dependents({&dio_var}, {&lumi});
  }

  var_t rpc_unc("RPC uncertainty", rpc_bkg*rpc_frac_unc);
  var_t rpc_beta("RPC beta", 0., -10., 10.);
  var_t rpc_var("RPC variation", 0., -1.*rpc_bkg, 5.);
  var_t rpc("RPC expectation", rpc_bkg, 0., std::max(1., rpc_bkg + 10.*sqrt(rpc_bkg)));
  if(useLogNormal_) {
    rpc_var.nom_ = 1.; rpc_var.val_ = 1.;
    rpc_var.set_dependents({&rpc_unc}, {}, {&rpc_beta});
    rpc.set_dependents({}, {&rpc_var, &lumi});
  } else {
    rpc_var.set_dependents({&rpc_unc}, {&rpc_beta});
    rpc.set_dependents({&rpc_var}, {&lumi});
  }

  var_t rpc_oot_unc("RPC OOT uncertainty", rpc_oot_bkg*rpc_oot_frac_unc);
  var_t rpc_oot_beta("RPC OOT beta", 0., -10., 10.);
  var_t rpc_oot_var("RPC OOT variation", 0., -1.*rpc_oot_bkg, 5.);
  var_t rpc_oot("RPC OOT expectation", rpc_oot_bkg, 0., 5.);
  if(useLogNormal_) {
    rpc_oot_var.nom_ = 1.; rpc_oot_var.val_ = 1.;
    rpc_oot_var.set_dependents({&rpc_oot_unc}, {}, {&rpc_oot_beta});
    rpc_oot.set_dependents({}, {&rpc_oot_var, &lumi});
  } else {
    rpc_oot_var.set_dependents({&rpc_oot_unc}, {&rpc_oot_beta});
    rpc_oot.set_dependents({&rpc_oot_var}, {&lumi});
  }

  var_t pbar_unc("Pbar uncertainty", pbar_bkg*pbar_frac_unc);
  var_t pbar_beta("Pbar beta", 0., -10., 10.);
  var_t pbar_var("Pbar variation", 0., -1.*pbar_bkg, 5.);
  var_t pbar("Pbar expectation", pbar_bkg, 0., std::max(1., 11.*pbar_bkg));
  if(useLogNormal_) {
    pbar_var.nom_ = 1.; pbar_var.val_ = 1.;
    pbar_var.set_dependents({&pbar_unc}, {}, {&pbar_beta});
    pbar.set_dependents({}, {&pbar_var, &lumi});
  } else {
    pbar_var.set_dependents({&pbar_unc}, {&pbar_beta});
    pbar.set_dependents({&pbar_var}, {&lumi});
  }

  var_t cr_light_yd("Cosmic-ray veto light yield", 0., 0., 1., "Flat"); //consider a flat light yield uncertainty, in arbitrary units
  var_t cr_light_coeff("Cosmic-ray light yield linear coeff", 2.*(cr_hi_ly - cr_lo_ly));
  var_t cr_light_offset("Cosmic-ray light yield linear offset", -1.*(cr_hi_ly - cr_lo_ly));
  pow_t cr_light_var("Cosmic-ray variation due to light yield"); //variation of the yield, linear fit to measured points
  var_t cr_light_bkg("Cosmic-ray background due to light yield", cr_lo_ly, cr_lo_ly, 10.*cr_hi_ly); //delta function at min, flat above
  cr_light_var.set_dependents({&one, &cr_light_yd}, {&cr_light_offset, &cr_light_coeff}, {0, 1});
  cr_light_bkg.set_dependents({&cr_light_var}, {});

  var_t cr_lo_unc("Cosmic-ray (lo) uncertainty", cr_lo_frac_unc);
  var_t cr_lo_beta("Cosmic-ray (lo) beta", 0., -10., 10.);
  var_t cr_lo_var("Cosmic-ray (lo) variation", 0., -1.*cr_lo_bkg, 5.);
  var_t cr_lo("Cosmic-ray (lo) expectation", 0., 0., std::max(1., cr_lo_bkg + 10.*sqrt(cr_lo_bkg)));
  if(useLogNormal_) {
    cr_lo_var.nom_ = 1.; cr_lo_var.val_ = 1.;
    cr_lo_var.set_dependents({&cr_lo_unc}, {}, {&cr_lo_beta});
    cr_lo.set_dependents({&cr_light_bkg}, {&cr_lo_var, &livetime_frac, &trigger_eff});
  } else {
    cr_lo_var.set_dependents({&cr_lo_unc}, {&cr_light_bkg, &cr_lo_beta});
    cr_lo.set_dependents({&cr_lo_var}, {&livetime_frac, &trigger_eff});
  }
  if(!useCRLightYieldUnc_) {
    cr_light_bkg.set_val(cr_lo_bkg);
    cr_light_bkg.set_constant();
    cr_light_yd.set_constant();
  } else {
    cr_lo_beta.set_constant();
  }

  var_t cr_hi_unc("Cosmic-ray (hi) uncertainty", cr_hi_bkg*cr_hi_frac_unc);
  var_t cr_hi_beta("Cosmic-ray (hi) beta", 0., -10., 10.);
  var_t cr_hi_var("Cosmic-ray (hi) variation", 0., -1.*cr_hi_bkg, 5.);
  var_t cr_hi("Cosmic-ray (hi) expectation", cr_hi_bkg, 0., std::max(1., cr_hi_bkg + 10.*sqrt(cr_hi_bkg)));
  if(useLogNormal_) {
    cr_hi_var.nom_ = 1.; cr_hi_var.val_ = 1.;
    cr_hi_var.set_dependents({&cr_hi_unc}, {}, {&cr_hi_beta});
    cr_hi.set_dependents({}, {&cr_hi_var, &livetime_frac, &trigger_eff});
  } else {
    cr_hi_var.set_dependents({&cr_hi_unc}, {&cr_hi_beta});
    cr_hi.set_dependents({&cr_hi_var}, {&livetime_frac, &trigger_eff});
  }

  var_t sig_mu("Nominal signal strength", 0., 0., 10./signal_acceptance);
  var_t sig_eff("Signal efficiency", signal_acceptance);
  var_t sig_unc("Signal uncertainty", sig_frac_unc);
  var_t sig_beta("Signal beta", 0., -10., 10.);
  var_t sig_var("Signal variation", 0., -1., 5.);
  var_t signal("Signal expectation", 0., 0., 20.);
  if(useLogNormal_) {
    sig_var.nom_ = 1.; sig_var.val_ = 1.;
    sig_var.set_dependents({&sig_unc}, {}, {&sig_beta});
    signal.set_dependents({&sig_mu}, {&sig_var, &sig_eff, &lumi});
  } else {
    sig_var.set_dependents({&sig_unc}, {&sig_beta, &sig_mu});
    signal.set_dependents({&sig_mu, &sig_var}, {&sig_eff, &lumi});
  }

  if(verbose_ > 0) {
    cout << "Nominal values:\n";
    lumi.print();
    lumi_unc.print();
    lumi_var.print();
    dio.print();
    dio_unc.print();
    dio_var.print();
    rpc.print();
    rpc_unc.print();
    rpc_var.print();
    rpc_oot.print();
    rpc_oot_unc.print();
    rpc_oot_var.print();
    pbar.print();
    pbar_unc.print();
    pbar_var.print();
    cr_lo.print();
    cr_lo_unc.print();
    cr_lo_var.print();
    cr_hi.print();
    cr_hi_unc.print();
    cr_hi_var.print();
    sig_eff.print();
    sig_unc.print();
    sig_var.print();
    signal.print();
    deadtime_unc.print();
    deadtime_var.print();
    deadtime.print();
    livetime_frac.print();
    intensity_cutoff_1b.print();
    intensity_cutoff_2b.print();
    intensity_frac_1b.print();
    intensity_frac_2b.print();
    intensity_frac.print();
    trigger_eff.print();
  }

  if(!doConstraints_) {
    cout << "Setting systematics to 0!\n";
    lumi_beta.set_constant();
    dio_beta.set_constant();
    rpc_beta.set_constant();
    rpc_oot_beta.set_constant();
    pbar_beta.set_constant();
    cr_light_yd.set_constant();
    cr_lo_beta.set_constant();
    cr_hi_beta.set_constant();
    sig_beta.set_constant();
    deadtime_beta.set_constant();
    intensity_cutoff_2b_beta.set_constant();
  }
  ///////////////////////////////////////////////////////
  // Initialize model
  ///////////////////////////////////////////////////////

  var_t nobs("Number observed", 0., 0., 20.);
  Poisson_t model("Counting model", nobs,
                  {&dio     , &rpc     , &rpc_oot     , &pbar     , &cr_lo     , &cr_hi     , &signal},
                  {&dio_beta, &rpc_beta, &rpc_oot_beta, &pbar_beta, &cr_lo_beta, &cr_hi_beta, &sig_beta,
                   &deadtime_beta, &intensity_cutoff_2b_beta, &lumi_beta, &cr_light_yd
                   },
                  {&signal, &sig_mu, &dio, &rpc, &rpc_oot, &pbar, &cr_lo, &cr_hi, &lumi, &deadtime, } //variables to plot values for
                 );
  model.ngen_ = (doConstraints_) ? scalePrecision_*1e5 : 1; //without constraints, it's just a Poisson PDF, no random shifting

  cout << "Model:\n";
  model.Print();
  double nexp_bkg = dio_bkg + rpc_bkg + rpc_oot_bkg + pbar_bkg + cr_lo_bkg + cr_hi_bkg;
  cout << "Nominal background mean = " << nexp_bkg*trig_eff << endl;
  cout << "Generating the NULL observable PDF\n";
  if(auto o = gDirectory->Get("null")) delete o;
  TH1D* hpdf = model.GeneratePDF(*rnd_);
  hpdf->SetName("null");
  cout << "Poisson model has a nominal mean of: " << model.GetNominalMean()
       << " and NULL PDF has a mean of: " << (hpdf->GetMean()) << endl;
  nominalBkg_ = model.GetNominalMean();
  nullBkg_ = hpdf->GetMean();

  TCanvas* c = new TCanvas();
  hpdf->SetLineColor(kRed);
  hpdf->SetLineWidth(2);
  hpdf->SetMarkerColor(kRed);
  hpdf->SetMarkerStyle(20);
  hpdf->SetMarkerSize(0.8);
  hpdf->Draw();
  hpdf->SetAxisRange(1.e-20, 10., "Y");
  hpdf->SetAxisRange(0, 19.9, "X");
  c->SetLogy();

  if(comparePDF_) {
    cout << "Randomly sampling the NULL model PDF\n";
    const int nentries = 1e5;
    TH1D* hexp = new TH1D("hexp", "", 20, 0, 20);
    for(int i = 0; i < nentries; ++i) {
      if(verbose_ && (i % (nentries/5) == 0)) {
        cout << "Samping " << i << ":\n";
        model.SetVerbose(10);
      } else model.SetVerbose(0);
      int n = model.RandomSample(*rnd_);
      if(i % (nentries/10) == 0) {
        signal.print();
        dio.print();
        dio_var.print();
        cout << endl;
      }

      hexp->Fill(n);
    }
    TCanvas* c2 = new TCanvas("c2", "c2", 1000, 700);
    TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1., 1.0); pad1->Draw();
    TPad* pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.3); pad2->Draw();
    pad1->cd();
    pad1->SetLogy();
    TH1D* hpoisson = new TH1D("hpoisson", "", 20, 0, 20);
    for(int i = 0; i < 20; ++i) hpoisson->Fill(i, ROOT::Math::poisson_pdf(i, nullBkg_));
    hpdf->Draw();
    hexp->Scale(1./nentries);
    hexp->Draw("hist E1 sames");
    hexp->SetLineWidth(2);
    hpoisson->Draw("hist sames");
    hpoisson->SetLineWidth(2);
    hpoisson->SetLineColor(kGreen);
    pad2->cd();
    TH1D* hexp_r = (TH1D*) hexp->Clone("hexp_r");
    TH1D* hpoisson_r = (TH1D*) hpoisson->Clone("hpoisson_r");
    for(int bin = 1; bin <= hexp_r->GetNbinsX(); ++bin) {
      double binval = hpdf->GetBinContent(bin);
      double binerr = hpdf->GetBinError(bin);
      hexp_r->SetBinContent(bin, (binval > 0.) ? hexp->GetBinContent(bin) / binval : -1.);
      hexp_r->SetBinError(bin, (binval > 0.) ? hexp->GetBinError(bin) / binval : 0.);
      hpoisson_r->SetBinContent(bin, (binval > 0.) ? hpoisson->GetBinContent(bin) / binval : -1.);
      hpoisson_r->SetBinError(bin, (binval > 0.) ? binerr / binval : 0.);
    }
    hexp_r->Draw("hist E1");
    hexp_r->GetYaxis()->SetRangeUser(0.,2.);
    hpoisson_r->Draw("hist E1 same");
    TLine* line = new TLine(0., 1., 20., 1.);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");
    pad2->Update();
    c2->Update();
    c->cd();
  }

  ///////////////////////////////////////////////////////
  // Initialize Feldman-Cousins calculator for the model
  ///////////////////////////////////////////////////////

  cout << "Initializing Feldman-Cousins calculator\n";
  FCCalculator fc(model, sig_mu, *rnd_, 0.90/*confidence level*/, 0/*verbosity*/);
  fc.res_ = 2.e-4/scalePrecision_; //roughly control the resolution of the calculation, as a fraction of POI range

  double mu_min, mu_max;
  int nseen;
  //get the median expected events for the null hypothesis PDF
  nseen = fc.GetMedian(hpdf);
  cout << "Using a single event sensitivity of: " << ses/signal_acceptance << endl;
  cout << "Performing Feldman-Cousins calculation for median nseen = " << nseen << endl;
  fc.verbose_ = verbose_;
  fc.CalculateInterval(nseen, mu_min, mu_max);
  printf("For %i seen, R_mue interval is: %.3e - %.3e (%.3f - %.3f mean events = mu * (nominal signal eff))\n",
         nseen, mu_min*ses, mu_max*ses,
         mu_min*sig_eff.nom_*scaleLuminosity_, mu_max*sig_eff.nom_*trig_eff*scaleLuminosity_);
  upperLimit_ = mu_max*ses;

  //get median discovery information
  int ndisc = fc.NSigmaThreshold(hpdf, 5.);
  cout << "N(discovery) for NULL model = " << ndisc << endl;
  double mu_disc = fc.FindForMedianN(ndisc);
  printf("For a median of %i, minimum R_mue is: %.3e (%.3f mean events = mu * (nominal signal eff))\n",
         ndisc, mu_disc*ses, mu_disc*sig_eff.nom_*trig_eff*scaleLuminosity_);

  discovery_ = mu_disc*ses;
  sig_mu.val_ = mu_max;
  TH1D* hUL = model.GeneratePDF(*rnd_);
  hUL->SetName("hUL");

  sig_mu.val_ = mu_disc;
  model.spectate_ = true; //make plots of some variables
  model.specBins_ = 200;
  TH1D* hDisc = model.GeneratePDF(*rnd_);
  hDisc->SetName("hDisc");

  hUL->SetLineColor(kBlue);
  hUL->SetLineWidth(2);
  hUL->Draw("same");

  hDisc->SetLineColor(kGreen+2);
  hDisc->SetLineWidth(2);
  hDisc->Draw("same");

  TLegend* leg = new TLegend();
  leg->AddEntry(hpdf, "Null PDF");
  leg->AddEntry(hUL, "Upper limit PDF");
  leg->AddEntry(hDisc, "Discovery PDF");
  leg->Draw();

  if(printSpectators_) {
    gSystem->Exec("[ ! -d figures ] && mkdir figures");
    for(auto spectator : model.spec_) {
      TCanvas* c = new TCanvas();
      model.specMap_[spectator]->Draw();
      c->Print(Form("figures/c_%s.png", spectator->name_.Data()));
      delete c;
    }
  }
}
