using namespace FCSys;

int model_builder(const char* file = "mu2eii/scripts/datacard.txt", bool print = false) {
  cout << "- Using file " << file << endl;
  ModelBuilder builder("builder", "builder");
  builder.SetVerbose(1);
  int status = builder.LoadModel(file);
  if(status) return status;
  Poisson_t* model = builder.GetModel();
  cout << "- Retrieved the Poisson model!\n";

  var_t& signal_mu = builder.GetPOI();
  TRandom3 rnd(90);
  model->verbose_ = 1;
  model->ngen_ = 1e5;
  model->nmax_ = std::max(20., 20*model->GetMean());

  //Generate the null hypothesis
  signal_mu.val_ = 0.;
  cout << "- Null model mean = " << model->GetMean() << endl;
  TH1D* hnull = model->GeneratePDF(rnd);
  hnull->SetName("hNullPdf");

  TCanvas* c = new TCanvas();
  hnull->SetLineWidth(2);
  hnull->SetFillStyle(3003);
  hnull->SetFillColor(kBlue);
  hnull->Draw("hist");


  cout << "- Initializing Feldman-Cousins calculator\n";
  FCCalculator fc(*model, signal_mu, rnd, 0.90/*confidence level*/, 0/*verbosity*/);
  fc.res_ = 1.e-3/fabs(signal_mu.max_ - signal_mu.min_); //roughly control the resolution of the calculation, as a fraction of POI range
  cout << "- Signal strength min - max = " << signal_mu.min_ << " - " << signal_mu.max_ << endl;
  double mu_min, mu_max;
  int nseen;
  //get the median expected events for the null hypothesis PDF
  nseen = fc.GetMedian(hnull);
  cout << "* Performing Feldman-Cousins calculation for median nseen = " << nseen << endl;
  fc.CalculateInterval(nseen, mu_min, mu_max);
  printf("** For %i seen, R_mue interval is: %.3e - %.3e **\n",
         nseen, mu_min, mu_max);

  //get median discovery information
  int ndisc = fc.NSigmaThreshold(hnull, 5.);
  cout << "* N(discovery) for NULL model = " << ndisc << endl;
  double mu_disc = fc.FindForMedianN(ndisc);
  printf("** For a median of %i, minimum R_mue is: %.3e **\n",
         ndisc, mu_disc);

  if(print) {
    //Generating a nominal signal model:
    model->spectate_ = true;
    signal_mu.val_ = 1.;
    model->GeneratePDF(rnd);
    model->PrintSpectators("model_figures");
  }

  return status;
}
