#include "Mu2e_model.C"
TGraph* gDisc_;
TGraph* gLim_;

void plot_sensitivity_vs_systematics(double scale_min = 0.2, double scale_max = 5., int nsteps = 10, double scalePrecision = 1.) {
  scalePrecision_ = scalePrecision;
  verbose_ = 0;
  double limits[nsteps+1], discoveries[nsteps+1], scales[nsteps+1];
  double tot_min(1.e10), tot_max(0.);
  for(int istep = 0; istep <= nsteps; ++istep) {
    if(nsteps < 50 || istep % 10 == 0) cout << "-- Beginning iteration " << istep << endl;
    scaleUncertainty_ = scale_min + (scale_max - scale_min)/nsteps * istep;
    Mu2e_model();
    limits[istep] = upperLimit_;
    discoveries[istep] = discovery_;
    scales[istep] = scaleUncertainty_;
    tot_min = min(min(discovery_, upperLimit_), tot_min);
    tot_max = max(max(discovery_, upperLimit_), tot_max);
  }
  TCanvas* c = new TCanvas("c_scale", "c_scale", 900, 600);
  gDisc_ = new TGraph(nsteps+1, scales, discoveries);
  gDisc_->SetName("gDiscovery");
  gLim_  = new TGraph(nsteps+1, scales, limits);
  gLim_->SetName("gLimit");

  gDisc_->SetTitle("Sensitivity vs Systematic Uncertainty Scale; Systematic uncertainty scale factor; Sensitivity");
  gDisc_->SetLineWidth(2);
  gDisc_->SetMarkerStyle(20);
  gDisc_->SetMarkerSize(0.8);
  gDisc_->SetMarkerColor(kBlue+1);
  gDisc_->SetLineColor(kBlue);
  gDisc_->Draw();

  gLim_->SetLineWidth(2);
  gLim_->SetMarkerStyle(20);
  gLim_->SetMarkerSize(0.8);
  gLim_->SetLineColor(kRed);
  gLim_->SetMarkerColor(kRed+1);
  gLim_->Draw("PL");

  TLegend* leg = new TLegend();
  leg->AddEntry(gDisc_, "Median discovery", "PL");
  leg->AddEntry(gLim_ , "Median 90% CL"   , "PL");
  leg->Draw();

  gDisc_->GetYaxis()->SetRangeUser(0.9*tot_min, 1.1*tot_max);
  c->Print("sensitivity_vs_systematics.pdf");
}
