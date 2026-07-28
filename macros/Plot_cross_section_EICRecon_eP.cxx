#include<iostream>
#include<fstream>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TF1.h"
#include"TGraphErrors.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TChain.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TAxis.h"
#include"TTree.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"
#include"TLorentzVector.h"
#include"TDatabasePDG.h"
#include"TParticlePDG.h"
#include"TRandom.h"
#include"TLine.h"

using namespace std;
//for tracking efficiency definition see D0 analysis note, page 24

const double M_D = 1.875612945;//deuteron mass in GeV

const float lumi_design_fb = 1.5; //expected lumi. in fb-1

//paper luminosity: sigma_tot = 4.5e-5 mb
//new BeAGLE production:
//cross section from PYTHIA/BeAGLE log files

//10x130 GeV
const float sigma_BeAGLE = 3.707e-4; //mb (eP)
const float sigma_BeAGLE_fb = 3.707e8; //mb (eP)

void Plot_cross_section_EICRecon_eP(const int plot_RC_flag = 1, const int plot_ACC_flag = 1, TString inFileName = "Input")
{
  //check that we have input file to analyze
  if( !inFileName.Contains(".root") )
  {
    cout<<"Missing input file name!"<<endl;
    cout<<"Please, provide input root file as an argument of run script."<<endl;
    
    return;  
  }
  
  //-------------------------------------------------------------

  TDatabasePDG *myPDGdatabase = new TDatabasePDG();

  TParticlePDG *e_PDG = myPDGdatabase->GetParticle(11);
  TParticlePDG *p_PDG = myPDGdatabase->GetParticle(2212);
  TParticlePDG *n_PDG = myPDGdatabase->GetParticle(2112);

  //aT2 evaluated at alpha_p = 1 (center of alpha_p bin)
  const float aT2 = n_PDG->Mass()*n_PDG->Mass() - M_D*M_D/4.;

  //---------------------------------------------------------------

  TRandom *my_random = new TRandom();//check if we need to set a specific seed

  //---------------------------------------------------------------

  float sys_err_lumi = 0.015; //1.5% sys. err. of luminosity
  float sys_err_theory = 0.02; //2% for theory - from old paper
  float sys_err_beam_eff = 0.015; //1.5% for beam effects (crossing angle, angular divergence, beam energy spread)
  float sys_err_neutron = 0.01; //1% on PV error for neutron reconstruction

  float sys_err_smear = sys_err_neutron;

  float sys_err_tot_corr = sqrt(sys_err_lumi*sys_err_lumi + sys_err_theory*sys_err_theory  + sys_err_beam_eff*sys_err_beam_eff + sys_err_neutron*sys_err_neutron);

  //---------------------------------------------------------------

  const int nx_bins = 11;
  double x_bins[nx_bins+1] = {0.001, 0.002, 0.004, 0.007, 0.01, 0.02, 0.04, 0.07, 0.09, 0.2, 0.3, 0.5};

  const int nQ2_bins = 4;
  float Q2_bins[nQ2_bins+1] = {1., 10., 20., 40., 100.};

  float alpha_p_min = 0.99;
  float alpha_p_max = 1.01;
  float Delta_apha_p = alpha_p_max - alpha_p_min;

  float two_times_2pi_cubed = 2.*(2.*TMath::Pi())*(2.*TMath::Pi())*(2.*TMath::Pi());

  //---------------------------------------------------------------


  TFile *inFile = new TFile(inFileName, "read");

  
  ofstream outputFile("./output/p_cross_section.txt");
  
  //maybe remove from here? Results produced using eN code
  ofstream outputFile_D("./output/D_cross_section_eP.txt");

  if (!outputFile.is_open()) {
    cout<< "Error: Could not open or create the file!" <<endl;
    return;
  }

  if (!outputFile_D.is_open())
  {
    cout<< "Error: Could not open or create the file!" <<endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //event stats
	TH1F *h_nEvents = (TH1F*)inFile->Get("h_nEvents");

	TH1F *h_nEvents_bins[nx_bins][nQ2_bins];

  //deuteron reduced crossection

	TH1F *h_eff_nucleon_func_bins[nx_bins][nQ2_bins];

	//----------------------------

	TH1F *h_pT2_bins[nx_bins][nQ2_bins];
  TH1F *h_sigma_D_red[nx_bins][nQ2_bins];
  TH1F *h_sigma_p_red[nx_bins][nQ2_bins];

  TH1F *h_N_vs_x[nQ2_bins];
  TH1F *h_sigma_D_red_vs_x[nQ2_bins];

  //this is for clone of TH1F *h_sigma_p_red[nx_bins][nQ2_bins]; for new systematic errors
  //Needs to be TH1D for Clone
  const int n_variations = 10; //number of variations for sys. err.

  //----------------------------

  TH1F *h_pT2_RC_bins[nx_bins][nQ2_bins];
  TH1F *h_sigma_D_red_RC[nx_bins][nQ2_bins];
  TH1F *h_sigma_p_red_RC[nx_bins][nQ2_bins];

  TH1F *h_N_vs_x_RC[nQ2_bins];
  TH1F *h_sigma_D_red_vs_x_RC[nQ2_bins];

  //----------------------------

  TH1F *h_pT2_MC_ACC_bins[nx_bins][nQ2_bins];
  TH1F *h_sigma_D_red_MC_ACC[nx_bins][nQ2_bins];
  TH1F *h_sigma_p_red_MC_ACC[nx_bins][nQ2_bins];

  TH1F *h_N_vs_x_MC_ACC[nQ2_bins];
  TH1F *h_sigma_D_red_vs_x_MC_ACC[nQ2_bins];

  //----------------------------

  TGraphErrors *sigma_free_p_vs_x[nQ2_bins];
  TGraphErrors *sigma_free_p_vs_x_RC_corr[nQ2_bins];

  //----------------------------

  TPaveText *text = new TPaveText(0.65, 0.7, 0.95, 0.89, "NDC"); //main results;
  text->SetTextFont(43);
  text->SetTextSize(30);
  text->AddText("ePIC simulation");
  text->AddText("eD 10x130 GeV^{2}");
  text->AddText("#gamma* + d #rightarrow X + n'");
  text->AddText("Tagged neutron");
  text->SetBorderSize(0);
  text->SetFillColorAlpha(0, 0.01);


  TCanvas *sigma_red_p_vs_x_one_can = new TCanvas("sigma_red_p_vs_x_one_can", "sigma_red_p_vs_x_one_can", 2400,2000);
  sigma_red_p_vs_x_one_can->Divide(2,2);

  TCanvas *sigma_red_p_vs_x_one_can_2 = new TCanvas("sigma_red_p_vs_x_one_can_2", "sigma_red_p_vs_x_one_can_2", 4800,1000);
  TPad *pad_0 = new TPad("pad_0", "pad_0", 0.0, 0.0, 0.28, 1.0);
  pad_0->Draw();
  TPad *pad_1 = new TPad("pad_1", "pad_1", 0.28, 0.0, 0.52, 1.0);
  pad_1->Draw();
  TPad *pad_2 = new TPad("pad_2", "pad_2", 0.52, 0.0, 0.76, 1.0);
  pad_2->Draw();
  TPad *pad_3 = new TPad("pad_3", "pad_3", 0.76, 0.0, 1.0, 1.0);
  pad_3->Draw();

  //-------------------

  TCanvas *sigma_red_p_vs_x_one_can_RC_corr = new TCanvas("sigma_red_p_vs_x_one_can_RC_corr", "sigma_red_p_vs_x_one_can_RC_corr", 2400,2000);
  sigma_red_p_vs_x_one_can_RC_corr->Divide(2,2);

  TCanvas *sigma_red_p_vs_x_one_can_RC_corr_2 = new TCanvas("sigma_red_p_vs_x_one_can_RC_corr_2", "sigma_red_p_vs_x_one_can_RC_corr_2", 4800,1000);
  TPad *pad_0_RC = new TPad("pad_0_RC", "pad_0_RC", 0.0, 0.0, 0.28, 1.0);
  pad_0_RC->Draw();
  TPad *pad_1_RC = new TPad("pad_1_RC", "pad_1_RC", 0.28, 0.0, 0.52, 1.0);
  pad_1_RC->Draw();
  TPad *pad_2_RC = new TPad("pad_2_RC", "pad_2_RC", 0.52, 0.0, 0.76, 1.0);
  pad_2_RC->Draw();
  TPad *pad_3_RC = new TPad("pad_3_RC", "pad_3_RC", 0.76, 0.0, 1.0, 1.0);
  pad_3_RC->Draw();

  //-----------------------------

  //convert to cross section
  //unit conversion: mb = 2.568 GeV^-2
  //sigma_tot = 4.5e-5 mb
  //float lumi_factor = 2.568*4.5e-5/iEvent;

  int iEvent_MC = h_nEvents->GetBinContent(2);
  int iEvent_RC = h_nEvents->GetBinContent(3);

  float lumi_factor_MC = 2.568*sigma_BeAGLE/iEvent_MC;
  float lumi_factor_RC = 2.568*sigma_BeAGLE/iEvent_RC;

  float lumi_fb = iEvent_MC/sigma_BeAGLE_fb;

  float statistics_scale = lumi_design_fb/lumi_fb/2.; //design lumi is for eN+eP, we have just 1/2

  for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
  {
    outputFile<<Q2_bins[Q2_bin]<<" < Q2 < "<<Q2_bins[Q2_bin+1]<<" GeV^2"<<endl;
    outputFile<<"x  sigma_p,red  stat_err  tot_err"<<endl;

    outputFile_D<<Q2_bins[Q2_bin]<<" < Q2 < "<<Q2_bins[Q2_bin+1]<<" GeV^2"<<endl;
    outputFile_D<<"x  sigma_D,red  stat_err  sys_err"<<endl;

    sigma_free_p_vs_x[Q2_bin] = new TGraphErrors(nx_bins);
    sigma_free_p_vs_x[Q2_bin]->SetNameTitle(Form("g_sigma_p_red_vs_x_Q2bin_%i", Q2_bin));

    sigma_free_p_vs_x_RC_corr[Q2_bin] = new TGraphErrors(nx_bins);
    sigma_free_p_vs_x_RC_corr[Q2_bin]->SetNameTitle(Form("g_sigma_p_red_vs_x_Q2bin_%i_RC_corr", Q2_bin));

    //------

    h_sigma_D_red_vs_x[Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_D_red_vs_x_%i", Q2_bin));
    h_sigma_D_red_vs_x_RC[Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_D_red_vs_x_RC_%i", Q2_bin));
    h_sigma_D_red_vs_x_MC_ACC[Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_D_red_vs_x_MC_ACC_%i", Q2_bin));

    //------

    h_N_vs_x[Q2_bin] = (TH1F*)inFile->Get(Form("h_N_vs_x_%i", Q2_bin));
    h_N_vs_x_RC[Q2_bin] = (TH1F*)inFile->Get(Form("h_N_vs_x_RC_%i", Q2_bin));
    h_N_vs_x_MC_ACC[Q2_bin] = (TH1F*)inFile->Get(Form("h_N_vs_x_MC_ACC_%i", Q2_bin));

    //------

    float Q_bin_center = (Q2_bins[Q2_bin+1] + Q2_bins[Q2_bin])/2.;
    float Delta_Q2 = Q2_bins[Q2_bin+1] - Q2_bins[Q2_bin];

    int n_fill = 0;

    for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
    {
      float x_bin_center = (x_bins[x_bin+1] + x_bins[x_bin])/2.;
      float Delta_xbj = x_bins[x_bin+1] - x_bins[x_bin];

      h_pT2_bins[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_pT2_bins_xbin_%i_Q2bin_%i", x_bin, Q2_bin));
      h_pT2_RC_bins[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_pT2_RC_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin));
	    h_pT2_MC_ACC_bins[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_pT2_MC_ACC_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin));

	    //---------------

      //first load cross section histogram to determine, if there's
      h_sigma_D_red[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_D_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin));
      h_sigma_D_red_RC[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_D_red_RC_xbin_%i_Q2bin_%i", x_bin, Q2_bin));
      h_sigma_D_red_MC_ACC[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_D_red_MC_ACC_xbin_%i_Q2bin_%i", x_bin, Q2_bin));

      //---------------

      h_sigma_p_red[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_p_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin));
      h_sigma_p_red_RC[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_p_red_RC_xbin_%i_Q2bin_%i", x_bin, Q2_bin));
      h_sigma_p_red_MC_ACC[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_sigma_p_red_MC_ACC_xbin_%i_Q2bin_%i", x_bin, Q2_bin));

      //-----------------

      float reduce_factor_2 = 1./( Delta_xbj*Delta_Q2 ); //integrated over phi_p and pT - removed factor of 2 (dpT2/2)?

      float old_bin_content_MC = h_sigma_D_red_vs_x[Q2_bin]->GetBinContent(x_bin+1);

      float relat_stat_err_x_MC = 0;
      if(old_bin_content_MC > 0) relat_stat_err_x_MC = h_sigma_D_red_vs_x[Q2_bin]->GetBinError(x_bin+1)/old_bin_content_MC;

      //---

      float old_bin_content_RC = h_sigma_D_red_vs_x_RC[Q2_bin]->GetBinContent(x_bin+1);

      float relat_stat_err_x_RC = 0;
      if(old_bin_content_RC > 0) relat_stat_err_x_RC = h_sigma_D_red_vs_x_RC[Q2_bin]->GetBinError(x_bin+1)/old_bin_content_RC;


      h_sigma_D_red_vs_x[Q2_bin]->SetBinContent(x_bin+1, old_bin_content_MC * reduce_factor_2 * lumi_factor_MC);
      h_sigma_D_red_vs_x[Q2_bin]->SetBinError(x_bin+1, relat_stat_err_x_MC * h_sigma_D_red_vs_x[Q2_bin]->GetBinContent(x_bin+1));

      //---

      h_sigma_D_red_vs_x_RC[Q2_bin]->SetBinContent(x_bin+1, old_bin_content_RC * reduce_factor_2 * lumi_factor_MC);
      h_sigma_D_red_vs_x_RC[Q2_bin]->SetBinError(x_bin+1, relat_stat_err_x_RC * h_sigma_D_red_vs_x_RC[Q2_bin]->GetBinContent(x_bin+1));

      //---

      //statistical error rescaling
      float old_err_D_RC = 1./sqrt(h_N_vs_x_RC[Q2_bin]->GetBinContent(x_bin+1));
      float new_err_D_RC = 1./sqrt(h_N_vs_x_RC[Q2_bin]->GetBinContent(x_bin+1)*statistics_scale);
      float scale_err_D_RC = new_err_D_RC/old_err_D_RC;
      //float scale_err_D_RC = 1.;

      //---

      float sigma_D_red_MC = h_sigma_D_red_vs_x[Q2_bin]->GetBinContent(x_bin+1);
      float sigma_D_red_RC_relat_err = h_sigma_D_red_vs_x_RC[Q2_bin]->GetBinError(x_bin+1)*scale_err_D_RC/h_sigma_D_red_vs_x_RC[Q2_bin]->GetBinContent(x_bin+1);

      float sigma_D_red_MC_err = sigma_D_red_MC*sigma_D_red_RC_relat_err;

      outputFile_D<<x_bin_center<<" "<<sigma_D_red_MC<<" "<<sigma_D_red_MC_err<<" "<<sigma_D_red_MC*sys_err_tot_corr<<endl;

      //-----------------

      if( h_sigma_p_red[x_bin][Q2_bin]->Integral() < 1 ) continue; //skip bins with too low statistics

      float Delta_pT2 = h_sigma_D_red[x_bin][Q2_bin]->GetBinWidth(1);

      float reduce_factor = two_times_2pi_cubed*2./( Delta_xbj*Delta_Q2*Delta_apha_p*Delta_pT2*2.*TMath::Pi() ); //integrated over phi_p

      n_fill++;

      //-------------------------------------

      //scale all errors to the expected luminosity
      for( unsigned int pT2_bin = 1; pT2_bin <= h_pT2_bins[x_bin][Q2_bin]->GetNbinsX(); pT2_bin++)
      {
        //true MC
        float old_err_MC = 1./sqrt(h_pT2_bins[x_bin][Q2_bin]->GetBinContent(pT2_bin));
        float new_err_MC = 1./sqrt(h_pT2_bins[x_bin][Q2_bin]->GetBinContent(pT2_bin)*statistics_scale);
        float scale_err_MC = new_err_MC/old_err_MC;

        h_sigma_D_red[x_bin][Q2_bin]->SetBinError(pT2_bin, h_sigma_D_red[x_bin][Q2_bin]->GetBinError(pT2_bin)*scale_err_MC);
        h_sigma_p_red[x_bin][Q2_bin]->SetBinError(pT2_bin, h_sigma_p_red[x_bin][Q2_bin]->GetBinError(pT2_bin)*scale_err_MC);

        //---

        //acc. only
        float old_err_MC_ACC = 1./sqrt(h_pT2_MC_ACC_bins[x_bin][Q2_bin]->GetBinContent(pT2_bin));
        float new_err_MC_ACC = 1./sqrt(h_pT2_MC_ACC_bins[x_bin][Q2_bin]->GetBinContent(pT2_bin)*statistics_scale);
        float scale_err_MC_ACC = new_err_MC_ACC/old_err_MC_ACC;

        h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->SetBinError(pT2_bin, h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->GetBinError(pT2_bin)*scale_err_MC_ACC);
        h_sigma_p_red_MC_ACC[x_bin][Q2_bin]->SetBinError(pT2_bin, h_sigma_p_red_MC_ACC[x_bin][Q2_bin]->GetBinError(pT2_bin)*scale_err_MC_ACC);

        //---

        //RC
        float old_err_RC = 1./sqrt(h_pT2_RC_bins[x_bin][Q2_bin]->GetBinContent(pT2_bin));
        float new_err_RC = 1./sqrt(h_pT2_RC_bins[x_bin][Q2_bin]->GetBinContent(pT2_bin)*statistics_scale);
        float scale_err_RC = new_err_RC/old_err_RC;

        h_sigma_D_red_RC[x_bin][Q2_bin]->SetBinError(pT2_bin, h_sigma_D_red_RC[x_bin][Q2_bin]->GetBinError(pT2_bin)*scale_err_RC);
        h_sigma_p_red_RC[x_bin][Q2_bin]->SetBinError(pT2_bin, h_sigma_p_red_RC[x_bin][Q2_bin]->GetBinError(pT2_bin)*scale_err_RC);

      }

      //----------------

      //number of events per Q2 and x bin
      //for correction of RC/MC events in each bing

      h_nEvents_bins[x_bin][Q2_bin] = (TH1F*)inFile->Get(Form("h_nEvents_bins_xbin_%i_Q2bin_%i", x_bin, Q2_bin));

      float nEventCorr = h_nEvents_bins[x_bin][Q2_bin]->GetBinContent(3)/h_nEvents_bins[x_bin][Q2_bin]->GetBinContent(2);

	    //---------------

      TCanvas *sigma_red_D_can = new TCanvas(Form("sigma_red_D_can_%i_%i",x_bin, Q2_bin), Form("sigma_red_D_can_%i_%i",x_bin, Q2_bin), 1200,1000);
      sigma_red_D_can->cd();

      gPad->SetLogy();

      h_sigma_D_red[x_bin][Q2_bin]->Scale( reduce_factor * lumi_factor_MC );
      h_sigma_D_red[x_bin][Q2_bin]->GetXaxis()->SetTitle("p^{2}_{T} (GeV/c)^{2}");
      h_sigma_D_red[x_bin][Q2_bin]->GetXaxis()->CenterTitle();
      h_sigma_D_red[x_bin][Q2_bin]->GetYaxis()->SetTitle("#bar{#sigma}_{red,d}");
      h_sigma_D_red[x_bin][Q2_bin]->GetYaxis()->CenterTitle();
      h_sigma_D_red[x_bin][Q2_bin]->GetYaxis()->SetRangeUser(h_sigma_D_red[x_bin][Q2_bin]->GetMinimum()/5,  5*h_sigma_D_red[x_bin][Q2_bin]->GetMaximum() );
      h_sigma_D_red[x_bin][Q2_bin]->SetMarkerStyle(25);
      h_sigma_D_red[x_bin][Q2_bin]->SetMarkerSize(1.5);
      h_sigma_D_red[x_bin][Q2_bin]->SetMarkerColor(1);
      h_sigma_D_red[x_bin][Q2_bin]->SetLineColor(1);
      h_sigma_D_red[x_bin][Q2_bin]->Draw("p e");

      h_sigma_D_red_RC[x_bin][Q2_bin]->Scale( reduce_factor * lumi_factor_MC );
      h_sigma_D_red_RC[x_bin][Q2_bin]->SetMarkerStyle(20);
      h_sigma_D_red_RC[x_bin][Q2_bin]->SetMarkerSize(1.5);
      h_sigma_D_red_RC[x_bin][Q2_bin]->SetMarkerColor(kRed);
      h_sigma_D_red_RC[x_bin][Q2_bin]->SetLineColor(kRed);
      if(plot_RC_flag > 0) h_sigma_D_red_RC[x_bin][Q2_bin]->Draw("p e same");

      h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->Scale( reduce_factor * lumi_factor_MC );
      h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->SetMarkerStyle(24);
      h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->SetMarkerSize(1.5);
      h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->SetMarkerColor(9);
      h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->SetLineColor(9);
      if(plot_ACC_flag == 1) h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->Draw("p e same");


      TLegend *legend = new TLegend(0.15, 0.78, 0.35, 0.89);
      legend->SetTextFont(43);
      legend->SetTextSize(30);
      legend->AddEntry(h_sigma_D_red[x_bin][Q2_bin], "True MC", "p e");
      if(plot_ACC_flag == 1) legend->AddEntry(h_sigma_D_red_MC_ACC[x_bin][Q2_bin], "MC, accept. only", "p e");
      if(plot_RC_flag > 0) legend->AddEntry(h_sigma_D_red_RC[x_bin][Q2_bin], "RC (uncorr.)", "p e");
      legend->SetBorderSize(0);
      legend->SetFillColorAlpha(0, 0.01);
      legend->Draw("same");


      TPaveText *text_kine = new TPaveText(0.36, 0.75, 0.66, 0.89, "NDC"); //main results;
      text_kine->SetTextFont(43);
      text_kine->SetTextSize(30);
      text_kine->AddText("0.99 < #alpha_{p} < 1.01");
      text_kine->AddText(Form("%.3f < x < %.3f", x_bins[x_bin], x_bins[x_bin+1]));
      text_kine->AddText(Form("%.0f < Q^{2} < %.0f GeV^{2}", Q2_bins[Q2_bin], Q2_bins[Q2_bin+1]));
      text_kine->SetBorderSize(0);
      text_kine->SetFillColorAlpha(0, 0.01);
      text_kine->Draw("same");

      text->Draw("same");


      sigma_red_D_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_D_Q2bin_%i_x_bin_%i.png", Q2_bin, x_bin));
      sigma_red_D_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_D_Q2bin_%i_x_bin_%i.pdf", Q2_bin, x_bin));

      //_______________________________________________________________________________


      TCanvas *sigma_red_p_can = new TCanvas(Form("sigma_red_p_can_%i_%i",x_bin, Q2_bin), Form("sigma_red_p_can_%i_%i",x_bin, Q2_bin), 1200,1000);
      sigma_red_p_can->cd();

      TH1F *default_hist = new TH1F(Form("default_hist_%i_%i", x_bin, Q2_bin), Form("default_hist_%i_%i", x_bin, Q2_bin), 100, -0.005, 0.01);
      default_hist->GetXaxis()->SetTitle("p^{2}_{T} (GeV/c)^{2}");
      default_hist->GetXaxis()->CenterTitle();
      default_hist->GetYaxis()->SetTitle("#sigma_{red,p}");
      default_hist->GetYaxis()->CenterTitle();
      default_hist->GetYaxis()->SetRangeUser(0 , 1.1);
      default_hist->SetLineColor(1);
      default_hist->Draw();

      h_sigma_p_red[x_bin][Q2_bin]->Scale( reduce_factor * lumi_factor_MC );
      h_sigma_p_red[x_bin][Q2_bin]->SetMarkerStyle(25);
      h_sigma_p_red[x_bin][Q2_bin]->SetMarkerSize(1.5);
      h_sigma_p_red[x_bin][Q2_bin]->SetMarkerColor(1);
      h_sigma_p_red[x_bin][Q2_bin]->SetLineColor(1);
      h_sigma_p_red[x_bin][Q2_bin]->Draw("p e same");

      //-------------

      h_sigma_p_red_RC[x_bin][Q2_bin]->Scale( reduce_factor * lumi_factor_MC );
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetMarkerStyle(20);
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetMarkerSize(1.5);
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetMarkerColor(kRed);
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetLineColor(kRed);
      if(plot_RC_flag == 1) h_sigma_p_red_RC[x_bin][Q2_bin]->Draw("p e same");

      TF1 *sigma_p_red_RC_fit = new TF1("sigma_p_red_RC_fit", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0, 0.01);
      sigma_p_red_RC_fit->SetParameters(0.1, 0.01, -0.1, 0.5);

      TF1 *sigma_p_red_RC_fit_2 = new TF1("sigma_p_red_RC_fit_2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 0.01);
      sigma_p_red_RC_fit_2->SetParameters(0.1, 0.01, -0.1, 0.5, 0.5);

      h_sigma_p_red_RC[x_bin][Q2_bin]->Fit(sigma_p_red_RC_fit, "0 R");
      h_sigma_p_red_RC[x_bin][Q2_bin]->Fit(sigma_p_red_RC_fit_2, "0 R");

      sigma_p_red_RC_fit->SetLineColor(1);
      sigma_p_red_RC_fit->SetLineStyle(8);
      if(plot_RC_flag == 1) sigma_p_red_RC_fit->Draw("same");

      sigma_p_red_RC_fit_2->SetLineColor(kRed);
      sigma_p_red_RC_fit_2->SetLineStyle(8);
      //if(plot_RC_flag == 1) sigma_p_red_RC_fit_2->Draw("same");

      //-------------

      TLine *zero_line = new TLine(0,0,0,1.1);
      zero_line->SetLineColor(1);
      zero_line->SetLineStyle(9);
      zero_line->Draw("same");

      TF1 *lin_fit = new TF1("lin_fit", "[0]+[1]*x", -0.005, 0.01);
      lin_fit->SetLineStyle(2);
      lin_fit->SetLineColor(1);
      lin_fit->SetParameters(0.2, 0.01);

      TFitResultPtr fit_ptr = h_sigma_p_red[x_bin][Q2_bin]->Fit(lin_fit, "0 R S");

      double sigma_p_at_pole[1] = { -aT2 };
      double sigma_p_at_pole_err[1];
      fit_ptr->GetConfidenceIntervals(1, 1, 1, sigma_p_at_pole, sigma_p_at_pole_err, 0.683, false); // 0.683 for 1-sigma

      lin_fit->Draw("same");


      TGraphErrors *sigma_free_n = new TGraphErrors(1);
      sigma_free_n->SetNameTitle("sigma_free_n", "sigma_free_n");
      sigma_free_n->SetPoint(0, -aT2, lin_fit->Eval(sigma_p_at_pole[0]));
      sigma_free_n->SetPointError(0, 0, sigma_p_at_pole_err[0]);
      sigma_free_n->SetMarkerStyle(45);
      sigma_free_n->SetMarkerColor(kRed);
      sigma_free_n->SetMarkerSize(2);
      sigma_free_n->SetLineColor(kRed);
      sigma_free_n->Draw("p e same");

      //-------

      sigma_free_p_vs_x[Q2_bin]->SetPoint(x_bin, x_bin_center, lin_fit->Eval(sigma_p_at_pole[0]));
      sigma_free_p_vs_x[Q2_bin]->SetPointError(x_bin, 0, sigma_p_at_pole_err[0]);

      //-------
      TLegend *legend_spectator = new TLegend(0.12, 0.7, 0.32, 0.89);
      legend_spectator->SetTextFont(43);
      legend_spectator->SetTextSize(30);
      legend_spectator->AddEntry(h_sigma_p_red[x_bin][Q2_bin], "True MC", "p e");
      legend_spectator->AddEntry(lin_fit, "Fit to MC", "l");
      if(plot_RC_flag == 1)
      {
        legend_spectator->AddEntry(h_sigma_p_red_RC[x_bin][Q2_bin], "RC (uncorr.)", "p e");
        legend_spectator->AddEntry(sigma_p_red_RC_fit, "Fit to RC (uncorr.)", "l");
      }
      legend_spectator->AddEntry(sigma_free_n, "Free p (MC)", "p e");
      legend_spectator->SetBorderSize(0);
      legend_spectator->SetFillColorAlpha(0, 0.01);
      legend_spectator->Draw("same");

      text_kine->Draw("same");

      text->Draw("same");

      sigma_red_p_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_p_Q2bin_%i_x_bin_%i.png", Q2_bin, x_bin));
      sigma_red_p_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_p_Q2bin_%i_x_bin_%i.pdf", Q2_bin, x_bin));

      //_______________________________________________________________________________

      //variation of MC points
      TGraphErrors *sigma_free_p_RC_var[n_variations];

      for(unsigned int var_i = 0; var_i < n_variations; var_i++)
      {
        TH1D *h_sigma_p_red_var = (TH1D*)h_sigma_p_red[x_bin][Q2_bin]->Clone(Form("h_sigma_p_red_xbin_%i_Q2bin_%i_var_%i", x_bin, Q2_bin, var_i));

        for(unsigned int i_point = 0; i_point < h_sigma_p_red_var->GetNbinsX(); i_point++)
        {
          float old_mean = h_sigma_p_red_var->GetBinContent(i_point+1);
          float stat_uncert = h_sigma_p_red_var->GetBinError(i_point+1);

          float new_mean = my_random->Gaus(old_mean, old_mean*sys_err_smear);

          h_sigma_p_red_var->SetBinContent(i_point+1, new_mean);
          h_sigma_p_red_var->SetBinError(i_point+1, stat_uncert);
        }

        TF1 *lin_fit_var = new TF1(Form("lin_fit_var_%i", var_i), "[0]+[1]*x", -0.005, 0.01);
        lin_fit_var->SetLineStyle(2);
        lin_fit_var->SetLineColor(1);
        lin_fit_var->SetParameters(0.2, 0.01);

        TFitResultPtr fit_ptr_var = h_sigma_p_red_var->Fit(lin_fit_var, "0 R S");


        TH1D *h_eff_nucleon_func_bins_var = (TH1D*)h_sigma_p_red_var->Clone(Form("h_eff_nucleon_func_bins_xbin_%i_Q2bin_%i_var_%i", x_bin, Q2_bin, var_i));

        //set eff. errors as bin-by-bin sys. error based on magnitude of the point shift
        for(unsigned int eff_bin = 1; eff_bin <= h_eff_nucleon_func_bins_var->GetNbinsX(); eff_bin++ )
        {
          float bin_center = h_eff_nucleon_func_bins_var->GetBinCenter(eff_bin);

          //The shift is calculated from uncorrected - we use the "total" efficiency correction, applied to uncorrected
          h_eff_nucleon_func_bins_var->SetBinContent(eff_bin, sigma_p_red_RC_fit->Eval(bin_center)/lin_fit_var->Eval(bin_center));
          h_eff_nucleon_func_bins_var->SetBinError(eff_bin, 0); //new version

        }


        TH1D *h_sigma_p_red_RC_var = (TH1D*)h_sigma_p_red_RC[x_bin][Q2_bin]->Clone(Form("h_sigma_p_red_RC_xbin_%i_Q2bin_%i_var_%i", x_bin, Q2_bin, var_i));
        h_sigma_p_red_RC_var->Divide(h_eff_nucleon_func_bins_var);


        TF1 *lin_fit_var_RC = new TF1(Form("lin_fit_var_RC_%i", var_i), "[0]+[1]*x", -0.005, 0.01);
        lin_fit_var_RC->SetLineStyle(2);
        lin_fit_var_RC->SetLineColor(1);
        lin_fit_var_RC->SetParameters(0.2, 0.01);

        TFitResultPtr fit_ptr_var_RC = h_sigma_p_red_RC_var->Fit(lin_fit_var_RC, "0 R S");

        double sigma_p_at_pole_var[1] = { -aT2 };
        double sigma_p_at_pole_err_var[1];
        fit_ptr_var_RC->GetConfidenceIntervals(1, 1, 1, sigma_p_at_pole_var, sigma_p_at_pole_err_var, 0.683, false); // 0.683 for 1-sigma

        sigma_free_p_RC_var[var_i] = new TGraphErrors(1);
        sigma_free_p_RC_var[var_i]->SetNameTitle(Form("sigma_free_p_var_%i", var_i), Form("sigma_free_p_var_%i", var_i));
        sigma_free_p_RC_var[var_i]->SetPoint(0, -aT2, lin_fit_var_RC->Eval(sigma_p_at_pole_var[0]));
        sigma_free_p_RC_var[var_i]->SetPointError(0, 0, sigma_p_at_pole_err_var[0]);

      }

      //_______________________________________________________________________________

      TCanvas *eff_nucleon_func_can = new TCanvas(Form("eff_nucleon_func_can%i_%i",x_bin, Q2_bin), Form("eff_nucleon_func_can%i_%i",x_bin, Q2_bin), 1200,1000);
      eff_nucleon_func_can->cd();

	    h_eff_nucleon_func_bins[x_bin][Q2_bin] = (TH1F*)h_sigma_p_red_RC[x_bin][Q2_bin]->Clone(Form("h_eff_nucleon_func_bins_xbin_%i_Q2bin_%i", x_bin, Q2_bin));
	    h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetXaxis()->SetTitle("p^{2}_{T} (GeV/c)^{2}");
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetXaxis()->CenterTitle();
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetYaxis()->SetTitle("RC fit/MC fit");
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetYaxis()->CenterTitle();
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->SetMarkerStyle(25);
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->SetMarkerSize(1.5);
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->SetMarkerColor(1);
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->SetLineColor(1);

      //set eff. errors as bin-by-bin sys. error based on magnitude of the point shift
      for(unsigned int eff_bin = 1; eff_bin <= h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetNbinsX(); eff_bin++ )
      {
        float bin_center = h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetBinCenter(eff_bin);

        //The shift is calculated from uncorrected - we use the "total" efficiency correction, applied to uncorrected
        h_eff_nucleon_func_bins[x_bin][Q2_bin]->SetBinContent(eff_bin, sigma_p_red_RC_fit->Eval(bin_center)/lin_fit->Eval(bin_center));
        h_eff_nucleon_func_bins[x_bin][Q2_bin]->SetBinError(eff_bin, 0); //new version

      }

      h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetYaxis()->SetRangeUser(0, h_eff_nucleon_func_bins[x_bin][Q2_bin]->GetMaximum()*3);
      h_eff_nucleon_func_bins[x_bin][Q2_bin]->Draw("p e");

      text_kine->Draw("same");

      text->Draw("same");

      eff_nucleon_func_can->SaveAs(Form("./plots/eP/eff/eff_nucleon_func_Q2bin_%i_x_bin_%i.png", Q2_bin, x_bin));
      eff_nucleon_func_can->SaveAs(Form("./plots/eP/eff/eff_nucleon_func_Q2bin_%i_x_bin_%i.pdf", Q2_bin, x_bin));


      //--------------------------------------------------------------------------------------------------------------------------

      TCanvas *sigma_red_p_corr_can = new TCanvas(Form("sigma_red_p_corr_can_%i_%i",x_bin, Q2_bin), Form("sigma_red_p_corr_can_%i_%i",x_bin, Q2_bin), 1200,1000);
      sigma_red_p_corr_can->cd();

      //gPad->SetLogy();

      default_hist->Draw();

      h_sigma_p_red_RC[x_bin][Q2_bin]->Divide(h_eff_nucleon_func_bins[x_bin][Q2_bin]);
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetMarkerStyle(20);
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetMarkerSize(1.5);
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetMarkerColor(kRed);
      h_sigma_p_red_RC[x_bin][Q2_bin]->SetLineColor(kRed);
      //h_sigma_p_red_RC[x_bin][Q2_bin]->Draw("p e same");


      TH1D *h_sigma_p_red_RC_stat_sys = (TH1D*)h_sigma_p_red_RC[x_bin][Q2_bin]->Clone(Form("h_sigma_p_red_RC_stat_sys_%i_%i", x_bin, Q2_bin));

      for(unsigned int i_bin = 1; i_bin <= h_sigma_p_red_RC_stat_sys->GetNbinsX(); i_bin++)
      {
        if( h_sigma_p_red_RC_stat_sys->GetBinContent(i_bin) == 0 ) continue;

        float stat_err = h_sigma_p_red_RC_stat_sys->GetBinError(i_bin)/h_sigma_p_red_RC_stat_sys->GetBinContent(i_bin);

        float tot_err = sqrt(stat_err*stat_err + sys_err_tot_corr*sys_err_tot_corr);

        h_sigma_p_red_RC_stat_sys->SetBinError(i_bin, tot_err*h_sigma_p_red_RC_stat_sys->GetBinContent(i_bin));

      }

      h_sigma_p_red_RC_stat_sys->Draw("p e same");


      //--------------------------

      zero_line->Draw("same");

      TF1 *lin_fit_RC_corr = new TF1("lin_fit_RC_corr", "[0]+[1]*x", -0.005, 0.01);
      lin_fit_RC_corr->SetLineStyle(2);
      lin_fit_RC_corr->SetLineColor(1);
      lin_fit_RC_corr->SetParameters(0.2, 0.01);

      TFitResultPtr fit_ptr_RC = h_sigma_p_red_RC[x_bin][Q2_bin]->Fit(lin_fit_RC_corr, "0 R S");

      double sigma_p_at_pole_RC_corr[1] = { -aT2 };
      double sigma_p_at_pole_err_RC_corr[1];
      fit_ptr_RC->GetConfidenceIntervals(1, 1, 1, sigma_p_at_pole_RC_corr, sigma_p_at_pole_err_RC_corr, 0.683, false); // 0.683 for 1-sigma

      //------

      TF1 *lin_fit_RC_corr_tot_err = new TF1("lin_fit_RC_corr_tot_err", "[0]+[1]*x", -0.005, 0.01);
      lin_fit_RC_corr_tot_err->SetLineStyle(2);
      lin_fit_RC_corr_tot_err->SetLineColor(1);
      lin_fit_RC_corr_tot_err->SetParameters(0.2, 0.01);

      TFitResultPtr fit_ptr_RC_tot_err = h_sigma_p_red_RC_stat_sys->Fit(lin_fit_RC_corr_tot_err, "0 R S");


      double sigma_p_at_pole_RC_corr_tot_err[1] = { -aT2 };
      double sigma_p_at_pole_err_RC_corr_tot_err[1];
      fit_ptr_RC_tot_err->GetConfidenceIntervals(1, 1, 1, sigma_p_at_pole_RC_corr_tot_err, sigma_p_at_pole_err_RC_corr_tot_err, 0.683, false); // 0.683 for 1-sigma

      //-------

      lin_fit_RC_corr_tot_err->Draw("same");

      TGraphErrors *sigma_free_p_RC_corr = new TGraphErrors(1);
      sigma_free_p_RC_corr->SetNameTitle("sigma_free_p_RC_corr", "sigma_free_p_RC_corr");
      sigma_free_p_RC_corr->SetPoint(0, -aT2, lin_fit_RC_corr_tot_err->Eval(sigma_p_at_pole_RC_corr[0]));
      sigma_free_p_RC_corr->SetPointError(0, 0, sigma_p_at_pole_err_RC_corr_tot_err[0]);
      sigma_free_p_RC_corr->SetMarkerStyle(47);
      sigma_free_p_RC_corr->SetMarkerColor(9);
      sigma_free_p_RC_corr->SetMarkerSize(2);
      sigma_free_p_RC_corr->SetLineColor(9);
      sigma_free_p_RC_corr->Draw("p e same"); //plot last to see point over sys. err.

      //-------------------------------------------------------------------------------

      sigma_free_p_vs_x_RC_corr[Q2_bin]->SetPoint(x_bin, x_bin_center*1.1, lin_fit_RC_corr_tot_err->Eval(sigma_p_at_pole_RC_corr[0]));
      sigma_free_p_vs_x_RC_corr[Q2_bin]->SetPointError(x_bin, 0, sigma_p_at_pole_err_RC_corr_tot_err[0]);

      //-------

      //outputFile<<x_bin_center<<" "<<lin_fit_RC_corr_tot_err->Eval(sigma_p_at_pole_RC_corr[0])<<" "<<sigma_p_at_pole_err_RC_corr[0]<<" "<<sigma_p_at_pole_err_RC_corr_tot_err[0]<<" "<<y_mean<<endl;
      outputFile<<x_bin_center<<" "<<lin_fit_RC_corr_tot_err->Eval(sigma_p_at_pole_RC_corr[0])<<" "<<sigma_p_at_pole_err_RC_corr[0]<<" "<<sigma_p_at_pole_err_RC_corr_tot_err[0]<<endl;

      //-------

      TLegend *legend_spectator_RC_corr = new TLegend(0.15, 0.72, 0.35, 0.89);
      legend_spectator_RC_corr->SetTextFont(43);
      legend_spectator_RC_corr->SetTextSize(30);
      legend_spectator_RC_corr->AddEntry(h_sigma_p_red_RC[x_bin][Q2_bin], "Corrected RC", "p e");
      legend_spectator_RC_corr->AddEntry(lin_fit_RC_corr, "Fit to RC", "l");
      legend_spectator_RC_corr->AddEntry(sigma_free_p_RC_corr, "Free p (RC)", "p e");
      legend_spectator_RC_corr->SetBorderSize(0);
      legend_spectator_RC_corr->SetFillColorAlpha(0, 0.01);
      legend_spectator_RC_corr->Draw("same");

      text_kine->Draw("same");

      text->Draw("same");

      sigma_red_p_corr_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_p_corr_Q2bin_%i_x_bin_%i.png", Q2_bin, x_bin));
      sigma_red_p_corr_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_p_corr_Q2bin_%i_x_bin_%i.pdf", Q2_bin, x_bin));

    }

    outputFile<<endl;
    outputFile_D<<endl;

    //-----------------------------

    TCanvas *sigma_red_D_vs_x_can = new TCanvas(Form("sigma_red_D_vs_x_can_%i", Q2_bin), Form("sigma_red_D_vs_x_can_%i", Q2_bin), 1200,1000);
    sigma_red_D_vs_x_can->cd();

    //gPad->SetLogy();
    gPad->SetLogx();

    h_sigma_D_red_vs_x[Q2_bin]->GetXaxis()->SetTitle("x");
    h_sigma_D_red_vs_x[Q2_bin]->GetXaxis()->CenterTitle();
    h_sigma_D_red_vs_x[Q2_bin]->GetYaxis()->SetTitle("#bar{#sigma}_{red,d}");
    h_sigma_D_red_vs_x[Q2_bin]->GetYaxis()->CenterTitle();
    h_sigma_D_red_vs_x[Q2_bin]->GetYaxis()->SetRangeUser(0,  1.2 );
    h_sigma_D_red_vs_x[Q2_bin]->SetMarkerStyle(25);
    h_sigma_D_red_vs_x[Q2_bin]->SetMarkerSize(1.5);
    h_sigma_D_red_vs_x[Q2_bin]->SetMarkerColor(1);
    h_sigma_D_red_vs_x[Q2_bin]->SetLineColor(1);
    h_sigma_D_red_vs_x[Q2_bin]->Draw("p e");

    h_sigma_D_red_vs_x_RC[Q2_bin]->SetMarkerStyle(20);
    h_sigma_D_red_vs_x_RC[Q2_bin]->SetMarkerSize(1.5);
    h_sigma_D_red_vs_x_RC[Q2_bin]->SetMarkerColor(kRed);
    h_sigma_D_red_vs_x_RC[Q2_bin]->SetLineColor(kRed);
    if(plot_RC_flag > 0) h_sigma_D_red_vs_x_RC[Q2_bin]->Draw("p e same");

    h_sigma_D_red_vs_x_MC_ACC[Q2_bin]->SetMarkerStyle(24);
    h_sigma_D_red_vs_x_MC_ACC[Q2_bin]->SetMarkerSize(1.5);
    h_sigma_D_red_vs_x_MC_ACC[Q2_bin]->SetMarkerColor(9);
    h_sigma_D_red_vs_x_MC_ACC[Q2_bin]->SetLineColor(9);

    TLegend *legend_vs_x = new TLegend(0.15, 0.78, 0.35, 0.89);
    legend_vs_x->SetTextFont(43);
    legend_vs_x->SetTextSize(30);
    legend_vs_x->AddEntry(h_sigma_D_red_vs_x[Q2_bin], "True MC", "p e");
    if(plot_ACC_flag == 1) legend_vs_x->AddEntry(h_sigma_D_red_vs_x_MC_ACC[Q2_bin], "MC, accept. only", "p e");
    if(plot_RC_flag > 0) legend_vs_x->AddEntry(h_sigma_D_red_vs_x_RC[Q2_bin], "RC (uncorr.)", "p e");
    legend_vs_x->SetBorderSize(0);
    legend_vs_x->SetFillColorAlpha(0, 0.01);
    legend_vs_x->Draw("same");


    TPaveText *text_kine_vs_x_D = new TPaveText(0.36, 0.75, 0.66, 0.89, "NDC"); //main results;
    text_kine_vs_x_D->SetTextFont(43);
    text_kine_vs_x_D->SetTextSize(30);
    text_kine_vs_x_D->AddText(Form("%.0f < Q^{2} < %.0f GeV^{2}", Q2_bins[Q2_bin], Q2_bins[Q2_bin+1]));
    text_kine_vs_x_D->SetBorderSize(0);
    text_kine_vs_x_D->SetFillColorAlpha(0, 0.01);
    text_kine_vs_x_D->Draw("same");

    TPaveText *text_vs_x = new TPaveText(0.65, 0.7, 0.95, 0.89, "NDC"); //main results;
    text_vs_x->SetTextFont(43);
    text_vs_x->SetTextSize(30);
    text_vs_x->AddText("ePIC simulation");
    text_vs_x->AddText("eD 10x130 GeV^{2}");
    text_vs_x->AddText("#gamma* + d #rightarrow X");
    text_vs_x->AddText("Inclusive");
    text_vs_x->SetBorderSize(0);
    text_vs_x->SetFillColorAlpha(0, 0.01);
    text_vs_x->Draw("same");


    sigma_red_D_vs_x_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_D_vs_x_Q2bin_%i.png", Q2_bin));
    sigma_red_D_vs_x_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_D_vs_X_Q2bin_%i.pdf", Q2_bin));



    //-----------------------------

    TCanvas *sigma_red_p_vs_x_can = new TCanvas(Form("sigma_red_p_vs_x_can_%i", Q2_bin), Form("sigma_red_p_vs_x_can_%i", Q2_bin), 1200,1000);
    sigma_red_p_vs_x_can->cd();

    gPad->SetLogx();

    TH1F *default_hist_vs_x = new TH1F(Form("default_hist_vs_x_%i", Q2_bin), Form("default_hist_vs_x_%i", Q2_bin), 10000, 0, 1);
    default_hist_vs_x->GetXaxis()->SetTitle("x_{bj}");
    default_hist_vs_x->GetXaxis()->CenterTitle();
    default_hist_vs_x->GetXaxis()->SetRangeUser(1e-3 , 0.8);
    default_hist_vs_x->GetYaxis()->SetTitle("#sigma_{red,p}");
    default_hist_vs_x->GetYaxis()->CenterTitle();
    default_hist_vs_x->GetYaxis()->SetRangeUser(0 , 1.0);

    default_hist_vs_x->SetLineColor(1);
    default_hist_vs_x->Draw();


    if(Q2_bin == 0) sigma_free_p_vs_x[Q2_bin]->RemovePoint(9);
    if(Q2_bin == 1 ) sigma_free_p_vs_x[Q2_bin]->RemovePoint(10);
    sigma_free_p_vs_x[Q2_bin]->RemovePoint(Q2_bin);
    if(Q2_bin == 3) sigma_free_p_vs_x[Q2_bin]->RemovePoint(3);
    sigma_free_p_vs_x[Q2_bin]->SetMarkerStyle(25);
    sigma_free_p_vs_x[Q2_bin]->SetMarkerSize(1.5);
    sigma_free_p_vs_x[Q2_bin]->SetMarkerColor(1);
    sigma_free_p_vs_x[Q2_bin]->SetLineColor(1);
    sigma_free_p_vs_x[Q2_bin]->Draw("p e same");

    if(Q2_bin == 0) sigma_free_p_vs_x_RC_corr[Q2_bin]->RemovePoint(9);
    if(Q2_bin == 1 ) sigma_free_p_vs_x_RC_corr[Q2_bin]->RemovePoint(10);
    sigma_free_p_vs_x_RC_corr[Q2_bin]->RemovePoint(Q2_bin);
    if(Q2_bin == 3) sigma_free_p_vs_x_RC_corr[Q2_bin]->RemovePoint(3);
    sigma_free_p_vs_x_RC_corr[Q2_bin]->SetMarkerStyle(20);
    sigma_free_p_vs_x_RC_corr[Q2_bin]->SetMarkerSize(1.5);
    sigma_free_p_vs_x_RC_corr[Q2_bin]->SetMarkerColor(kRed);
    sigma_free_p_vs_x_RC_corr[Q2_bin]->SetLineColor(kRed);
    sigma_free_p_vs_x_RC_corr[Q2_bin]->Draw("p e same");

    TPaveText *text_kine_vs_x_alpha = new TPaveText(0.36, 0.75, 0.66, 0.89, "NDC"); //main results;
    text_kine_vs_x_alpha->SetTextFont(43);
    text_kine_vs_x_alpha->SetTextSize(30);
    text_kine_vs_x_alpha->AddText("0.99 < #alpha_{p} < 1.01");
    text_kine_vs_x_alpha->AddText(Form("%.0f < Q^{2} < %.0f GeV^{2}", Q2_bins[Q2_bin], Q2_bins[Q2_bin+1]));
    text_kine_vs_x_alpha->SetBorderSize(0);
    text_kine_vs_x_alpha->SetFillColorAlpha(0, 0.01);
    text_kine_vs_x_alpha->Draw("same");

    text->Draw("same");

    TLegend *legend_sigma = new TLegend(0.15, 0.77, 0.35, 0.89);
    legend_sigma->SetTextFont(43);
    legend_sigma->SetTextSize(30);
    legend_sigma->AddEntry(sigma_free_p_vs_x[Q2_bin], "True MC");
    legend_sigma->AddEntry(sigma_free_p_vs_x_RC_corr[Q2_bin], "Corrected RC");
    legend_sigma->SetBorderSize(0);
    legend_sigma->SetFillColorAlpha(0, 0.01);
    legend_sigma->Draw("same");

    sigma_red_p_vs_x_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_p_vs_x_Q2bin_%i.png", Q2_bin));
    sigma_red_p_vs_x_can->SaveAs(Form("./plots/eP/cross_section/sigma_red_p_vs_x_Q2bin_%i.pdf", Q2_bin));

    //----------------------

    //same plot, but in 4-panel canvas
    sigma_red_p_vs_x_one_can->cd(Q2_bin+1);

    gPad->SetLogx();

    default_hist_vs_x->Draw();

    sigma_free_p_vs_x[Q2_bin]->Draw("p e same");
    sigma_free_p_vs_x_RC_corr[Q2_bin]->Draw("p e same");

    TPaveText *text_kine_vs_x = new TPaveText(0.36, 0.75, 0.66, 0.89, "NDC"); //main results;
    text_kine_vs_x->SetTextFont(43);
    text_kine_vs_x->SetTextSize(30);
    text_kine_vs_x->AddText("0.99 < #alpha_{p} < 1.01");
    text_kine_vs_x->AddText(Form("%.0f < Q^{2} < %.0f GeV^{2}", Q2_bins[Q2_bin], Q2_bins[Q2_bin+1]));
    text_kine_vs_x->SetBorderSize(0);
    text_kine_vs_x->SetFillColorAlpha(0, 0.01);
    text_kine_vs_x->Draw("same");

    if(Q2_bin == 1)
    {
      legend_sigma->Draw("same");

      text->Draw("same");
    }

    //----------------------------------------

    if(Q2_bin == 0)
    {
      pad_0->cd();

      gPad->SetRightMargin(0);
    }
    if(Q2_bin == 1)
    {
      pad_1->cd();

      gPad->SetRightMargin(0);
      gPad->SetLeftMargin(0);
    }
    if(Q2_bin == 2)
    {
      pad_2->cd();

      gPad->SetRightMargin(0);
      gPad->SetLeftMargin(0);
    }
    if(Q2_bin == 3)
    {
      pad_3->cd();

      gPad->SetLeftMargin(0);
    }

    gPad->SetLogx();


    default_hist_vs_x->Draw();

    sigma_free_p_vs_x[Q2_bin]->Draw("p e same");
    sigma_free_p_vs_x_RC_corr[Q2_bin]->Draw("p e same");

    text_kine_vs_x->Draw("same");

    if(Q2_bin == 1)
    {
      legend_sigma->Draw("same");

      text->Draw("same");
    }

    //--------------------------------------------------------------

  }

  //-------------------------------------

  sigma_red_p_vs_x_one_can->SaveAs("./plots/eP/cross_section/sigma_red_p_vs_x.png");
  sigma_red_p_vs_x_one_can->SaveAs("./plots/eP/cross_section/sigma_red_p_vs_x.pdf");

  sigma_red_p_vs_x_one_can_2->SaveAs("./plots/eP/cross_section/sigma_red_p_vs_x_2.png");
  sigma_red_p_vs_x_one_can_2->SaveAs("./plots/eP/cross_section/sigma_red_p_vs_x_2.pdf");


  //------------------------------------------------------------


  inFile->Close();
  outputFile.close();
  outputFile_D.close();

  return;
}
