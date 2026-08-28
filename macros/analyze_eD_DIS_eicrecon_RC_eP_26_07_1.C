//-------------------------
//
// Macro to run TDIS analysis on EICRecon output
//
// Author: Jan Vanek
// Based on macro by: Alex Jentsch
//
// Originally created on: 11/16/2025
//
//------------------------

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TMath.h>

#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include "TFitResultPtr.h"

#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include"TDatabasePDG.h"
#include"TParticlePDG.h"


using namespace std;


//ECAL cut range
const float E_over_p_cut_min = 0.8;
const float E_over_p_cut_max = 1.2;

//This is the same code as analyze_eD_DIS_eicrecon_RC_eP.C but is slightly modified to reflect changes to the ZDC cluster collection inside EICRecon.
void analyze_eD_DIS_eicrecon_RC_eP_26_07_1(TString fileList = "./fileLists/26.03.1/fileList_eD_10x130_DIS_eP_short.list", TString outputName = "./output/test_output_RC_eP")
{

	cout << "Input FileList: " << fileList << endl;
	TString fileType_ROOT = ".root";
	TString fileType_PDF = ".pdf";

	TString outputFileName = outputName + fileType_ROOT;
	string fileName;

	TTree * rootTree;
	cout << "Output file: " << outputFileName << endl;


	ifstream fileListStream;
	fileListStream.open(fileList);
	if(!fileListStream) { cout << "NO_LIST_FILE " << fileList << endl; return;}

	//--------------------------------------------------------------------------
	// important variables

	TDatabasePDG *myPDGdatabase = new TDatabasePDG();

  TParticlePDG *e_PDG = myPDGdatabase->GetParticle(11);
  TParticlePDG *p_PDG = myPDGdatabase->GetParticle(2212);
  TParticlePDG *n_PDG = myPDGdatabase->GetParticle(2112);

  float mDeuteron = 1.87561;  //GeV
  float mProton   = .93827;   //GeV
  float mNeutron  = .93957;   //GeV
	float mElectron = 0.000511; //GeV
	float mAvg      = 0.5*(mProton+mNeutron);

	float capGammaSquared = 0.007885; //0.0078;
  float deutBindEner = 0.0022;
  float aSquare = 0.002121; //deutBindEner*mAvg;

	float phaseFactor = 1.0/(2.0*TMath::Power(TMath::TwoPi(), 3));

	const float alpha_em = 1./137.;

	float two_times_2pi_cubed = 2.*(2.*TMath::Pi())*(2.*TMath::Pi())*(2.*TMath::Pi());

	//--------------------------------------------------------------------------
	//bins

  const int nx_bins = 11;
  const float x_bins[nx_bins+1] = {0.001, 0.002, 0.004, 0.007, 0.01, 0.02, 0.04, 0.07, 0.09, 0.2, 0.3, 0.5};

  const int nQ2_bins = 4;
  const float Q2_bins[nQ2_bins+1] = {1., 10., 20., 40., 100.};


  float alpha_p_min = 0.99;
  float alpha_p_max = 1.01;
  float Delta_apha_p = alpha_p_max - alpha_p_min;


	const int n_theta_bins = 4;
	float theta_bins[n_theta_bins+1] = {0.,1.,2.,3.,4.};

	const int n_ZDC_z_bins = 4;
	float ZDC_z_bins[n_ZDC_z_bins+1] = {35800,36000,36200,36400,37400};


	const int n_pT2_bins = 5;
	float pT2_bins[n_pT2_bins+1] = {0.,0.002,0.004,0.006,0.008,0.01};

	//--------------------------------------------------------------------------

	//histograms
	//event stats
	TH1F *h_nEvents = new TH1F("h_nEvents", "h_nEvents", 3, 0, 3);


  //scattered electron
  TH2F *h_x_vs_Q2 = new TH2F("h_x_vs_Q2", "h_x_vs_Q2", 200, 0, 2, 100, 0, 100);
  TH2F *h_x_vs_Q2_RC = new TH2F("h_x_vs_Q2_RC", "h_x_vs_Q2_RC", 200, 0, 2, 100, 0, 100);

  TH2F *h_x_nucleon_vs_Q2 = new TH2F("h_x_nucleon_vs_Q2", "h_x_nucleon_vs_Q2", 200, 0, 2, 100, 0, 100);
  TH2F *h_x_nucleon_vs_Q2_RC = new TH2F("h_x_nucleon_vs_Q2_RC", "h_x_nucleon_vs_Q2_RC", 200, 0, 2, 100, 0, 100);


  //MC vs. RC comparison
  TH2F *h_x_MC_vs_x_RC = new TH2F("h_x_MC_vs_x_RC", "h_x_MC_vs_x_RC", 200, 0, 2, 200, 0, 2);
  TH2F *h_x_nucleon_MC_vs_x_nucleon_RC = new TH2F("h_x_nucleon_MC_vs_x_nucleon_RC", "h_x_nucleon_MC_vs_x_nucleon_RC", 200, 0, 2, 200, 0, 2);

  TH2F *h_Q2_MC_vs_Q2_RC = new TH2F("h_Q2_MC_vs_Q2_RC", "h_Q2_MC_vs_Q2_RC", 100, 0, 100, 100, 0, 100);

  //---------
  //after "cuts" (especially for alpha_p)
  TH2F *h_x_vs_Q2_cuts = new TH2F("h_x_vs_Q2_cuts", "h_x_vs_Q2_cuts", 200, 0, 2, 100, 0, 100);
  TH2F *h_x_vs_Q2_RC_cuts = new TH2F("h_x_vs_Q2_RC_cuts", "h_x_vs_Q2_RC_cuts", 200, 0, 2, 100, 0, 100);

  TH2F *h_x_nucleon_vs_Q2_cuts = new TH2F("h_x_nucleon_vs_Q2_cuts", "h_x_nucleon_vs_Q2_cuts", 200, 0, 2, 100, 0, 100);
  TH2F *h_x_nucleon_vs_Q2_RC_cuts = new TH2F("h_x_nucleon_vs_Q2_RC_cuts", "h_x_nucleon_vs_Q2_RC_cuts", 200, 0, 2, 100, 0, 100);

  //MC vs. RC comparison
  TH2F *h_x_MC_vs_x_RC_cuts = new TH2F("h_x_MC_vs_x_RC_cuts", "h_x_MC_vs_x_RC_cuts", 200, 0, 2, 200, 0, 2);
  TH2F *h_x_nucleon_MC_vs_x_nucleon_RC_cuts = new TH2F("h_x_nucleon_MC_vs_x_nucleon_RC_cuts", "h_x_nucleon_MC_vs_x_nucleon_RC_cuts", 200, 0, 2, 200, 0, 2);

  TH2F *h_Q2_MC_vs_Q2_RC_cuts = new TH2F("h_Q2_MC_vs_Q2_RC_cuts", "h_Q2_MC_vs_Q2_RC_cuts", 100, 0, 100, 100, 0, 100);

  //---------
  //after "cuts" (especially for alpha_p)
  TH2F *h_x_vs_Q2_cuts_pT2_cut = new TH2F("h_x_vs_Q2_cuts_pT2_cut", "h_x_vs_Q2_cuts_pT2_cut", 200, 0, 2, 100, 0, 100);
  TH2F *h_x_vs_Q2_RC_cuts_pT2_cut = new TH2F("h_x_vs_Q2_RC_cuts_pT2_cut", "h_x_vs_Q2_RC_cuts_pT2_cut", 200, 0, 2, 100, 0, 100);

  TH2F *h_x_nucleon_vs_Q2_cuts_pT2_cut = new TH2F("h_x_nucleon_vs_Q2_cuts_pT2_cut", "h_x_nucleon_vs_Q2_cuts_pT2_cut", 200, 0, 2, 100, 0, 100);
  TH2F *h_x_nucleon_vs_Q2_RC_cuts_pT2_cut = new TH2F("h_x_nucleon_vs_Q2_RC_cuts_pT2_cut", "h_x_nucleon_vs_Q2_RC_cuts_pT2_cut", 200, 0, 2, 100, 0, 100);

  //MC vs. RC comparison
  TH2F *h_x_MC_vs_x_RC_cuts_pT2_cut = new TH2F("h_x_MC_vs_x_RC_cuts_pT2_cut", "h_x_MC_vs_x_RC_cuts_pT2_cut", 200, 0, 2, 200, 0, 2);
  TH2F *h_x_nucleon_MC_vs_x_nucleon_RC_cuts_pT2_cut = new TH2F("h_x_nucleon_MC_vs_x_nucleon_RC_cuts_pT2_cut", "h_x_nucleon_MC_vs_x_nucleon_RC_cuts_pT2_cut", 200, 0, 2, 200, 0, 2);

  TH2F *h_Q2_MC_vs_Q2_RC_cuts_pT2_cut = new TH2F("h_Q2_MC_vs_Q2_RC_cuts_pT2_cut", "h_Q2_MC_vs_Q2_RC_cuts_pT2_cut", 100, 0, 100, 100, 0, 100);

  TH2F *h_y_vs_x_vs_Q2_cuts_pT2_cut = new TH2F("h_y_vs_x_vs_Q2_cuts_pT2_cut", "h_y_vs_x_vs_Q2_cuts_pT2_cut", nx_bins, x_bins, nQ2_bins, Q2_bins);
  TH2F *h_y_vs_x_vs_Q2_cuts_pT2_cut_base = new TH2F("h_y_vs_x_vs_Q2_cuts_pT2_cut_base", "h_y_vs_x_vs_Q2_cuts_pT2_cut_base", nx_bins, x_bins, nQ2_bins, Q2_bins);
  
  TH2F *h_y_vs_x_vs_Q2_RC_cuts_pT2_cut = new TH2F("h_y_vs_x_vs_Q2_RC_cuts_pT2_cut", "h_y_vs_x_vs_Q2_RC_cuts_pT2_cut", nx_bins, x_bins, nQ2_bins, Q2_bins);
  TH2F *h_y_vs_x_vs_Q2_RC_cuts_pT2_cut_base = new TH2F("h_y_vs_x_vs_Q2_RC_cuts_pT2_cut_base", "h_y_vs_x_vs_Q2_RC_cuts_pT2_cut_base", nx_bins, x_bins, nQ2_bins, Q2_bins);

  //--------------------------

  TH1F *h_scat_e_pT = new TH1F("h_scat_e_pT", "h_scat_e_pT", 100, 0, 10);
  TH1F *h_scat_e_eta = new TH1F("h_scat_e_eta", "h_scat_e_eta", 100, -5, 5);
  TH1F *h_scat_e_phi = new TH1F("h_scat_e_phi", "h_scat_e_phi", 360, -TMath::Pi(), TMath::Pi());

  TH1F *h_track_eta_RC = new TH1F("h_track_eta_RC", "h_track_eta_RC", 100, -5, 5);
  TH1F *h_nECAL_cluster_eta_RC = new TH1F("h_nECAL_cluster_eta_RC", "h_nECAL_cluster_eta_RC", 100, -5, 5);
  TH1F *h_bECAL_cluster_eta_RC = new TH1F("h_bECAL_cluster_eta_RC", "h_bECAL_cluster_eta_RC", 100, -5, 5);

  TH1F *h_scat_e_pT_RC = new TH1F("h_scat_e_pT_RC", "h_scat_e_pT_RC", 100, 0, 10);
  TH1F *h_scat_e_eta_RC = new TH1F("h_scat_e_eta_RC", "h_scat_e_eta_RC", 100, -5, 5);
  TH1F *h_scat_e_phi_RC = new TH1F("h_scat_e_phi_RC", "h_scat_e_phi_RC", 360, -TMath::Pi(), TMath::Pi());

  //MC vs. RC comparison
  TH2F *h_scat_e_pT_MC_vs_e_pT_RC = new TH2F("h_scat_e_pT_MC_vs_e_pT_RC", "h_scat_e_pT_MC_vs_e_pT_RC", 100, 0, 10, 100, 0, 10);
  TH2F *h_scat_e_eta_MC_vs_e_eta_RC = new TH2F("h_scat_e_eta_MC_vs_e_eta_RC", "h_scat_e_eta_MC_vs_e_eta_RC", 100, -5, 5, 100, -5, 5);
  TH2F *h_scat_e_phi_MC_vs_e_phi_RC = new TH2F("h_scat_e_phi_MC_vs_e_phi_RC", "h_scat_e_phi_MC_vs_e_phi_RC", 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());

  //----
  //after "cuts" (especially for alpha_p)
  TH1F *h_scat_e_pT_cuts = new TH1F("h_scat_e_pT_cuts", "h_scat_e_pT_cuts", 100, 0, 10);
  TH1F *h_scat_e_eta_cuts = new TH1F("h_scat_e_eta_cuts", "h_scat_e_eta_cuts", 100, -5, 5);
  TH1F *h_scat_e_phi_cuts = new TH1F("h_scat_e_phi_cuts", "h_scat_e_phi_cuts", 360, -TMath::Pi(), TMath::Pi());

  TH1F *h_scat_e_pT_RC_cuts = new TH1F("h_scat_e_pT_RC_cuts", "h_scat_e_pT_RC_cuts", 100, 0, 10);
  TH1F *h_scat_e_eta_RC_cuts = new TH1F("h_scat_e_eta_RC_cuts", "h_scat_e_eta_RC_cuts", 100, -5, 5);
  TH1F *h_scat_e_phi_RC_cuts = new TH1F("h_scat_e_phi_RC_cuts", "h_scat_e_phi_RC_cuts", 360, -TMath::Pi(), TMath::Pi());

  //MC vs. RC comparison
  TH2F *h_scat_e_pT_MC_vs_e_pT_RC_cuts = new TH2F("h_scat_e_pT_MC_vs_e_pT_RC_cuts", "h_scat_e_pT_MC_vs_e_pT_RC_cuts", 100, 0, 10, 100, 0, 10);
  TH2F *h_scat_e_eta_MC_vs_e_eta_RC_cuts = new TH2F("h_scat_e_eta_MC_vs_e_eta_RC_cuts", "h_scat_e_eta_MC_vs_e_eta_RC_cuts", 100, -5, 5, 100, -5, 5);
  TH2F *h_scat_e_phi_MC_vs_e_phi_RC_cuts = new TH2F("h_scat_e_phi_MC_vs_e_phi_RC_cuts", "h_scat_e_phi_MC_vs_e_phi_RC_cuts", 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());

  //----
  //after "cuts" (especially for alpha_p)
  TH1F *h_scat_e_pT_cuts_pT2_cut = new TH1F("h_scat_e_pT_cuts_pT2_cut", "h_scat_e_pT_cuts_pT2_cut", 100, 0, 10);
  TH1F *h_scat_e_eta_cuts_pT2_cut = new TH1F("h_scat_e_eta_cuts_pT2_cut", "h_scat_e_eta_cuts_pT2_cut", 100, -5, 5);
  TH1F *h_scat_e_phi_cuts_pT2_cut = new TH1F("h_scat_e_phi_cuts_pT2_cut", "h_scat_e_phi_cuts_pT2_cut", 360, -TMath::Pi(), TMath::Pi());

  TH1F *h_scat_e_pT_RC_cuts_pT2_cut = new TH1F("h_scat_e_pT_RC_cuts_pT2_cut", "h_scat_e_pT_RC_cuts_pT2_cut", 100, 0, 10);
  TH1F *h_scat_e_eta_RC_cuts_pT2_cut = new TH1F("h_scat_e_eta_RC_cuts_pT2_cut", "h_scat_e_eta_RC_cuts_pT2_cut", 100, -5, 5);
  TH1F *h_scat_e_phi_RC_cuts_pT2_cut = new TH1F("h_scat_e_phi_RC_cuts_pT2_cut", "h_scat_e_phi_RC_cuts_pT2_cut", 360, -TMath::Pi(), TMath::Pi());

  //MC vs. RC comparison
  TH2F *h_scat_e_pT_MC_vs_e_pT_RC_cuts_pT2_cut = new TH2F("h_scat_e_pT_MC_vs_e_pT_RC_cuts_pT2_cut", "h_scat_e_pT_MC_vs_e_pT_RC_cuts_pT2_cut", 100, 0, 10, 100, 0, 10);
  TH2F *h_scat_e_eta_MC_vs_e_eta_RC_cuts_pT2_cut = new TH2F("h_scat_e_eta_MC_vs_e_eta_RC_cuts_pT2_cut", "h_scat_e_eta_MC_vs_e_eta_RC_cuts_pT2_cut", 100, -5, 5, 100, -5, 5);
  TH2F *h_scat_e_phi_MC_vs_e_phi_RC_cuts_pT2_cut = new TH2F("h_scat_e_phi_MC_vs_e_phi_RC_cuts_pT2_cut", "h_scat_e_phi_MC_vs_e_phi_RC_cuts_pT2_cut", 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());

  //--------------------------

  //spectator neutron
  TH1F *h_pT = new TH1F("h_pT","h_pT", 100, 0,1);
  TH1F *h_pT2 = new TH1F("h_pT2","h_pT2", 100, 0,1);

  TH1F *h_phi = new TH1F("h_phi", "h_phi", 360, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta = new TH1F("h_eta", "h_eta", 100, 0, 10);


  TH1F *h_pT_RC = new TH1F("h_pT_RC","h_pT_RC", 100, 0,1);
  TH1F *h_pT2_RC = new TH1F("h_pT2_RC","h_pT2_RC", 100, 0,1);

  TH1F *h_phi_RC = new TH1F("h_phi_RC", "h_phi_RC", 360, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta_RC = new TH1F("h_eta_RC", "h_eta_RC", 100, 0, 10);

  //----

  TH1F *h_alpha_p = new TH1F("h_alpha_p","h_alpha_p", 100, 0,2);
  TH1F *h_alpha_p_RC = new TH1F("h_alpha_p_RC","h_alpha_p_RC", 100, 0,2);

  TH1F *h_alpha_p_pT2_cut = new TH1F("h_alpha_p_pT2_cut","h_alpha_p_pT2_cut", 100, 0,2);
  TH1F *h_alpha_p_RC_pT2_cut = new TH1F("h_alpha_p_RC_pT2_cut","h_alpha_p_RC_pT2_cut", 100, 0,2);

  //-------------------

  //MC vs. RC comparison

  TH2F *h_pT_MC_vs_pT_RC = new TH2F("h_pT_MC_vs_pT_RC","h_pT_MC_vs_pT_RC", 100, 0, 1, 100, 0,1);
  TH2F *h_phi_MC_vs_phi_RC = new TH2F("h_phi_MC_vs_phi_RC", "h_phi_MC_vs_phi_RC", 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());
  TH2F *h_eta_MC_vs_eta_RC = new TH2F("h_eta_MC_vs_eta_RC", "h_eta_MC_vs_eta_RC", 100, 0, 10, 100, 0, 10);

  TH2F *h_E_MC_vs_E_RC = new TH2F("h_E_MC_vs_E_RC", "h_E_MC_vs_E_RC", 260, 0, 260, 260, 0, 260);

  //------
  //Resolution plots
  TH1F *h_theta_resolution = new TH1F("h_theta_resolution", "h_theta_resolution", 100, -5, 5);
  
  TH1F *h_theta_resolution_theta_1 = new TH1F("h_theta_resolution_theta_1", "h_theta_resolution_theta_1", 20, -1, 1); //theta < 1 mrad
  TH1F *h_theta_resolution_theta_2 = new TH1F("h_theta_resolution_theta_2", "h_theta_resolution_theta_2", 20, -1, 1); //theta > 1 mrad
  
  //----
  
  TH1F *h_theta_resolution_alpha = new TH1F("h_theta_resolution_alpha", "h_theta_resolution_alpha", 100, -5, 5);
  
  TH1F *h_theta_resolution_alpha_theta_1 = new TH1F("h_theta_resolution_alpha_theta_1", "h_theta_resolution_alpha_theta_1", 20, -1, 1); //theta < 1 mrad
  TH1F *h_theta_resolution_alpha_theta_2 = new TH1F("h_theta_resolution_alpha_theta_2", "h_theta_resolution_alpha_theta_2", 20, -1, 1); //theta > 1 mrad  
  
  //----
  
  TH1F *h_theta_resolution_cut = new TH1F("h_theta_resolution_cut", "h_theta_resolution_cut", 100, -5, 5);
  
  TH1F *h_theta_resolution_cut_theta_1 = new TH1F("h_theta_resolution_cut_theta_1", "h_theta_resolution_cut_theta_1", 20, -1, 1); //theta < 1 mrad
  TH1F *h_theta_resolution_cut_theta_2 = new TH1F("h_theta_resolution_cut_theta_2", "h_theta_resolution_cut_theta_2", 20, -1, 1); //theta > 1 mrad
  
  //TH1F *h_R_resolution = new TH1F("h_R_resolution", "h_R_resolution", 100, -100, 100);

  TH1F *h_theta_resolution_cut_pT2_bins[n_pT2_bins];
  TH1F *h_theta_resolution_cut_pT2_bins_MC[n_pT2_bins];
  TH1F *h_theta_resolution_cut_pT2_bins_beam_eff[n_pT2_bins];

  for(unsigned int i_pT2_bin = 0; i_pT2_bin < n_pT2_bins; i_pT2_bin++)
  {
    h_theta_resolution_cut_pT2_bins[i_pT2_bin] = new TH1F(Form("h_theta_resolution_cut_%i", i_pT2_bin), Form("h_theta_resolution_cut_%i", i_pT2_bin), 40, -1, 1);
    h_theta_resolution_cut_pT2_bins_MC[i_pT2_bin] = new TH1F(Form("h_theta_resolution_cut_MC_%i", i_pT2_bin), Form("h_theta_resolution_cut_MC_%i", i_pT2_bin), 40, -1, 1);
    h_theta_resolution_cut_pT2_bins_beam_eff[i_pT2_bin] = new TH1F(Form("h_theta_resolution_cut_beam_eff_%i", i_pT2_bin), Form("h_theta_resolution_cut_beam_eff_%i", i_pT2_bin), 40, -1, 1);

  }



  TH1F *h_E_resolution = new TH1F("h_E_resolution", "h_E_resolution", 100, -50, 50);
  TH1F *h_E_resolution_cut = new TH1F("h_E_resolution_cut", "h_E_resolution_cut", 100, -50, 50);

  //------

  TH2F *h_alpha_p_MC_vs_alpha_p_RC = new TH2F("h_alpha_p_MC_vs_alpha_p_RC","h_alpha_p_MC_vs_alpha_p_RC", 100, 0, 2, 100, 0, 2);
  TH2F *h_alpha_p_MC_vs_alpha_p_RC_pT2_cut = new TH2F("h_alpha_p_MC_vs_alpha_p_RC_pT2_cut","h_alpha_p_MC_vs_alpha_p_RC_pT2_cut", 100, 0, 2, 100, 0, 2);

  //----------------------------------

	//after "cuts" (especially for alpha_p)
	TH1F *h_pT_cuts = new TH1F("h_pT_cuts","h_pT_cuts", 100, 0,1);
  TH1F *h_phi_cuts = new TH1F("h_phi_cuts", "h_phi_cuts", 360, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta_cuts = new TH1F("h_eta_cuts", "h_eta_cuts", 100, 0, 10);


  TH1F *h_pT_RC_cuts = new TH1F("h_pT_RC_cuts","h_pT_RC_cuts", 100, 0,1);
  TH1F *h_phi_RC_cuts = new TH1F("h_phi_RC_cuts", "h_phi_RC_cuts", 360, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta_RC_cuts = new TH1F("h_eta_RC_cuts", "h_eta_RC_cuts", 100, 0, 10);

  //------------------------

  //MC vs. RC comparison

  TH2F *h_pT_MC_vs_pT_RC_cuts = new TH2F("h_pT_MC_vs_pT_RC_cuts","h_pT_MC_vs_pT_RC_cuts", 100, 0, 1, 100, 0,1);
  TH2F *h_phi_MC_vs_phi_RC_cuts = new TH2F("h_phi_MC_vs_phi_RC_cuts", "h_phi_MC_vs_phi_RC_cuts", 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());
  TH2F *h_eta_MC_vs_eta_RC_cuts = new TH2F("h_eta_MC_vs_eta_RC_cuts", "h_eta_MC_vs_eta_RC_cuts", 100, 0, 10, 100, 0, 10);

  TH2F *h_E_MC_vs_E_RC_cuts = new TH2F("h_E_MC_vs_E_RC_cuts", "h_E_MC_vs_E_RC_cuts", 260, 0, 260, 260, 0, 260);

  //----------------------------------

	//after "cuts" (especially for alpha_p)
	TH1F *h_pT_cuts_pT2_cut = new TH1F("h_pT_cuts_pT2_cut","h_pT_cuts_pT2_cut", 100, 0,1);
  TH1F *h_phi_cuts_pT2_cut = new TH1F("h_phi_cuts_pT2_cut", "h_phi_cuts_pT2_cut", 360, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta_cuts_pT2_cut = new TH1F("h_eta_cuts_pT2_cut", "h_eta_cuts_pT2_cut", 100, 0, 10);


  TH1F *h_pT_RC_cuts_pT2_cut = new TH1F("h_pT_RC_cuts_pT2_cut","h_pT_RC_cuts_pT2_cut", 100, 0,1);
  TH1F *h_phi_RC_cuts_pT2_cut = new TH1F("h_phi_RC_cuts_pT2_cut", "h_phi_RC_cuts_pT2_cut", 360, -TMath::Pi(), TMath::Pi());
  TH1F *h_eta_RC_cuts_pT2_cut = new TH1F("h_eta_RC_cuts_pT2_cut", "h_eta_RC_cuts_pT2_cut", 100, 0, 10);

  //------------------------

  //MC vs. RC comparison

  TH2F *h_pT_MC_vs_pT_RC_cuts_pT2_cut = new TH2F("h_pT_MC_vs_pT_RC_cuts_pT2_cut","h_pT_MC_vs_pT_RC_cuts_pT2_cut", 100, 0, 1, 100, 0,1);
  TH2F *h_phi_MC_vs_phi_RC_cuts_pT2_cut = new TH2F("h_phi_MC_vs_phi_RC_cuts_pT2_cut", "h_phi_MC_vs_phi_RC_cuts_pT2_cut", 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());
  TH2F *h_eta_MC_vs_eta_RC_cuts_pT2_cut = new TH2F("h_eta_MC_vs_eta_RC_cuts_pT2_cut", "h_eta_MC_vs_eta_RC_cuts_pT2_cut", 100, 0, 10, 100, 0, 10);


  //------------------------
  //ZDC QA hists

  TH2F *h_ZDC_cluster_x_vs_y = new TH2F("h_ZDC_cluster_x_vs_y", "h_ZDC_cluster_x_vs_y", 300, -1200, -600, 300, -300, 300);
  TH1F *h_ZDC_cluster_z = new TH1F("h_ZDC_cluster_z", "h_ZDC_cluster_z", 100, 35800, 37800);

  TH2F *h_ZDC_cluster_x_vs_y_no_X_angle = new TH2F("h_ZDC_cluster_x_vs_y_no_X_angle", "h_ZDC_cluster_x_vs_y_no_X_angle", 300, -300, 300, 300, -300, 300);


  TH2F *h_ZDC_cluster_x_vs_y_no_X_angle_z_bins[n_ZDC_z_bins];

  for(unsigned int ZDC_z_bin = 0; ZDC_z_bin < n_ZDC_z_bins; ZDC_z_bin++)
  {
    h_ZDC_cluster_x_vs_y_no_X_angle_z_bins[ZDC_z_bin] = new TH2F(Form("h_ZDC_cluster_x_vs_y_no_X_angle_z_bins_%i", ZDC_z_bin), Form("h_ZDC_cluster_x_vs_y_no_X_angle_z_bins_%i", ZDC_z_bin), 300, -300, 300, 300, -300, 300);
  }


  //const int n_theta_bins = 4;
  TH1F *h_phi_RC_theta_bins[n_theta_bins]; //QA
  TH1F *h_phi_RC_theta_bins_alpha[n_theta_bins]; //QA
  TH1F *h_phi_RC_theta_bins_cuts[n_theta_bins]; //QA

  for(unsigned int theta_bin = 0; theta_bin < n_theta_bins; theta_bin++)
  {
    h_phi_RC_theta_bins[theta_bin] = new TH1F(Form("h_phi_RC_theta_bin_%i", theta_bin), Form("h_phi_RC_theta_bin_%i", theta_bin), 360, -TMath::Pi(), TMath::Pi());
    h_phi_RC_theta_bins_alpha[theta_bin] = new TH1F(Form("h_phi_RC_theta_alpha_bin_%i", theta_bin), Form("h_phi_RC_theta_aplha_bin_%i", theta_bin), 360, -TMath::Pi(), TMath::Pi());
    h_phi_RC_theta_bins_cuts[theta_bin] = new TH1F(Form("h_phi_RC_theta_cuts_bin_%i", theta_bin), Form("h_phi_RC_theta_cuts_bin_%i", theta_bin), 360, -TMath::Pi(), TMath::Pi());
  }
  


	//------------------
	//for efficiency
	TH1F *h_nEvents_bins[nx_bins][nQ2_bins];

	TH1F *h_pT2_bins[nx_bins][nQ2_bins];
	TH1F *h_pT2_MCPart_bins[nx_bins][nQ2_bins];
	TH1F *h_pT2_RC_bins[nx_bins][nQ2_bins];
	TH1F *h_pT2_MC_ACC_bins[nx_bins][nQ2_bins];

	//true MC reduced cross section
	TH1F *h_N_vs_x[nQ2_bins];
	TH1F *h_sigma_D_red_vs_x[nQ2_bins];

	TH1F *h_sigma_D_red[nx_bins][nQ2_bins];
  TH1F *h_sigma_p_red[nx_bins][nQ2_bins];
  
  //----

  TH1F *h_sigma_D_red_MCPart[nx_bins][nQ2_bins];
  TH1F *h_sigma_p_red_MCPart[nx_bins][nQ2_bins];
  
  //----
  
  TH1F *h_N_vs_x_RC[nQ2_bins];
  TH1F *h_sigma_D_red_vs_x_RC[nQ2_bins];

  TH1F *h_sigma_D_red_RC[nx_bins][nQ2_bins];
  TH1F *h_sigma_p_red_RC[nx_bins][nQ2_bins];
  
  //----
  
  TH1F *h_N_vs_x_MC_ACC[nQ2_bins];
  TH1F *h_sigma_D_red_vs_x_MC_ACC[nQ2_bins];

  TH1F *h_sigma_D_red_MC_ACC[nx_bins][nQ2_bins];
  TH1F *h_sigma_p_red_MC_ACC[nx_bins][nQ2_bins];

  for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
  {
    for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
    {
      h_nEvents_bins[x_bin][Q2_bin] = new TH1F(Form("h_nEvents_bins_xbin_%i_Q2bin_%i", x_bin, Q2_bin), Form("h_nEvents_bins_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 3, 0, 3);

	    h_pT2_bins[x_bin][Q2_bin] = new TH1F(Form("h_pT2_bins_xbin_%i_Q2bin_%i", x_bin, Q2_bin), Form("h_pT2_bins_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
	    h_pT2_MCPart_bins[x_bin][Q2_bin] = new TH1F(Form("h_pT2_MCPart_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), Form("h_pT2_MCPart_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
	    h_pT2_RC_bins[x_bin][Q2_bin] = new TH1F(Form("h_pT2_RC_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), Form("h_pT2_RC_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
	    h_pT2_MC_ACC_bins[x_bin][Q2_bin] = new TH1F(Form("h_pT2_MC_ACC_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), Form("h_pT2_MC_ACC_bins_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);

      //---------------------------------

      h_sigma_D_red[x_bin][Q2_bin] = new TH1F(Form("h_sigma_D_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_D_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
      h_sigma_p_red[x_bin][Q2_bin] = new TH1F(Form("h_sigma_p_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_p_red_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
      
      h_sigma_D_red_MCPart[x_bin][Q2_bin] = new TH1F(Form("h_sigma_D_red_MCPart_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_D_red_MCPart_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
      h_sigma_p_red_MCPart[x_bin][Q2_bin] = new TH1F(Form("h_sigma_p_red_MCPart_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_p_red_MCPart_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);

      h_sigma_D_red_RC[x_bin][Q2_bin] = new TH1F(Form("h_sigma_D_red_RC_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_D_red_RC_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
      h_sigma_p_red_RC[x_bin][Q2_bin] = new TH1F(Form("h_sigma_p_red_RC_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_p_red_RC_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);

      h_sigma_D_red_MC_ACC[x_bin][Q2_bin] = new TH1F(Form("h_sigma_D_red_MC_ACC_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_D_red_MC_ACC_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
      h_sigma_p_red_MC_ACC[x_bin][Q2_bin] = new TH1F(Form("h_sigma_p_red_MC_ACC_xbin_%i_Q2bin_%i", x_bin, Q2_bin),Form("h_sigma_p_red_MC_ACC_xbin_%i_Q2bin_%i", x_bin, Q2_bin), 20, 0, 0.01);
      
      //----
      
      if(x_bin == 0)
      {
        h_sigma_D_red_vs_x[Q2_bin] = new TH1F(Form("h_sigma_D_red_vs_x_%i", Q2_bin), Form("h_sigma_D_red_vs_x_%i", Q2_bin), nx_bins, x_bins);
        h_sigma_D_red_vs_x_RC[Q2_bin] = new TH1F(Form("h_sigma_D_red_vs_x_RC_%i", Q2_bin), Form("h_sigma_D_red_vs_x_RC_%i", Q2_bin), nx_bins, x_bins);
        h_sigma_D_red_vs_x_MC_ACC[Q2_bin] = new TH1F(Form("h_sigma_D_red_vs_x_MC_ACC_%i", Q2_bin), Form("h_sigma_D_red_vs_x_MC_ACC_%i", Q2_bin), nx_bins, x_bins);      
        
        h_N_vs_x[Q2_bin] = new TH1F(Form("h_N_vs_x_%i", Q2_bin), Form("h_N_vs_x_%i", Q2_bin), nx_bins, x_bins);
        h_N_vs_x_RC[Q2_bin] = new TH1F(Form("h_N_vs_x_RC_%i", Q2_bin), Form("h_N_vs_x_RC_%i", Q2_bin), nx_bins, x_bins);
        h_N_vs_x_MC_ACC[Q2_bin] = new TH1F(Form("h_N_vs_x_MC_ACC_%i", Q2_bin), Form("h_N_vs_x_MC_ACC_%i", Q2_bin), nx_bins, x_bins); 
      }
    }
  }


	//--------------------------------------------------------------------------


	int fileCounter = 0;

	int iEvent = 0;
	int iEvent_MC = 0;
	int iEvent_RC = 0;

	while(getline(fileListStream, fileName) )
	{

	  TString tmp = fileName;

	  cout << "Input file " << fileCounter << ": " << fileName << endl;

	  //inputRootFile = new TFile(tmp, "read");
	  auto inputRootFile = std::unique_ptr<TFile>{TFile::Open(tmp)};

		if(inputRootFile->IsZombie()){ cout << "MISSING_ROOT_FILE"<< fileName << endl; continue;}

		fileCounter++;

		TTree * evtTree = (TTree*)inputRootFile->Get("events");

    //------------------------------------------------------------------

		int numEvents = evtTree->GetEntries();

  	TTreeReader tree_reader(evtTree);       // !the tree reader

  	//Verticies
    TTreeReaderArray<float> vertex_x = {tree_reader, "CentralTrackVertices.position.x"};
    TTreeReaderArray<float> vertex_y = {tree_reader, "CentralTrackVertices.position.y"};
    TTreeReaderArray<float> vertex_z = {tree_reader, "CentralTrackVertices.position.z"};

		//MCParticles (only for cross check with head-on frame)
  	TTreeReaderArray<double> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
  	TTreeReaderArray<double> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
  	TTreeReaderArray<double> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
  	TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
  	TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};
		TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
		TTreeReaderArray<unsigned int> mc_parents_begin = {tree_reader, "MCParticles.parents_begin"};
    TTreeReaderArray<unsigned int> mc_parents_end = {tree_reader, "MCParticles.parents_end"};

    TTreeReaderArray<int> mc_parents = {tree_reader, "_MCParticles_parents.index"}; //ID of parents in MCParticlesHeadOnFrameNoBeamFX array
    
    //-----

		TTreeReaderArray<double> truth_px_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.x"};
    TTreeReaderArray<double> truth_py_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.y"};
    TTreeReaderArray<double> truth_pz_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.z"};
    TTreeReaderArray<double> truth_mass_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.mass"};
    TTreeReaderArray<int> truth_pdg_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.PDG"};
    TTreeReaderArray<int> truth_genStatus_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.generatorStatus"};
    TTreeReaderArray<unsigned int> truth_parents_begin = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.parents_begin"};
    TTreeReaderArray<unsigned int> truth_parents_end = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.parents_end"};

    TTreeReaderArray<int> truth_parents = {tree_reader, "_MCParticlesHeadOnFrameNoBeamFX_parents.index"}; //ID of parents in MCParticlesHeadOnFrameNoBeamFX array

    //-----
	
    //reco ZDC neutral perticles
		TTreeReaderArray<float> reco_ZDC_px = {tree_reader, "ReconstructedHcalFarForwardZDCNeutrals.momentum.x"};
    TTreeReaderArray<float> reco_ZDC_py = {tree_reader, "ReconstructedHcalFarForwardZDCNeutrals.momentum.y"};
    TTreeReaderArray<float> reco_ZDC_pz = {tree_reader, "ReconstructedHcalFarForwardZDCNeutrals.momentum.z"};
    TTreeReaderArray<float> reco_ZDC_E = {tree_reader, "ReconstructedHcalFarForwardZDCNeutrals.energy"};

    TTreeReaderArray<float> reco_ZDC_x = {tree_reader, "ReconstructedHcalFarForwardZDCNeutrals.referencePoint.x"};
    TTreeReaderArray<float> reco_ZDC_y = {tree_reader, "ReconstructedHcalFarForwardZDCNeutrals.referencePoint.y"};
    TTreeReaderArray<float> reco_ZDC_z = {tree_reader, "ReconstructedHcalFarForwardZDCNeutrals.referencePoint.z"};

    //-----

		//reco tracks
		TTreeReaderArray<float> reco_track_x = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  	TTreeReaderArray<float> reco_track_y = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  	TTreeReaderArray<float> reco_track_z = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  	TTreeReaderArray<float> reco_track_charge = {tree_reader, "ReconstructedChargedParticles.charge"};
		TTreeReaderArray<int> reco_track_PDG = {tree_reader, "ReconstructedChargedParticles.PDG"};

		//MC RC charged track associations
  	TTreeReaderArray<int> trk_simID = {tree_reader, "_ReconstructedChargedParticleAssociations_sim.index"};
  	TTreeReaderArray<int> trk_recID = {tree_reader, "_ReconstructedChargedParticleAssociations_rec.index"};

  	//---------------------------

		//forward EMCAL endcap
		TTreeReaderArray<float> nECAL_cluster_x = {tree_reader, "EcalEndcapNClusters.position.x"};
    TTreeReaderArray<float> nECAL_cluster_y = {tree_reader, "EcalEndcapNClusters.position.y"};
    TTreeReaderArray<float> nECAL_cluster_z = {tree_reader, "EcalEndcapNClusters.position.z"};
    TTreeReaderArray<float> nECAL_cluster_energy = {tree_reader, "EcalEndcapNClusters.energy"};

    //MC RC nECAL associations
  	TTreeReaderArray<int> nECAL_simID = {tree_reader, "_EcalEndcapNClusterAssociations_sim.index"};
  	TTreeReaderArray<int> nECAL_recID = {tree_reader, "_EcalEndcapNClusterAssociations_rec.index"};

  	//cluster to charged track association
  	TTreeReaderArray<int> nECAL_clust_ID = {tree_reader, "_EcalEndcapNTrackClusterMatches_cluster.index"};
  	TTreeReaderArray<int> nECAL_clust_recTrackID = {tree_reader, "_EcalEndcapNTrackClusterMatches_track.index"};

    //---------------------------

    //barrel ECAL clusters
    TTreeReaderArray<float> bECAL_cluster_x = {tree_reader, "EcalBarrelClusters.position.x"};
    TTreeReaderArray<float> bECAL_cluster_y = {tree_reader, "EcalBarrelClusters.position.y"};
    TTreeReaderArray<float> bECAL_cluster_z = {tree_reader, "EcalBarrelClusters.position.z"};
    TTreeReaderArray<float> bECAL_cluster_energy = {tree_reader, "EcalBarrelClusters.energy"};

    TTreeReaderArray<unsigned int> bECAL_cluster_begin = {tree_reader, "EcalBarrelClusters.clusters_begin"};
    TTreeReaderArray<unsigned int> bECAL_cluster_end = {tree_reader, "EcalBarrelClusters.clusters_end"};

    //cluster to charged track association
  	TTreeReaderArray<int> bECAL_clust_ID = {tree_reader, "_EcalBarrelTrackClusterMatches_cluster.index"};
  	TTreeReaderArray<int> bECAL_clust_recTrackID = {tree_reader, "_EcalBarrelTrackClusterMatches_track.index"};

  	//cluster IDs for asociation (for global, SciFi, imaging?)
  	//accessed by clusters_begin and clusters_end
  	//This should then be index in association arrays above
  	TTreeReaderArray<int> bECAL_clust_clusterID = {tree_reader, "_EcalBarrelClusters_clusters.index"};

  	//RC - MC matching for bECAL
  	TTreeReaderArray<int> bECAL_simID = {tree_reader, "_EcalBarrelClusterAssociations_sim.index"};
  	TTreeReaderArray<int> bECAL_recID = {tree_reader, "_EcalBarrelClusterAssociations_rec.index"};

  	//---------------------------

		cout << "file has " << evtTree->GetEntries() <<  " events..." << endl;

		tree_reader.SetEntriesRange(0, evtTree->GetEntries());
		while (tree_reader.Next())
		{

			if( iEvent % 10000 == 0 ) cout<<"Working on event "<<iEvent<<endl;

	    //MCParticles
	    //finding the far-forward neutron;

			TVector3 mctrk;
			TVector3 rctrk;
			TVector3 final_mc_track;

			TLorentzVector scatteredElectron_MC(-9999, -9999, -9999, -9999);
			TLorentzVector scatteredElectron_TRUTH(-9999, -9999, -9999, -9999);
			TLorentzVector scatteredElectron_RECO(-9999, -9999, -9999, -9999);

			//---

			TLorentzVector photonVect_MC(-9999, -9999, -9999, -9999);
			TLorentzVector photonVect_TRUTH(-9999, -9999, -9999, -9999);
			TLorentzVector photonVect_RECO(-9999, -9999, -9999, -9999);

			//---

			TLorentzVector spectator_MC(-9999, -9999, -9999, -9999);
			TLorentzVector spectator_MC_2(-9999, -9999, -9999, -9999);			
			
      TLorentzVector spectator_TRUTH(-9999, -9999, -9999, -9999);
      TLorentzVector spectator_TRUTH_2(-9999, -9999, -9999, -9999);

      TLorentzVector spectator_RECO(-9999, -9999, -9999, -9999);
      TLorentzVector spectator_RECO_2(-9999, -9999, -9999, -9999);


			//electron and deuteron beam vectors and boost to IRF
			TLorentzVector deuteron_beam(0.0, 0.0, 260.0, TMath::Sqrt(260*260 + mDeuteron*mDeuteron));
			TLorentzVector electron_beam(0.0, 0.0, -10.0, TMath::Sqrt(10*10 + mElectron*mElectron));


			TVector3 boost_IRF = deuteron_beam.BoostVector();

      float calculatedQ2_MC = -9999;
      int Q2_bin_event_MC = -1;

			float calculatedX_bjorken_MC = -9999;
			float calculatedX_bjorken_nucleon_MC = -9999;
			int x_bin_event_MC = -1;
			
			float calculated_y_MC = -9999;
			
			float Flux_MC = -9999;
			
			float one_minus_epsilon_MC;

      //--------

			float calculatedQ2_TRUTH = -9999;
			int Q2_bin_event_TRUTH = -1;

			float calculatedX_bjorken_TRUTH = -9999;
			float calculatedX_bjorken_nucleon_TRUTH = -9999;
			int x_bin_event_TRUTH = -1;

			float calculated_y_TRUTH = -9999;

			float Flux_TRUTH = -9999;

			float one_minus_epsilon_TRUTH;

			//--------

			float calculatedQ2_RC = -9999;
			int Q2_bin_event_RC = -1;

			float calculatedX_bjorken_RC = -9999;
			float calculatedX_bjorken_nucleon_RC = -9999;
			int x_bin_event_RC = -1;

			float calculated_y_RC = -9999;

			float Flux_RC = -9999;

			float one_minus_epsilon_RC;

			//--------

			float alphaSpect_MC = -9999;
			float alphaSpect_TRUTH = -9999;
			float alphaSpect_RECO = -9999;

      //for MC to RC matching
      int mcElectronIdx = -99; //MC track index
      int rcElectronIdx_trk = -99; //RC charged track index
      int rcElectronIdx_nECAL = -99; //RC cluster in backward ECAL

			int mcNeutronIdx = -99;


      int numOfElectrons_MC = 0;
			int numOfNeutrons_MC = 0;

			int numOfElectrons_TRUTH = 0;
			int numOfNeutrons_TRUTH = 0;

			int numOfElectrons_RC = 0;
			int numOfNeutrons_RC = 0;
			
			int beamParentElectronFound = 0;
			int beamParentElectronFound_MCParticle = 0;
			
			
			//loop over truth MC information (without beam effects and crossing angle)
			for(int imc = 0; imc < mc_px_array.GetSize(); imc++)
			{
				mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);

				//find the scattered electron
				if(mc_pdg_array[imc] == 11 && mc_genStatus_array[imc] == 1 && beamParentElectronFound_MCParticle == 0 )
				{
          //mctrk.RotateY(0.025);

          float_t E_electron_MC = TMath::Sqrt(mctrk.Mag2() + mElectron*mElectron);
          
          //find beam electron as parent
          for( unsigned int parent_arr_index = mc_parents_begin[imc]; parent_arr_index <= mc_parents_end[imc]; parent_arr_index++ )
          {
            int parent_MC_index = mc_parents[parent_arr_index];

            if( mc_pdg_array[parent_MC_index] == 11 &&  mc_genStatus_array[parent_MC_index] == 4 ) //Is this correct? - Taking gen. status from PYTHIA.
            {
              beamParentElectronFound_MCParticle = 1;

            }

          }

          if(beamParentElectronFound_MCParticle != 1) continue;

					scatteredElectron_MC.SetPxPyPzE(mctrk.Px(), mctrk.Py(), mctrk.Pz(), E_electron_MC);

          photonVect_MC = electron_beam - scatteredElectron_MC;

          calculatedQ2_MC        = -photonVect_MC.Mag2();   //Jan

          calculatedX_bjorken_MC = calculatedQ2_MC/(deuteron_beam*photonVect_MC);       //Jan - updated, removed 1/2, see Eq. (5) in paper

          calculated_y_MC = ( deuteron_beam*photonVect_MC )/( deuteron_beam*electron_beam ); //Jan

          //---------------------

          one_minus_epsilon_MC = calculated_y_MC*calculated_y_MC/( 1. + (1.-calculated_y_MC)*(1.-calculated_y_MC) );

          //------------------------------------

          //find Q2 bin
          for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
          {
            if( calculatedQ2_MC > Q2_bins[Q2_bin] && calculatedQ2_MC < Q2_bins[Q2_bin+1] )
            {
              Q2_bin_event_MC = Q2_bin;
              break;
            }

          }

          numOfElectrons_MC++;

		    }

		    //------------------------------
		    
		    mctrk.RotateY(0.025);

		    //spectator neutron
				if(mc_pdg_array[imc] == 2112 && mc_genStatus_array[imc] == 1 /*&& numOfNeutrons == 0 */)
				{
				  //simple selection based on spectator kinematics
					if(mctrk.Eta() < 5.5){ continue; }

					TLorentzVector neutron_4_vect(mctrk.Px(), mctrk.Py(), mctrk.Pz(), TMath::Sqrt(mctrk.Mag2() + mNeutron*mNeutron));


					if( numOfNeutrons_MC == 0 ) spectator_MC = neutron_4_vect;

					if( numOfNeutrons_MC == 1 )
          {
            if( neutron_4_vect.Eta() > spectator_MC.Eta() )
            {
              spectator_MC_2 = spectator_MC;

              spectator_MC = neutron_4_vect;
            }
            else
            {
              spectator_MC_2 = neutron_4_vect;
            }
          }


					alphaSpect_MC = 2*(spectator_MC.E() + spectator_MC.Pz()) / ( deuteron_beam.E() + deuteron_beam.Pz() ); //in head-on lab frame - Jan

					calculatedX_bjorken_nucleon_MC = calculatedX_bjorken_MC/(2.-alphaSpect_MC);

					//find x bin
          for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
          {
            if( calculatedX_bjorken_nucleon_MC > x_bins[x_bin] && calculatedX_bjorken_nucleon_MC < x_bins[x_bin+1] )
            {
              x_bin_event_MC = x_bin;
              break;
            }

          }

					Flux_MC = 2.*TMath::Pi()*alpha_em*alpha_em*calculated_y_MC*calculated_y_MC / ( calculatedQ2_MC*calculatedQ2_MC*one_minus_epsilon_MC*calculatedX_bjorken_MC );

					//mcNeutronIdx = imc;

					numOfNeutrons_MC++;
				} //neutrons


      }//end loop over MCParticles


			//loop over truth MC information (without beam effects and crossing angle)
			for(int imc = 0; imc < truth_px_array.GetSize(); imc++)
			{
				mctrk.SetXYZ(truth_px_array[imc], truth_py_array[imc], truth_pz_array[imc]);

				//find the scattered electron
				if(truth_pdg_array[imc] == 11 && truth_genStatus_array[imc] == 1 && beamParentElectronFound == 0 )
				{
          //mctrk.RotateY(0.025);

          float_t E_electron_TRUTH = TMath::Sqrt(mctrk.Mag2() + mElectron*mElectron);

          //find beam electron as parent
          for( unsigned int parent_arr_index = truth_parents_begin[imc]; parent_arr_index <= truth_parents_end[imc]; parent_arr_index++ )
          {
            int parent_MC_index = truth_parents[parent_arr_index];

            if( truth_pdg_array[parent_MC_index] == 11 &&  truth_genStatus_array[parent_MC_index] == 4 ) //Is this correct? - Taking gen. status from PYTHIA.
            {
              beamParentElectronFound = 1;

            }

          }

          if(beamParentElectronFound != 1) continue;

					scatteredElectron_TRUTH.SetPxPyPzE(mctrk.Px(), mctrk.Py(), mctrk.Pz(), E_electron_TRUTH);

          photonVect_TRUTH = electron_beam - scatteredElectron_TRUTH;

          calculatedQ2_TRUTH        = -photonVect_TRUTH.Mag2();   //Jan
          
          calculatedX_bjorken_TRUTH = calculatedQ2_TRUTH/(deuteron_beam*photonVect_TRUTH);       //Jan - updated, removed 1/2, see Eq. (5) in paper


          calculated_y_TRUTH = ( deuteron_beam*photonVect_TRUTH )/( deuteron_beam*electron_beam ); //Jan

          //---------------------

          one_minus_epsilon_TRUTH = calculated_y_TRUTH*calculated_y_TRUTH/( 1. + (1.-calculated_y_TRUTH)*(1.-calculated_y_TRUTH) );

          //------------------------------------

          //find Q2 bin
          for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
          {
            if( calculatedQ2_TRUTH > Q2_bins[Q2_bin] && calculatedQ2_TRUTH < Q2_bins[Q2_bin+1] )
            {
              Q2_bin_event_TRUTH = Q2_bin;
              break;
            }

          }

          mcElectronIdx = imc;

          numOfElectrons_TRUTH++;

		    }

		    //------------------------------

		    //spectator neutron
				if(truth_pdg_array[imc] == 2112 && truth_genStatus_array[imc] == 1 /*&& numOfNeutrons == 0 */)
				{
				  //simple selection based on spectator kinematics
					if(mctrk.Eta() < 5.5){ continue; }

					TLorentzVector neutron_4_vect(mctrk.Px(), mctrk.Py(), mctrk.Pz(), TMath::Sqrt(mctrk.Mag2() + mNeutron*mNeutron));

					if( numOfNeutrons_TRUTH == 0 ) spectator_TRUTH = neutron_4_vect;

					if( numOfNeutrons_TRUTH == 1 )
          {
            if( neutron_4_vect.Eta() > spectator_TRUTH.Eta() )
            {
              spectator_TRUTH_2 = spectator_TRUTH;

              spectator_TRUTH = neutron_4_vect;
            }
            else
            {
              spectator_TRUTH_2 = neutron_4_vect;
            }
          }


					alphaSpect_TRUTH = 2*(spectator_TRUTH.E() + spectator_TRUTH.Pz()) / ( deuteron_beam.E() + deuteron_beam.Pz() ); //in head-on lab frame - Jan

					calculatedX_bjorken_nucleon_TRUTH = calculatedX_bjorken_TRUTH/(2.-alphaSpect_TRUTH);

					//find x bin
          for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
          {
            if( calculatedX_bjorken_nucleon_TRUTH > x_bins[x_bin] && calculatedX_bjorken_nucleon_TRUTH < x_bins[x_bin+1] )
            {
              x_bin_event_TRUTH = x_bin;
              break;
            }

          }

					Flux_TRUTH = 2.*TMath::Pi()*alpha_em*alpha_em*calculated_y_TRUTH*calculated_y_TRUTH / ( calculatedQ2_TRUTH*calculatedQ2_TRUTH*one_minus_epsilon_TRUTH*calculatedX_bjorken_TRUTH );

					mcNeutronIdx = imc;

					numOfNeutrons_TRUTH++;
				} //neutrons


      }//end loop over truth MC

      //----------------------------------------------------------------

      //RC scattered electron

      //First, find matching RC scattered electron ID in charged particle array
      //double-check that this is working
      for(unsigned int assoc_i = 0; assoc_i < trk_simID.GetSize(); assoc_i++)
			{
			  if( trk_simID[assoc_i] == mcElectronIdx )
			  {
			    rcElectronIdx_trk = trk_recID[assoc_i];

			    break;
			  }

			}

			for(int rc_trk_i = 0; rc_trk_i < reco_track_x.GetSize(); rc_trk_i++)
      {
        TVector3 rc_trk(reco_track_x[rc_trk_i], reco_track_y[rc_trk_i], reco_track_z[rc_trk_i]);

        if( reco_track_charge[rc_trk_i] == -1 ) h_track_eta_RC->Fill(rc_trk.Eta());

      }


      for(int clust_i = 0; clust_i < bECAL_cluster_energy.GetSize(); clust_i++)
      {

        TVector3 bECAL_hit_vector(bECAL_cluster_x[clust_i], bECAL_cluster_y[clust_i], bECAL_cluster_z[clust_i]);

        h_bECAL_cluster_eta_RC->Fill(bECAL_hit_vector.Eta());
      }


			//---------------------------
			//Loop over negative ecal clusters
			//Find matching RC track
			//Accept only negatively charged tracks with 0.8 < E/p < 1.2
			//Find track with max. pT (from Win). Maybe optimize to max. p or E?

      float max_trk = 0; //maximum energy/momentum/pT (pick the optimal version)

			for(int clust_i = 0; clust_i < nECAL_cluster_energy.GetSize(); clust_i++)
      {

        TVector3 nECAL_hit_vector(nECAL_cluster_x[clust_i], nECAL_cluster_y[clust_i], nECAL_cluster_z[clust_i]);

        h_nECAL_cluster_eta_RC->Fill(nECAL_hit_vector.Eta());


        int rcTrkClustMatchID = -1;

        //first find matching RC track to the cluster
        for(int assoc_clust_i = 0; assoc_clust_i < nECAL_clust_ID.GetSize(); assoc_clust_i++)
        {
          if( nECAL_clust_ID[assoc_clust_i] == clust_i )
          {
            rcTrkClustMatchID = nECAL_clust_recTrackID[assoc_clust_i];

            break;
          }

        }

        if( rcTrkClustMatchID < 0 ) continue; //skip, if no matching RC track is found

        if( reco_track_charge[rcTrkClustMatchID] != -1 ) continue; //accept only negatively charged tracks

        rctrk.SetXYZ(reco_track_x[rcTrkClustMatchID], reco_track_y[rcTrkClustMatchID], reco_track_z[rcTrkClustMatchID]);

        float_t E_electron_RECO = nECAL_cluster_energy[clust_i];


        if( E_electron_RECO/rctrk.Mag() > E_over_p_cut_min && E_electron_RECO/rctrk.Mag() < E_over_p_cut_max )
        {
          //get momentum vector with size given by energy in calo and direction from tracking
          rctrk.SetMag(TMath::Sqrt(E_electron_RECO*E_electron_RECO - mElectron*mElectron));

          scatteredElectron_RECO.SetPxPyPzE(rctrk.Px(), rctrk.Py(), rctrk.Pz(), E_electron_RECO);

          if( rctrk.Mag() > max_trk )
          {
            max_trk = rctrk.Mag();
          }
          else
          {
            continue;
          }

          photonVect_RECO = electron_beam - scatteredElectron_RECO;

          calculatedQ2_RC        = -photonVect_RECO.Mag2();   //Jan

          calculatedX_bjorken_RC = calculatedQ2_RC/(deuteron_beam*photonVect_RECO);       //Jan - updated, removed 1/2, see Eq. (5) in paper

          calculated_y_RC = ( deuteron_beam*photonVect_RECO )/( deuteron_beam*electron_beam ); //Jan

          //---------------------

          one_minus_epsilon_RC = calculated_y_RC*calculated_y_RC/( 1. + (1.-calculated_y_RC)*(1.-calculated_y_RC) );

          //------------------------------------

          //find Q2 bin
          for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
          {
            if( calculatedQ2_RC > Q2_bins[Q2_bin] && calculatedQ2_RC < Q2_bins[Q2_bin+1] )
            {
              Q2_bin_event_RC = Q2_bin;
              break;
            }

          }

        } //end if E/p cut

      }//end loop over N ecal clusers

      if(max_trk > 0) numOfElectrons_RC++;

      //------------------------------------------------------------------------------------

      //If no good scattered electron is found in nECAL, look into bECAL
      if( numOfElectrons_RC == 0 )
      {
        max_trk = 0; //reset this value to be safe


        for(int clust_i = 0; clust_i < bECAL_cluster_energy.GetSize(); clust_i++)
        {
          int clusterID = bECAL_clust_clusterID[bECAL_cluster_begin[clust_i]];

          int rcTrkClustMatchID = -1;
          int mcTrkClustMatchID = -1;

          for(unsigned int match_clust_i = 0; match_clust_i < bECAL_recID.GetSize(); match_clust_i++)
          {
            if( bECAL_recID[match_clust_i] == clust_i )
            {
              mcTrkClustMatchID = bECAL_simID[match_clust_i];

              break;
            }
          }

          if(mcTrkClustMatchID < 0 ) continue;

          for(unsigned int match_track_i = 0; match_track_i < trk_simID.GetSize(); match_track_i++)
          {
            if( trk_simID[match_track_i] == mcTrkClustMatchID )
            {
              rcTrkClustMatchID = trk_recID[match_track_i];

              break;
            }
          }

          if( rcTrkClustMatchID < 0 ) continue; //skip, if no matching RC track is found

          if( reco_track_charge[rcTrkClustMatchID] != -1 ) continue; //accept only negatively charged tracks

          rctrk.SetXYZ(reco_track_x[rcTrkClustMatchID], reco_track_y[rcTrkClustMatchID], reco_track_z[rcTrkClustMatchID]);

          float_t E_electron_RECO = bECAL_cluster_energy[clust_i];

          //cout<<"bECAL electron"<<endl;

          if( E_electron_RECO/rctrk.Mag() > E_over_p_cut_min && E_electron_RECO/rctrk.Mag() < E_over_p_cut_max )
          {

            //cout<<"bECAL electron"<<endl;

            //get momentum vector with size given by energy in calo and direction from tracking
            rctrk.SetMag(TMath::Sqrt(E_electron_RECO*E_electron_RECO - mElectron*mElectron));

            scatteredElectron_RECO.SetPxPyPzE(rctrk.Px(), rctrk.Py(), rctrk.Pz(), E_electron_RECO);

            //if( rctrk.Pt() > max_trk )
            if( rctrk.Mag() > max_trk )
            {
              //max_trk = rctrk.Pt();
              max_trk = rctrk.Mag();
            }
            else
            {
              continue;
            }

            photonVect_RECO = electron_beam - scatteredElectron_RECO;

            calculatedQ2_RC        = -photonVect_RECO.Mag2();   //Jan

            calculatedX_bjorken_RC = calculatedQ2_RC/(deuteron_beam*photonVect_RECO);       //Jan - updated, removed 1/2, see Eq. (5) in paper

            calculated_y_RC = ( deuteron_beam*photonVect_RECO )/( deuteron_beam*electron_beam ); //Jan

            //---------------------

            one_minus_epsilon_RC = calculated_y_RC*calculated_y_RC/( 1. + (1.-calculated_y_RC)*(1.-calculated_y_RC) );

            //------------------------------------

            //find Q2 bin
            for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
            {
              if( calculatedQ2_RC > Q2_bins[Q2_bin] && calculatedQ2_RC < Q2_bins[Q2_bin+1] )
              {
                Q2_bin_event_RC = Q2_bin;
                break;
              }

            }


          } //end if E/p cut

        }//end loop over N ecal clusers

        if(max_trk > 0) numOfElectrons_RC++;

      } //end if numOfElectrons_RC == 0

      //---------------------------------------------------------------------------------------------------------------------------------------
 
      //loop over ZDC reconstructed neutrals (probably don't want to use this)
			for(int ZDC_neutral_part = 0; ZDC_neutral_part < reco_ZDC_px.GetSize(); ZDC_neutral_part++)
			{
			  TVector3 neutral_hit(reco_ZDC_x[ZDC_neutral_part], reco_ZDC_y[ZDC_neutral_part], reco_ZDC_z[ZDC_neutral_part]);

			  neutral_hit.RotateY(0.025);

			  //------

			  h_ZDC_cluster_x_vs_y->Fill(reco_ZDC_x[ZDC_neutral_part], reco_ZDC_y[ZDC_neutral_part]);
			  h_ZDC_cluster_z->Fill(reco_ZDC_z[ZDC_neutral_part]);

			  h_ZDC_cluster_x_vs_y_no_X_angle->Fill(neutral_hit.X(), neutral_hit.Y());

			  int ZDC_z_bin_event = -1;

			  for(unsigned int ZDC_z_bin = 0; ZDC_z_bin < n_ZDC_z_bins; ZDC_z_bin++)
			  {
			    if( reco_ZDC_z[ZDC_neutral_part] > ZDC_z_bins[ZDC_z_bin] && reco_ZDC_z[ZDC_neutral_part] < ZDC_z_bins[ZDC_z_bin+1] )
			    {
			      ZDC_z_bin_event = ZDC_z_bin;
			    }
			  }

			  if(ZDC_z_bin_event != -1) h_ZDC_cluster_x_vs_y_no_X_angle_z_bins[ZDC_z_bin_event]->Fill(neutral_hit.X(), neutral_hit.Y());


			  //using RC ZDC energy
			  if(reco_ZDC_E[ZDC_neutral_part] < 0) continue;

			  TLorentzVector neutron_4_vect;
			  neutron_4_vect.SetPxPyPzE(reco_ZDC_px[ZDC_neutral_part], reco_ZDC_py[ZDC_neutral_part], reco_ZDC_pz[ZDC_neutral_part], reco_ZDC_E[ZDC_neutral_part]);

			  neutron_4_vect.RotateY(0.025);


				if( numOfNeutrons_RC == 0 ) spectator_RECO = neutron_4_vect;

			  if( numOfNeutrons_RC == 1 )
        {
          if( neutron_4_vect.Eta() > spectator_RECO.Eta() )
          {
            spectator_RECO_2 = spectator_RECO;

            spectator_RECO = neutron_4_vect;
          }
          else
          {
            spectator_RECO_2 = neutron_4_vect;
          }

        }

        alphaSpect_RECO = 2*(spectator_RECO.E() + spectator_RECO.Pz()) / ( deuteron_beam.E() + deuteron_beam.Pz() ); //in head-on lab frame - Jan

				calculatedX_bjorken_nucleon_RC = calculatedX_bjorken_RC/(2.-alphaSpect_RECO);

				//find x bin
        for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
        {
          if( calculatedX_bjorken_nucleon_RC > x_bins[x_bin] && calculatedX_bjorken_nucleon_RC < x_bins[x_bin+1] )
          {
            x_bin_event_RC = x_bin;
            break;
          }

        }

				Flux_RC = 2.*TMath::Pi()*alpha_em*alpha_em*calculated_y_RC*calculated_y_RC / ( calculatedQ2_RC*calculatedQ2_RC*one_minus_epsilon_RC*calculatedX_bjorken_nucleon_RC );

				numOfNeutrons_RC++;

			}


      //----------------------------------------------------------------




			//Use the TRUTH information to calculate all of the proper kinematic variables and factors

			//Apply simple x and Q2 cut here for initial analysis -- will bin and revise later

			iEvent++;
			
			
			//-----------------------
			
			//inclusive deuteron reduced cross section
			//no dependence on spectator, just scattered electron
			//just check that we have a good Q2 bin
			
			if( Q2_bin_event_TRUTH != -1) 			
			{
			  h_N_vs_x[Q2_bin_event_TRUTH]->Fill(calculatedX_bjorken_nucleon_TRUTH);
			  
			  h_sigma_D_red_vs_x[Q2_bin_event_TRUTH]->Fill(calculatedX_bjorken_nucleon_TRUTH, 1./Flux_TRUTH);
			}
			


			if( numOfNeutrons_TRUTH > 0 && numOfNeutrons_TRUTH < 3 ) //throw away events with neutrons != 1, we should have one max. energz
      {
        //---------------

        //scattered electron hists

        h_x_vs_Q2->Fill(calculatedX_bjorken_TRUTH, calculatedQ2_TRUTH);
        h_x_nucleon_vs_Q2->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedQ2_TRUTH);

        h_scat_e_pT->Fill(scatteredElectron_TRUTH.Pt());
        h_scat_e_eta->Fill(scatteredElectron_TRUTH.Eta());
        h_scat_e_phi->Fill(scatteredElectron_TRUTH.Phi());

        //----------------
        //spectator hists

        h_pT->Fill(spectator_TRUTH.Pt());
        h_pT2->Fill(spectator_TRUTH.Pt()*spectator_TRUTH.Pt());

        h_phi->Fill(spectator_TRUTH.Phi());
        h_eta->Fill(spectator_TRUTH.Eta());

        h_alpha_p->Fill(alphaSpect_TRUTH);
        if(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() < 0.01) h_alpha_p_pT2_cut->Fill(alphaSpect_TRUTH);

        //----------------

        if( x_bin_event_TRUTH != -1 && Q2_bin_event_TRUTH != -1)
        {
          if( calculatedX_bjorken_nucleon_TRUTH > x_bins[x_bin_event_TRUTH] && calculatedX_bjorken_nucleon_TRUTH < x_bins[x_bin_event_TRUTH+1] && calculatedQ2_TRUTH > Q2_bins[Q2_bin_event_TRUTH] && calculatedQ2_TRUTH < Q2_bins[Q2_bin_event_TRUTH+1] && alphaSpect_TRUTH > alpha_p_min && alphaSpect_TRUTH < alpha_p_max)
          {
            h_x_vs_Q2_cuts->Fill(calculatedX_bjorken_TRUTH, calculatedQ2_TRUTH);
            h_x_nucleon_vs_Q2_cuts->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedQ2_TRUTH);

            h_scat_e_pT_cuts->Fill(scatteredElectron_TRUTH.Pt());
            h_scat_e_eta_cuts->Fill(scatteredElectron_TRUTH.Eta());
            h_scat_e_phi_cuts->Fill(scatteredElectron_TRUTH.Phi());

            //----------------
            //spectator hists

            h_pT_cuts->Fill(spectator_TRUTH.Pt());
            h_phi_cuts->Fill(spectator_TRUTH.Phi());
            h_eta_cuts->Fill(spectator_TRUTH.Eta());

            //-------------

            if(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() < 0.01)
            {
              h_x_vs_Q2_cuts_pT2_cut->Fill(calculatedX_bjorken_TRUTH, calculatedQ2_TRUTH);
              h_x_nucleon_vs_Q2_cuts_pT2_cut->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedQ2_TRUTH);

              h_scat_e_pT_cuts_pT2_cut->Fill(scatteredElectron_TRUTH.Pt());
              h_scat_e_eta_cuts_pT2_cut->Fill(scatteredElectron_TRUTH.Eta());
              h_scat_e_phi_cuts_pT2_cut->Fill(scatteredElectron_TRUTH.Phi());

              //----------------
              //spectator hists

              h_pT_cuts_pT2_cut->Fill(spectator_TRUTH.Pt());
              h_phi_cuts_pT2_cut->Fill(spectator_TRUTH.Phi());
              h_eta_cuts_pT2_cut->Fill(spectator_TRUTH.Eta());
              
              //-------
              
              h_y_vs_x_vs_Q2_cuts_pT2_cut->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedQ2_TRUTH, calculated_y_TRUTH);
              h_y_vs_x_vs_Q2_cuts_pT2_cut_base->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedQ2_TRUTH);
            }

            //-----------------------------------------------------------------------------------------------

            h_nEvents_bins[x_bin_event_TRUTH][Q2_bin_event_TRUTH]->Fill(1.5); //nEvents_MC is bin 1-2

            //-------------

            h_pT2_bins[x_bin_event_TRUTH][Q2_bin_event_TRUTH]->Fill(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() );

            h_sigma_D_red[x_bin_event_TRUTH][Q2_bin_event_TRUTH]->Fill(spectator_TRUTH.Pt()*spectator_TRUTH.Pt(), alphaSpect_TRUTH/Flux_TRUTH );

            //p reduced cross section
            float R = alphaSpect_TRUTH*alphaSpect_TRUTH*mAvg*capGammaSquared*(2.-alphaSpect_TRUTH);
            float aT2 = mAvg*mAvg - alphaSpect_TRUTH*(2.-alphaSpect_TRUTH)*mDeuteron*mDeuteron/4.;

            float Sd = R/(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() + aT2)/(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() + aT2);

            h_sigma_p_red[x_bin_event_TRUTH][Q2_bin_event_TRUTH]->Fill( spectator_TRUTH.Pt()*spectator_TRUTH.Pt(), alphaSpect_TRUTH/Flux_TRUTH /( two_times_2pi_cubed ) / Sd );

          }
        }

			  iEvent_MC++;

			} //end MC if


			//----------------------------------------------------------------------------------------------------------------------------
			
			if( numOfNeutrons_MC > 0 && numOfNeutrons_MC < 3 ) //throw away events with neutrons != 1, we should have one max. energz
      {
        
        //----------------

        if( x_bin_event_MC != -1 && Q2_bin_event_MC != -1)
        {
          if( calculatedX_bjorken_nucleon_MC > x_bins[x_bin_event_MC] && calculatedX_bjorken_nucleon_MC < x_bins[x_bin_event_MC+1] && calculatedQ2_MC > Q2_bins[Q2_bin_event_MC] && calculatedQ2_MC < Q2_bins[Q2_bin_event_MC+1] && alphaSpect_MC > alpha_p_min && alphaSpect_MC < alpha_p_max)
          {
           

            //-------------

            h_pT2_MCPart_bins[x_bin_event_MC][Q2_bin_event_MC]->Fill(spectator_MC.Pt()*spectator_MC.Pt() );

            h_sigma_D_red_MCPart[x_bin_event_MC][Q2_bin_event_MC]->Fill(spectator_MC.Pt()*spectator_MC.Pt(), alphaSpect_MC/Flux_MC );

            //p reduced cross section
            float R = alphaSpect_MC*alphaSpect_MC*mAvg*capGammaSquared*(2.-alphaSpect_MC);
            float aT2 = mAvg*mAvg - alphaSpect_MC*(2.-alphaSpect_MC)*mDeuteron*mDeuteron/4.;

            float Sd = R/(spectator_MC.Pt()*spectator_MC.Pt() + aT2)/(spectator_MC.Pt()*spectator_MC.Pt() + aT2);

            h_sigma_p_red_MCPart[x_bin_event_MC][Q2_bin_event_MC]->Fill( spectator_MC.Pt()*spectator_MC.Pt(), alphaSpect_MC/Flux_MC /( two_times_2pi_cubed ) / Sd );

          }
        }

			  //iEvent_MC++;
			  //cout << endl;
			} //end MC if


			//----------------------------------------------------------------------------------------------------------------------------
			
			
			//RC analysis
			
			
			//-------------------------------------------------------------------------------------------
			
			//inclusive deuteron reduced cross section
			//no dependence on spectator, just scattered electron
			//just check that we have a good Q2 bin
			
			if( numOfElectrons_RC == 1  )
			{
			  if( Q2_bin_event_RC != -1 ) 
			  {
			    h_N_vs_x_RC[Q2_bin_event_RC]->Fill(calculatedX_bjorken_nucleon_RC);
			    
			    h_sigma_D_red_vs_x_RC[Q2_bin_event_RC]->Fill(calculatedX_bjorken_nucleon_RC, 1./Flux_RC );
			  }
			  
			  //--------
			  
			  if( Q2_bin_event_TRUTH != -1 )
			  {
			    h_N_vs_x_MC_ACC[Q2_bin_event_TRUTH]->Fill(calculatedX_bjorken_nucleon_TRUTH);
			    
			    h_sigma_D_red_vs_x_MC_ACC[Q2_bin_event_TRUTH]->Fill(calculatedX_bjorken_nucleon_TRUTH, 1./Flux_TRUTH );
			  }
			}

			if( numOfNeutrons_RC > 0 && numOfNeutrons_RC < 3 && numOfElectrons_RC == 1 )//throw away events with neutrons != 1, we should have one max. energy
			//if( numOfNeutrons_RC > 0 && numOfNeutrons_RC < 3 && numOfElectrons_RC == 1 && vertex_z.GetSize() == 1  && fabs(vertex_z[0]) < 5 )//test vertex cut
			{
			  //---------------
			  
			  //this is just for theta resolution
			  float pT2_MC = spectator_TRUTH.Pt()*spectator_TRUTH.Pt();
			  
			  int pT2_bin_event = -1;
			  
			  //find pT2 bin for resolution
        for(unsigned int pT2_bin = 0; pT2_bin < n_pT2_bins; pT2_bin++)
        {
          if( pT2_MC > pT2_bins[pT2_bin] && pT2_MC < pT2_bins[pT2_bin+1] )
          {
            pT2_bin_event = pT2_bin;
            break;
          }

        }
        
        //-----------------------------

        //scattered electron hists

        h_x_vs_Q2_RC->Fill(calculatedX_bjorken_RC, calculatedQ2_RC);

        h_x_nucleon_vs_Q2_RC->Fill(calculatedX_bjorken_nucleon_RC, calculatedQ2_RC);

        h_scat_e_pT_RC->Fill(scatteredElectron_RECO.Pt());
        h_scat_e_eta_RC->Fill(scatteredElectron_RECO.Eta());
        h_scat_e_phi_RC->Fill(scatteredElectron_RECO.Phi());

        //MC vs. RC comparison
        h_x_MC_vs_x_RC->Fill(calculatedX_bjorken_TRUTH, calculatedX_bjorken_RC);
        h_x_nucleon_MC_vs_x_nucleon_RC->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedX_bjorken_nucleon_RC);

        h_Q2_MC_vs_Q2_RC->Fill(calculatedQ2_TRUTH, calculatedQ2_RC);

        h_scat_e_pT_MC_vs_e_pT_RC->Fill(scatteredElectron_TRUTH.Pt(), scatteredElectron_RECO.Pt());
        h_scat_e_eta_MC_vs_e_eta_RC->Fill(scatteredElectron_TRUTH.Eta(), scatteredElectron_RECO.Eta());
        h_scat_e_phi_MC_vs_e_phi_RC->Fill(scatteredElectron_TRUTH.Phi(), scatteredElectron_RECO.Phi());


        //----------------
        //spectator hists

        h_pT_RC->Fill(spectator_RECO.Pt());
        h_pT2_RC->Fill(spectator_RECO.Pt()*spectator_RECO.Pt());

        h_phi_RC->Fill(spectator_RECO.Phi());

        //ZDC QA histograms
        int theta_bin_event = -1;
        for( unsigned int theta_bin = 0; theta_bin < n_theta_bins; theta_bin++ )
        {
          float theta_mrad = spectator_TRUTH.Theta()*1000.;

          if( theta_mrad > theta_bins[theta_bin] && theta_mrad < theta_bins[theta_bin+1] )
          {
            theta_bin_event = theta_bin;
          }

        }

        if( theta_bin_event != -1 )
        {
          h_phi_RC_theta_bins[theta_bin_event]->Fill(spectator_RECO.Phi());
        }



        h_eta_RC->Fill(spectator_RECO.Eta());

        h_alpha_p_RC->Fill(alphaSpect_RECO);
        if(spectator_RECO.Pt()*spectator_RECO.Pt() < 0.01) h_alpha_p_RC_pT2_cut->Fill(alphaSpect_RECO);

        //MC vs. RC comparison
        h_pT_MC_vs_pT_RC->Fill(spectator_TRUTH.Pt(), spectator_RECO.Pt());
        h_phi_MC_vs_phi_RC->Fill(spectator_TRUTH.Phi(), spectator_RECO.Phi());
        h_eta_MC_vs_eta_RC->Fill(spectator_TRUTH.Eta(), spectator_RECO.Eta());

        h_E_MC_vs_E_RC->Fill(spectator_TRUTH.E(), spectator_RECO.E());

        h_alpha_p_MC_vs_alpha_p_RC->Fill(alphaSpect_TRUTH, alphaSpect_RECO);
        if(spectator_RECO.Pt()*spectator_RECO.Pt() < 0.01) h_alpha_p_MC_vs_alpha_p_RC_pT2_cut->Fill(alphaSpect_TRUTH, alphaSpect_RECO);

        //----------------
        //Resolution plost
        h_theta_resolution->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.); //converted to mrad
        
        if( spectator_TRUTH.Theta()*1000 < 1 && spectator_RECO.Theta()*1000 < 1 ) h_theta_resolution_theta_1->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.);
        if( spectator_TRUTH.Theta()*1000 > 1 && spectator_RECO.Theta()*1000 > 1 ) h_theta_resolution_theta_2->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.);
        
        //h_R_resolution

        h_E_resolution->Fill(spectator_TRUTH.E() - spectator_RECO.E());

        //----------------

        if( alphaSpect_RECO > alpha_p_min && alphaSpect_RECO < alpha_p_max  )
        {
          if( theta_bin_event != -1 )
          {
            h_phi_RC_theta_bins_alpha[theta_bin_event]->Fill(spectator_RECO.Phi());
          }
          
          //h_theta_resolution_alpha->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.); //converted to mrad
          h_theta_resolution_alpha->Fill((spectator_MC.Theta() - spectator_RECO.Theta())*1000.); //converted to mrad
          
          if( spectator_TRUTH.Theta()*1000 < 1 && spectator_RECO.Theta()*1000 < 1 ) h_theta_resolution_alpha_theta_1->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.);
          if( spectator_TRUTH.Theta()*1000 > 1 && spectator_RECO.Theta()*1000 > 1 ) h_theta_resolution_alpha_theta_2->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.);
        
          if( x_bin_event_RC != -1 && Q2_bin_event_RC != -1)
          {
            if( calculatedX_bjorken_nucleon_RC > x_bins[x_bin_event_RC] && calculatedX_bjorken_nucleon_RC < x_bins[x_bin_event_RC+1] && calculatedQ2_RC > Q2_bins[Q2_bin_event_RC] && calculatedQ2_RC < Q2_bins[Q2_bin_event_RC+1] )
            {

              //scattered electron hists

              h_x_vs_Q2_RC_cuts->Fill(calculatedX_bjorken_RC, calculatedQ2_RC);

              h_x_nucleon_vs_Q2_RC_cuts->Fill(calculatedX_bjorken_nucleon_RC, calculatedQ2_RC);

              h_scat_e_pT_RC_cuts->Fill(scatteredElectron_RECO.Pt());
              h_scat_e_eta_RC_cuts->Fill(scatteredElectron_RECO.Eta());
              h_scat_e_phi_RC_cuts->Fill(scatteredElectron_RECO.Phi());

              //MC vs. RC comparison
              h_x_MC_vs_x_RC_cuts->Fill(calculatedX_bjorken_TRUTH, calculatedX_bjorken_RC);
              h_x_nucleon_MC_vs_x_nucleon_RC_cuts->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedX_bjorken_nucleon_RC);

              h_Q2_MC_vs_Q2_RC_cuts->Fill(calculatedQ2_TRUTH, calculatedQ2_RC);

              h_scat_e_pT_MC_vs_e_pT_RC_cuts->Fill(scatteredElectron_TRUTH.Pt(), scatteredElectron_RECO.Pt());
              h_scat_e_eta_MC_vs_e_eta_RC_cuts->Fill(scatteredElectron_TRUTH.Eta(), scatteredElectron_RECO.Eta());
              h_scat_e_phi_MC_vs_e_phi_RC_cuts->Fill(scatteredElectron_TRUTH.Phi(), scatteredElectron_RECO.Phi());

              h_E_MC_vs_E_RC_cuts->Fill(spectator_TRUTH.E(), spectator_RECO.E());

              //----------------
              //spectator hists

              h_pT_RC_cuts->Fill(spectator_RECO.Pt());
              h_phi_RC_cuts->Fill(spectator_RECO.Phi());
              h_eta_RC_cuts->Fill(spectator_RECO.Eta());

              //MC vs. RC comparison
              h_pT_MC_vs_pT_RC_cuts->Fill(spectator_TRUTH.Pt(), spectator_RECO.Pt());
              h_phi_MC_vs_phi_RC_cuts->Fill(spectator_TRUTH.Phi(), spectator_RECO.Phi());
              h_eta_MC_vs_eta_RC_cuts->Fill(spectator_TRUTH.Eta(), spectator_RECO.Eta());

              if(spectator_RECO.Pt()*spectator_RECO.Pt() < 0.01)
              {
                //scattered electron hists

                h_x_vs_Q2_RC_cuts_pT2_cut->Fill(calculatedX_bjorken_RC, calculatedQ2_RC);

                h_x_nucleon_vs_Q2_RC_cuts_pT2_cut->Fill(calculatedX_bjorken_nucleon_RC, calculatedQ2_RC);

                h_scat_e_pT_RC_cuts_pT2_cut->Fill(scatteredElectron_RECO.Pt());
                h_scat_e_eta_RC_cuts_pT2_cut->Fill(scatteredElectron_RECO.Eta());
                h_scat_e_phi_RC_cuts_pT2_cut->Fill(scatteredElectron_RECO.Phi());

                //MC vs. RC comparison
                h_x_MC_vs_x_RC_cuts_pT2_cut->Fill(calculatedX_bjorken_TRUTH, calculatedX_bjorken_RC);
                h_x_nucleon_MC_vs_x_nucleon_RC_cuts_pT2_cut->Fill(calculatedX_bjorken_nucleon_TRUTH, calculatedX_bjorken_nucleon_RC);

                h_Q2_MC_vs_Q2_RC_cuts_pT2_cut->Fill(calculatedQ2_TRUTH, calculatedQ2_RC);

                h_scat_e_pT_MC_vs_e_pT_RC_cuts_pT2_cut->Fill(scatteredElectron_TRUTH.Pt(), scatteredElectron_RECO.Pt());
                h_scat_e_eta_MC_vs_e_eta_RC_cuts_pT2_cut->Fill(scatteredElectron_TRUTH.Eta(), scatteredElectron_RECO.Eta());
                h_scat_e_phi_MC_vs_e_phi_RC_cuts_pT2_cut->Fill(scatteredElectron_TRUTH.Phi(), scatteredElectron_RECO.Phi());

                //----------------
                //spectator hists
                h_pT_RC_cuts_pT2_cut->Fill(spectator_RECO.Pt());
                h_phi_RC_cuts_pT2_cut->Fill(spectator_RECO.Phi());
                h_eta_RC_cuts_pT2_cut->Fill(spectator_RECO.Eta());

                //MC vs. RC comparison
                h_pT_MC_vs_pT_RC_cuts_pT2_cut->Fill(spectator_TRUTH.Pt(), spectator_RECO.Pt());
                h_phi_MC_vs_phi_RC_cuts_pT2_cut->Fill(spectator_TRUTH.Phi(), spectator_RECO.Phi());
                h_eta_MC_vs_eta_RC_cuts_pT2_cut->Fill(spectator_TRUTH.Eta(), spectator_RECO.Eta());

                //----------------
                
                h_y_vs_x_vs_Q2_RC_cuts_pT2_cut->Fill(calculatedX_bjorken_nucleon_RC, calculatedQ2_RC, calculated_y_RC);
                h_y_vs_x_vs_Q2_RC_cuts_pT2_cut_base->Fill(calculatedX_bjorken_nucleon_RC, calculatedQ2_RC);
                
                //----------------
                
                //Resolution plost
                h_theta_resolution_cut->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.); //converted to mrad
                //h_R_resolution
                
                if( spectator_TRUTH.Theta()*1000 < 1 && spectator_RECO.Theta()*1000 < 1 ) h_theta_resolution_cut_theta_1->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.);
                if( spectator_TRUTH.Theta()*1000 > 1 && spectator_RECO.Theta()*1000 > 1 ) h_theta_resolution_cut_theta_2->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.);
                
                if(pT2_bin_event != -1)
                {
                  h_theta_resolution_cut_pT2_bins[pT2_bin_event]->Fill((spectator_TRUTH.Theta() - spectator_RECO.Theta())*1000.); //converted to mrad
                  h_theta_resolution_cut_pT2_bins_MC[pT2_bin_event]->Fill((spectator_MC.Theta() - spectator_RECO.Theta())*1000.); //converted to mrad
                  h_theta_resolution_cut_pT2_bins_beam_eff[pT2_bin_event]->Fill((spectator_MC.Theta() - spectator_TRUTH.Theta())*1000.); //converted to mrad
                }
                
                
                h_E_resolution_cut->Fill(spectator_TRUTH.E() - spectator_RECO.E());

                //----------------
                
                if( theta_bin_event != -1 )
                {
                  h_phi_RC_theta_bins_cuts[theta_bin_event]->Fill(spectator_RECO.Phi());
                }

              }

              //---------------------------------------------------------------------------------------------------------------------------------------

              h_nEvents_bins[x_bin_event_RC][Q2_bin_event_RC]->Fill(2.5); //nEvents_RC is bin 2-3

              h_pT2_RC_bins[x_bin_event_RC][Q2_bin_event_RC]->Fill(spectator_RECO.Pt()*spectator_RECO.Pt() );

              h_sigma_D_red_RC[x_bin_event_RC][Q2_bin_event_RC]->Fill(spectator_RECO.Pt()*spectator_RECO.Pt(), alphaSpect_RECO/Flux_RC );

              //n reduced cross section
              float R = alphaSpect_RECO*alphaSpect_RECO*mAvg*capGammaSquared*(2.-alphaSpect_RECO); //Eq. (45) in paper (no factor of 2)
              float aT2 = mAvg*mAvg - alphaSpect_RECO*(2.-alphaSpect_RECO)*mDeuteron*mDeuteron/4.; //this version is from Alex's slides

              float Sd = R/(spectator_RECO.Pt()*spectator_RECO.Pt() + aT2)/(spectator_RECO.Pt()*spectator_RECO.Pt() + aT2); //Eq. (44) in paper

              h_sigma_p_red_RC[x_bin_event_RC][Q2_bin_event_RC]->Fill( spectator_RECO.Pt()*spectator_RECO.Pt(), alphaSpect_RECO/Flux_RC /( two_times_2pi_cubed ) / Sd );

            }
          }

          //--------------------------------------------------------

          //this is for acceptance only histograms
          //contain MC info, but are only for RC accepted events
          if( x_bin_event_TRUTH != -1 && Q2_bin_event_TRUTH != -1 && numOfNeutrons_TRUTH > 0 && numOfNeutrons_TRUTH < 3)
          {
            if( calculatedX_bjorken_nucleon_TRUTH > x_bins[x_bin_event_TRUTH] && calculatedX_bjorken_nucleon_TRUTH < x_bins[x_bin_event_TRUTH+1] && calculatedQ2_TRUTH > Q2_bins[Q2_bin_event_TRUTH] && calculatedQ2_TRUTH < Q2_bins[Q2_bin_event_TRUTH+1] )
            //if( calculatedX_bjorken_nucleon > x_bins[x_bin_event_TRUTH] && calculatedX_bjorken_nucleon < x_bins[x_bin_event+1] && calculatedQ2_TRUTH > Q2_bins[Q2_bin_event_TRUTH] && calculatedQ2_TRUTH < Q2_bins[Q2_bin_event+1] && alphaSpect_TRUTH > alpha_p_min && alphaSpect_TRUTH < alpha_p_max )
            {

              h_pT2_MC_ACC_bins[x_bin_event_TRUTH][Q2_bin_event_TRUTH]->Fill(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() );

              h_sigma_D_red_MC_ACC[x_bin_event_TRUTH][Q2_bin_event_TRUTH]->Fill(spectator_TRUTH.Pt()*spectator_TRUTH.Pt(), alphaSpect_TRUTH/Flux_TRUTH );

              //n reduced cross section
              float R = alphaSpect_TRUTH*alphaSpect_TRUTH*mAvg*capGammaSquared*(2.-alphaSpect_TRUTH); //Eq. (45) in paper (no factor of 2)
              float aT2 = mAvg*mAvg - alphaSpect_TRUTH*(2.-alphaSpect_TRUTH)*mDeuteron*mDeuteron/4.; //this version is from Alex's slides
              //float aT2 = mAvg*mAvg - mDeuteron*mDeuteron/4.; //this version is from Alex's slides

              float Sd = R/(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() + aT2)/(spectator_TRUTH.Pt()*spectator_TRUTH.Pt() + aT2); //Eq. (44) in paper

              h_sigma_p_red_MC_ACC[x_bin_event_TRUTH][Q2_bin_event_TRUTH]->Fill( spectator_TRUTH.Pt()*spectator_TRUTH.Pt(), alphaSpect_TRUTH/Flux_TRUTH /( two_times_2pi_cubed ) / Sd );

            }

          }

        }//end if alpha RC cut


			  iEvent_RC++;
			  //cout << endl;
			}


			//reset variables
	    numOfElectrons_TRUTH = 0;
	    numOfNeutrons_TRUTH  = 0;
	    
	    numOfElectrons_MC = 0;
	    numOfNeutrons_MC  = 0;

	    numOfElectrons_RC = 0;
		  numOfNeutrons_RC   = 0;




		}// event loop

		inputRootFile->Close();

	}// input file loop

	//store and print number of events
	cout<<"Total events = "<<iEvent<<endl;
	cout<<"Total accepted MC events = "<<iEvent_MC<<endl;
	cout<<"Total accepted RC events = "<<iEvent_RC<<endl;

	h_nEvents->SetBinContent(1, iEvent);
	h_nEvents->SetBinContent(2, iEvent_MC);
	h_nEvents->SetBinContent(3, iEvent_RC);

  //------------------------------------------------------------------

  //save histograms to a file

	TFile * outputFile = new TFile(outputFileName, "RECREATE");

  h_nEvents->Write();

  for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
  {
    for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
    {
      h_nEvents_bins[x_bin][Q2_bin]->Write();
    }
  }

  //scattered electron
  h_x_vs_Q2->Write();
  h_x_nucleon_vs_Q2->Write();

  h_scat_e_pT->Write();
  h_scat_e_eta->Write();
  h_scat_e_phi->Write();



  //spectator neutron
  h_pT->Write();
  h_pT2->Write();

  h_eta->Write();
  h_phi->Write();

  h_alpha_p->Write();
  h_alpha_p_pT2_cut->Write();

  //------------------

	h_x_vs_Q2_cuts->Write();
  h_x_nucleon_vs_Q2_cuts->Write();

  h_scat_e_pT_cuts->Write();
  h_scat_e_eta_cuts->Write();
  h_scat_e_phi_cuts->Write();

  //spectator proton
  h_pT_cuts->Write();
  h_eta_cuts->Write();
  h_phi_cuts->Write();

  //------------------

	h_x_vs_Q2_cuts_pT2_cut->Write();
  h_x_nucleon_vs_Q2_cuts_pT2_cut->Write();

  h_scat_e_pT_cuts_pT2_cut->Write();
  h_scat_e_eta_cuts_pT2_cut->Write();
  h_scat_e_phi_cuts_pT2_cut->Write();

  //spectator proton
  h_pT_cuts_pT2_cut->Write();
  h_eta_cuts_pT2_cut->Write();
  h_phi_cuts_pT2_cut->Write();
  
  
  //---
  
  h_y_vs_x_vs_Q2_cuts_pT2_cut->Divide(h_y_vs_x_vs_Q2_cuts_pT2_cut_base);
  h_y_vs_x_vs_Q2_cuts_pT2_cut->Write();


	//------------------
	//true MC reduced cross section


  for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
  {
    for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
    {
      h_pT2_bins[x_bin][Q2_bin]->Write();
      
      if(x_bin == 0) 
      {
        h_N_vs_x[Q2_bin]->Write();
        
        h_sigma_D_red_vs_x[Q2_bin]->Write();      
      }

      h_sigma_D_red[x_bin][Q2_bin]->Write();
      h_sigma_p_red[x_bin][Q2_bin]->Write();
      
      
    }
  }

  //--------------------------------------------------------------------
  
  //------------------
	//MC reduced cross section


  for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
  {
    for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
    {
      h_pT2_MCPart_bins[x_bin][Q2_bin]->Write();
          
      h_sigma_D_red_MCPart[x_bin][Q2_bin]->Write();
      h_sigma_p_red_MCPart[x_bin][Q2_bin]->Write();      
      
    }
  }

  //--------------------------------------------------------------------

  //RC information
  //scattered electron
  h_x_vs_Q2_RC->Write();
  h_x_nucleon_vs_Q2_RC->Write();

  h_track_eta_RC->Write();
  h_nECAL_cluster_eta_RC->Write();
  h_bECAL_cluster_eta_RC->Write();

  h_scat_e_pT_RC->Write();
  h_scat_e_eta_RC->Write();
  h_scat_e_phi_RC->Write();

  //MC vs. RC comparison
  h_x_MC_vs_x_RC->Write();
  h_x_nucleon_MC_vs_x_nucleon_RC->Write();

  h_Q2_MC_vs_Q2_RC->Write();

  h_scat_e_pT_MC_vs_e_pT_RC->Write();
  h_scat_e_eta_MC_vs_e_eta_RC->Write();
  h_scat_e_phi_MC_vs_e_phi_RC->Write();

  //spectator neutron
  h_pT_RC->Write();
  h_pT2_RC->Write();

  h_eta_RC->Write();
  h_phi_RC->Write();

  h_alpha_p_RC->Write();
  h_alpha_p_RC_pT2_cut->Write();

  //MC vs. RC comparison
  h_pT_MC_vs_pT_RC->Write();
  h_eta_MC_vs_eta_RC->Write();
  h_phi_MC_vs_phi_RC->Write();

  h_E_MC_vs_E_RC->Write();

  h_alpha_p_MC_vs_alpha_p_RC->Write();
  h_alpha_p_MC_vs_alpha_p_RC_pT2_cut->Write();

  //------------------

	h_x_vs_Q2_RC_cuts->Write();
  h_x_nucleon_vs_Q2_RC_cuts->Write();

  h_scat_e_pT_RC_cuts->Write();
  h_scat_e_eta_RC_cuts->Write();
  h_scat_e_phi_RC_cuts->Write();

  //MC vs. RC comparison
  h_x_MC_vs_x_RC_cuts->Write();
  h_x_nucleon_MC_vs_x_nucleon_RC_cuts->Write();

  h_Q2_MC_vs_Q2_RC_cuts->Write();

  h_scat_e_pT_MC_vs_e_pT_RC_cuts->Write();
  h_scat_e_eta_MC_vs_e_eta_RC_cuts->Write();
  h_scat_e_phi_MC_vs_e_phi_RC_cuts->Write();

  //spectator proton
  h_pT_RC_cuts->Write();
  h_eta_RC_cuts->Write();
  h_phi_RC_cuts->Write();

  //MC vs. RC comparison
  h_pT_MC_vs_pT_RC_cuts->Write();
  h_eta_MC_vs_eta_RC_cuts->Write();
  h_phi_MC_vs_phi_RC_cuts->Write();

  h_E_MC_vs_E_RC_cuts->Write();

  //------------------

	h_x_vs_Q2_RC_cuts_pT2_cut->Write();
  h_x_nucleon_vs_Q2_RC_cuts_pT2_cut->Write();

  h_scat_e_pT_RC_cuts_pT2_cut->Write();
  h_scat_e_eta_RC_cuts_pT2_cut->Write();
  h_scat_e_phi_RC_cuts_pT2_cut->Write();

  //MC vs. RC comparison
  h_x_MC_vs_x_RC_cuts_pT2_cut->Write();
  h_x_nucleon_MC_vs_x_nucleon_RC_cuts_pT2_cut->Write();

  h_Q2_MC_vs_Q2_RC_cuts_pT2_cut->Write();

  h_scat_e_pT_MC_vs_e_pT_RC_cuts_pT2_cut->Write();
  h_scat_e_eta_MC_vs_e_eta_RC_cuts_pT2_cut->Write();
  h_scat_e_phi_MC_vs_e_phi_RC_cuts_pT2_cut->Write();

  //spectator proton
  h_pT_RC_cuts_pT2_cut->Write();
  h_eta_RC_cuts_pT2_cut->Write();
  h_phi_RC_cuts_pT2_cut->Write();
  
  //---
  
  h_y_vs_x_vs_Q2_RC_cuts_pT2_cut->Divide(h_y_vs_x_vs_Q2_RC_cuts_pT2_cut_base);
  h_y_vs_x_vs_Q2_RC_cuts_pT2_cut->Write();

  //MC vs. RC comparison
  h_pT_MC_vs_pT_RC_cuts_pT2_cut->Write();
  h_eta_MC_vs_eta_RC_cuts_pT2_cut->Write();
  h_phi_MC_vs_phi_RC_cuts_pT2_cut->Write();

	//------------------

	h_theta_resolution->Write();
	h_theta_resolution_theta_1->Write();
	h_theta_resolution_theta_2->Write();
	
	h_theta_resolution_alpha->Write();
	h_theta_resolution_alpha_theta_1->Write();
	h_theta_resolution_alpha_theta_2->Write();
	
	h_theta_resolution_cut->Write();
	h_theta_resolution_cut_theta_1->Write();
	h_theta_resolution_cut_theta_2->Write();

	for(unsigned int i_pT2_bin = 0; i_pT2_bin < n_pT2_bins; i_pT2_bin++)
  {
    h_theta_resolution_cut_pT2_bins[i_pT2_bin]->Write();
    
    h_theta_resolution_cut_pT2_bins_MC[i_pT2_bin]->Write();
    
    h_theta_resolution_cut_pT2_bins_beam_eff[i_pT2_bin]->Write();
  }

	h_E_resolution->Write();
	h_E_resolution_cut->Write();

	//------------------

	h_ZDC_cluster_x_vs_y->Write();
	h_ZDC_cluster_z->Write();

	h_ZDC_cluster_x_vs_y_no_X_angle->Write();

	for(unsigned int ZDC_z_bin = 0; ZDC_z_bin < n_ZDC_z_bins; ZDC_z_bin++)
  {
    h_ZDC_cluster_x_vs_y_no_X_angle_z_bins[ZDC_z_bin]->Write();
  }

	for(unsigned int theta_bin = 0; theta_bin < n_theta_bins; theta_bin++)
  {
    h_phi_RC_theta_bins[theta_bin]->Write();
    h_phi_RC_theta_bins_alpha[theta_bin]->Write();
    h_phi_RC_theta_bins_cuts[theta_bin]->Write();
  }

	//MC reduced cross section


  for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
  {
    for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
    {
      h_pT2_RC_bins[x_bin][Q2_bin]->Write();
      
      if(x_bin == 0)
      {
        h_N_vs_x_RC[Q2_bin]->Write();
        
        h_sigma_D_red_vs_x_RC[Q2_bin]->Write();
      }

      h_sigma_D_red_RC[x_bin][Q2_bin]->Write();
      h_sigma_p_red_RC[x_bin][Q2_bin]->Write();
    }
  }

  //------------------------------------------------------
  //acceptance only histograms
  for(unsigned int x_bin = 0; x_bin < nx_bins; x_bin++)
  {
    for(unsigned int Q2_bin = 0; Q2_bin < nQ2_bins; Q2_bin++)
    {
      h_pT2_MC_ACC_bins[x_bin][Q2_bin]->Write();
      
      if(x_bin == 0)
      {
        h_N_vs_x_MC_ACC[Q2_bin]->Write();
        
        h_sigma_D_red_vs_x_MC_ACC[Q2_bin]->Write();
      }

      h_sigma_D_red_MC_ACC[x_bin][Q2_bin]->Write();
      h_sigma_p_red_MC_ACC[x_bin][Q2_bin]->Write();
    }
  }

	outputFile->Close();

	//--------------------------------------------------------------------

  return;

}


