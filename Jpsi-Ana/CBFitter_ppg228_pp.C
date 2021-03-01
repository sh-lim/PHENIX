// To run this code, enter bin number "1 - 24" and binwidth "1" for all the 0.25 GeV/c bins 0.0 - 6.0
// Then enter bin "13" and binwidth "2" for 0.5 pt bin 6.0 - 6.5 GeV/c
// again, enter  bin "14" and binwidth "2" for 0.5 pt bin 6.5 - 7.0 GeV/c

// For p + p

#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <Math/MinimizerOptions.h>
#include <fstream>
#include <iostream>
#include <TAttText.h>

using namespace std;








TH1 *recomass;
TH1 *recomassFG;
TH1 *recomassBG;
TH1 *recomassBG1;
TH1 *recomassBG2;
TH1 *recomassBG3;
TF1 *combg_fit;
char nameN[500];
char nameS[500];
char nameN2[500];
char nameS2[500];


Double_t CBcalc_jpsi(Double_t *xx, Double_t *par)  // CB fucntion for J/psi 
{
  double f;
  double x = xx[0];

  // The four CB parameters (alpha, n, x_mean, sigma) plus normalization (N) are:

  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];

  double N_2nd = par[5];
  double X_mean_2nd = x_mean;
  double Sigma_2nd = sigma;


  int s=0;

  const double sigma_step = 0.001;
  double fixed = 0.0 + s*sigma_step;


  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_jpsi = N_2nd  * exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));
 
  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_jpsi;
    }
  else
    {
      f = N  * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_jpsi;
    }

  if(x == 3.20)
    cout << "check 1: " << N_2nd << " , " << X_mean_2nd << " , " << Sigma_2nd << " , " << endl;

  return f;
}

Double_t CBcalc_psi2s(Double_t *xx, Double_t *par)  // CB fucntion for J/psi 
{
  double f;
  double x = xx[0];

  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];

  double n_2nd = par[5];
  double x_mean_2nd = x_mean;
  double sigma_2nd = sigma;

  int s=0;

  const double sigma_step = 0.001;
  double fixed = 0.0 + s*sigma_step;

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_psi2s = n_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2)));

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_psi2s;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_psi2s;
    }

  if(x == 3.20)
    cout << "check 2: " << n_2nd << " , " << x_mean_2nd << " , " << sigma_2nd << " , " << endl;


  return f;
}


Double_t CBcalc_jpsi_G(Double_t *xx, Double_t *par)  // CB fucntion for J/psi 
{
  double f;
  double x = xx[0];

  // The four CB parameters (alpha, n, x_mean, sigma) plus normalization (N) are:

  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];

  double N_2nd = par[5];
  double X_mean_2nd = par[6];
  double Sigma_2nd = par[7];

  int s=0;

  const double sigma_step = 0.001;
  double fixed = 0.0 + s*sigma_step;


  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_jpsi = N_2nd * exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));
 
  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_jpsi;
    }
  else
    {
      f = N  * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_jpsi;
    }


 if(x == 2.60)
   cout << "check 3: " << N_2nd << " , " << X_mean_2nd << " , " << Sigma_2nd << " , " << endl;

  return f;
}


Double_t CBcalc_psi2s_G(Double_t *xx, Double_t *par)  // CB fucntion for psi 2s
{
  double f;
  double x = xx[0];

  // The four CB parameters (alpha, n, x_mean, sigma) plus normalization (N) are:

  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];

  double n_2nd = par[5];
  double x_mean_2nd = par[6];
  double sigma_2nd = par[7];

  int s=0;

  const double sigma_step = 0.001;
  double fixed = 0.0 + s*sigma_step;


  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_psi2s = n_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2)));

  if( (x-x_mean)/sigma > -alpha)
    {
      f = N  * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_psi2s;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_psi2s;
    }


 if(x == 2.60)
   cout << "check 4: " << n_2nd << " , " << x_mean_2nd << " , " << sigma_2nd << " , " << endl;

  return f;
}

///the second gaussian alone to draw fit for J/psi and psi2s
Double_t Gauss2(Double_t *xx, Double_t *par) 
{
  double f;
  double x = xx[0];

  double N_2nd = par[0];
  double x_mean_2nd = par[1];
  double sigma_2nd = par[2];

  int s=0;

  const double sigma_step = 0.001;
  double fixed = 0.0 + s*sigma_step;


  f = N_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2))); 


 if(x == 2.60)
   cout << "check 5: " << N_2nd << " , " << x_mean_2nd << " , " << sigma_2nd << " , " << endl;

  return f;

}

Double_t CBcalc_LL(Double_t *xx, Double_t *par) // total fit function
{
  double f;
  double x = xx[0];
 

  int s=0;

  const double sigma_step = 0.001;
  double fixed = 0.0 + s*sigma_step;


  // Jpsi CB parameters
 
  double alpha = par[0];
  double n = par[1];
  double x_mean = par[2];
  double sigma = par[3];
  double N = par[4];
  
 // the second gaussian for J/psi peak
  double N_2nd = par[5];
  double X_mean_2nd = x_mean;
  double Sigma_2nd = fixed;
  
  // psi2s CB parameters (alpha and n equal J/psi values)
  double mass_offset = par[6];
  double width_factor = par[7];
  double N_psi2s = par[8];
  
  // YH correlated bg parameters
  double par_a = par[9];
  double par_b = par[10];
  double par_c = par[11];
  double par_d = par[12];
  double par_e = par[13];

// the second gaussian for psi2s peak
  double n_2nd = N_2nd * N_psi2s / N;
  double x_mean_2nd = x_mean + mass_offset;// + 0.589*x_mean/3.0969;
   double sigma_2nd = fixed;
  
  // The other psi2s parameters will be constrained by the jpsi parameters, only the normalization will vary
  double x_mean_psi2s = x_mean + mass_offset;
  double sigma_psi2s = sigma*width_factor;

  double A = pow( (n/TMath::Abs(alpha)),n) * exp(-pow(alpha,2)/2.0);
  double B = n/TMath::Abs(alpha) - TMath::Abs(alpha);

  double gauss2_jpsi = N_2nd * exp( -pow( x - X_mean_2nd, 2)/ (2.0*pow(fixed,2)));
  double gauss2_psi2s = n_2nd * exp( -pow( x - x_mean_2nd, 2)/ (2.0*pow(fixed,2)));


  // for J/psi CB function
  if( (x-x_mean)/sigma > -alpha)
    {
      f = N * exp( -pow(x-x_mean,2) / (2.0*pow(sigma,2))) + gauss2_jpsi;
    }
  else
    {
      f = N * A * pow(B - (x-x_mean)/sigma, -n) + gauss2_jpsi;
    }
  // for psi(2S) CB function
  if( (x-x_mean_psi2s)/sigma_psi2s > -alpha)
    {
      f += N_psi2s * exp( -pow(x-x_mean_psi2s,2) / (2.0*pow(sigma_psi2s,2))) + gauss2_psi2s;
    }
  else
    {
      f += N_psi2s * A * pow(B - (x-x_mean_psi2s)/sigma_psi2s, -n) + gauss2_psi2s;
    }
  // add correlated backgroudn to CB functions
  f += par_c /( pow( exp(-par_a*x - par_b*x*x) + x/par_d, par_e)); // NEW corrbg formula

  //Like sign combinatoric background fit
  double bg12 = combg_fit->Eval(x);

  // Add total fit with like sign background fit
  f += bg12;  // comb bg fit results are added to jpsi CB, psi(2s)CB and corr bg 

   if(x == 2.60)
    {
      cout << "==================== " << endl;
      cout << "J/psi 2nd gaussian parameters:" << endl;
      cout << "==================== " << endl;
      cout <<" N_2nd : " << N_2nd << endl;
      cout << " X_mean_2nd : " << X_mean_2nd << endl;
      cout << " Sigma 2nd : " << sigma << endl;
      cout << " N CB : " << N << endl;
      cout << "==================== " << endl;
      cout << "psi2S 2nd gaussian parameters:" << endl;
      cout << "==================== " << endl;
      cout <<" n_2nd : " << n_2nd << endl;
      cout << " x_mean_2nd : " << x_mean_2nd << endl;
      cout << " sigma 2nd : " << sigma_psi2s << endl;
      cout << " n CB : " << N_psi2s << endl;
      cout << "==================== " << endl;

      if( N_2nd/N == n_2nd/N_psi2s)
	cout << " --> 1. normalization ratios are equal " << endl;
      if(  N_2nd/N != n_2nd/N_psi2s)
	cout << " --> 1. f: " << N_2nd / N << " and f': " << n_2nd / N_psi2s << " " << endl;

      if(Sigma_2nd == sigma_2nd)
	cout << " --> 2. gaussian widths both fixed to " << fixed << " " << endl;

      if(X_mean_2nd == x_mean)
	cout << " --> 3. 2nd gaussian mean equal to Jpsi CB mean " << endl;

      if(x_mean_2nd == x_mean_psi2s)
	cout << " --> 4. 2nd gaussian mean equal to psi2S CB mean " << endl;

    }

  return f;
}

void CBFitter_ppg228_pp()
{


 // for 500 mev function descriptions of parameters a,b,d,e
#include "fit_coeff_S_array_twelve_point.C" 
#include "fit_coeff_N_array_twelve_point.C" 



double mass_array_corrbg[100] = {0};
double mass_array_combg[100] = {0};
double mass_array_combg_err[100] = {0};
double mass_array_jpsi[100] = {0};
double mass_array_jpsi_err[100] = {0};
double mass_array_psi2s[100] = {0};
double mass_array_psi2s_err[100] = {0};
double mass_array_total[100] = {0};
double mass_array_total_err[100] = {0};

double x_array[100] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95};

 double bg_array[60] = {2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95};

  bool fixed_corr_bg = false; // set to true for bin 26 on North Arm
  bool calculated_fit = false; // set both to false for min bias fit using parameter f(x) from yuehang_pt_integrated macro

  bool fx_pt_int = true;  // set to true to fit pT/rapidity integrated (bin == 49 in north arm and bin 29 in South arm)
  bool width_write = false; 

  bool write_pAl_syst_err = false;

  ////////////////////////////////////////////
  // bool north_arm = true;
  bool north_arm = false;
  ////////////////////////////////////////

  int s=0;

  

  const double sigma_step = 0.001;
  double fixed = 0.0 + s*sigma_step;


  bool chisquare_write = false; 

  bool sngtrk = false;
  bool fvtx = true;
  bool fvtx_sngtrk = false;

  bool print_screen = false;

  bool fixed_width = false;
  bool fixed_center_of_peak = false;
 
  bool test_write = false; // write for uncertainty study 
  bool test = false;

  bool pAu = false; // for 500 MeV (4-4.5 & 4.5-5) and 2 GeV (5-7) binwidths in RAA plots
  bool HeAu = false;
  bool pAl = false;

  bool psi_2s = false;

  int bin;
  int testbin;
  double chisquare;
  double ndf;

  // cout << "Enter the bin number " << endl;
  // cin >> bin;
  //bin = 49;
  // bin = 13;
  //bin = 48;
  bin = 1;

 int bg_fit;
 // cout << "BG fit? (yes = 1 or no = 0) " << endl;
 //  cin >> bg_fit;
 bg_fit = 1;

  int special = 0;
  //int special = 1; // use for 5-7 pt bin only

  int case_type;

  //cout << "Enter Case A, F, G or H ('0', '1', '2' or '3')" << endl;
  // cin >> case_type;
  case_type = 0;
  double pt_bin = 0;
  //double pt_bin = 2; // use for 5-7 pt bin only

   int bin_width = 1;
  //  int bin_width = 2; // use for 500 MeV or  2 GeV binwidths
  // int bin_width = 4;
   //  int bin_width = 8;

  int fit = 1;

  int com_bg;
  cout << "Enter combinatoric BG initial parameter set (0-13) " << endl;
  cin >> com_bg;
  // if(north_arm)
  //   com_bg = 7;
  //   else
  //    com_bg = 11;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  cout << "start" << endl;
  bool LL_fit = true;

  double maxmass = 5.0;

  TFile *file1S;  


  /// comparison for PPG188 and PPG228 
  
  if(north_arm)
	{     
	  if(special == 0)
	    {
	      if(bin_width == 4) 
		{
		  //  file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_PPG188.root");   // pT int is bin 13 in this file
			//	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_PPG188_pro108.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_nofvtx_dimu_02.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_sngtrk_dimu_02.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_mass_fvtx_dimu02_noprob_int.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_fvtx_sngtrk_dimu_02.root");
		  //  file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_PPG188_binwidth_48.root");
		  file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_105intersect108_mixed_local.root");
		}
	      else if(bin_width == 8)
		{
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_PPG188_2.root"); 
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_PPG188_pro108_2.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_nofvtx_dimu_02_2.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_sngtrk_dimu_02_2.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_mass_fvtx_dimu02_noprob_int.root");
		  //  file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_fvtx_sngtrk_dimu_02_2.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_PPG188_binwidth_48.root");
		  file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_105intersect108_mixed_local.root");
		}
	      else 
		{
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_PPG188_int.root"); 
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_PPG188_pro108_int.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_nofvtx_dimu_02_int.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_sngtrk_dimu_02_int.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_mass_fvtx_dimu02_noprob_int.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pp_N_rebinned_mass_pro108_fvtx_sngtrk_dimu_02_int.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_PPG188_binwidth_48.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pAu_N_rebinned_105intersect108_mixed_local.root");

		  //  file1S = TFile::Open("psi2S_root_files/fvtx_newlib_notrig_root5_N.root");
		  //  file1S = TFile::Open("psi2S_root_files/nofvtx_newlib_N.root");
		  //  file1S = TFile::Open("psi2S_root_files/mutr_newlib_root5_local_N.root");
		  ////////////////////////////////////////////////////////////////////////////////////////
		  file1S = TFile::Open("Run15pp_DiMuons_rebinned_N.root");
		  ////////////////////////////////////////////////////////////////////////////////////////
		}
	    } // special == 0
	  else
	    {
	      if(pt_bin == 2)
		file1S = TFile::Open("Run15pp_N_rebinned_PPG_binwidth_special_2.root"); 
	    } // special == 1
	  cout << "You are running North arm data " << endl;
	} // north arm 
      else
	{
	  if(special == 0)
	    {
	      if(bin_width == 4)
		{
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_PPG188.root"); 
			//	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_PPG188_pro108.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_nofvtx_dimu_02.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_sngtrk_dimu_02.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_mass_fvtx_dimu02_noprob_int.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_PPG188_binwidth_48.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_fvtx_sngtrk_dimu_02.root");
		  file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_105intersect108_mixed_local.root");
		}
	      else if(bin_width == 8)
		{
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_PPG188_2.root"); 
			//	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_PPG188_pro108_2.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_nofvtx_dimu_02_2.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_sngtrk_dimu_02_2.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_mass_fvtx_dimu02_noprob_int.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_fvtx_sngtrk_dimu_02_2.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_PPG188_binwidth_48.root");
		  file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_105intersect108_mixed_local.root");
		}
	      else      
		{
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_PPG188_int.root"); 
			//	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_PPG188_pro108_int.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_nofvtx_dimu_02_int.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_sngtrk_dimu_02_int.root");
		  //	file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_mass_fvtx_dimu02_noprob_int.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pp_S_rebinned_mass_pro108_fvtx_sngtrk_dimu_02_int.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_PPG188_binwidth_48.root");
		  // file1S = TFile::Open("psi2S_root_files/Run15pAu_S_rebinned_105intersect108_mixed_local.root");

		  //	 file1S = TFile::Open("psi2S_root_files/fvtx_newlib_notrig_root5_S.root");
		 //  file1S = TFile::Open("psi2S_root_files/nofvtx_newlib_S.root");
		  //  file1S = TFile::Open("psi2S_root_files/mutr_newlib_root5_local_S.root");
		  ////////////////////////////////////////////////////////////////////////////////////////
		  file1S = TFile::Open("Run15pp_DiMuons_rebinned_S.root");
		  ////////////////////////////////////////////////////////////////////////////////////////
		}
	    } // special == 0
	  else
	    {
	     if(pt_bin == 2)
		file1S = TFile::Open("Run15pp_S_rebinned_PPG_binwidth_special_2.root"); 
	    } // special == 1
	  cout << "You are running South arm data " << endl;
	}
  

  if(!file1S)
    {
      cout << " Failed to open input root file" << endl;
      return;
    }
  

  TCanvas *cups = new TCanvas("cups","cups",5,5,800,600);
  cups->SetLogy();
  cups->SetTicky();
  cups->SetTickx();

  cups->cd();


  int nrebin = 1;  

  char ulname[500];
  char mmname[500];
  char ppname[500];
  char mixname[500];
  
  sprintf(ulname,"ul_%i",bin);
  cout << "Open " << ulname << endl;
  recomassFG = (TH1 *)file1S->Get(ulname);
  if(!recomassFG)
    cout << "Failed to get " << ulname << endl; 
 
  sprintf(mmname,"mm_%i",bin);
  cout << "Open " << mmname << endl;
  recomassBG1 = (TH1 *)file1S->Get(mmname);

  if(!recomassBG1)
    cout << "Failed to get " << mmname << endl; 

  sprintf(ppname,"pp_%i",bin);
  cout << "Open " << ppname << endl;
  recomassBG2 = (TH1 *)file1S->Get(ppname);

  if(!recomassBG1)
    cout << "Failed to get " << ppname << endl;

  sprintf(mixname,"mix_ul_%i",bin);
  cout << "Open " << mixname << endl;
  recomassBG3 = (TH1 *)file1S->Get(mixname);

  if(!recomassBG3)
    cout << "Failed to get " << mixname << endl;
  
  double dNpp, dNmm, dNmix;
  double Npp = recomassBG2->IntegralAndError(recomassBG2->GetXaxis()->FindBin(2), recomassBG2->GetXaxis()->FindBin(5),dNpp,"");
  double Nmm = recomassBG1->IntegralAndError(recomassBG1->GetXaxis()->FindBin(2), recomassBG1->GetXaxis()->FindBin(5),dNmm,"");
  double bgnorm = 2*sqrt(Npp*Nmm)/(Npp+Nmm);
  cout << "bgnorm " << bgnorm << " Npp " << Npp << " Nmm " << Nmm << endl;

  recomassBG = (TH1D*)recomassBG1->Clone("recomassBG");  
  recomassBG->SetName("recomassBG"); 
  recomassBG->SetTitle("recomassBG"); 
  recomassBG->SetLineColor(kRed);
  for(int i = 0; i < recomassBG->GetNbinsX(); i++)
    {
      recomassBG->SetBinContent(i+1,(recomassBG2->GetBinContent(i+1) + recomassBG1->GetBinContent(i+1))); 
      recomassBG->SetBinError(i+1,(recomassBG2->GetBinError(i+1) + recomassBG1->GetBinError(i+1)));
    }
  cout <<  recomassBG->GetXaxis()->GetNbins() << endl;


  /////////////////////////////////////////////////////////////////////////////////////////
  double Nmix = recomassBG3->IntegralAndError(recomassBG3->GetXaxis()->FindBin(2), recomassBG3->GetXaxis()->FindBin(5),dNmix,"");
  double mixnorm = 2*sqrt(Npp*Nmm)/Nmix;
  recomassBG->Scale(bgnorm);
  recomassBG3->Scale(mixnorm);
  /////////////////////////////////////////////////////////////////////////////////////////
   

  double m0,m1,m2,m3,m4;
   
  for(int i = 0; i < 100; i++)
    {
      mass_array_combg[i] = recomassBG->GetBinContent(i+1);
      mass_array_combg_err[i] = recomassBG->GetBinError(i+1);
    }

  bool CaseF = false;
  bool CaseG = false;
  bool CaseH = false;
  bool CaseA = false;
  


  if(case_type == 0)
    CaseA = true;
  if(case_type == 1)
    CaseF = true;
  if(case_type == 2)
    CaseG = true;
  if(case_type == 3)
    CaseH =true;

  // combinatoric BG initital parameters for pT dependent fitting - USE SAME SET FOR ALL SYSTEMS AND ALL DEPENDENCES

  if(com_bg == 0)
    {
      m0 = 0.1;
      m1 = 0.01;
      m2 = 10; 
      m3 = 10; 
      m4 = 10;
    }
  if(com_bg == 1)
    {
      m0 = -0.1;
      m1 = 0.01;      
      m2 = 10;
      m3 = 10;
      m4 = 10;
    }
     if(com_bg == 2)   
    {
      m0 = 0.01;
      m1 = 0.001;
      m2 = 100;
      m3 = 10;
      m4 = 10;
    }
     if(com_bg == 3)   
    {
      m0 = 0.1;
      m1 = 0.1;      
      m2 = 10;
      m3 = 10;
      m4 = 10;
    }
     if(com_bg == 4)    
    {
      m0 = 0.1;
      m1 = 0.01;      
      m2 = 0.01;
      m3 = 10;
      m4 = 10;
    }
   if(com_bg == 5)  
    {
      m0 = -0.001;
      m1 = -0.01;
      m2 = 10; 
      m3 = 1; 
      m4 = 1;
    }
  
  if(com_bg == 6)
    {
      m0 = 0.01;
      m1 = 0.1;
      m2 = 1; 
      m3 = 10; 
      m4 = 10;
    }
 
 if(com_bg == 7)
    {
      m0 = -0.01;
      m1 = 0.01;
      m2 = 10; 
      m3 = 10; 
      m4 = 1;
    }

 if(com_bg == 8)
    {
      m0 = 0.01;
      m1 = 0.01;
      m2 = 10; 
      m3 = 10; 
      m4 = 1;
    }
  if(com_bg == 9)
    {
      m0 = -0.1;
      m1 = 0.01;      
      m2 = 100;
      m3 = 10;
      m4 = 10;
    }
  if(com_bg == 10)
    {
      m0 = 0.1;
      m1 = -0.01;      
      m2 = 1;
      m3 = 10;
      m4 = 10;
    }
 if(com_bg == 11)
    {
      m0 = 10;
      m1 = 10;
      m2 = 10;
      m3 = 100; 
      m4 = 0.01;
    }
 if(com_bg == 12)
    {
      m0 = 0.01;
      m1 = 10;
      m2 = 10;
      m3 = 100; 
      m4 = 0.01;
    }
  if(com_bg == 13)
    {
      m0 = -0.1;
      m1 = 0.01;      
      m2 = 1;
      m3 = 10;
      m4 = 10;
    }

if(com_bg == 14)
    {
       m0 = -0.01;
      m1 = 0.01;
      m2 = 10; 
      m3 = 10; 
      m4 = 10;
    }



  combg_fit = new TF1("combg_fit"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",2,5);
  combg_fit->SetParameter(0, m0);  
  combg_fit->SetParameter(1, m1);
  combg_fit->SetParameter(2, m2); 
  combg_fit->SetParameter(3, m3);
  combg_fit->SetParameter(4, m4);

  combg_fit->SetLineColor(kRed);
  combg_fit->SetLineStyle(2);

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000); 
 
  if(bg_fit == 1)
    {
      //////////////////////////////////////////////////////////
      // recomassBG->Fit(combg_fit,"LL","",2,5); 
      recomassBG3->Fit(combg_fit,"LL","",2,5); 
      //////////////////////////////////////////////////////////
    }
  else
    combg_fit->SetParameter(2,0); 

  TF1 *total_fit;

  if(LL_fit)
    { 
      total_fit = new TF1("total_fitLL",CBcalc_LL,2.0,maxmass,14);
    }
  // else
  //  { 
  //    total_fit = new TF1("total_fit",CBcalc_exp,2.0,maxmass,10);    // not using anymore (simple exponential bg)
  // }

  cout << recomassFG->GetXaxis()->GetBinWidth(1) << endl;

  double fit_combg[100] = {0};
  double fit_combg_err[100] = {0};
  double mass_combg_err[100] = {0};
  double temp[100] = {0};

////////////////////////////////////
  for(int i = 0; i < 100; i++)
    {
      
      mass_array_combg[i] = recomassBG->GetBinContent(i+1);
      mass_combg_err[i] = recomassBG->GetBinError(i+1);
      //  mass_array_corrbg[i] = corrbg_fit->Eval(x_array[i]);
       mass_array_jpsi[i] = recomassFG->GetBinContent(i+1);
       mass_array_jpsi_err[i] = recomassFG->GetBinError(i+1);
      // mass_array_total[i] = total_fit->Eval(x_array[i]);
       // mass_array_psi2s[i] = psi2s_fit->Eval(x_array[i]);	 
       //cout << combg_fit->Eval(bg_array[i]) << "for bin " <<  i+1 << endl;
       //cout << mass_array_combg[i] << "for bin " <<  i+1 << endl;
    	 
    }
  for(int i = 0; i < 60; i++)
    {
      fit_combg[i] = combg_fit->Eval(bg_array[i]);
      fit_combg_err[i] = sqrt(combg_fit->Eval( bg_array[i] ));
    }

  int Bins = 60;
  cout << "fit_LS_N[" << Bins << "] = {";
      for(int i = 0; i < Bins-1; i++)
 	{
	  cout << mass_array_combg[i] << ", ";
 	}
       cout << mass_array_combg[Bins-1] << "};" << endl;
      

  // Get the starting estimate for N from the integral of the data from 2.8-3.4
  double binw = recomassFG->GetBinWidth(1);
  double renorm = 1.0/binw;   // (1 / (bin_width of data in GeV) )
  cout << "total_fit created " << endl;

  double pt_initial = bin*0.25 - 0.25;
  double pt_final = bin*0.25;
  
  double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
  Double_t NJpsi,err_NJpsi;
 Double_t Npsi2s,err_psi2s;

  int binl = recomassFG->FindBin(2.6); // 2.8 //2.6
  int binh = recomassFG->FindBin(3.6); // 3.4 // 3.6

  double coeff;
  
  

  double N_start = 0;

  bool write = false; 
  //bool write = true; 

   if(north_arm)
    {
      coeff = 1.5;
      N_start = coeff*recomassFG->Integral(binl,binh) / renorm;  // 1 // 1.5
       f0 = 1;
       f1 = 10;
      f2 = 3.1; // Jpsi x_mean
      f3 = 0.08;  // Jpsi sigma
      f4 = N_start;  // Jpsi norm

      f5 = 0; //  norm 2nd gaussian
      f6 = 0.589;  //  x_mean ppsi2s
      f7 = 1.15;   // width factor psi2s
 
      f8 = 0.1*N_start; // norm psi2s

      if(fit == 1)
 	{
 	  f9 = 0.01;
 	  f10 = 0.01;
 	  f11 = 0.01;
 	  f12 = 10;
 	  f13 = 1;
 	}
    
    }
 
  else
    {
      coeff = 1.5;
      N_start = coeff*recomassFG->Integral(binl,binh) / renorm;  // 1 // 1.5
      f0 = 1;  // alpha
      f1 = 10.0; // n
      f2 = 3.1; // Jpsi x_mean
      f3 = 0.08;  // Jpsi sigma
      f4 = N_start;  // Jpsi norm
      f5 = 0; //  norm 2nd gaussian for fvtx
      f6 = 0.589;  //  x_mean ppsi2s
      f7 = 1.15;   // width factor psi2s
 
      f8 = 0.1*N_start; // norm psi2s

      if(fit == 1)
 	{
 	  f9 = 0.01;
 	  f10 = 0.01;
 	  f11 = 0.01;
 	  f12 = 10;
 	  f13 = 1;
 	}
      
    }
   // //////////////////////////////////////////////////////////////
 
  //  if(north_arm)
  //    {
  //      total_fit->SetParameter(0, f0);
  //      total_fit->SetParLimits(0, 0, 2);
  //      total_fit->SetParameter(1, f1);  
  //      total_fit->SetParLimits(1, 0, 5); 
  //    }
  //  else
  //    {
  //      total_fit->SetParameter(0, f0);
  //      total_fit->SetParLimits(0, 0, 2);
  //      total_fit->SetParameter(1, f1);  
  //      total_fit->SetParLimits(1, 0, 5); 
  //    }
	 
  // total_fit->SetParameter(0, f0);
  // total_fit->SetParameter(1, f1);  
  // total_fit->SetParameter(2, f2); // Jpsi CB mean
  // total_fit->SetParameter(3, f3);  // Jpsi CB sigma


 //////////////////////////////////////////////////////////////////////  

   if(north_arm)  // from mixed evnts floater stretch in matt_centrality directory
     {
       f2 = 3.17482;
       f3 = 0.082614;
     }
   else
     {
       f2 = 3.12959;
       f3 = 0.0925465;
     }


   if(north_arm)
     {
       total_fit->SetParameter(0, f0);
       total_fit->SetParLimits(0, 0, 2);
       total_fit->SetParameter(1, f1);
       total_fit->SetParLimits(1, 0, 100);

       // total_fit->SetParameter(0, f0);
       // total_fit->SetParLimits(0, 1.1, 1.3);
       // total_fit->SetParameter(1, f1);
       // total_fit->SetParLimits(1, 85, 95);
     }
   else
     {
       total_fit->SetParameter(0, f0);
       total_fit->SetParLimits(0, 0, 2);
       total_fit->SetParameter(1, f1);
       total_fit->SetParLimits(1, 0, 100);

       // total_fit->SetParameter(0, f0);
       // total_fit->SetParLimits(0, 1.1, 1.3);
       // total_fit->SetParameter(1, f1);
       // total_fit->SetParLimits(1, 85, 95);
     }
  
   //////////////////////////////////////////////////////////////////////  
	 
   total_fit->SetParameter(0, f0);
   total_fit->SetParameter(1, f1);  

   //////////////////////////////////////////////////////////////////////  
   // total_fit->FixParameter(2, f2); // Jpsi CB mean 
   // total_fit->FixParameter(3, f3);  // Jpsi CB sigma
   ////////////////////////////////////////////////////////////////////
     total_fit->SetParameter(2, f2); // Jpsi CB mean
    total_fit->SetParameter(3, f3);  // Jpsi CB sigma
   ////////////////////////////////////////////////////////////////////

   total_fit->SetParameter(4, f4);    // N J/psi
   total_fit->SetParameter(5, f5);    //norm 2nd G Jpsi
   total_fit->FixParameter(6, f6);    // mass offset psi2s
   total_fit->FixParameter(7, f7);    // width factor psi2s
   total_fit->SetParameter(8, f8);  //  N psi2s
   total_fit->SetParameter(9, f9);  // a
   total_fit->SetParameter(10, f10);  // b
   total_fit->SetParameter(11, f11);  // c
   total_fit->SetParameter(12, f12);  // d
   total_fit->SetParameter(13, f13);  // e
   
  double a,b,c,d,e,a_err,b_err,c_err,d_err,e_err;
  char name100[500];

  if(fixed_corr_bg == true)
    {
      if( (north_arm) && (bin < 3) )
	sprintf(name100,"tony_bestfit_parameters/yuehang_bestfit_200mev/bestfit_parameters_low_pt_N_%i.dat", bin);
      if( (north_arm) && ( (bin == 3) || (bin == 4) || (bin == 5)) )
	sprintf(name100,"tony_bestfit_parameters/yuehang_bestfit_200mev/bestfit_parameters_low_pt_N_%i.dat", bin+1);
      if( (north_arm == false) && (bin < 3) )
	sprintf(name100,"tony_bestfit_parameters/yuehang_bestfit_200mev/bestfit_parameters_low_pt_S_%i.dat", bin);
      if( (north_arm == false) && ( (bin == 3) || (bin == 4) || (bin == 5)) )
	sprintf(name100,"tony_bestfit_parameters/yuehang_bestfit_200mev/bestfit_parameters_low_pt_S_%i.dat", bin+1);
      if( north_arm && bin == 14 && bin_width == 2)
	sprintf(name100,"tony_bestfit_parameters/yuehang_bestfit_500mev/bestfit_parameters_N_14.dat");
      if( bin == 14 && bin_width == 2)
	sprintf(name100,"tony_bestfit_parameters/yuehang_bestfit_500mev/bestfit_parameters_N_14.dat");
      
      ifstream Run15pp_testfit(name100);
      if(Run15pp_testfit)
	{
	  do 
	    {
	      double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9;	  
	      Run15pp_testfit >> f0 >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9;
	      a = f0; a_err = f1; b = f2; b_err = f3; c = f4; c_err = f5; d = f6; d_err = f7; e = f8; e_err = f9;
	      cout << "c: " << c << endl;
	    } while(Run15pp_testfit.good());
	}
      else
	cout << "test read file does not exist" << endl;
      
      cout << "Open " << name100 << endl;
           

      double c_array[4][5] = { // NORTH ARM

	1*pow(10,-17),1*pow(10,-13),1*pow(10,-9),1*pow(10,-6),1*pow(10,-13), //   Case A - 19

	//1*pow(10,-19),1*pow(10,-15),1*pow(10,-8),1*pow(10,-5),1*pow(10,-13), //   Case A for pAl sys study
      

      1*pow(10,-18),1*pow(10,-18),1*pow(10,-16),1*pow(10,-34),1*pow(10,-13), //   Case F

      1*pow(10,-19),1*pow(10,-18),1*pow(10,-11),1*pow(10,-28),1*pow(10,-13), //  Case G

      1*pow(10,-19),1*pow(10,-16),1*pow(10,-16),1*pow(10,-16),1*pow(10,-13) }; //   Case H
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      // south arm
      double c_array_S[4][5] = { // SOUTH ARM

	1*pow(10,-17),1*pow(10,-17),1*pow(10,-7),1*pow(10,-2),1*pow(10,-1), // CaseA
	//1*pow(10,-19),1*pow(10,-17),1*pow(10,-9),1*pow(10,-5),1*pow(10,-1), // for pAl sys study

	1*pow(10,-17),1*pow(10,-18),1*pow(10,-2),1*pow(10,-5),1*pow(10,-3), //  CaseF   
	1*pow(10,-29),1*pow(10,-17),1*pow(10,-3),1*pow(10,-5),1*pow(10,-5), //  CaseG
	1*pow(10,-13),1*pow(10,-15),1*pow(10,-3),1*pow(10,-5),1*pow(10,-2)}; // CaseH 
      
	      
      c = c_array[case_type][bin-1];
      
      if(!north_arm)
	c = c_array_S[case_type][bin-1];
      
      // set correlated bg parameters to the bestfit values from Yue hang corr bg fits
      if(CaseA)
	{
	  total_fit->FixParameter(9,a);
	  if( bin == 4)
	    b = 0;
	  total_fit->FixParameter(10,b);
	  total_fit->SetParameter(11,c);
	  if(bin == 14)
	    //c = 0.001;
	    c = 0.1;
	  total_fit->FixParameter(12,d);
	  total_fit->SetParameter(13,e);
	}
      
      if(CaseF)
	{
	  total_fit->SetParameter(9,a);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,c);
	  total_fit->SetParameter(12,d);
	  total_fit->FixParameter(13,e);
	}
      if(CaseG)
	{
	  total_fit->FixParameter(9,a);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,c);
	  total_fit->SetParameter(12,d);
	  total_fit->SetParameter(13,e);
	}
      if(CaseH)
	{
	  total_fit->SetParameter(9,a);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,c);
	  total_fit->FixParameter(12,d);
	  total_fit->SetParameter(13,e);
	  
	}
      
      
      cout <<"par c: " << c << endl;
      
    } // end fixed bg
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  double par_a[10];
  double par_b[10];
  double par_d[10];
  double par_e[10];
  
  double par_a_fx;
  double par_b_fx;
  double par_d_fx;
  double par_e_fx;
  
  double x;
  

  double pt_center[26] = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,6.25,6.75}; 

  double par_c_array_North[4][28] = {
     {pow(10,-5),pow(10,-5),pow(10,-1),pow(10,-1),0.1,0.1,10,10,10,10,10,10,1,0.1,0.1,0.1,0.01,0.01,0.1,1,1,1,0.1,0.1,pow(10,-11),1,100,100}, // Case A
     //{pow(10,-5),pow(10,-5),pow(10,-1),pow(10,-1),0.01,0.1,10,100,10,10,10,10,1,0.1,0.1,0.1,0.01,1,0.1,1,0.01,0.01,0.1,pow(10,-10),pow(10,-11),10,100,100}, // testing for unstable fit funtion here
    // {pow(10,-5),pow(10,-5),pow(10,-1),pow(10,-1),0.1,0.01,0.001,1,10,10,10,10,1,0.1,0.1,0.01,1,1,1,1,0.01,0.01,0.1,0.001,0.1,10,100,100}, // for HeAu Centrality
    // {pow(10,-5),pow(10,-5),pow(10,-2),pow(10,-1),0.1,0.1,0.1,10,0.01,0.01,0.001,10,1,0.1,0.1,0.1,0.01,0.01,0.1,1,1,1,0.1,0.1,pow(10,-11),1,100,100}, // pAl bins 17,18,19 ---> 9,10,3
     

    {pow(10,-7),pow(10,-6),pow(10,-4),pow(10,-3),0.1,0.01,0.1,0.01,0.01,0.1,0.1,0.1,0.1,0.01,0.1,0.1,0.1,10,0.01,0.01,1,0.1,0.1,0.1,10,10,100,100}, // Case F
    {pow(10,-5),pow(10,-3),pow(10,-2),pow(10,-3),0.1,0.01,0.1,0.1,1,1,0.1,0.01,0.1,0.01,0.1,0.1,1,1,1,0.1,0.0001,0.01,1,10,1,10,100,100}, // Case G
    {pow(10,-5),pow(10,-3),pow(10,-2),pow(10,-1),0.1,0.1,0.1,0.1,10,10,1,0.1,0.1,0.01,0.1,0.1,0.001,0.01,0.001,0.01,0.01,0.1,0.1,0.1,10,10,100,100}}; // Case H

  double par_c_array[4][28] = {
     {pow(10,-1),pow(10,-1),pow(10,-1),pow(10,-1),0.001,10,10,1,10,10,1,10,0.1,0.01,0.0001,0.1,0.01,1,0.1,0.01,0.001,0.1,0.001,0.1,0.1,0.01,100,100}, // Case A redone after binning issue  (25 = 0.1)
     //  {pow(10,-1),pow(10,-1),pow(10,-2),pow(10,-1),0.001,0.001,0.001,0.01,0.001,0.001,1,10,0.1,0.01,0.0001,0.1,0.01,1,0.1,0.0001,0.001,0.1,0.001,0.001,0.1,0.01,100,100}, // pAl bins
   


   {1,pow(10,0),pow(10,-1),pow(10,-1),0.01,1,0.1,10,10,10,0.001,10,0.001,0.1,0.01,0.1,1,1,pow(10,-10),0.01,0.01,0.1,0.01,0.001,0.1,0.01,1,1}, // Case F
   {pow(10,0),pow(10,0),pow(10,0),pow(10,-13),1,0.1,1,0.001,100,100,10,0.1,1,10,0.1,0.1,1,1,pow(10,-6),0.01,0.001,0.01,0.01,0.001,0.01,0.01,100,100}, // Case G
    {pow(10,0),pow(10,0),pow(10,-1),pow(10,-1),0.01,0.1,0.1,10,10,10,0.1,0.1,0.1,10,1,0.1,0.01,0.1,1,1,0.01,0.01,0.01,0.001,1,1,100,100}}; // Case H
  
  if(calculated_fit == true)
    {
      
      if(bin_width == 1)
	x = pt_center[bin-1];
      if(bin == 13 && bin_width == 2)
	{
	  x = pt_center[24];
	  f13 = par_c_array[case_type][24];
	}
      if(bin == 14 && bin_width == 2)
	{
	  x = pt_center[25];
	  f13 = par_c_array[case_type][25];
	}
      
      if(pAu == true) // centrality binning
	{
	  if(bin == 8 && bin_width == 2)
	    {
	      x = 4.25;
	      f13 = pow(10,-3);
	    }
	  if(bin == 9 && bin_width == 2)
	    {
	      x = 4.75;
	      f13 = pow(10,-1);
	    }
	  if(special == 1 && bin == 3)
	    {
	      x = 6.00;
	      f13 = pow(10,-4);
	      if(!north_arm)
		f13 = pow(10,-4);
	    }
	}
      
      if(HeAu == true) // centrality binning here
	{
	  if(bin == 6 && bin_width == 2)
	    {
	      x = 2.75;
	      f13 = pow(10,0);
	    }
	  if(bin == 7 && bin_width == 2)
	    {
	      x = 3.25;
	      f13 = pow(10,-4);
	    }
	  if(bin == 8 && bin_width == 2)
	    {
	      x = 3.75;
	      f13 = pow(10,-4);
	    }
	  if(bin == 5 && bin_width == 4)
	    {
	      x = 4.5;
	      f13 = pow(10,-2);
	    }
	}

      if(pAl == true) // 0-100 binning
	{
	  if(bin == 6 && bin_width == 4)
	    {
	      x = 5.5;
	    }
	  if(bin == 7 && bin_width == 4)
	    {
	      x = 6.5;
	    }
	  if(bin == 9 && bin_width == 2)
	    {
	      x = 4.25;
	    }
	  if(bin == 10 && bin_width == 2)
	    {
	      x = 4.75;
	    }
	}
      
      cout << "x = pT center: " << x << endl;
    

      for(int par = 0; par < 10; par++)
	{	 
	  if(north_arm == true)
	    {
	      par_a[par] = coeff_par_N[0][par];
	      par_b[par] = coeff_par_N[1][par];
	      par_d[par] = coeff_par_N[2][par];
	      par_e[par] = coeff_par_N[3][par];
	      par_c_array[case_type][bin-1] = par_c_array_North[case_type][bin-1];
	      if(bin_width == 2 && bin == 13)
		par_c_array[case_type][24] = par_c_array_North[case_type][24];
	      if(bin_width == 2 && bin == 14)
		par_c_array[case_type][25] = par_c_array_North[case_type][25];
	     
	    }
	  else
	    {
	      par_a[par] = coeff_par_S[0][par];
	      par_b[par] = coeff_par_S[1][par];
	      par_d[par] = coeff_par_S[2][par];
	      par_e[par] = coeff_par_S[3][par];
	    }
	}
      
      par_a_fx = par_a[0] + par_a[1]*x + par_a[2]*x*x + par_a[3]*x*x*x + par_a[4]*x*x*x*x + par_a[5]*x*x*x*x*x + par_a[6]*x*x*x*x*x*x + par_a[7]*x*x*x*x*x*x*x + par_a[8]*x*x*x*x*x*x*x*x + par_a[9]*x*x*x*x*x*x*x*x*x;
      par_b_fx = par_b[0] + par_b[1]*x + par_b[2]*x*x + par_b[3]*x*x*x + par_b[4]*x*x*x*x + par_b[5]*x*x*x*x*x + par_b[6]*x*x*x*x*x*x + par_b[7]*x*x*x*x*x*x*x + par_b[8]*x*x*x*x*x*x*x*x + par_b[9]*x*x*x*x*x*x*x*x*x;
      par_d_fx = par_d[0] + par_d[1]*x + par_d[2]*x*x + par_d[3]*x*x*x + par_d[4]*x*x*x*x + par_d[5]*x*x*x*x*x + par_d[6]*x*x*x*x*x*x + par_d[7]*x*x*x*x*x*x*x + par_d[8]*x*x*x*x*x*x*x*x + par_d[9]*x*x*x*x*x*x*x*x*x;
      par_e_fx = par_e[0] + par_e[1]*x + par_e[2]*x*x + par_e[3]*x*x*x + par_e[4]*x*x*x*x + par_e[5]*x*x*x*x*x + par_e[6]*x*x*x*x*x*x + par_e[7]*x*x*x*x*x*x*x + par_e[8]*x*x*x*x*x*x*x*x + par_e[9]*x*x*x*x*x*x*x*x*x;
      
      cout << "a0: " << par_a[0] << ", a1: " << par_a[1] << ", a2: " << par_a[2] << ", a3: " << par_a[3] << ", a4: " << par_a[4] << ", a5: " << par_a[5] << ", a6: " << par_a[6] << ", a7: " << par_a[7] << ", a8: " << par_a[8] << ", a9: " << par_a[9] << endl;
      
      cout << "b0: " << par_b[0] << ", b1: " << par_b[1] << ", b2: " << par_b[2] << ", b3: " << par_b[3] << ", b4: " << par_b[4] << ", b5: " << par_b[5] << ", b6: " << par_b[6] << ", b7: " << par_b[7] << ", b8: " << par_b[8] << ", b9: " << par_b[9] << endl;
      
      cout << "d0: " << par_d[0] << ", d1: " << par_d[1] << ", d2: " << par_d[2] << ", d3: " << par_d[3] << ", d4: " << par_d[4] << ", d5: " << par_d[5] << ", d6: " << par_d[6] << ", d7: " << par_d[7] << ", d8: " << par_d[8] << ", d9: " << par_d[9] << endl;
      
      cout << "e0: " << par_e[0] << ", e1: " << par_e[1] << ", e2: " << par_e[2] << ", e3: " << par_e[3] << ", e4: " << par_e[4] << ", e5: " << par_e[5] << ", e6: " << par_e[6] << ", e7: " << par_e[7] << ", e8: " << par_e[8] << ", e9: " << par_e[9] << endl;
      
      if(!pAu && !HeAu)
	f13 = par_c_array[case_type][bin-1];
      
      if(CaseA)
	{
	  total_fit->FixParameter(9, par_a_fx);
	  if( (bin == 21 || bin == 7) && !north_arm)
	    par_b_fx = 0;
	  total_fit->FixParameter(10, par_b_fx);
	  total_fit->SetParameter(11, f13);
	  total_fit->FixParameter(12, par_d_fx);
	  total_fit->SetParameter(13, par_e_fx);
	}

      if(CaseF)
	{
	  total_fit->SetParameter(9,par_a_fx);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,f13);
	  total_fit->SetParameter(12,par_d_fx);
	  total_fit->FixParameter(13,par_e_fx);
	}
      if(CaseG)
	{
	  total_fit->FixParameter(9,par_a_fx);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,f13);
	  total_fit->SetParameter(12,par_d_fx);
	  total_fit->SetParameter(13,par_e_fx);
	}
      if(CaseH)
	{
	  total_fit->SetParameter(9,par_a_fx);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,f13);
	  total_fit->FixParameter(12,par_d_fx);
	  total_fit->SetParameter(13,par_e_fx);
	}
      
      
    } // pt calculated
  
  
  if(fx_pt_int)
    {
      //  March 9, 2019:
      // bestfit parameters from f(x) pT integrated NORTH method (yuehang_pt_integrated.C)
      if(north_arm == true)
	{
	  a = 1.01851e-02;
	  b = 3.15392e-01;
	  c = par_c_array_North[case_type][bin-1];
	  if(bin == 49)
	    {
	      if(CaseA)
		c = 10;
	      if(CaseF)
		c = 0.1;
	      if(CaseG)
		c = 1000;
	      if(CaseH)
		c = 1;
	    }
	  if(bin !=49)
	    c = 0.001;

	      d = 2.78569e+00;
	  e =  3.65850e+00;
	}
      else
	{
	  a =  -1.96190e-01;
	  b =  3.90264e-01;
	  c = par_c_array[case_type][bin-1];
	  if(bin == 49)
	    {
	      if(CaseA)
		c = 100;
	      if(CaseF)
		c = 1;
	      if(CaseG)
		c = 1000;
	      if(CaseH)
		c = 10;
	    }
	  //if(bin !=49)
	      d =  2.40817e+00;
	      e = 3.74928e+00;
	}
      
      if(CaseF)
	{
	  total_fit->SetParameter(9,a);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,c);
	  total_fit->SetParameter(12,d);
	  total_fit->FixParameter(13,e);
	}
      if(CaseG)
	{
	  total_fit->FixParameter(9,a);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,c);
	  total_fit->SetParameter(12,d);
	  total_fit->SetParameter(13,e);
	}
      if(CaseH)
	{
	  total_fit->SetParameter(9,a);
	  total_fit->FixParameter(10,0);
	  total_fit->SetParameter(11,c);
	  total_fit->FixParameter(12,d);
	  total_fit->SetParameter(13,e);
	}
      
      if(CaseA)
	{
	  
	  total_fit->FixParameter(9,a);
	  total_fit->FixParameter(10,b);
	  total_fit->SetParameter(11,c);
	  total_fit->FixParameter(12,d);
	  total_fit->SetParameter(13,e);
	}
      
    }// fx_pt_int
  // *******************************ABOVE FOR NEW YUE HANG STUDY *******************************************************
  
    
  total_fit->SetLineColor(kBlue);
  total_fit->SetMarkerStyle(42);
  total_fit->SetLineWidth(3);
  total_fit->SetLineStyle(kDashed);
  
  TF1 *fresult= 0;

 
  /////////////////////////////////////////////////////////////BEGIN ORIGINAL CODE
  


  
  if(LL_fit)
    {
      // log likelihood fit with integer data option is "L"
      TFitResultPtr r = recomassFG->Fit(total_fit,"LS");////// does the fit
      fresult = recomassFG->GetFunction("total_fitLL");
      
      TMatrixDSym cov = r->GetCovarianceMatrix();
      r->Print("V"); //verbose setting
      Double_t *fullmat = 0;
      fullmat = cov.GetMatrixArray();

      double Jpsi_array[36] = {0};

      for(Int_t i = 0;i < 6; i++)
	{
	  for(Int_t j = 0;j < 6; j++)
	    {
	      // Jpsi_array[5*i+j] = fullmat[10*i+j];
	      Jpsi_array[6*i+j] = fullmat[14*i+j];
	      // cout << "Jpsi array " << i << " " << j << " " << Jpsi_array[8*i+j] << endl;
	    }  
	}
     

   
      

      TF1 *jpsi_fit = new TF1("jpsi_fit",CBcalc_jpsi,2.0,5.0,6); // use CBcalc_jpsi instead
      for(int i=0;i<6;i++)
	{
	  jpsi_fit->SetParameter(i,total_fit->GetParameter(i));
	  jpsi_fit->SetParError(i,total_fit->GetParError(i));
	}
       
      jpsi_fit->SetLineColor(kBlack);
      
            
      Double_t width = recomassFG->GetBinWidth(1);
      NJpsi = jpsi_fit->Integral(2.0, 5.0)/width;
           
      Double_t Jpsipar[6] = {0};
      
      for(int i = 0;i < 6; i++)
	{
	  Jpsipar[i] = total_fit->GetParameter(i); 
	}
   
  
      err_NJpsi = jpsi_fit->IntegralError(2.0,5.0,Jpsipar,Jpsi_array)/width;
      recomassFG->GetYaxis()->SetTitle("N_{#mu#mu} / (50 MeV/c^{2})");
      recomassFG->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (Gev/c^{2})");
      cups->SetLogy();
      recomassFG->DrawCopy();

      double percent_error;
      percent_error = err_NJpsi/NJpsi*100; 
      if(north_arm)
	cout << "North p + p jpsi counts " << NJpsi << " +/- " << err_NJpsi << " , sqrt error " << sqrt(NJpsi) << " and % error " << percent_error << endl;
      //cout << "North p + p jpsi counts " << NJpsi << " , sqrt error " << sqrt(NJpsi) << " , error " << err_NJpsi << " and % error " << percent_error << endl;
      if(!north_arm)
	cout << "South p + p jpsi counts " << NJpsi << " +/- " << err_NJpsi << " , sqrt error " << sqrt(NJpsi) << " and % error " << percent_error << endl;
      //	cout << "South p + p jpsi counts " << NJpsi << " , sqrt error " << sqrt(NJpsi) << " , error " << err_NJpsi << " and % error " << percent_error << endl;
      cout <<" \n\nIntegral error: " << err_NJpsi*width << endl;
      cout <<"Error upper bound: " << (err_NJpsi*width)*2*(0.01) << endl;
   
      percent_error = 0;




      //////// psi(2s)/////////////////////////
   
     
      double psi2s_array[36] = {0};

      //////////////////////////////////////
      // par 0 row (alpha)
      psi2s_array[0] = fullmat[0]; // alpha
      psi2s_array[1] = fullmat[1]; // n
      psi2s_array[2] = fullmat[2]; // x_mean J/psi
      psi2s_array[3] = fullmat[3]; // sigma Jpsi
      psi2s_array[4] = fullmat[8]; // N psi2S
      psi2s_array[5] = fullmat[5]; // N_2nd gaussian
      
      // par 1 row (n)
      psi2s_array[6] = fullmat[14];
      psi2s_array[7] = fullmat[15];
      psi2s_array[8] = fullmat[16];
      psi2s_array[9] = fullmat[17];
      psi2s_array[10] = fullmat[22];
      psi2s_array[11] = fullmat[19];

      // par 2 row (x_mean)
      psi2s_array[12] = fullmat[28];
      psi2s_array[13] = fullmat[29];
      psi2s_array[14] = fullmat[30];
      psi2s_array[15] = fullmat[31];
      psi2s_array[16] = fullmat[35];
      psi2s_array[17] = fullmat[32];

      // par 3 row (sigma)
      psi2s_array[18] = fullmat[42];
      psi2s_array[19] = fullmat[43];
      psi2s_array[20] = fullmat[44];
      psi2s_array[21] = fullmat[45];
      psi2s_array[22] = fullmat[50];
      psi2s_array[23] = fullmat[47];
     
      // par 8 row (N_psi2s)
      psi2s_array[24] = fullmat[112];
      psi2s_array[25] = fullmat[113];
      psi2s_array[26] = fullmat[114];
      psi2s_array[27] = fullmat[115];
      psi2s_array[28] = fullmat[120];
      psi2s_array[29] = fullmat[117];
   
      // par 5 row (N_2nd)
      psi2s_array[30] = fullmat[70];
      psi2s_array[31] = fullmat[71];
      psi2s_array[32] = fullmat[72];
      psi2s_array[33] = fullmat[73];
      psi2s_array[34] = fullmat[78];
      psi2s_array[34] = fullmat[75];
   
      //////////////////////////////////////

     
   
      TF1 *psi2s_fit = new TF1("psi2s_fit",CBcalc_psi2s,2.0, 5.0, 6);  // use CBcalc_psi2s instead
      psi2s_fit->SetLineColor(kBlue);
    
      psi2s_fit->SetParameter(0,total_fit->GetParameter(0));
      psi2s_fit->SetParError(0,total_fit->GetParError(0));
      psi2s_fit->SetParameter(1,total_fit->GetParameter(1));
      psi2s_fit->SetParError(1,total_fit->GetParError(1));
      psi2s_fit->SetParameter(2,total_fit->GetParameter(2) + total_fit->GetParameter(6));
      psi2s_fit->SetParError(2,total_fit->GetParError(2));  // par 6 is fixed
      psi2s_fit->SetParameter(3,total_fit->GetParameter(3)*total_fit->GetParameter(7));
      psi2s_fit->SetParError(3,total_fit->GetParError(3)*total_fit->GetParameter(7)); 
      psi2s_fit->SetParameter(4,total_fit->GetParameter(8));
      psi2s_fit->SetParError(4,total_fit->GetParError(8));
      psi2s_fit->SetParameter(5,total_fit->GetParameter(8)*total_fit->GetParameter(5)/total_fit->GetParameter(4));
      psi2s_fit->SetParError(5,total_fit->GetParError(8)*total_fit->GetParError(5)/total_fit->GetParError(4));
        
      
      
      Double_t psi2s_par[6] = {0}; // 5 CB parameters 

      psi2s_par[0] =  total_fit->GetParameter(0);     
      psi2s_par[1] =  total_fit->GetParameter(1);     
      psi2s_par[2] =  total_fit->GetParameter(2) + total_fit->GetParameter(6);         
      psi2s_par[3] =  total_fit->GetParameter(3)*total_fit->GetParameter(7);       
      psi2s_par[4] =  total_fit->GetParameter(8);          
      psi2s_par[5] =  total_fit->GetParameter(8)*total_fit->GetParameter(5)/total_fit->GetParameter(4);
  
      cout << "Check 6: " << endl;
      cout << "psi2s par 0: " <<  psi2s_par[0] << ", psi2s par 1: " <<  psi2s_par[1] << ", psi2s par 2: " <<  psi2s_par[2] << ", psi2s par 3: " <<  psi2s_par[3] << ", psi2s par 4: " <<  psi2s_par[4] << ", psi2s par 5: " <<  psi2s_par[5]  << endl;


      Npsi2s =  psi2s_fit->Integral(2.0, 5.0)/width; 
      err_psi2s = psi2s_fit->IntegralError(2.0,5.0,psi2s_par,psi2s_array)/width;
 
    
      recomassFG->SetTitle("");
      recomassFG->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (GeV/c^{2})    ");
      recomassFG->GetXaxis()->SetLabelSize(0.05);
      recomassFG->GetXaxis()->SetNdivisions(5,5,0);
      recomassFG->GetXaxis()->SetTitleOffset(1.11);  
      recomassFG->GetYaxis()->SetLabelSize(0.05);
      recomassFG->GetYaxis()->SetTitleOffset(1.3);  
      recomassFG->GetYaxis()->SetTitle("raw counts/(50 MeV/c^{2})");
      recomassFG->DrawCopy();

    }
  else
    {
      recomassFG->Fit(total_fit); 
      recomassFG->DrawCopy();
      fresult = recomassFG->GetFunction("CBcalc_exp");
    }
  
  if(fresult)
    {
      cout << endl << "From Fit:  Fit chisquare = " <<  fresult->GetChisquare()
	   << " NDF  = " <<   fresult->GetNDF() 
	   << " chisquare/NDF = " << fresult->GetChisquare() / fresult->GetNDF() << endl;
    }
  
  chisquare = fresult->GetChisquare();
  ndf =  fresult->GetNDF();

  
  double percent_error = err_psi2s/Npsi2s*100; 
  if(north_arm)
    cout << "North p + p jpsi counts " << Npsi2s << " +/- " << err_psi2s << " , sqrt error " << sqrt(Npsi2s) << " and % error " << percent_error << endl;
  //cout << "North p + p jpsi counts " << NJpsi << " , sqrt error " << sqrt(NJpsi) << " , error " << err_NJpsi << " and % error " << percent_error << endl;
  if(!north_arm)
    cout << "South p + p jpsi counts " << Npsi2s << " +/- " << err_psi2s << " , sqrt error " << sqrt(Npsi2s) << " and % error " << percent_error << endl;
    
  ///////////////////////////////////////////////////////////////////////////////END ORIGINAL CODE

 
  
  double alpha =  total_fit->GetParameter(0);
  double n =  total_fit->GetParameter(1);
  double CB_mean =  total_fit->GetParameter(2);
  double sigma =  total_fit->GetParameter(3);
  double N =  total_fit->GetParameter(4);


  ///////////////////////////////////////
  // char name3500[500]; 

  // if(chisquare_write)
  //   {
    
  //     if(north_arm)
  // 	sprintf(name3500,"CB_parameter_scan_vs_sigma_2nd/north/N_gaussian_10%i.txt", s);	 
  //     else
  // 	sprintf(name3500,"CB_parameter_scan_vs_sigma_2nd/south/S_gaussian_10%i.txt", s);	  	
     
  //     // whether north arm or south arm, write to unique filename  
  //     std::fstream Run15pp_width(name3500,std::ofstream::out); 
  //     Run15pp_width >>  NJpsi >> " " >>  err_NJpsi >> " " >> Npsi2s >> " " >> err_psi2s >> " " >> Npsi2s/NJpsi >> " " >> chisquare >> " " >> chisquare/ndf >> " " >> fixed >> " " >> CB_mean >> " " >> alpha >> " " >> n >> " " >> sigma >> " " >> N << " " << endl;
  //     Run15pp_width.close();
   

  //   } // if width_
// write
  ///////////////////////////////////////

  ///////////////////////////////////////
  char name200[500];
  //TFile *file2S; 

  if(width_write)
    {
      if(bin_width == 1)
	{
	  if(north_arm)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/jpsi_widths/Run15pp_N_width_%i.dat", bin);	 
	  else
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/jpsi_widths/Run15pp_S_width_%i.dat", bin);	  
	}
      if(bin_width == 2)
	{
	  if(north_arm && bin == 13)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/jpsi_widths/Run15pp_N_width_25.dat");
	  if(north_arm && bin == 14)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/jpsi_widths/Run15pp_N_width_26.dat");
	  if(!north_arm && bin == 13)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/jpsi_widths/Run15pp_S_width_25.dat");
	  if(!north_arm && bin == 14)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/jpsi_widths/Run15pp_S_width_26.dat");
	}
      
      // whether north arm or south arm, capture same information  
      a = total_fit->GetParameter(3);
      b = total_fit->GetParError(3);
    
      // whether north arm or south arm, write to unique filename  
      std::fstream Run15pp_width(name200,std::ofstream::out); 
      Run15pp_width <<  a << " " << b << " " << endl;
      Run15pp_width.close();
   

    } // if width_write
  ///////////////////////////////////////

  if(write_pAl_syst_err)
    {
      if(north_arm)
	{
	  if(bin < 25)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/pAl_syst_err_yields/Run15pp_N_yield_%i.dat", bin);
	  if(bin == 13 && bin_width == 2)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/pAl_syst_err_yields/Run15pp_N_yield_25.dat");
	  if(bin == 14 && bin_width == 2)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/pAl_syst_err_yields/Run15pp_N_yield_26.dat");
	}
      else
	{
	  if(bin < 25)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/pAl_syst_err_yields/Run15pp_S_yield_%i.dat", bin);	
	  if(bin == 13 && bin_width == 2)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/pAl_syst_err_yields/Run15pp_S_yield_25.dat");
	  if(bin == 14 && bin_width == 2)
	    sprintf(name200,"tony_bestfit_parameters/sanghoon_corrbg/pAl_syst_err_yields/Run15pp_S_yield_26.dat");
	}  
    

      // whether north arm or south arm, write to unique filename  
      std::fstream Run15pp_yield(name200,std::ofstream::out); 
      Run15pp_yield <<  NJpsi << " " << err_NJpsi << " " << endl;
      Run15pp_yield.close();
   

    } // if pAl_write
  ///////////////////////////////////////
  


  TF1 *jpsi_fit = new TF1("jpsi_fit",CBcalc_jpsi_G,2.0,5.0,8); // only total fit fucntion has correct normalzation for 2nd gaussian
  for(int i=0;i<6;i++)
    {
      jpsi_fit->SetParameter(i,total_fit->GetParameter(i));
      jpsi_fit->SetParError(i,total_fit->GetParError(i));
    }

  jpsi_fit->FixParameter(6,total_fit->GetParameter(2));
  jpsi_fit->FixParameter(7,fixed);
  
  jpsi_fit->SetLineColor(kMagenta);
  jpsi_fit->Draw("same");


  TF1 *gauss_fit = new TF1("gauss_fit",Gauss2,2.0,5.0,3); // gauss2  J/psi
  gauss_fit->SetParameter(0,total_fit->GetParameter(5));
  gauss_fit->SetParError(0,total_fit->GetParError(5));
  gauss_fit->FixParameter(1,total_fit->GetParameter(2));
  gauss_fit->FixParameter(2,fixed);
  
  gauss_fit->SetLineColor(kBlack);
  gauss_fit->SetLineStyle(10);
  gauss_fit->Draw("same");


  TF1 *psi2s_fit = new TF1("psi2s_fit",CBcalc_psi2s_G,2.0,5.0,8);  // CBcalc_psi2s
  psi2s_fit->SetParameter(0,total_fit->GetParameter(0));
  psi2s_fit->SetParError(0,total_fit->GetParError(0));
  psi2s_fit->SetParameter(1,total_fit->GetParameter(1));
  psi2s_fit->SetParError(1,total_fit->GetParError(1));
  psi2s_fit->SetParameter(2,total_fit->GetParameter(2) + total_fit->GetParameter(6));
  psi2s_fit->SetParError(2,total_fit->GetParError(2));  // par 6 is fixed
  psi2s_fit->SetParameter(3,total_fit->GetParameter(3)*total_fit->GetParameter(7));
  psi2s_fit->SetParError(3,total_fit->GetParError(3)*total_fit->GetParameter(7)); 
  psi2s_fit->SetParameter(4,total_fit->GetParameter(8));
  psi2s_fit->SetParError(4,total_fit->GetParError(8));
  psi2s_fit->SetParameter(5,total_fit->GetParameter(8)*total_fit->GetParameter(5)/total_fit->GetParameter(4));
  psi2s_fit->SetParError(5,total_fit->GetParError(8)*total_fit->GetParError(5)/total_fit->GetParError(4));
      
  psi2s_fit->FixParameter(6,total_fit->GetParameter(2) + total_fit->GetParameter(6));
  psi2s_fit->FixParameter(7,fixed);

  psi2s_fit->SetLineColor(kMagenta);
  psi2s_fit->Draw("same");

 cout << "Check 7: " << endl;
 cout << "psi2s par 0: " <<  psi2s_fit->GetParameter(0) << ", psi2s par 1: " <<  psi2s_fit->GetParameter(1) << ", psi2s par 2: " << psi2s_fit->GetParameter(2)<< ", psi2s par 3: " <<  psi2s_fit->GetParameter(3) << ", psi2s par 4: " << psi2s_fit->GetParameter(4) << ", psi2s par 5: " <<  psi2s_fit->GetParameter(5)  << ", psi2s par 6: " <<  psi2s_fit->GetParameter(6)  << ", psi2s par 7: " <<  psi2s_fit->GetParameter(7)  << endl;


  TF1 *gauss2_fit = new TF1("gauss2_fit",Gauss2,2.0,5.0,3); // gauss2_psi2s
  gauss2_fit->SetParameter(0,total_fit->GetParameter(8)*total_fit->GetParameter(5)/total_fit->GetParameter(4));
  gauss2_fit->SetParError(0,total_fit->GetParError(8)*total_fit->GetParError(5)/total_fit->GetParError(4));
  gauss2_fit->FixParameter(1,total_fit->GetParameter(2) + total_fit->GetParameter(6));
  gauss2_fit->FixParameter(2,fixed);
			   
  gauss2_fit->SetLineColor(kViolet+2);
  gauss2_fit->SetLineStyle(10);
  gauss2_fit->Draw("same");

cout << "Check 8: " << endl;
 cout << "gauss2 par 0: " <<  gauss2_fit->GetParameter(0) << ", gauss2 par 1: " <<  gauss2_fit->GetParameter(1) << ", gauss2 par 2: " << gauss2_fit->GetParameter(2)<< ", psi2s par 3: " <<  endl;

  
  TF1 *corrbg_fit = new TF1("corrbg_fit","[2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4])",2.0,5.0);
  corrbg_fit->SetParameter(0,total_fit->GetParameter(9));
  corrbg_fit->SetParameter(1,total_fit->GetParameter(10));
  corrbg_fit->SetParameter(2,total_fit->GetParameter(11));
  corrbg_fit->SetParameter(3,total_fit->GetParameter(12));
  corrbg_fit->SetParameter(4,total_fit->GetParameter(13));
 
  double corr_bg = corrbg_fit->Eval(3.1);
  double comb_bg = combg_fit->Eval(3.1);

  for(int i = 0; i < 100; i++)
    {
      mass_array_corrbg[i] = corrbg_fit->Eval(x_array[i]);
    }

  
  corrbg_fit->SetLineColor(kGreen);
  corrbg_fit->SetLineWidth(3);
  corrbg_fit->SetLineStyle(kDashed);
  corrbg_fit->Draw("same");

////////////////////////////////////////////////////////////
  //  LS background
  if(LL_fit)
    {
      //recomassBG->SetLineColor(kRed); 
      // recomassBG->DrawCopy("same");
    }
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  //  mixed event background
  if(LL_fit)
    {
      recomassBG3->SetLineColor(kRed); 
      recomassBG3->DrawCopy("same");
    }
  ////////////////////////////////////////////////////////////

TLegend *leg2 = new TLegend(0.65, 0.52, 0.84, 0.74);  //(start x, start y, end x, end y)
  leg2->SetFillColor(19); 
  leg2->SetFillStyle(3003); 
  leg2->SetLineWidth(1);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.03); 
  //leg2->SetTextFont(102);

  leg2->AddEntry(corrbg_fit,"YueHang Corr BG fit", "l");
  leg2->AddEntry(combg_fit,"Mixed events BG fit", "l");
  leg2->AddEntry(total_fit,"Total fit", "l");
  leg2->AddEntry(psi2s_fit,"Crystal Ball Fit", "l");
  leg2->AddEntry(gauss2_fit,"2nd Guassian Fit", "l");
  leg2->Draw("same");



  if(write == true)
    {
      if(north_arm)
	{
	  if(bin_width == 1 && special == 0)
	    {
	      if(CaseF == true)
		sprintf(nameN,"tony_bestfit_parameters/CaseF/Run15pp_N_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(CaseG == true)
		sprintf(nameN,"tony_bestfit_parameters/CaseG/Run15pp_N_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(CaseH == true)
		sprintf(nameN,"tony_bestfit_parameters/CaseH/Run15pp_N_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(CaseA == true)
		{
		  sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/pp_yields/CaseA/Run15pp_N_bestfit_parameters_tony_1_%i.dat",bin); 
		  sprintf(nameN2,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_N_bestfit_parameters_tony_1_%i.dat",bin); 
		  // if(bin < 17)
		  // sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/pAl_yields/Run15pp_N_bestfit_parameters_tony_1_%i.dat",bin); 
		  // if(bin < 17)
		  // sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_N_bestfit_parameters_tony_1_%i.dat",bin); 
		}
	    }
	  if(bin_width == 2  && special == 0)
	    {
	      if(CaseF == true)
		sprintf(nameN,"tony_bestfit_parameters/CaseF/Run15pp_N_bestfit_parameters_tony_2_%i.dat",bin); 
	      if(CaseG == true)
		sprintf(nameN,"tony_bestfit_parameters/CaseG/Run15pp_N_bestfit_parameters_tony_2_%i.dat",bin); 
	      if(CaseH == true)
		sprintf(nameN,"tony_bestfit_parameters/CaseH/Run15pp_N_bestfit_parameters_tony_2_%i.dat",bin); 
	      if(CaseA == true)
		{
		  if(bin == 13)
		    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/pp_yields/CaseA/Run15pp_N_bestfit_parameters_tony_1_25.dat"); 
		  if(bin == 14)
		    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/pp_yields/CaseA/Run15pp_N_bestfit_parameters_tony_1_26.dat"); 
		  // for centrality binning
		  // if(bin == 6 && HeAu)
		  //   sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_N_bestfit_parameters_tony_1_11.dat"); 
		  // if(bin == 7 && HeAu)
		  //   sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_N_bestfit_parameters_tony_1_12.dat"); 
		  // if(bin == 8 && HeAu)
		  //   sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_N_bestfit_parameters_tony_1_13.dat"); 
		  if(bin == 9 && pAu)
		    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_N_bestfit_parameters_tony_1_17.dat"); 
		  if(bin == 10 && pAu)
		    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_N_bestfit_parameters_tony_1_18.dat"); 
		}
	    }
	  if(bin_width == 4 && pAl)
	    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/0100_pAl_yields/Run15pp_N_bestfit_parameters_tony_4_%i.dat",bin); 
	  if(bin == 3 && special == 1 && pAu)
	    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_N_bestfit_parameters_tony_1_19.dat"); 
	  if(bin == 3 && special == 1 && HeAu)
	    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/0100_HeAu_yields/Run15pp_N_bestfit_parameters_tony_1_19.dat"); 
	 if(bin == 5 && bin_width == 4 && HeAu)
	    sprintf(nameN,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_N_bestfit_parameters_tony_1_14.dat"); 
	  
	  std::fstream Run15pp_N_bestfit_parameters_tony(nameN,std::ofstream::out); 
	  Run15pp_N_bestfit_parameters_tony <<  f0  <<  " " << f1 << " " <<  f2 <<  " " << f3 << " " << f4 << " " <<  f5 << " " << f6 << " "  <<  f7 << " " << f8 << " " << f9 << " " << NJpsi  << " " << err_NJpsi << " " << chisquare << " " << ndf << " " << endl;
	  Run15pp_N_bestfit_parameters_tony.close();
	}
      else
	{
	  if(bin_width == 1  && special == 0)
	    {
	      if(CaseF == true)
		sprintf(nameS,"tony_bestfit_parameters/CaseF/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(CaseG == true)
		sprintf(nameS,"tony_bestfit_parameters/CaseG/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(CaseH == true)
		sprintf(nameS,"tony_bestfit_parameters/CaseH/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(CaseA == true)
		{
		  sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/pp_yields/CaseA/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
		  sprintf(nameS2,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
		  // if(bin < 17)
		  //  sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/pAl_yields/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
		  // if(bin < 17)
		  //  sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
		}
	      if(bin == 29 && CaseF)
		sprintf(nameS,"tony_bestfit_parameters/pt_int/CaseF/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(bin == 29 && CaseG)
		sprintf(nameS,"tony_bestfit_parameters/pt_int/CaseG/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
	      if(bin == 29 && CaseH)
		sprintf(nameS,"tony_bestfit_parameters/pt_int/CaseH/Run15pp_S_bestfit_parameters_tony_1_%i.dat",bin); 
	    }
	  
	  if(bin_width == 2  && special == 0)
	    {
	      if(CaseF == true)
		sprintf(nameS,"tony_bestfit_parameters/CaseF/Run15pp_S_bestfit_parameters_tony_2_%i.dat",bin); 
	      if(CaseG == true)
		sprintf(nameS,"tony_bestfit_parameters/CaseG/Run15pp_S_bestfit_parameters_tony_2_%i.dat",bin); 
	      if(CaseH == true)
		sprintf(nameS,"tony_bestfit_parameters/CaseH/Run15pp_S_bestfit_parameters_tony_2_%i.dat",bin); 
	      if(CaseA == true)
		{
		  if(bin == 13)
		    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/pp_yields/CaseA/Run15pp_S_bestfit_parameters_tony_1_25.dat"); 
		  if(bin == 14)
		    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/pp_yields/CaseA/Run15pp_S_bestfit_parameters_tony_1_26.dat"); 
		  // for centrality binning
		  // if(bin == 6 && HeAu)
		  //   sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_S_bestfit_parameters_tony_1_11.dat"); 
		  // if(bin == 7 && HeAu)
		  //   sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_S_bestfit_parameters_tony_1_12.dat"); 
		  // if(bin == 8 && HeAu)
		  //   sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_S_bestfit_parameters_tony_1_13.dat"); 
		  if(bin == 9 && pAu)
		    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_S_bestfit_parameters_tony_1_17.dat"); 
		  if(bin == 10 && pAu)
		    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_S_bestfit_parameters_tony_1_18.dat"); 
		}
	    }
	  
	  if(bin_width == 4)
	    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/0100_pAl_yields/Run15pp_S_bestfit_parameters_tony_4_%i.dat",bin); 
	  if(bin == 3 && special == 1 && pAu)
	    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/pAu_yields/Run15pp_S_bestfit_parameters_tony_1_19.dat"); 
	  if(bin == 3 && special == 1 && HeAu)
	    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/0100_HeAu_yields/Run15pp_S_bestfit_parameters_tony_1_19.dat"); 
	  if(bin == 5 && bin_width == 4 && HeAu)
	    sprintf(nameS,"tony_bestfit_parameters/sanghoon_corrbg/HeAu_yields/Run15pp_S_bestfit_parameters_tony_1_14.dat"); 

	  std::fstream Run15pp_S_bestfit_parameters_tony(nameS,std::ofstream::out); 
	  Run15pp_S_bestfit_parameters_tony <<  f0  <<  " " << f1 << " " <<  f2 <<  " " << f3 << " " << f4 << " " <<  f5 << " " << f6 << " "  <<  f7 << " " << f8 << " " << f9 << " " << NJpsi  << " " << err_NJpsi << " " << chisquare << " " << ndf << " " << endl;  //NJpsi is covariant
	  Run15pp_S_bestfit_parameters_tony.close();
	}
    }
  
  TLatex l2;
  l2.SetTextSize(0.04);
  l2.SetTextAlign(13);
  l2.SetTextColor(1);
  
  char text3[100];
  // if(bin_width == 1)
  //   {
  //     if(north_arm == true)
  // 	sprintf(text3,"North bin_%i,  %.2f - %.2f GeV/c", bin, pt_initial, pt_final);
  //     if(north_arm == false)
  // 	sprintf(text3,"South bin_%i,  %.2f - %.2f GeV/c", bin, pt_initial, pt_final);
  //   }
  if(bin_width == 2 && bin == 13)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_25, 6.0 - 6.5 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_25, 6.0 - 6.5 GeV");
    }
  if(bin_width == 2 && bin == 14)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_26, 6.5 - 7.0 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_26, 6.5 - 7.0 GeV");
    }
  if(bin_width == 4 && bin == 8)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_27, 7.0 - 8.0 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_27, 7.0 - 8.0 GeV");
    }
  
  if(bin_width == 8 && bin == 5)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_28, 8.0 - 10.0 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_28, 8.0 - 10.0 GeV");
    }
  ////////////////////////////
  
  if(bin_width == 4 && bin == 1)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_1, 0.0 - 1.0 GeV/c");
      if(north_arm == false)
	sprintf(text3,"South bin_1, 0.0 - 1.0 GeV/c");
    }
  if(bin_width == 4 && bin == 2)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_2, 1.0 - 2.0 GeV/c");
      if(north_arm == false)
	sprintf(text3,"South bin_2, 1.0 - 2.0 GeV/c");
    }
  if(bin_width == 4 && bin == 3)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_3, 2.0 - 3.0 GeV/c");
      if(north_arm == false)
	sprintf(text3,"South bin_3, 2.0 - 3.0 GeV/c");
    }
  if(bin_width == 4 && bin == 4)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_4, 3.0 - 4.0 GeV/c");
      if(north_arm == false)
	sprintf(text3,"South bin_4, 3.0 - 4.0 GeV/c");
    }
 if(bin_width == 8 && bin == 3)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_5, 4.0 - 6.0 GeV/c");
      if(north_arm == false)
	sprintf(text3,"South bin_5, 4.0 - 6.0 GeV/c");
    }
 if(bin_width == 8 && bin == 5)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_5, 4.0 - 6.0 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_5, 4.0 - 6.0 GeV");
    }
  if(bin_width == 2 && bin == 10)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_18, 4.5 - 5.0 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_18, 4.5 - 5.0 GeV");
    }
 if(pt_bin == 2 && bin == 3)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_19, 5.0 - 7.0 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_19, 5.0 - 7.0 GeV");
    }
 if(pt_bin == 2 && bin == 4)
    {
      if(north_arm == true)
	sprintf(text3,"North bin_21, 7.0 - 9.0 GeV");
      if(north_arm == false)
	sprintf(text3,"South bin_21, 7.0 - 9.0 GeV");
    }
 if(bin == 48)
   {
     if(north_arm == true)
       sprintf(text3,"North p_{T} Int ");
     if(north_arm == false)
       sprintf(text3,"South p_{T} Int ");
    }



 if(bin == 49)
    {
      if(north_arm == true)
	//	sprintf(text3,"North p_{T} Int. (Local), FVTX");
      sprintf(text3,"North p_{T} Integrated");
      if(north_arm == false)
	//	sprintf(text3,"South p_{T} Int. (Local), FVTX");
      sprintf(text3,"South p_{T} Integrated ");
    }
 if(bin == 1)
    {
      if(north_arm == true)
	//	sprintf(text3,"North p_{T} Int. (Local), FVTX");
      sprintf(text3,"North p_{T} Int ");
      if(north_arm == false)
	//	sprintf(text3,"South p_{T} Int. (Local), FVTX");
      sprintf(text3,"South p_{T} Int ");
    }


 l2.SetTextSize(0.035);
  l2.SetTextAlign(12);
  // l2.DrawLatexNDC(0.70, 0.48, text3); //4.4,150   //0.65,0.84




  double ratio = Npsi2s/NJpsi;
  double free=fixed;


  TLatex l3;
  l3.SetTextSize(0.06);
  l3.SetTextAlign(13);
  l3.SetTextColor(4);

  char text4[100];
  if(fvtx)
    sprintf(text4,"p + Au, mutr+mutr");
  if(fvtx_sngtrk)
    sprintf(text4,"p + Au, dbl+sngtrk");
  if(sngtrk)
    sprintf(text4,"p + Au, fvtx+mutr");
  l3.SetTextAlign(12);
  l3.DrawLatexNDC(0.56, 0.93, text4); //4.4,150
 

  char text5[100];
  l3.SetTextColor(1);
  l3.SetTextSize(0.035);
  if(north_arm)
    sprintf(text5,"north p_{T} Int");
  else
    sprintf(text5,"South p_{T} Int");
  l3.SetTextAlign(12);
   l3.DrawLatexNDC(0.13, 0.79, text5); //4.4,150
  
  char text6[100];
  l3.SetTextColor(1);
  l3.SetTextSize(0.035);
  sprintf(text6,"ratio =  %.4f ",ratio);
  l3.SetTextAlign(12);
  l3.DrawLatexNDC(0.13, 0.84, text6);

char text7[100];
  l3.SetTextColor(1);
  l3.SetTextSize(0.035);
  sprintf(text7,"Inclusive mix evt");
  l3.SetTextAlign(12);
  //  l3.DrawLatexNDC(0.13, 0.84, text7);
 

  Char_t message[200];
  Char_t message2[200];
  Char_t message3[200];
  Char_t message4[200];
 

  sprintf(message,"#chi^{2}/NDF = %.1f / %.d",fresult->GetChisquare(),fresult->GetNDF());
   sprintf(message2,"Psi(2S) Counts =  %.1f +/- %.2f",Npsi2s,err_psi2s);
  sprintf(message3,"Prob =  %.6f",fresult->GetProb());
  sprintf(message4,"J/Psi Counts =  %.1f +/- %.2f",NJpsi,err_NJpsi);
   
  TPaveText *mytext = new TPaveText(0.47,0.88,0.89,0.76,"NDC"); // x0,y0,x1,y1
  mytext->SetTextSize(0.035);
  mytext->SetFillColor(0);
  mytext->SetTextAlign(12);
  mytext->AddText(message);
   mytext->AddText(message2);
  // mytext->AddText(message3);
  mytext->AddText(message4);
  mytext->Draw();


  cout << "f13: " << f13 << endl;

  if(CaseF)
    cout << "   ------->   Case F" << endl;
  if(CaseG)
    cout << "   ------->   Case G" << endl;
  if(CaseH)
    cout << "   ------->   Case H" << endl;
  if(CaseA)
    cout << "   ------->   Case A" << endl;



   int bins = 100;

   // for corrbg
 if(print_screen == true)
    {
      cout << "mass_array_corrbg[" << bins << "] = {";
      for(int i = 0;i < bins-1; i++)
      	{
      	  cout << mass_array_corrbg[i] << ", ";
      	}
      cout << mass_array_corrbg[bins-1] << "};" << endl;
      cout << "_______________________________ " << endl;
     
    } // print out

 

  
 TH1D *pp_MB_S_combg = 0;
 TH1D *pp_MB_S_corrbg = 0;
 TH1D *pp_MB_S_jpsi = 0;
 TH1D *pp_MB_S_psi2s = 0;
 TH1D *pp_MB_S_total = 0;

 TH1F *pp_MB_S_combg_fit = 0;
 TH1F *pp_MB_S_corrbg_fit = 0;
 TH1F *pp_MB_S_jpsi_fit = 0;
 TH1F *pp_MB_S_psi2s_fit = 0;
 TH1F *pp_MB_S_total_fit = 0;

 // Rpp 0-100 and centralities
 pp_MB_S_combg = new TH1D("pp_pT_S_MB_combg","pp_pT_S_MB_combg",100,0,5); 
 pp_MB_S_combg->GetXaxis()->SetRangeUser(2,5);
 pp_MB_S_corrbg = new TH1D("pp_pT_S_MB_corrbg","pp_pT_S_MB_corrbg",bins,2,5);
 pp_MB_S_corrbg->GetXaxis()->SetRangeUser(2,5);
 pp_MB_S_jpsi = new TH1D("pp_pT_S_MB_jpsi","pp_pT_S_MB_jpsi",60,2,5);
 pp_MB_S_jpsi->GetXaxis()->SetRangeUser(2,5);
 pp_MB_S_psi2s = new TH1D("pp_pT_S_MB_psi2s","pp_pT_S_MB_psi2s",bins,2,5);
 pp_MB_S_psi2s->GetXaxis()->SetRangeUser(2,5);
 pp_MB_S_total = new TH1D("pp_pT_S_MB_total","pp_pT_S_MB_total",bins,2,5);
pp_MB_S_total->GetXaxis()->SetRangeUser(2,5);

 pp_MB_S_combg_fit = new TH1F("pp_pT_S_MB_combg_fit"," [2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4]",bins,2,5);
 pp_MB_S_corrbg_fit = new TH1F("pp_pT_S_MB_corrbg_fit","[2]/ pow( ((exp(-[0]*x - [1]*x*x)) + x/[3]), [4]",bins,2,5);
 pp_MB_S_jpsi_fit = new TH1F("pp_pT_S_MB_jpsi_fit","CB Function",bins,2,5);
 pp_MB_S_psi2s_fit = new TH1F("pp_pT_S_MB_psi2s_fit","CB Function",bins,2,5);
 pp_MB_S_total_fit = new TH1F("pp_pT_S_MB_total_fit","Total",bins,2,5);
 
  // set bin content for 0-100 plots
  for(int i = 0; i < 100;i++)
    { 
      pp_MB_S_combg->SetBinContent(i+1,mass_array_combg[i]); 
      pp_MB_S_combg->SetBinError(i+1,mass_array_combg_err[i]); 
      pp_MB_S_jpsi->SetBinContent(i+1,mass_array_jpsi[i]);
      pp_MB_S_jpsi->SetBinError(i+1,mass_array_jpsi_err[i]);
      pp_MB_S_psi2s->SetBinContent(i+1,mass_array_psi2s[i]);
      pp_MB_S_psi2s->SetBinError(i+1,mass_array_psi2s_err[i]);
      pp_MB_S_total->SetBinContent(i+1,mass_array_total[i]);
      pp_MB_S_total->SetBinError(i+1,mass_array_total_err[i]);
    }

  
 // write histograms with fits to root file

  TFile *h = new TFile("pp_MB_S_fit_hist_fx.root", "RECREATE");
  h->cd();
  
  pp_MB_S_combg->Write();
  pp_MB_S_jpsi->Write();
 
  combg_fit->Write();
  corrbg_fit->Write();
  jpsi_fit->Write();
  psi2s_fit->Write();
  total_fit->Write();
 
  

  if(sqrt(NJpsi) < err_NJpsi)
    cout << " --> 5. J/psi error consistent" << endl;



  if(sqrt(Npsi2s) < err_psi2s)
    cout << " --> 6. psi2S error consistent" << endl;





  }




