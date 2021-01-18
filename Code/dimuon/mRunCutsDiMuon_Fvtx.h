class TTree;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class TFile;
class TChain;
class TProfile;

#include <set>

class mRunCutsDiMuon_Fvtx
{ 
 public:

  //! default constructor
  mRunCutsDiMuon_Fvtx();
  
  //! destructor
  virtual 
  ~mRunCutsDiMuon_Fvtx();

  //! global initialization
  int Init(const char *listname);
  
  //! event method
  int process_event();
  
  //! global termination
  int End();

	//! set dataset
	void set_dataset(const char *dataset);

	//! set filetype
	void set_filetype(const char *filetype);

	//! set triggertype
	void set_triggertype(const char *triggertype);

	//! set basic cuts 
	void set_basic_cuts();

	//! set pT array 
	void set_pT_array();

	//! book histograms
	void book_histos();

	//! load analysis cuts
	void load_cuts();

	//! load femtoDSTs 
	void load_files(const char *listname);

	//! check run quality 
	bool check_run(int arm);

	//! check event quality 
	bool check_event();

	//! check dimuon quality 
	bool check_dimuon();

	//! check fvtx quality 
	bool check_fvtx(int tr_index);

	//! check mutr fiducial
	bool check_mutr_fiducial(int index);

	//! mass cut 
	void set_mass_cut(float min, float max);

	//! set swap 
	void use_fvtx(bool use);

	//! fill cut histograms 
	//void fill_cut_histos(int index, float weight);

	void get_run_range(std::string dataset="Run12pp200", int *runrange=NULL){

		if ( dataset=="Run9pp500" ){
			runrange[0] = 276000;
			runrange[1] = 281000;
		}else if ( dataset=="Run9pp200" ){
			runrange[0] = 281000;
			runrange[1] = 292000;
		}else if ( dataset=="Run12pp510" ){
			runrange[0] = 364500;
			runrange[1] = 368800;
		}else if ( dataset=="Run14AuAu200" ){
			runrange[0] = 405800;
			runrange[1] = 415000;
		}else if ( dataset=="Run15pp200" ){
			runrange[0] = 422000;
			runrange[1] = 432100;
		}else if ( dataset=="Run15pAl200" ){
			runrange[0] = 436700;
			runrange[1] = 438500;
		}else if ( dataset=="Run15pAu200" ){
			runrange[0] = 432600;
			runrange[1] = 436100;
		}else if ( dataset=="Run14HeAu200" ){
			runrange[0] = 415500;
			runrange[1] = 417000;
		}

		std::cout << "DATASET : " << dataset << ", RANGE " << runrange[0] << " - " << runrange[1] << std::endl;

		return;
	}

 private:

	//bool first;
	bool _is_sim;
	bool _use_fvtx;

	TFile *file_out;

	TChain *fChain;

	std::string _dataset;
	std::string _filetype;
	std::string _triggertype;

	std::vector<float> vec_dca_r;
	std::set<float> set_dca_r;
	int _prev_runnumber;
	float _prev_bbcZ;

	//Jpsi info.
	int _ndimuon, _same_event, _charge, _arm;
	float _mass, _mass_fvtx, _mass_fvtxmutr, _rapidity, _rapidity_fvtxmutr;
	float _px, _py, _pz, _pT;
	float _pT_fvtxmutr, _pz_fvtxmutr;
	float _phi_dca, _dca_r, _dca_z;
	float _evt_vtxchi2, _evt_vtxchi2_fvtxmutr;

	//dimu info.
	int _tr_ntrhits[2], _tr_nidhits[2], _tr_lastgap[2];
	int _tr_trhits[2], _tr_idhits[2];
	float _tr_trchi2[2], _tr_idchi2[2], _tr_DG0[2], _tr_DDG0[2];
	float _tr_px[2], _tr_py[2], _tr_pz[2], _tr_pT[2], _tr_rapidity[2];
	float _tr_st1px[2], _tr_st1py[2], _tr_st1pz[2];
	float _tr_xst1[2], _tr_xst2[2], _tr_xst3[2];
	float _tr_yst1[2], _tr_yst2[2], _tr_yst3[2];
	float _tr_dca_z[2], _tr_dca_r[2], _tr_dca_phi[2];
	float _tr_dca_r_fvtx[2], _tr_dca_r_smear[2], _tr_dca_r_sim[2];
	float _tr_dca_r_refit[2];
	int _tr_mc_hits_mutr_true[2], _tr_mc_hits_muid_true[2];

	//fvtx info.
	int _tr_nhits_fvtx[2], _tr_nhits_fvtx_charge[2];
	int _tr_hit_pattern[2];
	float _tr_dphi_fvtx[2], _tr_dtheta_fvtx[2], _tr_dr_fvtx[2];
	float _tr_chi2_fvtx[2], _tr_chi2_fvtxmutr[2];
	float _tr_x0_fvtxmutr[2], _tr_y0_fvtxmutr[2], _tr_z0_fvtxmutr[2];
	float _tr_px_fvtxmutr[2], _tr_py_fvtxmutr[2], _tr_pz_fvtxmutr[2];
	float _tr_x0_fvtx[2], _tr_y0_fvtx[2], _tr_z0_fvtx[2];
	float _tr_px_fvtx[2], _tr_py_fvtx[2], _tr_pz_fvtx[2];
	float _tr_charge_fvtx[2][4];

	//etc.
	float _tr_phist1[2], _tr_phist2[2], _tr_phist3[2];
	float _tr_radst1[2], _tr_radst2[2], _tr_radst3[2];
	int _tr_halfoctst1[2], _tr_halfoctst2[2], _tr_halfoctst3[2];
	float _tr_pdtheta[2];
	int _tr_vtx_index[2];

	//event info.
	float _bbcZ, _bbcqn, _bbcqs;
	float _fvtxX, _fvtxY, _fvtxZ;
	float _fvtxX2, _fvtxY2, _fvtxZ2;
	float _fvtxX_err, _fvtxY_err, _fvtxZ_err;
	int _mult_fvtxS, _mult_fvtxN, _mult_svx;
	int _mult_fvtx_dimu;
	int _centrality,  _runnumber;
	int _mult_fvtx_prim_cut;

	float _vtx_res_x, _vtx_res_y, _vtx_res_z;
	float _doubleint_frac;

	//trigger info.
	int _scaled_trigbit;
	int _live_trigbit;
	int _raw_trigbit;

	int _trig_emul_2D_S, _trig_emul_2D_N;

	//mc info.
	int _tr_mis_match[2], _tr_mc_g_pid[2];
	int _mc_pid, _mc_g_pid;
	float _mc_g_px, _mc_g_py, _mc_g_pT, _mc_g_pz;
	float _mc_g_z, _mc_z;
	float _simZ;

	int _tr_mc_hits_fvtx[2], _tr_mc_hits_svx[2];
	int _tr_mc_hits_fvtx_true[2], _tr_mc_hits_svx_true[2];

	int _nhepmc;
	int _hepmc_g_pid[10];
	float _hepmc_g_pT[10], _hepmc_g_rapidity[10];

	//MuTR location
	float _zst1[2];

	//constants
	int _nptbin;
	int _nphibin;
	int _ncentbin;
	enum { _narm = 2};
	enum { _ncharge = 3};
	enum { _nmuidgap = 5};
	enum { _nmuidpanel = 6};
	enum { _nmutrsta = 3};
	double *_pT_array;
	double *_p_array;

	//analysis cuts
	float cut_eta[2];
	float cut_z[_narm][2];
	float cut_pz[_narm];
	float cut_mass[2];
	float cut_chi2_mutr;
	float cut_chi2_muid;
	float cut_chi2_evtvtx;
	float cut_pdtheta;
	float cut_ptot;

	float dca_shift_sim_vtx_2fvtx[_narm];
	float dca_shift_sim_vtx_3fvtx[_narm];
	float dca_shift_sim_vtx_4fvtx[_narm];
	float dca_shift_sim_fvtx[_narm];

	float dca_shift_data_vtx[_narm];
	float dca_shift_data_fvtx[_narm];

	TF1* fcut_dr_fvtx[_narm];
	TF1* fcut_dphi_fvtx[_narm];
	TF1* fcut_dtheta_fvtx[_narm];
	TF1* fcut_dg0[_narm];
	TF1* fcut_ddg0[_narm];
	TF1* fcut_dca_phi;
	TFile *cut_file;

	TFile *badrun_file;
	TH1F *hbadrun[_narm];
	int _prev_run;
	bool _is_goodrun[2];

	TFile *rate_file;
	TH1F *hBBC_rate;

	//cut histograms
	
	//QA histograms
	TH1F *N_RUN_LOOSE_SAME[_narm][_ncharge];
	TH1F *QA_MUTR_CHI2[_narm][_ncharge];
	TH1F *QA_MUID_CHI2[_narm][_ncharge];
	TH1F *QA_EVTVTX_CHI2[_narm][_ncharge];
	TH1F *QA_FVTXMUTR_CHI2[_narm][_ncharge];
	TH1F *QA_MUTR_NHIT[_narm][_ncharge];
	TH1F *QA_DG0[_narm][_ncharge];
	TH1F *QA_DDG0[_narm][_ncharge];

	TH1F *QA_MUID_DIFF_GAP[_narm];
	TH1F *QA_MUID_DIFF_NHIT[_narm];
	TH2F *QA_MUID_DEPTH[_narm];
	TH2F *QA_MUID_NHIT[_narm];
	TH1F *QA_MUID_NHIT_1D[_narm];

	TH3F *MASS_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_2MATCH_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_1MATCH_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_2MATCH_MISMATCH_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_1MATCH_MISMATCH_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_FVTX_2MATCH_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_FVTX_1MATCH_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_FVTX_2MATCH_MISMATCH_Y_PT_SAME[_narm][_ncharge];
	TH3F *MASS_FVTX_1MATCH_MISMATCH_Y_PT_SAME[_narm][_ncharge];

	TH3F *MASS_FVTX_2MATCH_Y_PT_SAME_CENT[10][_narm][_ncharge];
	TH3F *MASS_FVTX_1MATCH_Y_PT_SAME_CENT[10][_narm][_ncharge];
	TH3F *MASS_FVTX_2MATCH_MISMATCH_Y_PT_SAME_CENT[10][_narm][_ncharge];
	TH3F *MASS_FVTX_1MATCH_MISMATCH_Y_PT_SAME_CENT[10][_narm][_ncharge];

	//
	TH2F *MASS_PT_SAME[_narm][_ncharge];
	TH2F *MASS_PT_MIXED[_narm][_ncharge];

	TH2F *MASS_FVTXMUTR_PT_SAME[_narm][_ncharge];
	TH2F *MASS_FVTXMUTR_PT_MIXED[_narm][_ncharge];

	TH2F *MASS_HIGH_PT_SAME[_narm][_ncharge];
	TH2F *MASS_HIGH_PT_MIXED[_narm][_ncharge];

	TH2F *MASS_Y_SAME[_narm][_ncharge];
	TH2F *MASS_Y_MIXED[_narm][_ncharge];

	TH2F *MASS_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Z_MIXED[_narm][_ncharge];

	TH2F *MASS_Y_AND_SAME[_narm][_ncharge];
	TH2F *MASS_Y_AND_MIXED[_narm][_ncharge];
	TH2F *MASS_Y_OR_SAME[_narm][_ncharge];
	TH2F *MASS_Y_OR_MIXED[_narm][_ncharge];

	TH2F *MASS_Z_AND_SAME[_narm][_ncharge];
	TH2F *MASS_Z_AND_MIXED[_narm][_ncharge];
	TH2F *MASS_Z_OR_SAME[_narm][_ncharge];
	TH2F *MASS_Z_OR_MIXED[_narm][_ncharge];

	TH2F *MASS_CENTS_SAME[_narm][_ncharge];
	TH2F *MASS_CENTS_MIXED[_narm][_ncharge];
	TH2F *MASS_CENTN_SAME[_narm][_ncharge];
	TH2F *MASS_CENTN_MIXED[_narm][_ncharge];

	TH2F *MASS_YY_SAME[_narm][_ncharge];
	TH2F *MASS_YY_MIXED[_narm][_ncharge];
	TH2F *MASS_YY_SAME_CENTBIN[4][_narm][_ncharge];

	TH2F *MASS_Y_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Y_AND_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Y_OR_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Z_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Z_AND_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Z_OR_W_Z_SAME[_narm][_ncharge];

	TH2F *MASS_Y_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_Y_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_YY_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_YY_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_YY_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_Y_TRIG_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Y_TRIG_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_Y_TRIG_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_YY_TRIG_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_YY_TRIG_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_YY_TRIG_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_Y_W_PT_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_Y_W_PT_pp_SAME[_narm][_ncharge];

	TH2F *MASS_Y_W_Y_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_Y_W_Y_pp_SAME[_narm][_ncharge];


	//Run15pAl200 half MuTr dead
	TH2F *MASS_Y_HALF_SAME[_narm][_ncharge];
	TH2F *MASS_YY_HALF_SAME[_narm][_ncharge];

	TH2F *MASS_Y_HALF_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Y_HALF_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_Y_HALF_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_YY_HALF_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_YY_HALF_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_YY_HALF_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_Y_HALF_TRIG_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_Y_HALF_TRIG_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_Y_HALF_TRIG_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_YY_HALF_TRIG_W_Z_SAME[_narm][_ncharge];
	TH2F *MASS_YY_HALF_TRIG_W_dAu_SAME[_narm][_ncharge];
	TH2F *MASS_YY_HALF_TRIG_W_pp_SAME[_narm][_ncharge];

	TH2F *MASS_Y_A_SAME[_narm][_ncharge];
	TH2F *MASS_Y_B_SAME[_narm][_ncharge];
	TH2F *MASS_Y_C_SAME[_narm][_ncharge];

	TH2F *MASS_Y_A_MIXED[_narm][_ncharge];
	TH2F *MASS_Y_B_MIXED[_narm][_ncharge];
	TH2F *MASS_Y_C_MIXED[_narm][_ncharge];

	TH2F *MASS_YY_A_SAME[_narm][_ncharge];
	TH2F *MASS_YY_B_SAME[_narm][_ncharge];
	TH2F *MASS_YY_C_SAME[_narm][_ncharge];

	TH2F *MASS_YY_A_MIXED[_narm][_ncharge];
	TH2F *MASS_YY_B_MIXED[_narm][_ncharge];
	TH2F *MASS_YY_C_MIXED[_narm][_ncharge];


	//mass & pT distribution for acc x eff
	TH2F *N_PT_ETA_SAME_CUT0[_narm][_ncharge];
	TH2F *N_PT_ETA_SAME_CUT1[_narm][_ncharge];
	TH2F *N_PT_ETA_SAME_CUT2[_narm][_ncharge];
	TH2F *N_PT_ETA_SAME_CUT3[_narm][_ncharge];
	TH2F *N_PT_ETA_SAME_CUT4[_narm][_ncharge];
	TH2F *N_PT_ETA_SAME_CUT5[_narm][_ncharge];

	//mass & pT distribution
	TH2F *N_PT_ETA_LOOSE_SAME[_narm][_ncharge];
	TH2F *N_PT_ETA_SAME_AND[_narm][_ncharge];
	TH2F *N_PT_ETA_SAME_OR[_narm][_ncharge];
	TH2F *N_PT_ETA_MIXED_AND[_narm][_ncharge];
	TH2F *N_PT_ETA_MIXED_OR[_narm][_ncharge];

	TH2F *N_PT_Y_MU_SAME[_narm][_ncharge];
	TH2F *N_PT_Y_MU_W_dAu_SAME[_narm][_ncharge];
	TH2F *N_PT_Y_MU_W_pp_SAME[_narm][_ncharge];

	//TH1F *N_RAP_SAME_OR[_narm][_ncharge];

	TH2F *MASS_PT_CUT1_SAME[_narm][_ncharge];
	TH2F *MASS_PT_CUT1_MIXED[_narm][_ncharge];
	TH2F *MASS_PT_CUT2_SAME[_narm][_ncharge];
	TH2F *MASS_PT_CUT2_MIXED[_narm][_ncharge];
	TH2F *MASS_PT_CUT3_SAME[_narm][_ncharge];
	TH2F *MASS_PT_CUT3_MIXED[_narm][_ncharge];
	TH2F *MASS_PT_CUT4_SAME[_narm][_ncharge];
	TH2F *MASS_PT_CUT4_MIXED[_narm][_ncharge];

	TH2F *FVTX_PROB_PT[_narm];
	TH2F *FVTX_PROB_PT_MISMATCH[_narm];

	TH2F *FVTX_Q_PROB_PT[_narm];
	TH2F *FVTX_Q_PROB_PT_MISMATCH[_narm];

	//DCA QA histograms
	//picodst variable
	TH1F *VTX_CHI2[_narm];
	TH1F *FVTX_CHI2[_narm];
	TH2F *MUTR_RPHI[_narm];

	//smearing function
	TFile *file_smear;
	TRandom3 *rand3;
	TF1 *fsmear_fvtx[_narm];
	TF1 *fsmear_vtx[_narm];
	TF1 *fsmear_gaus;
	TF1 *fsmear_f;
	TF1 *fsmear_v;
	TH1F *hzvtx_smear;

	//weight function
	TFile *file_weight0;
	TFile *file_weight1;
	TF1 *fwt_pT_dAu_fwd;
	TF1 *fwt_pT_dAu_bwd;
	TF1 *fwt_pT_pp;
	TF1 *fwt_y_dAu;
	TF1 *fwt_y_pp;
	TF1 *fwt_z;
	TF1 *fwt_z_S;
	TF1 *fwt_z_N;

	TFile *file_weight2;
	TF1 *ftrig_wt[2];

	TFile *file_trigeff;
	TF1 *ftrigeff[2][6];

	std::set<float> set_dca_swap;
	std::pair<std::set<float>::iterator,bool> pr_dca_swap;


};
