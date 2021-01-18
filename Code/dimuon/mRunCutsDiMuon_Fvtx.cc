#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <set>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TVector3.h>

#include "mRunCutsDiMuon_Fvtx.h"

using namespace std;

//___________________________________________________________________
mRunCutsDiMuon_Fvtx::mRunCutsDiMuon_Fvtx()
{

	_dataset = "Run12pp200";
	_nptbin = 16;
	_ncentbin = 1;
	_filetype = "DATA";
	_triggertype = "NONE";
	_use_fvtx = true;

	_zst1[0] = -188.9619;
	_zst1[1] = +188.84249;

	cut_mass[0] = 2.0;
	cut_mass[1] = 14.0;

	_prev_runnumber = 0;
	_prev_bbcZ = 0.0;

	_prev_run = 0;
	_is_goodrun[0] = _is_goodrun[1] = true;

	return;
}

//___________________________________________________________________
mRunCutsDiMuon_Fvtx::~mRunCutsDiMuon_Fvtx()
{

	return;
}

//___________________________________________________________________
int mRunCutsDiMuon_Fvtx::Init(const char *listname = "file.lst")
{  

	cout << "mRunCutsDiMuon_Fvtx::Init" << endl;

	//first = true;
	
	set_pT_array();
	set_basic_cuts();

	load_files(listname);

	char tmp_filename[100];
	sprintf(tmp_filename,"Runcuts_DIMU_FVTXMUTR_%s_%s_%s.root",_dataset.data(),_filetype.data(),_triggertype.data());
	file_out = new TFile(tmp_filename,"RECREATE");
	std::cout << "opening TFile " << tmp_filename << std::endl;

	book_histos();

	process_event();

	End();

  return 0;
}

//___________________________________________________________________
void mRunCutsDiMuon_Fvtx::set_dataset(const char *dataset)
{

	cout << "mRunCutsDiMuon_Fvtx::set_dataset : " << dataset << endl; 
	_dataset = dataset;

	return;
}

//___________________________________________________________________
void mRunCutsDiMuon_Fvtx::set_mass_cut(float min, float max)
{

	cout << "mRunCuts_femto::set_mass_cut : " << min << " - " << max << endl; 
	cut_mass[0] = min;
	cut_mass[1] = max;

	return;
}

//___________________________________________________________________
void mRunCutsDiMuon_Fvtx::use_fvtx(bool use)
{

	cout << "mRunCuts_femto::use_fvtx : " << use << endl; 
	_use_fvtx = use;

	return;
}

//___________________________________________________________________
void mRunCutsDiMuon_Fvtx::set_filetype(const char *filetype)
{

	_filetype = filetype;

	/*
	if ( _filetype=="DATA" ) _is_sim = false;
	else _is_sim = true;
	*/
	if ( _filetype.find("SIM") != std::string::npos ){
		_is_sim = true;
	}else{
		_is_sim = false;
	}

	cout << "mRunCutsDiMuon_Fvtx::set_filetype : " << filetype << endl; 
	cout << "mRunCutsDiMuon_Fvtx::is_sim : " << _is_sim << endl; 

	return;
}

//___________________________________________________________________
void mRunCutsDiMuon_Fvtx::set_triggertype(const char *triggertype)
{

	cout << "mRunCutsDiMuon_Fvtx::set_triggertype : " << triggertype << endl; 
	_triggertype = triggertype;

	return;
}

//______________________________________________________
void mRunCutsDiMuon_Fvtx::load_files(const char *listname)
{


	string tmp_name;
	fChain = new TChain("ana_tree");

	//char tmp_name[300];
	int count = 0;
	//int tmp_runnum;

	if ( strstr(listname,"root") ){
		cout << "OPEN  : " << tmp_name << " #" << count << endl;
		fChain->Add(listname);
	}else{
		ifstream flist;
		flist.open(listname);
		while ( getline(flist,tmp_name) ){
			//while ( flist >> tmp_runnum ){
			//getline(flist, tmp_name);
			//sprintf(tmp_name,"/direct/phenix+hhj/shlim/work/09.run12_pdst/14.pp510/Jpsi/Run12pp510MU_dimuon_femtoDST/Run12pp510MU_femtoDST_dimu_%d.root",tmp_runnum);
			if ( count<10 ){
				cout << "OPEN  : " << tmp_name << " #" << count << endl;
			}
			fChain->Add(tmp_name.data());
			count++;
		}

		flist.close();
	}

	fChain->SetBranchAddress("_runnumber",&_runnumber);
	fChain->SetBranchAddress("_bbcZ",&_bbcZ);
	fChain->SetBranchAddress("_bbcqn",&_bbcqn);
	fChain->SetBranchAddress("_bbcqs",&_bbcqs);
	if ( _dataset!="Run15pp200" && _dataset!="Run9pp500" && _dataset!="Run9pp200" ){
		fChain->SetBranchAddress("_centrality",&_centrality);
	}
	fChain->SetBranchAddress("_scaled_trigbit",&_scaled_trigbit);

	fChain->SetBranchAddress("_ndimuon",&_ndimuon);
	fChain->SetBranchAddress("_same_event",&_same_event);
	fChain->SetBranchAddress("_charge",&_charge);
	fChain->SetBranchAddress("_mass",&_mass);
	fChain->SetBranchAddress("_pT",&_pT);
	fChain->SetBranchAddress("_pz",&_pz);
	fChain->SetBranchAddress("_rapidity",&_rapidity);
	fChain->SetBranchAddress("_evt_vtxchi2",&_evt_vtxchi2);

	fChain->SetBranchAddress("_tr_pT",_tr_pT);
	fChain->SetBranchAddress("_tr_pz",_tr_pz);
	fChain->SetBranchAddress("_tr_ntrhits",_tr_ntrhits);
	fChain->SetBranchAddress("_tr_nidhits",_tr_nidhits);
	fChain->SetBranchAddress("_tr_trhits",_tr_trhits);
	fChain->SetBranchAddress("_tr_idhits",_tr_idhits);
	fChain->SetBranchAddress("_tr_lastgap",_tr_lastgap);
	fChain->SetBranchAddress("_tr_trchi2",_tr_trchi2);
	fChain->SetBranchAddress("_tr_idchi2",_tr_idchi2);
	fChain->SetBranchAddress("_tr_DG0",_tr_DG0);
	fChain->SetBranchAddress("_tr_DDG0",_tr_DDG0);
	fChain->SetBranchAddress("_tr_xst1",_tr_xst1);
	fChain->SetBranchAddress("_tr_xst2",_tr_xst2);
	fChain->SetBranchAddress("_tr_xst3",_tr_xst3);
	fChain->SetBranchAddress("_tr_yst1",_tr_yst1);
	fChain->SetBranchAddress("_tr_yst2",_tr_yst2);
	fChain->SetBranchAddress("_tr_yst3",_tr_yst3);

	fChain->SetBranchAddress("_tr_pdtheta",_tr_pdtheta);
	fChain->SetBranchAddress("_tr_rapidity",_tr_rapidity);

	if ( _use_fvtx ){
		fChain->SetBranchAddress("_fvtxX",&_fvtxX);
		fChain->SetBranchAddress("_fvtxY",&_fvtxY);
		fChain->SetBranchAddress("_fvtxZ",&_fvtxZ);

		fChain->SetBranchAddress("_fvtxX_err",&_fvtxX_err);
		fChain->SetBranchAddress("_fvtxY_err",&_fvtxY_err);
		fChain->SetBranchAddress("_fvtxZ_err",&_fvtxZ_err);

		fChain->SetBranchAddress("_mult_fvtxS",&_mult_fvtxS);
		fChain->SetBranchAddress("_mult_fvtxN",&_mult_fvtxN);
		fChain->SetBranchAddress("_mult_fvtx_dimu",&_mult_fvtx_dimu);
		fChain->SetBranchAddress("_mult_svx",&_mult_svx);
		fChain->SetBranchAddress("_mass_fvtx",&_mass_fvtx);
		fChain->SetBranchAddress("_mass_fvtxmutr",&_mass_fvtxmutr);

		fChain->SetBranchAddress("_pT_fvtxmutr",&_pT_fvtxmutr);
		fChain->SetBranchAddress("_pz_fvtxmutr",&_pz_fvtxmutr);
		fChain->SetBranchAddress("_rapidity_fvtxmutr",&_rapidity_fvtxmutr);
		fChain->SetBranchAddress("_evt_vtxchi2_fvtxmutr",&_evt_vtxchi2_fvtxmutr);

		fChain->SetBranchAddress("_tr_dca_z",_tr_dca_z);
		fChain->SetBranchAddress("_tr_dca_r",_tr_dca_r);
		fChain->SetBranchAddress("_tr_dca_phi",_tr_dca_phi);
		//fChain->SetBranchAddress("_tr_dca_r_fvtx",_tr_dca_r_fvtx);
		//
		fChain->SetBranchAddress("_tr_nhits_fvtx",_tr_nhits_fvtx);
		fChain->SetBranchAddress("_tr_nhits_fvtx_charge",_tr_nhits_fvtx_charge);
		fChain->SetBranchAddress("_tr_dphi_fvtx",_tr_dphi_fvtx);
		fChain->SetBranchAddress("_tr_dtheta_fvtx",_tr_dtheta_fvtx);
		fChain->SetBranchAddress("_tr_dr_fvtx",_tr_dr_fvtx);
		fChain->SetBranchAddress("_tr_chi2_fvtx",_tr_chi2_fvtx);

		fChain->SetBranchAddress("_tr_chi2_fvtxmutr",_tr_chi2_fvtxmutr);
		fChain->SetBranchAddress("_tr_vtx_index",_tr_vtx_index);
		fChain->SetBranchAddress("_tr_charge_fvtx",_tr_charge_fvtx);

		fChain->SetBranchAddress("_tr_px_fvtxmutr",_tr_px_fvtxmutr);
		fChain->SetBranchAddress("_tr_py_fvtxmutr",_tr_py_fvtxmutr);
		fChain->SetBranchAddress("_tr_pz_fvtxmutr",_tr_pz_fvtxmutr);

		fChain->SetBranchAddress("_tr_px_fvtx",_tr_px_fvtx);
		fChain->SetBranchAddress("_tr_py_fvtx",_tr_py_fvtx);
		fChain->SetBranchAddress("_tr_pz_fvtx",_tr_pz_fvtx);

		fChain->SetBranchAddress("_tr_x0_fvtxmutr",_tr_x0_fvtxmutr);
		fChain->SetBranchAddress("_tr_y0_fvtxmutr",_tr_y0_fvtxmutr);
		fChain->SetBranchAddress("_tr_z0_fvtxmutr",_tr_z0_fvtxmutr);

		fChain->SetBranchAddress("_tr_x0_fvtx",_tr_x0_fvtx);
		fChain->SetBranchAddress("_tr_y0_fvtx",_tr_y0_fvtx);
		fChain->SetBranchAddress("_tr_z0_fvtx",_tr_z0_fvtx);

		fChain->SetBranchAddress("_tr_hit_pattern",_tr_hit_pattern);
	}

	if ( _is_sim ){
		fChain->SetBranchAddress("_simZ",&_simZ);
		fChain->SetBranchAddress("_nhepmc",&_nhepmc);
		fChain->SetBranchAddress("_hepmc_g_pid",_hepmc_g_pid);
		fChain->SetBranchAddress("_hepmc_g_pT",_hepmc_g_pT);
		fChain->SetBranchAddress("_hepmc_g_rapidity",_hepmc_g_rapidity);
		//fChain->SetBranchAddress("_tr_mis_match",_tr_mis_match);
		fChain->SetBranchAddress("_tr_dca_r_sim",_tr_dca_r_sim);
		fChain->SetBranchAddress("_tr_dca_r_smear",_tr_dca_r_smear);
		fChain->SetBranchAddress("_tr_mc_g_pid",_tr_mc_g_pid);

		fChain->SetBranchAddress("_tr_mc_hits_fvtx",_tr_mc_hits_fvtx);
		fChain->SetBranchAddress("_tr_mc_hits_svx",_tr_mc_hits_svx);

		fChain->SetBranchAddress("_tr_mc_hits_fvtx_true",_tr_mc_hits_fvtx_true);
		fChain->SetBranchAddress("_tr_mc_hits_svx_true",_tr_mc_hits_svx_true);

		fChain->SetBranchAddress("_tr_mc_hits_mutr_true",_tr_mc_hits_mutr_true);
		fChain->SetBranchAddress("_tr_mc_hits_muid_true",_tr_mc_hits_muid_true);

		fChain->SetBranchAddress("_trig_emul_2D_S",&_trig_emul_2D_S);
		fChain->SetBranchAddress("_trig_emul_2D_N",&_trig_emul_2D_N);
	}else{
		fChain->SetBranchAddress("_doubleint_frac",&_doubleint_frac);
	}


	//fChain->SetBranchAddress("_mult_fvtx_prim_cut",&_mult_fvtx_prim_cut);


	return;

}

//______________________________________________________
int mRunCutsDiMuon_Fvtx::process_event()
{

	cout << "mRunCutsDiMuon_Fvtx::process_event" << endl; 

	Long64_t nentries = fChain->GetEntries();

	cout << nentries << " tracks to be analyzed!!" << endl;

	for (Long64_t ien=0; ien<nentries; ien++){

		fChain->GetEntry(ien);

		if ( (ien%(nentries/10))==0 ){
			cout << setw(10) << int(ien*100.0/nentries + 0.1) << " percent completed..." << endl;
		}

		if ( _rapidity>0 ) _arm = 1;
		else _arm = 0;

		//CHECK GOOD RUN
		//if ( !check_run(_arm) && _filetype=="DATA" ) continue;

		int charge = 2;
		if ( _charge==0 ) charge = 0;
		else if ( _charge==-2 ) charge = 1;
		else if ( _charge==2 ) charge = 2;

		/*
		float good_mass = _mass;

		if ( !isnan(_mass_fvtxmutr) && _mass_fvtxmutr>cut_mass[0] && _mass_fvtxmutr<cut_mass[1] ){
			good_mass = _mass_fvtxmutr;
		}
		*/

		//QA histogram
		if ( _mass>2.7 && _mass<3.5 ){
			if ( _tr_ntrhits[0]>=11 && _tr_nidhits[0]>=6 ){
				QA_MUTR_NHIT[_arm][charge]->Fill(_tr_ntrhits[0]);
				QA_MUTR_CHI2[_arm][charge]->Fill(_tr_trchi2[0]);
				QA_MUID_CHI2[_arm][charge]->Fill(_tr_idchi2[0]);
				QA_DG0[_arm][charge]->Fill(_tr_DG0[0]);
				QA_DDG0[_arm][charge]->Fill(_tr_DDG0[0]);

				if ( fabs(_bbcZ)<10.0 && _tr_nhits_fvtx[0]>3.5 ){
					QA_FVTXMUTR_CHI2[_arm][charge]->Fill(_tr_chi2_fvtxmutr[0]);
				}
			}
			if ( _tr_ntrhits[1]>=11 && _tr_nidhits[1]>=6 ){
				QA_MUTR_NHIT[_arm][charge]->Fill(_tr_ntrhits[1]);
				QA_MUTR_CHI2[_arm][charge]->Fill(_tr_trchi2[1]);
				QA_MUID_CHI2[_arm][charge]->Fill(_tr_idchi2[1]);
				QA_DG0[_arm][charge]->Fill(_tr_DG0[1]);
				QA_DDG0[_arm][charge]->Fill(_tr_DDG0[1]);

				if ( fabs(_bbcZ)<10.0 && _tr_nhits_fvtx[1]>3.5 ){
					QA_FVTXMUTR_CHI2[_arm][charge]->Fill(_tr_chi2_fvtxmutr[1]);
				}
			}

			QA_EVTVTX_CHI2[_arm][charge]->Fill(_evt_vtxchi2);

			int muid_hit_pattern[10] = {0};
			for (int itrk=0; itrk<2; itrk++){
				if ( _tr_lastgap[itrk]<4 ) continue;

				for (int ip=0; ip<10; ip++)
				{    
					muid_hit_pattern[ip] = 0.;
					int iBit = ip/2 + 5*(ip%2);
					if( (_tr_idhits[itrk]&(1<<iBit)) ) muid_hit_pattern[ip] = 1; 
				}//ip 

				int depthH = 0, nhitsH = 0;
				int depthV = 0, nhitsV = 0;
				for (int igap=0; igap<5; igap++)
				{
					if ( muid_hit_pattern[igap*2]>0 ){
						depthH = igap;
						nhitsH++;
					}
					if ( muid_hit_pattern[igap*2+1]>0 ){
						depthV = igap;
						nhitsV++;
					}
				}

				QA_MUID_DIFF_GAP[_arm]->Fill(abs(depthH-depthV));
				QA_MUID_DIFF_NHIT[_arm]->Fill(abs(nhitsH-nhitsV));
				QA_MUID_DEPTH[_arm]->Fill(depthH, depthV);
				QA_MUID_NHIT[_arm]->Fill(nhitsH, nhitsV);
				QA_MUID_NHIT_1D[_arm]->Fill(muid_hit_pattern[0]+muid_hit_pattern[5]);
			}//itrk
		}//QA


		/*
		float weight_y_dAu = 1, weight_pT_dAu = 1;
		float weight_y_pp = 1, weight_pT_pp = 1;

		if ( _is_sim ){
			weight_z = fwt_z->Eval(_simZ);

			if ( fabs(_hepmc_g_rapidity[0])>0.6 && fabs(_hepmc_g_rapidity[0])<3.0 ){

				if ( _hepmc_g_rapidity[0]<0 ){
					weight_y_dAu = fwt_y_dAu->Eval(_hepmc_g_rapidity[0]);
					weight_pT_dAu = fwt_pT_dAu_bwd->Eval(_hepmc_g_pT[0]);
				}else{
					weight_y_dAu = fwt_y_dAu->Eval(_hepmc_g_rapidity[0]);
					weight_pT_dAu = fwt_pT_dAu_fwd->Eval(_hepmc_g_pT[0]);
				}

				weight_y_pp = fwt_y_pp->Eval(fabs(_hepmc_g_rapidity[0]));
				weight_pT_pp = fwt_pT_pp->Eval(_hepmc_g_pT[0]);

			}//
		}//_is_sim
		*/

		/*
		int centbin = 0; 
		int centbin_s = 0, centbin_n = 0;
		if ( _dataset=="Run15pAl200" ){
			if ( _centrality<=20 ) centbin = 0;
			else if ( _centrality<=40 ) centbin = 1;
			else if ( _centrality<=72 ) centbin = 2;
			else centbin = -1;
		}else if ( _dataset=="Run15pAu200" ){
			if ( _centrality<=20 ) centbin = 0;
			else if ( _centrality<=40 ) centbin = 1;
			else if ( _centrality<=60 ) centbin = 2;
			else if ( _centrality<=84 ) centbin = 3;
			else centbin = -1;
		}else if ( _dataset=="Run14HeAu200" ){
			if ( _centrality<=20 ) centbin = 0;
			else if ( _centrality<=40 ) centbin = 1;
			else if ( _centrality<=60 ) centbin = 2;
			else if ( _centrality<=88 ) centbin = 3;
			else centbin = -1;
		}else if ( _dataset=="Run14AuAu200" ){
			if ( _centrality<=20 ) { centbin = centbin_s = centbin_n = 0; }
			else if ( _centrality<=40 ) { centbin = centbin_s = centbin_n = 1; }
			else if ( _centrality<=60 ) { centbin = centbin_s = centbin_n = 2; }
			else if ( _centrality<=94 ) { centbin = centbin_s = centbin_n = 3; }
			else { centbin = centbin_s = centbin_n = -1; }
		}else if ( _dataset=="Run15pp200" ){
			if ( _bbcqs>6.5 ) centbin_s = 0;
			else if ( _bbcqs>3.7 ) centbin_s = 1;
			else if ( _bbcqs>1.85 ) centbin_s = 2;
			else centbin_s = 3;

			if ( _bbcqn>6.5 ) centbin_n = 0;
			else if ( _bbcqn>3.7 ) centbin_n = 1;
			else if ( _bbcqn>1.85 ) centbin_n = 2;
			else centbin_n = 3;
		}
		*/

		float weight_z = 1;

		if ( _is_sim ){
			weight_z = fwt_z->Eval(_simZ);
		}

		if ( check_event() && check_dimuon() ){

			float pvalue_tr0 = 0.0, pvalue_tr1 = 0.0;

			if ( _tr_nhits_fvtx[0]>0 ){
				pvalue_tr0 = TMath::Prob(_tr_chi2_fvtx[0]*(2*_tr_nhits_fvtx[0]-5), 2*_tr_nhits_fvtx[0]-5);
			}
			if ( _tr_nhits_fvtx[1]>0 ){
				pvalue_tr1 = TMath::Prob(_tr_chi2_fvtx[1]*(2*_tr_nhits_fvtx[1]-5), 2*_tr_nhits_fvtx[1]-5);
			}

			float purity_tr0 = 0.0, purity_tr1 = 0.0;
			if ( _tr_nhits_fvtx[0]>0 ){
				purity_tr0 = float(_tr_mc_hits_fvtx_true[0]+_tr_mc_hits_svx_true[0])/float(_tr_mc_hits_fvtx[0]+_tr_mc_hits_svx[0]);
			}
			if ( _tr_nhits_fvtx[1]>0 ){
				purity_tr1 = float(_tr_mc_hits_fvtx_true[1]+_tr_mc_hits_svx_true[1])/float(_tr_mc_hits_fvtx[1]+_tr_mc_hits_svx[1]);
			}

			if ( _tr_nhits_fvtx[0]>0 ){
				FVTX_PROB_PT[_arm]->Fill(pvalue_tr0, _tr_pT[0]);
				if ( purity_tr0<0.5 ){
					FVTX_PROB_PT_MISMATCH[_arm]->Fill(pvalue_tr0, _tr_pT[0]);
				}

				if ( check_fvtx(0) ){
					FVTX_Q_PROB_PT[_arm]->Fill(pvalue_tr0, _tr_pT[0]);
					if ( purity_tr0<0.5 ){
						FVTX_Q_PROB_PT_MISMATCH[_arm]->Fill(pvalue_tr0, _tr_pT[0]);
					}
				}
			}//

			if ( _tr_nhits_fvtx[1]>0 ){
				FVTX_PROB_PT[_arm]->Fill(pvalue_tr1, _tr_pT[1]);
				if ( purity_tr1<0.5 ){
					FVTX_PROB_PT_MISMATCH[_arm]->Fill(pvalue_tr1, _tr_pT[1]);
				}

				if ( check_fvtx(1) ){
					FVTX_Q_PROB_PT[_arm]->Fill(pvalue_tr1, _tr_pT[1]);
					if ( purity_tr1<0.5 ){
						FVTX_Q_PROB_PT_MISMATCH[_arm]->Fill(pvalue_tr1, _tr_pT[1]);
					}
				}
			}

			if ( _same_event ){
				MASS_Y_SAME[_arm][charge]->Fill(_mass, fabs(_rapidity));
				MASS_Z_SAME[_arm][charge]->Fill(_mass, _bbcZ);

				MASS_Y_W_Z_SAME[_arm][charge]->Fill(_mass, fabs(_rapidity), weight_z);
				MASS_Z_W_Z_SAME[_arm][charge]->Fill(_mass, _bbcZ, weight_z);
			}else{
				MASS_Y_MIXED[_arm][charge]->Fill(_mass, fabs(_rapidity));
				MASS_Z_MIXED[_arm][charge]->Fill(_mass, _bbcZ);
			}

			if ( check_fvtx(0) && check_fvtx(1) ){
				if ( _same_event ){
					MASS_Y_AND_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity));
					MASS_Z_AND_SAME[_arm][charge]->Fill(_mass_fvtxmutr, _bbcZ);

					MASS_Y_AND_W_Z_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), weight_z);
					MASS_Z_AND_W_Z_SAME[_arm][charge]->Fill(_mass_fvtxmutr, _bbcZ, weight_z);
				}else{
					MASS_Y_AND_MIXED[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity));
					MASS_Z_AND_MIXED[_arm][charge]->Fill(_mass_fvtxmutr, _bbcZ);
				}
			}//

			if ( check_fvtx(0) || check_fvtx(1) ){
				if ( _same_event ){
					MASS_Y_OR_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity));
					MASS_Z_OR_SAME[_arm][charge]->Fill(_mass_fvtxmutr, _bbcZ);

					MASS_Y_OR_W_Z_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), weight_z);
					MASS_Z_OR_W_Z_SAME[_arm][charge]->Fill(_mass_fvtxmutr, _bbcZ, weight_z);
				}else{
					MASS_Y_OR_MIXED[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity));
					MASS_Z_OR_MIXED[_arm][charge]->Fill(_mass_fvtxmutr, _bbcZ);
				}
			}//

			/*
			if ( _same_event ){
				MASS_Y_PT_SAME[_arm][charge]->Fill(_mass, fabs(_rapidity), _pT);

				if ( pvalue_tr0>0.05 && pvalue_tr1>0.05 ){
				//if ( _tr_nhits_fvtx[0]>0 && _tr_nhits_fvtx[1]>0 ){
					MASS_2MATCH_Y_PT_SAME[_arm][charge]->Fill(_mass, fabs(_rapidity), _pT);
					MASS_FVTX_2MATCH_Y_PT_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);
					MASS_FVTX_2MATCH_Y_PT_SAME_CENT[centbin][_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);

					if ( purity_tr0<0.6 || purity_tr1<0.6 ){
						MASS_2MATCH_MISMATCH_Y_PT_SAME[_arm][charge]->Fill(_mass, fabs(_rapidity), _pT);
						MASS_FVTX_2MATCH_MISMATCH_Y_PT_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);
						MASS_FVTX_2MATCH_MISMATCH_Y_PT_SAME_CENT[centbin][_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);
					}

				}else if ( pvalue_tr0<0.05 && pvalue_tr1<0.05 ){
				//}else if ( _tr_nhits_fvtx[0]==0 && _tr_nhits_fvtx[1]==0 ){
					//do nothing
				}else{
					MASS_1MATCH_Y_PT_SAME[_arm][charge]->Fill(_mass, fabs(_rapidity), _pT);
					MASS_FVTX_1MATCH_Y_PT_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);
					MASS_FVTX_1MATCH_Y_PT_SAME_CENT[centbin][_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);

					if ( (pvalue_tr0>0.05 && purity_tr0<0.6) || (pvalue_tr1>0.05 && purity_tr1<0.6) ){
						MASS_1MATCH_MISMATCH_Y_PT_SAME[_arm][charge]->Fill(_mass, fabs(_rapidity), _pT);
						MASS_FVTX_1MATCH_MISMATCH_Y_PT_SAME[_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);
						MASS_FVTX_1MATCH_MISMATCH_Y_PT_SAME_CENT[centbin][_arm][charge]->Fill(_mass_fvtxmutr, fabs(_rapidity), _pT);
					}

				}


			}//_same_event
			*/

		}//check_event && check_dimuon

/*
		if ( 0 ){
		//if ( check_event() && check_dimuon() && check_fvtx() ){

			 float phi0 = atan2(_tr_yst2[0],_tr_xst2[0]);
			 float phi1 = atan2(_tr_yst2[1],_tr_xst2[1]);

			if ( _same_event ){
				MASS_PT_SAME[_arm][charge]->Fill(_mass,_pT);
				MASS_Y_SAME[_arm][charge]->Fill(good_mass,fabs(_rapidity));

				MASS_CENTS_SAME[_arm][charge]->Fill(good_mass,centbin_s);
				MASS_CENTN_SAME[_arm][charge]->Fill(good_mass,centbin_n);

				if ( phi0>-1.2 && phi0<2.0 && phi1>-1.2 && phi1<2.0 ){
					MASS_Y_HALF_SAME[_arm][charge]->Fill(_mass,fabs(_rapidity));
				}

				//if ( (fabs(_tr_rapidity[0])<2.2+0.2*_arm) && (fabs(_tr_rapidity[1])<2.2+0.2*_arm) && fabs(_tr_rapidity[0])>1.2 && fabs(_tr_rapidity[1])>1.2 ){
				if ( 1 ){
					MASS_YY_SAME[_arm][charge]->Fill(_mass,fabs(_rapidity));
					MASS_YY_SAME_CENTBIN[centbin][_arm][charge]->Fill(_mass,fabs(_rapidity));
					if ( phi0>-1.2 && phi0<2.0 && phi1>-1.2 && phi1<2.0 ){
						MASS_YY_HALF_SAME[_arm][charge]->Fill(_mass,fabs(_rapidity));
					}
				}
				MASS_HIGH_PT_SAME[_arm][charge]->Fill(_mass,_pT);

				N_PT_Y_MU_SAME[_arm][charge]->Fill(_tr_pT[0],fabs(_tr_rapidity[0]));
				N_PT_Y_MU_SAME[_arm][charge]->Fill(_tr_pT[1],fabs(_tr_rapidity[1]));

				if ( charge==0 && _mass>2.7 && _mass<3.5  ){
					if ( _is_sim ){
						if ( (_arm==0 && _trig_emul_2D_S>0) || (_arm==1 && _trig_emul_2D_N>0) ){
							MUTR_RPHI[_arm]->Fill(atan2(_tr_yst2[0],_tr_xst2[0]),sqrt(_tr_xst2[0]*_tr_xst2[0]+_tr_yst2[0]*_tr_yst2[0]));
							MUTR_RPHI[_arm]->Fill(atan2(_tr_yst2[1],_tr_xst2[1]),sqrt(_tr_xst2[1]*_tr_xst2[1]+_tr_yst2[1]*_tr_yst2[1]));
						}
					}else{
						MUTR_RPHI[_arm]->Fill(atan2(_tr_yst2[0],_tr_xst2[0]),sqrt(_tr_xst2[0]*_tr_xst2[0]+_tr_yst2[0]*_tr_yst2[0]));
						MUTR_RPHI[_arm]->Fill(atan2(_tr_yst2[1],_tr_xst2[1]),sqrt(_tr_xst2[1]*_tr_xst2[1]+_tr_yst2[1]*_tr_yst2[1]));
					}
				}//J/psi candidate

			}else{
				MASS_PT_MIXED[_arm][charge]->Fill(_mass,_pT);
				MASS_Y_MIXED[_arm][charge]->Fill(good_mass,fabs(_rapidity));

				MASS_CENTS_MIXED[_arm][charge]->Fill(good_mass,centbin_s);
				MASS_CENTN_MIXED[_arm][charge]->Fill(good_mass,centbin_n);
				if ( (fabs(_tr_rapidity[0])<2.2+0.2*_arm) && (fabs(_tr_rapidity[1])<2.2+0.2*_arm) && fabs(_tr_rapidity[0])>1.2 && fabs(_tr_rapidity[1])>1.2 ){
					MASS_YY_MIXED[_arm][charge]->Fill(_mass,fabs(_rapidity));
				}

			}//same_event
		}//check_event && check_dimuon
			*/

	}//nentries

	return 0;
}

//______________________________________________________
int mRunCutsDiMuon_Fvtx::End()
{
  cout << "mRunCutsDiMuon_Fvtx::End" << endl;

	file_out->cd();
	file_out->Write();
	//file_out->Close();

	delete file_out;

  return 0;
}

//______________________________________________________
void mRunCutsDiMuon_Fvtx::book_histos()
{

  cout << "mRunCutsDiMuon_Fvtx::book_histos" << endl;

	//char hname[300];

	int run_range[2] = {0};
	get_run_range(_dataset,run_range);

	for (int iarm=0; iarm<_narm; iarm++){
		QA_MUID_DIFF_GAP[iarm] = new TH1F(Form("QA_MUID_DIFF_GAP_arm%d",iarm),"",10,0,10);
		QA_MUID_DIFF_NHIT[iarm] = new TH1F(Form("QA_MUID_DIFF_NHIT_arm%d",iarm),"",10,0,10);
		QA_MUID_DEPTH[iarm] = new TH2F(Form("QA_MUID_DEPTH_arm%d",iarm),"",10,0,10,10,0,10);
		QA_MUID_NHIT[iarm] = new TH2F(Form("QA_MUID_NHIT_arm%d",iarm),"",10,0,10,10,0,10);
		QA_MUID_NHIT_1D[iarm] = new TH1F(Form("QA_MUID_NHIT_1D_arm%d",iarm),"",12,0,12);

		for (int ichg=0; ichg<_ncharge; ichg++){

			QA_MUTR_NHIT[iarm][ichg] = new TH1F(Form("QA_MUTR_NHIT_arm%d_chg%d",iarm,ichg),"",10,6.5,16.5);
			QA_MUTR_CHI2[iarm][ichg] = new TH1F(Form("QA_MUTR_CHI2_arm%d_chg%d",iarm,ichg),"",100,0,20);
			QA_MUID_CHI2[iarm][ichg] = new TH1F(Form("QA_MUID_CHI2_arm%d_chg%d",iarm,ichg),"",100,0,10);
			QA_EVTVTX_CHI2[iarm][ichg] = new TH1F(Form("QA_EVTVTX_CHI2_arm%d_chg%d",iarm,ichg),"",100,0,10);
			QA_FVTXMUTR_CHI2[iarm][ichg] = new TH1F(Form("QA_FVTXMUTR_CHI2_arm%d_chg%d",iarm,ichg),"",100,0,20);

			QA_DG0[iarm][ichg] = new TH1F(Form("QA_DG0_arm%d_chg%d",iarm,ichg),"",100,0,60);
			QA_DDG0[iarm][ichg] = new TH1F(Form("QA_DDG0_arm%d_chg%d",iarm,ichg),"",100,0,40);

			MASS_Y_SAME[iarm][ichg] = new TH2F(Form("mass_y_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
			MASS_Y_MIXED[iarm][ichg] = new TH2F(Form("mass_y_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

			MASS_Z_SAME[iarm][ichg] = new TH2F(Form("mass_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);
			MASS_Z_MIXED[iarm][ichg] = new TH2F(Form("mass_z_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);

			MASS_Y_AND_SAME[iarm][ichg] = new TH2F(Form("mass_y_and_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
			MASS_Y_AND_MIXED[iarm][ichg] = new TH2F(Form("mass_y_and_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
			MASS_Y_OR_SAME[iarm][ichg] = new TH2F(Form("mass_y_or_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
			MASS_Y_OR_MIXED[iarm][ichg] = new TH2F(Form("mass_y_or_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

			MASS_Z_AND_SAME[iarm][ichg] = new TH2F(Form("mass_z_and_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);
			MASS_Z_AND_MIXED[iarm][ichg] = new TH2F(Form("mass_z_and_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);
			MASS_Z_OR_SAME[iarm][ichg] = new TH2F(Form("mass_z_or_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);
			MASS_Z_OR_MIXED[iarm][ichg] = new TH2F(Form("mass_z_or_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);

			MASS_Y_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_y_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
			MASS_Y_AND_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_y_and_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
			MASS_Y_OR_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_y_or_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

			MASS_Z_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_z_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);
			MASS_Z_AND_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_z_and_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);
			MASS_Z_OR_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_z_or_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,60,-30,30);

			MASS_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_2MATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_2match_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_1MATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_1match_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_2MATCH_MISMATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_2match_mismatch_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_1MATCH_MISMATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_1match_mismatch_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_FVTX_2MATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_fvtx_2match_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_FVTX_1MATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_fvtx_1match_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_FVTX_2MATCH_MISMATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_fvtx_2match_mismatch_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			MASS_FVTX_1MATCH_MISMATCH_Y_PT_SAME[iarm][ichg] = new TH3F(Form("mass_fvtx_1match_mismatch_y_pT_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);

			for (int icent=0; icent<_ncentbin; icent++){
				MASS_FVTX_2MATCH_Y_PT_SAME_CENT[icent][iarm][ichg] = new TH3F(Form("mass_fvtx_2match_y_pT_same_cent%d_arm%d_chg%d",icent,iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
				MASS_FVTX_1MATCH_Y_PT_SAME_CENT[icent][iarm][ichg] = new TH3F(Form("mass_fvtx_1match_y_pT_same_cent%d_arm%d_chg%d",icent,iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
				MASS_FVTX_2MATCH_MISMATCH_Y_PT_SAME_CENT[icent][iarm][ichg] = new TH3F(Form("mass_fvtx_2match_mismatch_y_pT_same_cent%d_arm%d_chg%d",icent,iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
				MASS_FVTX_1MATCH_MISMATCH_Y_PT_SAME_CENT[icent][iarm][ichg] = new TH3F(Form("mass_fvtx_1match_mismatch_y_pT_same_cent%d_arm%d_chg%d",icent,iarm,ichg),"",90,1.5,6.0,4,1.2,2.2,10,0,10);
			}//icent


			if ( 0 ){

				MASS_PT_SAME[iarm][ichg] = new TH2F(Form("mass_pT_same_arm%d_chg%d",iarm,ichg),"",80,2.0,6.0,10,0,10);
				MASS_PT_MIXED[iarm][ichg] = new TH2F(Form("mass_pT_mixed_arm%d_chg%d",iarm,ichg),"",80,2.0,6.0,10,0,10);

				MASS_FVTXMUTR_PT_SAME[iarm][ichg] = new TH2F(Form("mass_fvtxmutr_pT_same_arm%d_chg%d",iarm,ichg),"",80,2.0,6.0,10,0,10);
				MASS_FVTXMUTR_PT_MIXED[iarm][ichg] = new TH2F(Form("mass_fvtxmutr_pT_mixed_arm%d_chg%d",iarm,ichg),"",80,2.0,6.0,10,0,10);

				MASS_HIGH_PT_SAME[iarm][ichg] = new TH2F(Form("mass_high_pT_same_arm%d_chg%d",iarm,ichg),"",16,6.0,14.0,10,0,10);
				MASS_HIGH_PT_MIXED[iarm][ichg] = new TH2F(Form("mass_high_pT_mixed_arm%d_chg%d",iarm,ichg),"",16,6.0,14.0,10,0,10);

				MASS_Y_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_y_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_y_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_YY_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_yy_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_yy_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_yy_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_Y_TRIG_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_y_trig_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_TRIG_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_y_trig_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_TRIG_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_y_trig_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_YY_TRIG_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_yy_trig_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_TRIG_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_yy_trig_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_TRIG_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_yy_trig_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_Y_W_PT_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_y_w_pT_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_W_PT_pp_SAME[iarm][ichg] = new TH2F(Form("mass_y_w_pT_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_Y_W_Y_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_y_w_y_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_W_Y_pp_SAME[iarm][ichg] = new TH2F(Form("mass_y_w_y_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_Y_HALF_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_y_HALF_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_HALF_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_y_HALF_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_HALF_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_y_HALF_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_YY_HALF_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_yy_HALF_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_HALF_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_yy_HALF_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_HALF_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_yy_HALF_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_Y_HALF_TRIG_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_y_HALF_trig_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_HALF_TRIG_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_y_HALF_trig_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_Y_HALF_TRIG_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_y_HALF_trig_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_YY_HALF_TRIG_W_Z_SAME[iarm][ichg] = new TH2F(Form("mass_yy_HALF_trig_w_z_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_HALF_TRIG_W_dAu_SAME[iarm][ichg] = new TH2F(Form("mass_yy_HALF_trig_w_dAu_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_HALF_TRIG_W_pp_SAME[iarm][ichg] = new TH2F(Form("mass_yy_HALF_trig_w_pp_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_YY_SAME[iarm][ichg] = new TH2F(Form("mass_yy_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_MIXED[iarm][ichg] = new TH2F(Form("mass_yy_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				MASS_CENTS_SAME[iarm][ichg] = new TH2F(Form("mass_cents_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,0,4);
				MASS_CENTS_MIXED[iarm][ichg] = new TH2F(Form("mass_cents_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,0,4);
				MASS_CENTN_SAME[iarm][ichg] = new TH2F(Form("mass_centn_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,0,4);
				MASS_CENTN_MIXED[iarm][ichg] = new TH2F(Form("mass_centn_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,0,4);

				MASS_Y_HALF_SAME[iarm][ichg] = new TH2F(Form("mass_y_HALF_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				MASS_YY_HALF_SAME[iarm][ichg] = new TH2F(Form("mass_yy_HALF_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

				for (int icent=0; icent<4; icent++){
					MASS_YY_SAME_CENTBIN[icent][iarm][ichg] = new TH2F(Form("mass_yy_same_centbin%d_arm%d_chg%d",icent,iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				}

				if ( !_is_sim ){
					MASS_Y_A_SAME[iarm][ichg] = new TH2F(Form("mass_y_A_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_Y_B_SAME[iarm][ichg] = new TH2F(Form("mass_y_B_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_Y_C_SAME[iarm][ichg] = new TH2F(Form("mass_y_C_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

					MASS_YY_A_SAME[iarm][ichg] = new TH2F(Form("mass_yy_A_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_YY_B_SAME[iarm][ichg] = new TH2F(Form("mass_yy_B_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_YY_C_SAME[iarm][ichg] = new TH2F(Form("mass_yy_C_same_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

					MASS_Y_A_MIXED[iarm][ichg] = new TH2F(Form("mass_y_A_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_Y_B_MIXED[iarm][ichg] = new TH2F(Form("mass_y_B_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_Y_C_MIXED[iarm][ichg] = new TH2F(Form("mass_y_C_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);

					MASS_YY_A_MIXED[iarm][ichg] = new TH2F(Form("mass_yy_A_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_YY_B_MIXED[iarm][ichg] = new TH2F(Form("mass_yy_B_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
					MASS_YY_C_MIXED[iarm][ichg] = new TH2F(Form("mass_yy_C_mixed_arm%d_chg%d",iarm,ichg),"",90,1.5,6.0,4,1.2,2.2);
				}

				N_PT_Y_MU_SAME[iarm][ichg] = new TH2F(Form("n_pT_y_mu_same_arm%d_chg%d",iarm,ichg),"",50,0,10,10,1.2,2.2);
				N_PT_Y_MU_W_dAu_SAME[iarm][ichg] = new TH2F(Form("n_pT_y_mu_w_dAu_same_arm%d_chg%d",iarm,ichg),"",50,0,10,10,1.2,2.2);
				N_PT_Y_MU_W_pp_SAME[iarm][ichg] = new TH2F(Form("n_pT_y_mu_w_pp_same_arm%d_chg%d",iarm,ichg),"",50,0,10,10,1.2,2.2);
			}//1


		}//ichg
	}//iarm

	for (int iarm=0; iarm<_narm; iarm++){
		MUTR_RPHI[iarm] = new TH2F(Form("MUTR_RPHI_%d",iarm),"",128,-TMath::Pi(),TMath::Pi(),250,0,250);

		FVTX_PROB_PT[iarm] = new TH2F(Form("FVTX_PROB_PT_arm%d",iarm),"",100,0,1,20,0,10);
		FVTX_PROB_PT_MISMATCH[iarm] = new TH2F(Form("FVTX_PROB_PT_MISMATCH_arm%d",iarm),"",100,0,1,20,0,10);

		FVTX_Q_PROB_PT[iarm] = new TH2F(Form("FVTX_Q_PROB_PT_arm%d",iarm),"",100,0,1,20,0,10);
		FVTX_Q_PROB_PT_MISMATCH[iarm] = new TH2F(Form("FVTX_Q_PROB_PT_MISMATCH_arm%d",iarm),"",100,0,1,20,0,10);

	}

	return;

}

//______________________________________________________
void mRunCutsDiMuon_Fvtx::load_cuts()
{

}

//______________________________________________________
void mRunCutsDiMuon_Fvtx::set_pT_array()
{
  cout << "mRunCutsDiMuon_Fvtx::set_pT_array" << endl;

	//cout << "DATASET : " << _dataset << ", NPTBIN : " << _nptbin << endl;

	if ( _dataset=="Run12pp510" ){
		_nptbin = 9;
		float tmp_pT_array[10]  = {0.00, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 7.00, 10.00};
		float tmp_p_array[10] = {3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 10.00, 12.00, 15.00, 20.00};
		_pT_array = new double[10];
		_p_array = new double[10];

		for (int ipt=0; ipt<_nptbin+1; ipt++){
			_pT_array[ipt] = tmp_pT_array[ipt];
			_p_array[ipt] = tmp_p_array[ipt];
		}
	}else if ( _dataset=="Run15pp200" ){
		_nptbin = 7;
		float tmp_pT_array[8]  = {0.00, 0.50, 1.00, 1.50, 2.00, 3.00, 5.00, 10.00};
		float tmp_p_array[8] = {3.00, 4.00, 5.00, 6.00, 8.00, 10.00, 15.00, 20.00};
		_pT_array = new double[8];
		_p_array = new double[8];

		for (int ipt=0; ipt<_nptbin+1; ipt++){
			_pT_array[ipt] = tmp_pT_array[ipt];
			_p_array[ipt] = tmp_p_array[ipt];
		}
	}else{
		_nptbin = 7;
		float tmp_pT_array[8]  = {0.00, 0.50, 1.00, 1.50, 2.00, 3.00, 5.00, 10.00};
		float tmp_p_array[8] = {3.00, 4.00, 5.00, 6.00, 8.00, 10.00, 15.00, 20.00};
		_pT_array = new double[8];
		_p_array = new double[8];

		for (int ipt=0; ipt<_nptbin+1; ipt++){
			_pT_array[ipt] = tmp_pT_array[ipt];
			_p_array[ipt] = tmp_p_array[ipt];
		}
	}

	cout << "DATASET : " << _dataset << ", NPTBIN : " << _nptbin << endl;

	return;
}

//______________________________________________________
void mRunCutsDiMuon_Fvtx::set_basic_cuts()
{

  cout << "mRunCutsDiMuon_Fvtx::set_basic_cuts" << endl;
	//cout << "BASIC CUTS FOR " << dataset << endl;

	cut_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/RunCuts/cut_files/Run15pp200_dimuon_tight_cut.root","READ");

	for (int iarm=0; iarm<_narm; iarm++){
		fcut_dr_fvtx[iarm] = (TF1*)cut_file->Get(Form("fcut_dr_fvtx_arm%d_gap4",iarm));
		fcut_dphi_fvtx[iarm] = new TF1(Form("fcut_dphi_fvtx%d",iarm),"[0]+[1]/x",3,50);
		fcut_dphi_fvtx[iarm]->SetParameters(0.03,0.35);
		fcut_dtheta_fvtx[iarm] = new TF1(Form("fcut_dtheta_fvtx%d",iarm),"[0]+[1]/x",3,50);
		fcut_dtheta_fvtx[iarm]->SetParameters(0.02,0.18);
		fcut_dg0[iarm] = (TF1*)cut_file->Get(Form("fcut_dg0_arm%d_gap4",iarm));
		fcut_ddg0[iarm] = (TF1*)cut_file->Get(Form("fcut_ddg0_arm%d_gap4",iarm));
	}


	if ( _dataset=="Run15pp200" ){
		cut_eta[0] = 1.2;
		cut_eta[1] = 2.2;
		cut_pz[0] = 2.5; 
		cut_pz[1] = 2.5;
		cut_ptot = 3.0;
		cut_chi2_muid = 6.;
		cut_chi2_mutr = 15.;
		cut_chi2_evtvtx = 5.;
		cut_pdtheta = 0.3;

		//badrun_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/RunCuts/cut_files/Run15pp200_badrun_for_BJpsi_20161030.root","READ");
		badrun_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/02.datafiles/01.simfiles/mysetup/RunCuts/cut_files/Run15pp200_badrun_for_InclusiveJpsi_20161027.root","READ");
		for (int iarm=0; iarm<_narm; iarm++){
			hbadrun[iarm] = (TH1F*)badrun_file->Get(Form("hBADRUN_arm%d",iarm));
		}

		file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive.root","read");
		if ( file_weight0->IsOpen() ){
			cout << "OPEN: " << file_weight0->GetName() << endl;
			fwt_pT_dAu_fwd = (TF1*)file_weight0->Get("fwt_pT_dAu_fwd");
			fwt_pT_dAu_bwd = (TF1*)file_weight0->Get("fwt_pT_dAu_bwd");
			fwt_y_dAu = (TF1*)file_weight0->Get("fwt_y_dAu");
			if ( !fwt_pT_dAu_fwd || !fwt_pT_dAu_bwd || !fwt_y_dAu ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_pT_dAu_fwd->Print();
				fwt_pT_dAu_bwd->Print();
				fwt_y_dAu->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/05.Jpsi-FVTX/outfile_jpsi_zweight_Run15pp200.root","read");
		if ( file_weight1->IsOpen() ){
			cout << "OPEN: " << file_weight1->GetName() << endl;
			fwt_z = (TF1*)file_weight1->Get("fwt_z");
			if ( !fwt_z ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_z->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

	}else if ( _dataset=="Run15pAl200" ){
		//for South
		cut_z[0][0] = -30.0;
		cut_z[0][1] = 30.0;
		//for North
		cut_z[1][0] = -30.0;
		cut_z[1][1] = 30.0;
		cut_eta[0] = 1.2;
		cut_eta[1] = 2.2;
		cut_pz[0] = 2.0; 
		cut_pz[1] = 2.0;
		cut_ptot = 3.0;
		cut_chi2_muid = 6.;
		cut_chi2_mutr = 23.;
		//cut_chi2_mutr = 30.;
		cut_chi2_evtvtx = 5.;
		cut_pdtheta = 0.3;

		badrun_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/10.runQA/03.yield_run/Run15pAl200_badrun_for_dimuon_20170810.root","READ");
		for (int iarm=0; iarm<_narm; iarm++){
			hbadrun[iarm] = (TH1F*)badrun_file->Get(Form("hBADRUN_arm%d",iarm));
		}

		file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive.root","read");
		if ( file_weight0->IsOpen() ){
			fwt_pT_dAu_fwd = (TF1*)file_weight0->Get("fwt_pT_dAu_fwd");
			fwt_pT_dAu_bwd = (TF1*)file_weight0->Get("fwt_pT_dAu_bwd");
			fwt_pT_pp = (TF1*)file_weight0->Get("fwt_pT_pp");
			fwt_y_dAu = (TF1*)file_weight0->Get("fwt_y_dAu");
			fwt_y_pp = (TF1*)file_weight0->Get("fwt_y_pp");
			if ( !fwt_pT_dAu_fwd || !fwt_pT_dAu_bwd || !fwt_y_dAu || !fwt_pT_pp || !fwt_y_pp ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_pT_dAu_fwd->Print();
				fwt_pT_dAu_bwd->Print();
				fwt_pT_pp->Print();
				fwt_y_dAu->Print();
				fwt_y_pp->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run15pAl200.root","read");
		if ( file_weight1->IsOpen() ){
			fwt_z = (TF1*)file_weight1->Get("fwt_z");
			if ( !fwt_z ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_z->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		rate_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/07.genericQA/02.BBC_rates/Run15pAl200_BBC_rate.root","read");
		hBBC_rate = (TH1F*)rate_file->Get("hBBC_rate");

	}else if ( _dataset=="Run15pAu200" ){
		//for South
		cut_z[0][0] = -30.0;
		cut_z[0][1] = 30.0;
		//for North
		cut_z[1][0] = -30.0;
		cut_z[1][1] = 30.0;
		cut_eta[0] = 1.2;
		cut_eta[1] = 2.2;
		cut_pz[0] = 2.0; 
		cut_pz[1] = 2.0;
		cut_ptot = 3.0;
		cut_chi2_muid = 6.;
		cut_chi2_mutr = 23.;
		cut_chi2_evtvtx = 5.;
		cut_pdtheta = 0.3;

		badrun_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/10.runQA/03.yield_run/Run15pAu200_badrun_for_dimuon_20170906.root","READ");
		for (int iarm=0; iarm<_narm; iarm++){
			hbadrun[iarm] = (TH1F*)badrun_file->Get(Form("hBADRUN_arm%d",iarm));
		}

		file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive.root","read");
		if ( file_weight0->IsOpen() ){
			cout << "OPEN: " << file_weight0->GetName() << endl;
			fwt_pT_dAu_fwd = (TF1*)file_weight0->Get("fwt_pT_dAu_fwd");
			fwt_pT_dAu_bwd = (TF1*)file_weight0->Get("fwt_pT_dAu_bwd");
			fwt_pT_pp = (TF1*)file_weight0->Get("fwt_pT_pp");
			fwt_y_dAu = (TF1*)file_weight0->Get("fwt_y_dAu");
			fwt_y_pp = (TF1*)file_weight0->Get("fwt_y_pp");
			if ( !fwt_pT_dAu_fwd || !fwt_pT_dAu_bwd || !fwt_y_dAu || !fwt_pT_pp || !fwt_y_pp ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_pT_dAu_fwd->Print();
				fwt_pT_dAu_bwd->Print();
				fwt_pT_pp->Print();
				fwt_y_dAu->Print();
				fwt_y_pp->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run15pAu200.root","read");
		if ( file_weight1->IsOpen() ){
			cout << "OPEN: " << file_weight1->GetName() << endl;
			fwt_z = (TF1*)file_weight1->Get("fwt_z");
			fwt_z_S = (TF1*)file_weight1->Get("fwt_z_S");
			fwt_z_N = (TF1*)file_weight1->Get("fwt_z_N");
			if ( !fwt_z ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_z->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		rate_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/07.genericQA/02.BBC_rates/Run15pAu200_BBC_rate.root","read");
		hBBC_rate = (TH1F*)rate_file->Get("hBBC_rate");

	}else if ( _dataset=="Run14HeAu200" ){
		//for South
		cut_z[0][0] = -30.0;
		cut_z[0][1] = 30.0;
		//for North
		cut_z[1][0] = -30.0;
		cut_z[1][1] = 30.0;
		cut_eta[0] = 1.2;
		cut_eta[1] = 2.2;
		cut_pz[0] = 2.0; 
		cut_pz[1] = 2.0;
		cut_ptot = 3.0;
		cut_chi2_muid = 6.;
		cut_chi2_mutr = 23.;
		cut_chi2_evtvtx = 5.;
		cut_pdtheta = 0.3;
		_ncentbin = 4;

		badrun_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/10.runQA/03.yield_run/Run14HeAu200_badrun_for_dimuon_20171017.root","READ");
		for (int iarm=0; iarm<_narm; iarm++){
			hbadrun[iarm] = (TH1F*)badrun_file->Get(Form("hBADRUN_arm%d",iarm));
		}

		file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive.root","read");
		if ( file_weight0->IsOpen() ){
			fwt_pT_dAu_fwd = (TF1*)file_weight0->Get("fwt_pT_dAu_fwd");
			fwt_pT_dAu_bwd = (TF1*)file_weight0->Get("fwt_pT_dAu_bwd");
			fwt_pT_pp = (TF1*)file_weight0->Get("fwt_pT_pp");
			fwt_y_dAu = (TF1*)file_weight0->Get("fwt_y_dAu");
			fwt_y_pp = (TF1*)file_weight0->Get("fwt_y_pp");
			if ( !fwt_pT_dAu_fwd || !fwt_pT_dAu_bwd || !fwt_y_dAu || !fwt_pT_pp || !fwt_y_pp ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_pT_dAu_fwd->Print();
				fwt_pT_dAu_bwd->Print();
				fwt_pT_pp->Print();
				fwt_y_dAu->Print();
				fwt_y_pp->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_zweight_Run14HeAu200.root","read");
		if ( file_weight1->IsOpen() ){
			fwt_z = (TF1*)file_weight1->Get("fwt_z");
			if ( !fwt_z ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_z->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		rate_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/07.genericQA/02.BBC_rates/Run14HeAu200_BBC_rate.root","read");
		hBBC_rate = (TH1F*)rate_file->Get("hBBC_rate");

		file_trigeff = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/10.runQA/01.trig_eff/Run14HeAu200_trigeff_wt.root","read");
		for (int iarm=0; iarm<_narm; iarm++){
			for (int icent=0; icent<_ncentbin; icent++){
				ftrigeff[iarm][icent] = (TF1*)file_trigeff->Get(Form("fMUID1D_wt_eta_centbin%d_arm%d",icent,iarm));
				if ( !ftrigeff[iarm][icent] ){
					cout << "CAN NOT GET TRIG weight functions!!" << endl;
					exit(1);
				}
			}
		}

	}else if ( _dataset=="Run9pp500" ){

		//for South
		cut_z[0][0] = -30.0;
		cut_z[0][1] = 30.0;
		//for North
		cut_z[1][0] = -30.0;
		cut_z[1][1] = 30.0;
		cut_eta[0] = 1.2;
		cut_eta[1] = 2.2;
		cut_pz[0] = 2.0; 
		cut_pz[1] = 2.0;
		cut_ptot = 3.0;
		cut_chi2_muid = 6.;
		cut_chi2_mutr = 23.;
		cut_chi2_evtvtx = 5.;
		cut_pdtheta = 0.3;

		badrun_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/04.run9_pdst/12.dimuon/runQA/Run9pp500_badrun_for_dimuon_20180128.root","READ");
		for (int iarm=0; iarm<_narm; iarm++){
			hbadrun[iarm] = (TH1F*)badrun_file->Get(Form("hBADRUN_arm%d",iarm));
		}

		file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive.root","read");
		if ( file_weight0->IsOpen() ){
			fwt_pT_dAu_fwd = (TF1*)file_weight0->Get("fwt_pT_dAu_fwd");
			fwt_pT_dAu_bwd = (TF1*)file_weight0->Get("fwt_pT_dAu_bwd");
			fwt_pT_pp = (TF1*)file_weight0->Get("fwt_pT_pp");
			fwt_y_dAu = (TF1*)file_weight0->Get("fwt_y_dAu");
			fwt_y_pp = (TF1*)file_weight0->Get("fwt_y_pp");
			if ( !fwt_pT_dAu_fwd || !fwt_pT_dAu_bwd || !fwt_y_dAu || !fwt_pT_pp || !fwt_y_pp ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_pT_dAu_fwd->Print();
				fwt_pT_dAu_bwd->Print();
				fwt_pT_pp->Print();
				fwt_y_dAu->Print();
				fwt_y_pp->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/04.run9_pdst/12.dimuon/acceff/outfile_jpsi_zweight_Run9pp500.root","read");
		if ( file_weight1->IsOpen() ){
			fwt_z = (TF1*)file_weight1->Get("fwt_z");
			if ( !fwt_z ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_z->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

	}else if ( _dataset=="Run9pp200" ){

		//for South
		cut_z[0][0] = -30.0;
		cut_z[0][1] = 30.0;
		//for North
		cut_z[1][0] = -30.0;
		cut_z[1][1] = 30.0;
		cut_eta[0] = 1.2;
		cut_eta[1] = 2.2;
		cut_pz[0] = 2.0; 
		cut_pz[1] = 2.0;
		cut_ptot = 3.0;
		cut_chi2_muid = 6.;
		cut_chi2_mutr = 15.;
		cut_chi2_evtvtx = 6.;
		cut_pdtheta = 0.3;

		badrun_file = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/04.run9_pdst/12.dimuon/runQA/Run9pp200_badrun_for_dimuon_20180128.root","READ");
		for (int iarm=0; iarm<_narm; iarm++){
			hbadrun[iarm] = (TH1F*)badrun_file->Get(Form("hBADRUN_arm%d",iarm));
			cout << "ARM: " << iarm << ", BAD RUN: " << hbadrun[iarm]->Integral() << endl;
		}

		file_weight0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/14.acc_by_eff/01.Inclusive_Jpsi/outfile_jpsi_weight_for_inclusive.root","read");
		if ( file_weight0->IsOpen() ){
			fwt_pT_dAu_fwd = (TF1*)file_weight0->Get("fwt_pT_dAu_fwd");
			fwt_pT_dAu_bwd = (TF1*)file_weight0->Get("fwt_pT_dAu_bwd");
			fwt_pT_pp = (TF1*)file_weight0->Get("fwt_pT_pp");
			fwt_y_dAu = (TF1*)file_weight0->Get("fwt_y_dAu");
			fwt_y_pp = (TF1*)file_weight0->Get("fwt_y_pp");
			if ( !fwt_pT_dAu_fwd || !fwt_pT_dAu_bwd || !fwt_y_dAu || !fwt_pT_pp || !fwt_y_pp ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_pT_dAu_fwd->Print();
				fwt_pT_dAu_bwd->Print();
				fwt_pT_pp->Print();
				fwt_y_dAu->Print();
				fwt_y_pp->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

		file_weight1 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/04.run9_pdst/12.dimuon/acceff/outfile_jpsi_zweight_Run9pp200.root","read");
		if ( file_weight1->IsOpen() ){
			fwt_z = (TF1*)file_weight1->Get("fwt_z");
			if ( !fwt_z ){
				cout << "CAN NOT GET weight functions!!" << endl;
				exit(1);
			}else{
				fwt_z->Print();
			}
		}else{
			cout << "CAN NOT OPEN FILE containing weight functions!!" << endl;
		}

	}else if ( _dataset=="Run14AuAu200" ){
		cut_eta[0] = 1.2;
		cut_eta[1] = 2.2;
		cut_pz[0] = 2.0; 
		cut_pz[1] = 2.0;
		cut_chi2_muid = 6.;
		cut_chi2_mutr = 12.;
		cut_chi2_evtvtx = 6.;
		cut_pdtheta = 0.3;
	}//Run-dependent 

	cout << "ETA CUT : " << cut_eta[0] << " < |eta| <  " << cut_eta[1] << endl;

	return;
}

//______________________________________________________
bool mRunCutsDiMuon_Fvtx::check_event()
{

  //cout << "mRunCutsDiMuon_Fvtx::check_event" << endl;

	if ( _is_sim ){
		if ( _arm==0 ){
			if ( _trig_emul_2D_S==0 ) return false;
		}else{
			if ( _trig_emul_2D_N==0 ) return false;
		}
		return true;
	}

	if ( fabs(_bbcZ)>30.0 ) return false;
	if ( _dataset=="Run14AuAu200" ){
		if ( isnan(_fvtxZ) || fabs(_fvtxZ)>30.0 ) return false;
	}

	//if ( !_same_event ) return true;

	//don't use for ratio study
	//trigger cut
	/*
	if ( _dataset=="Run15pp200" ){
		if ( _triggertype=="MUID2D" && !(_scaled_trigbit&0x00100000) && !(_scaled_trigbit&0x00200000) ) return false;
		if ( _triggertype=="MUID1D" && !(_scaled_trigbit&0x00400000) && !(_scaled_trigbit&0x00800000) ) return false;
		if ( _triggertype=="SG3_MUID1DH" && !(_scaled_trigbit&0x01000000) && !(_scaled_trigbit&0x02000000) ) return false;
		if ( _triggertype=="COMBINED" && !(_scaled_trigbit&0x01000000) && !(_scaled_trigbit&0x02000000) && !(_scaled_trigbit&0x00100000) && !(_scaled_trigbit&0x00200000) && !(_scaled_trigbit&0x00400000) && !(_scaled_trigbit&0x00800000) ) return false;
	}else if ( _dataset=="Run14AuAu200" ){
		if ( _centrality<=40 || _centrality>94 ) return false;
	}else if ( _dataset=="Run15pAu200" || _dataset=="Run15pAl200" ){
		if ( _triggertype=="MUID2D" ){
			if ( _arm==1 && !(_scaled_trigbit&0x00100000) ) return false;
			if ( _arm==0 && !(_scaled_trigbit&0x00200000) ) return false;
		}
	}else if ( _dataset=="Run14HeAu200" ){
		if ( _triggertype=="MUID2D" ){
			if ( _arm==1 && !(_scaled_trigbit&0x00040000) ) return false;
			if ( _arm==0 && !(_scaled_trigbit&0x00080000) ) return false;
		}
	}else if ( _dataset=="Run9pp500" || _dataset=="Run9pp200" ){
		if ( _triggertype=="MUID2D" ){
			if ( !(_scaled_trigbit&0x00080000) ) return false;
		}
	}
	*/

  //cout << "mRunCutsDiMuon_Fvtx::**" << endl;

	return true;
}


//______________________________________________________
bool mRunCutsDiMuon_Fvtx::check_dimuon()
{

	if ( isnan(_mass) || isnan(_pT) ) return false;
	if ( _pT<_pT_array[0] || _pT>_pT_array[_nptbin] ) return false;
	if ( _mass<cut_mass[0] || _mass>cut_mass[1] ) return false;

	if ( fabs(_rapidity)<cut_eta[0] || fabs(_rapidity)>cut_eta[1] ) return false;

	if ( (fabs(_tr_rapidity[0])>2.2+0.2*_arm) || (fabs(_tr_rapidity[1])>2.2+0.2*_arm) || fabs(_tr_rapidity[0])<1.2 || fabs(_tr_rapidity[1])<1.2 ) return false;

	if ( _tr_trchi2[0]>cut_chi2_mutr || _tr_trchi2[1]>cut_chi2_mutr ) return false; 
	if ( fabs(_tr_pz[0])<cut_pz[_arm] || fabs(_tr_pz[1])<cut_pz[_arm] ) return false; 
	if ( _tr_pz[0]*_tr_pz[1]<0 ) return false;

	if ( _dataset=="Run15pp200" || _dataset=="Run15pAu200" || _dataset=="Run15pAl200" || _dataset=="Run14HeAu200" ){
		if ( _tr_lastgap[0]<2.5 || _tr_lastgap[1]<2.5 ) return false;
		if ( _tr_idhits[0]<14 || _tr_idhits[1]<14 ) return false;
		if ( !check_mutr_fiducial(0) || !check_mutr_fiducial(1) ) return false; //don't use for ratio study
	}else if ( _dataset=="Run14AuAu200" ){
		if ( _tr_lastgap[0]<3.5 && _tr_lastgap[1]<3.5 ) return false;
		if ( sqrt(_tr_pz[0]*_tr_pz[0]+_tr_pT[0]*_tr_pT[0])<3.0 || sqrt(_tr_pz[1]*_tr_pz[1]+_tr_pT[1]*_tr_pT[1])<3.0 ) return false; 
		if ( _evt_vtxchi2>cut_chi2_evtvtx ) return false;
		if ( _tr_ntrhits[0]<11.5 || _tr_ntrhits[1]<11.5 ) return false;
	}

	//if ( _tr_idchi2[0]>cut_chi2_muid && _tr_idchi2[1]>cut_chi2_muid ) return false; //AND for dimuon

	for (int ii=0; ii<2; ii++){
		float ptot = sqrt(_tr_pT[ii]*_tr_pT[ii] + _tr_pz[ii]*_tr_pz[ii]);
		float cut_dg0 = 80./ptot;
		float cut_ddg0 = 40./ptot;
		cut_dg0 = fcut_dg0[_arm]->Eval(ptot);
		cut_ddg0 = fcut_ddg0[_arm]->Eval(ptot);

		if ( _dataset=="Run15pAu200" ){
			if ( _rapidity<0 ){
				//South
				cut_dg0 = 9.71955+381.68/ptot;
				cut_ddg0 = 5.26039+130.321/ptot;
			}else{
				//North
				cut_dg0 = 8.29053+113.65/ptot;
				cut_ddg0 = 4.42153+149.754/ptot;
			}
		}else if ( _dataset=="Run15pAl200" ){
			if ( _rapidity<0 ){
				//South
				cut_dg0 = 8.27264+396.456/ptot;
				cut_ddg0 = 4.41174+136.571/ptot;
			}else{
				//North
				cut_dg0 = 7.87927+119.362/ptot;
				cut_ddg0 = 4.23478+143.84/ptot;
			}
		}else if ( _dataset=="Run14HeAu200" ){
			if ( _rapidity<0 ){
				//South
				cut_dg0 = 7.97212+426.77/ptot;
				cut_ddg0 = 4.78783+118.559/ptot;
			}else{
				//North
				cut_dg0 = 8.34527+105.503/ptot;
				cut_ddg0 = 4.04693+139.679/ptot;
			}
		}else if ( _dataset=="Run15pp200" ){
			if ( _rapidity<0 ){
				//South
				cut_dg0 = 8.85659+388.943/ptot;
				cut_ddg0 = 4.44408+132.448/ptot;
			}else{
				//North
				cut_dg0 = 8.68861+111.857/ptot;
				cut_ddg0 = 4.3318+154.743/ptot;
			}
		}

		if ( _tr_DG0[ii]>cut_dg0 ) return false;
		if ( _tr_DDG0[ii]>cut_ddg0 ) return false;
	}//dg0/ddg0 cuts

	if ( _dataset=="Run15pp200" && _dataset=="Run15pAl200" && _dataset=="Run15pAu200" && _dataset=="Run14HeAu200" && _dataset=="Run14AuAu200" ){
		if ( abs((((int((atan2(_tr_xst1[0],_tr_yst1[0])+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8) - (((int((atan2(_tr_xst1[1],_tr_yst1[1])+(4*atan(1)))/((4*atan(1))/8))+1)/2)%8))==0 ) return false;
	}


	if ( _is_sim ){
		//if ( (abs(_tr_mc_g_pid[0])!=13 || abs(_tr_mc_g_pid[1])!=13) ) return false;
		if ( _filetype.find("JPSI") != std::string::npos ){
			if ( _nhepmc<1 ) return false;
		}
		if ( _tr_mc_hits_mutr_true[0]*1.0/_tr_ntrhits[0]<0.6 ) return false;
		if ( _tr_mc_hits_mutr_true[1]*1.0/_tr_ntrhits[1]<0.6 ) return false;
		if ( _tr_mc_hits_muid_true[0]*1.0/_tr_nidhits[0]<0.6 ) return false;
		if ( _tr_mc_hits_muid_true[1]*1.0/_tr_nidhits[1]<0.6 ) return false;
	}//_is_sim

	return true;
}

//______________________________________________________
bool mRunCutsDiMuon_Fvtx::check_fvtx(int tr_index)
{

	//return true;
	if ( _mass_fvtxmutr<cut_mass[0] || _mass_fvtxmutr>cut_mass[1] ) return false;

	float ptot = sqrt(_tr_pT[tr_index]*_tr_pT[tr_index] + _tr_pz[tr_index]*_tr_pz[tr_index]);
	float cut_dr_fvtx = (_rapidity<0) ? 1.60504+31.1032/ptot/ptot : 1.56119+28.5809/ptot/ptot; 
	float cut_dtheta_fvtx = (_rapidity<0) ? 0.0617706+0.616019/ptot/ptot : 0.0522562+0.621762/ptot/ptot;
	float cut_dphi_fvtx = (_rapidity<0) ? 0.11861+1.84085/ptot/ptot : 0.127857+1.54661/ptot/ptot; 

	if ( !(_tr_chi2_fvtx[tr_index]>0 && _tr_chi2_fvtx[tr_index]<10.0) ) return false;
	if ( !(_tr_chi2_fvtxmutr[tr_index]>0 && _tr_chi2_fvtxmutr[tr_index]<10.0) ) return false;
	if ( fabs(_tr_dr_fvtx[tr_index])>cut_dr_fvtx ) return false;
	if ( fabs(_tr_dtheta_fvtx[tr_index])>cut_dtheta_fvtx ) return false;
	if ( fabs(_tr_dphi_fvtx[tr_index])>cut_dphi_fvtx ) return false;

	float pvalue = TMath::Prob(_tr_chi2_fvtx[tr_index]*(2*_tr_nhits_fvtx[tr_index]-5), 2*_tr_nhits_fvtx[tr_index]-5);
	if ( pvalue<0.01 ) return false;

	for (int ii=0; ii<2; ii++){

		//float pvalue = TMath::Prob(_tr_chi2_fvtx[ii]*(2*_tr_nhits_fvtx[ii]-5), 2*_tr_nhits_fvtx[ii]-5);
		//if ( pvalue<0.01 ) return false;

		/*
		if ( _tr_hit_pattern[ii]>256 ){
			if ( _rapidity<0 ){
				if ( _tr_nhits_fvtx[ii]<3.5 ) return false;
			}else{
				if ( _tr_nhits_fvtx[ii]<3.5 ) return false;
			}
		}else{
			if ( _tr_nhits_fvtx[ii]<2.5 ) return false;
		}
		*/
	}

	//if ( _tr_nhits_fvtx[0]<3.5 || _tr_nhits_fvtx[1]<3.5 ) return false;

	return true;
}

//______________________________________________________
bool mRunCutsDiMuon_Fvtx::check_mutr_fiducial(int index){

	//good for south arm
	//if ( _tr_rapidity[index]<0 ) return true;

	if ( _dataset=="Run15pAu200" || _dataset=="Run15pAl200" || _dataset=="Run14HeAu200" ){
		if ( _tr_rapidity[index]<0 ){
			//South
			if ( atan2(_tr_yst1[index],_tr_xst1[index])>0.00 && atan2(_tr_yst1[index],_tr_xst1[index])<0.40 ) return false;
			if ( atan2(_tr_yst2[index],_tr_xst2[index])>0.00 && atan2(_tr_yst2[index],_tr_xst2[index])<0.40 ) return false;
			if ( atan2(_tr_yst3[index],_tr_xst3[index])>0.00 && atan2(_tr_yst3[index],_tr_xst3[index])<0.40 ) return false;
			if ( atan2(_tr_yst1[index],_tr_xst1[index])>1.55 && atan2(_tr_yst1[index],_tr_xst1[index])<1.90 ) return false;
			if ( atan2(_tr_yst2[index],_tr_xst2[index])>1.55 && atan2(_tr_yst2[index],_tr_xst2[index])<1.90 ) return false;
			if ( atan2(_tr_yst3[index],_tr_xst3[index])>1.55 && atan2(_tr_yst3[index],_tr_xst3[index])<1.90 ) return false;
			if ( atan2(_tr_yst3[index],_tr_xst3[index])>-0.60 && atan2(_tr_yst3[index],_tr_xst3[index])<-0.40 ) return false;

			//if ( _dataset=="Run14HeAu200" && atan2(_tr_yst1[index],_tr_xst1[index])>-2.35 && atan2(_tr_yst1[index],_tr_xst1[index])<-1.95 ) return false;
			//if ( _dataset=="Run14HeAu200" && atan2(_tr_yst2[index],_tr_xst2[index])>-2.35 && atan2(_tr_yst2[index],_tr_xst2[index])<-1.95 ) return false;
			//if ( _dataset=="Run14HeAu200" && atan2(_tr_yst3[index],_tr_xst3[index])>-2.35 && atan2(_tr_yst3[index],_tr_xst3[index])<-1.95 ) return false;
			if ( atan2(_tr_yst1[index],_tr_xst1[index])>-2.35 && atan2(_tr_yst1[index],_tr_xst1[index])<-1.95 ) return false;
			if ( atan2(_tr_yst2[index],_tr_xst2[index])>-2.35 && atan2(_tr_yst2[index],_tr_xst2[index])<-1.95 ) return false;
			if ( atan2(_tr_yst3[index],_tr_xst3[index])>-2.35 && atan2(_tr_yst3[index],_tr_xst3[index])<-1.95 ) return false;

		}else{
			//North
			if ( atan2(_tr_yst3[index],_tr_xst3[index])>0.35 && atan2(_tr_yst3[index],_tr_xst3[index])<0.77 ) return false;
			if ( atan2(_tr_yst3[index],_tr_xst3[index])>0.77 && atan2(_tr_yst3[index],_tr_xst3[index])<1.15 
					&& sqrt(_tr_yst3[index]*_tr_yst3[index]+_tr_xst3[index]*_tr_xst3[index])>160 && sqrt(_tr_yst3[index]*_tr_yst3[index]+_tr_xst3[index]*_tr_xst3[index])<275 ) return false;
			if ( atan2(_tr_yst3[index],_tr_xst3[index])>-2.75 && atan2(_tr_yst3[index],_tr_xst3[index])<-1.92 
					&& sqrt(_tr_yst3[index]*_tr_yst3[index]+_tr_xst3[index]*_tr_xst3[index])>237 && sqrt(_tr_yst3[index]*_tr_yst3[index]+_tr_xst3[index]*_tr_xst3[index])<275 ) return false;
			if ( atan2(_tr_yst2[index],_tr_xst2[index])>-2.75 && atan2(_tr_yst2[index],_tr_xst2[index])<-1.92 
					&& sqrt(_tr_yst2[index]*_tr_yst2[index]+_tr_xst2[index]*_tr_xst2[index])>140 && sqrt(_tr_yst2[index]*_tr_yst2[index]+_tr_xst2[index]*_tr_xst2[index])<163 ) return false;
			if ( atan2(_tr_yst1[index],_tr_xst1[index])>-2.75 && atan2(_tr_yst1[index],_tr_xst1[index])<-2 
					&& sqrt(_tr_yst1[index]*_tr_yst1[index]+_tr_xst1[index]*_tr_xst1[index])>76 && sqrt(_tr_yst1[index]*_tr_yst1[index]+_tr_xst1[index]*_tr_xst1[index])<80 ) return false;
			if ( atan2(_tr_yst1[index],_tr_xst1[index])>0.52 && atan2(_tr_yst1[index],_tr_xst1[index])<1.11 
					&& sqrt(_tr_yst1[index]*_tr_yst1[index]+_tr_xst1[index]*_tr_xst1[index])>60 && sqrt(_tr_yst1[index]*_tr_yst1[index]+_tr_xst1[index]*_tr_xst1[index])<75 ) return false;
			if ( atan2(_tr_yst2[index],_tr_xst2[index])>-1.20 && atan2(_tr_yst2[index],_tr_xst2[index])<-0.90 ) return false;
		}
	}else if ( _dataset=="Run9pp200" ){
		if ( _tr_rapidity[index]>0 ){
			if ( atan2(_tr_yst2[index],_tr_xst2[index])>0.4 && atan2(_tr_yst2[index],_tr_xst2[index])<1.2
					&& sqrt(_tr_yst2[index]*_tr_yst2[index]+_tr_xst2[index]*_tr_xst2[index])<100 ) return false;

			if ( atan2(_tr_yst2[index],_tr_xst2[index])>-1.3 && atan2(_tr_yst2[index],_tr_xst2[index])<-0.4
					&& sqrt(_tr_yst2[index]*_tr_yst2[index]+_tr_xst2[index]*_tr_xst2[index])>195 ) return false;
		}
	}else{
		return true;
	}

	return true;
}

//______________________________________________________
bool mRunCutsDiMuon_Fvtx::check_run(int arm){

	if ( _is_sim ) return true;
	if ( _dataset!="Run15pp200" && _dataset!="Run15pAu200" && _dataset!="Run14HeAu200" && _dataset!="Run15pAl200" && _dataset!="Run9pp500" && _dataset!="Run9pp200" ) return true;

	if ( _prev_run!=_runnumber ){

		for (int iarm=0; iarm<_narm; iarm++){
			float val = hbadrun[iarm]->GetBinContent(hbadrun[iarm]->FindBin(_runnumber));
			if ( val>0.5 ){
				cout << "BADRUN ARM: " << iarm << ", RUN #: " << _runnumber << endl;
				_is_goodrun[iarm] = false;
			}else{
				_is_goodrun[iarm] = true;
			}
		}//iarm

		_prev_run = _runnumber;
	}

	return _is_goodrun[arm];
}
