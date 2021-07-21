#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TChain.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TVector3.h>
#include <iostream>
#include <cstring>
#include <TMath.h>
#include <fstream>

#include "/phenix/u/shlim/RunRange.h"

using namespace std;

void filter_hadron_fvtx_pp(string dataset="Run15pp200", string trigger="MB", bool bwide_eta=true){

	const int narm = 2;
	const int ngap = 3;
	const int nchg = 3;

	const int nptbin = 13; 
	double pT_array[nptbin+1] = {1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.0};

	const int netabin = 6; 
	//double eta_array[netabin+1] = {1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

	//For Gap3
	float pz_cut[narm][netabin] = {{4.5,4,4,4,4,4}, {4.5,4,4,4,4,4}};

	float cut_ddg0[narm] = {18, 16};
	float cut_dg0[narm] = {25, 20};

	//float cut_pdtheta[nptbin] = {0.20, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36};
	//float cut_vtxchi2[nptbin] = {5.00, 4.50, 4.50, 4.50, 4.50, 4.50, 4.50, 4.00, 3.50, 3.00};
	//float cut_pdtheta[narm] = {0.2, 0.2}; //cut_pdtheta[iarm]+0.2*ptbin
	//float cut_vtxchi2[narm] = {3, 3}; //cut_vtxchi2[iarm]-0.2*ptbin
	
	//bad run list
	TFile *infile_badrun = new TFile(Form("%s_badrun.root",dataset.c_str()),"read");
	TH1F *hbadrun[narm];
	if ( infile_badrun->IsOpen() ){
		cout << "OPEN BAD RUN FILE: " << infile_badrun->GetName() << endl;
		for (int iarm=0; iarm<narm; iarm++){
			hbadrun[iarm] = (TH1F*)infile_badrun->Get(Form("hBADRUN_arm%d",iarm));
			cout << "Number of bad runs: " << hbadrun[iarm]->GetEntries() << endl;
		}
	}else{
		cout << "CAN NOT OPEN BADRUN FILE!" << endl;
		return;
	}

	//trigger efficiency
	TFile *infile_trig = new TFile("Run15pp200_hadron_trig_eff_func.root","read");
	TF1 *ftrig_pT[narm][ngap];
	TH1F *htrig_eta[narm][ngap];

	for (int iarm=0; iarm<narm; iarm++){
		for (int igap=0; igap<ngap; igap++){
			ftrig_pT[iarm][igap] = (TF1*)infile_trig->Get(Form("ftrig_eff_pT_SG3_MUID1DH_arm%d_gap%d",iarm,igap+2));
			if ( !ftrig_pT[iarm][igap] ){
				cout << "CAN NOT FIND ftrig_eff_pT, ARM" << iarm << ", GAP" << igap+2 << endl;
				exit(1);
			}
		}
		for (int igap=0; igap<ngap; igap++){
			htrig_eta[iarm][igap] = (TH1F*)infile_trig->Get(Form("htrig_eff_eta_SG3_MUID1DH_arm%d_gap%d",iarm,igap+2));
			if ( !htrig_eta[iarm][igap] ){
				cout << "CAN NOT FIND htrig_eff_eta, ARM" << iarm << ", GAP" << igap+2 << endl;
				exit(1);
			}
		}
	}

	int run_range[2];
	get_run_range(dataset, run_range);

	//return;

	TH2F *hPT_ETA[narm][ngap][nchg];
	TH2F *hPT_ETA_cor_pT[narm][ngap][nchg];
	TH2F *hPT_ETA_cor_eta[narm][ngap][nchg];
	TH2F *hPT_ETA_cor_all[narm][ngap][nchg];
	TH2F *hPT_ETA_tight[narm][ngap][nchg];
	TH2F *hPT_ETA_loose[narm][ngap][nchg];
	TProfile *hPT_mean[narm][ngap][nchg];

	TH2F *hPT_ETA_fine[narm][ngap][nchg];
	TH2F *hPT_ETA_fine_cor_eta[narm][ngap][nchg];
	TH2F *hPT_ETA_fine_cor_pT[narm][ngap][nchg];
	TH2F *hPT_ETA_fine_cor_all[narm][ngap][nchg];

	TH2F *hRPHI_MUTR[narm];
	TH2F *hRPHI_FVTX[narm];

	TH1F *hRUN[narm][ngap];

	//cout << "check" << endl;

	for (int iarm=0; iarm<narm; iarm++){
		for (int igap=0; igap<ngap; igap++){
			for (int ichg=0; ichg<nchg; ichg++){
				hPT_ETA[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,12,1.2,2.4);
				hPT_ETA[iarm][igap][ichg]->Sumw2();

				hPT_ETA_cor_pT[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_cor_pT_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,12,1.2,2.4);
				hPT_ETA_cor_pT[iarm][igap][ichg]->Sumw2();

				hPT_ETA_cor_eta[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_cor_eta_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,12,1.2,2.4);
				hPT_ETA_cor_eta[iarm][igap][ichg]->Sumw2();

				hPT_ETA_cor_all[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_cor_all_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,12,1.2,2.4);
				hPT_ETA_cor_all[iarm][igap][ichg]->Sumw2();

				hPT_ETA_tight[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_tight_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,12,1.2,2.4);
				hPT_ETA_tight[iarm][igap][ichg]->Sumw2();

				hPT_ETA_loose[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_loose_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,12,1.2,2.4);
				hPT_ETA_loose[iarm][igap][ichg]->Sumw2();

				hPT_ETA_fine[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_fine_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,24,1.2,2.4);
				hPT_ETA_fine[iarm][igap][ichg]->Sumw2();
				hPT_ETA_fine_cor_eta[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_fine_cor_eta_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,24,1.2,2.4);
				hPT_ETA_fine_cor_eta[iarm][igap][ichg]->Sumw2();
				hPT_ETA_fine_cor_pT[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_fine_cor_pT_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,24,1.2,2.4);
				hPT_ETA_fine_cor_pT[iarm][igap][ichg]->Sumw2();
				hPT_ETA_fine_cor_all[iarm][igap][ichg] = new TH2F(Form("hPT_ETA_fine_cor_all_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array,24,1.2,2.4);
				hPT_ETA_fine_cor_all[iarm][igap][ichg]->Sumw2();

				hPT_mean[iarm][igap][ichg] = new TProfile(Form("hPT_mean_arm%d_gap%d_chg%d",iarm,igap+2,ichg),"",nptbin,pT_array);
				hPT_mean[iarm][igap][ichg]->Sumw2();
			}

			hRUN[iarm][igap] = new TH1F(Form("hRUN_arm%d_gap%d",iarm,igap+2),"",run_range[1]-run_range[0],run_range[0],run_range[1]);
		}

		hRPHI_MUTR[iarm] = new TH2F(Form("hRPHI_MUTR_arm%d",iarm),"",96,-TMath::Pi(),TMath::Pi(),200,0,200);
		hRPHI_FVTX[iarm] = new TH2F(Form("hRPHI_FVTX_arm%d",iarm),"",96,-TMath::Pi(),TMath::Pi(),200,0,20);
	}

	char fname[300];
	sprintf(fname,"%s%s_femtoDST.lst",dataset.c_str(),trigger.c_str());

	ifstream flist;
	flist.open(fname);

	TChain *fChain = new TChain("ana_tree");

	while (flist >> fname){
		cout << "OPEN : " << fname << endl;
		fChain->AddFile(fname);
	}
	flist.close();

	long nentries = fChain->GetEntries();
	cout << "NENTREIS : " << nentries << endl;

	int _arm, _runnumber;
	float _pT, _pz, _eta;
	float _DG0, _DDG0, _rf_vtxchi2pdf, _pdtheta, _rf_slope;
	float _trchi2, _idchi2;
	int _ntrhits, _nidhits, _charge, _lastgap;
	float _bbcZ, _fvtxZ;
	unsigned int _scaled_trigbit;
	int _nhits_fvtx, _nhits_fvtx_charge, _hit_pattern;
	float _chi2_fvtxmutr, _chi2_fvtx;
	float _px_fvtxmutr, _py_fvtxmutr, _pz_fvtxmutr;
	float _dr_fvtx, _dphi_fvtx, _dtheta_fvtx;
	int _mult_fvtxN, _mult_fvtxS;

	float _xst1, _xst2, _xst3;
	float _yst1, _yst2, _yst3;

	float _projx_fvtxmutr, _projy_fvtxmutr;

	fChain->SetBranchAddress("_runnumber",&_runnumber);
	fChain->SetBranchAddress("_arm",&_arm);
	fChain->SetBranchAddress("_pT",&_pT);
	fChain->SetBranchAddress("_pz",&_pz);
	fChain->SetBranchAddress("_eta",&_eta);
	fChain->SetBranchAddress("_DG0",&_DG0);
	fChain->SetBranchAddress("_DDG0",&_DDG0);
	fChain->SetBranchAddress("_trchi2",&_trchi2);
	fChain->SetBranchAddress("_idchi2",&_idchi2);
	fChain->SetBranchAddress("_chi2_fvtxmutr",&_chi2_fvtxmutr);
	fChain->SetBranchAddress("_ntrhits",&_ntrhits);
	fChain->SetBranchAddress("_nidhits",&_nidhits);
	fChain->SetBranchAddress("_charge",&_charge);
	fChain->SetBranchAddress("_lastgap",&_lastgap);
	fChain->SetBranchAddress("_bbcZ",&_bbcZ);
	fChain->SetBranchAddress("_fvtxZ",&_fvtxZ);

	fChain->SetBranchAddress("_scaled_trigbit",&_scaled_trigbit);
	fChain->SetBranchAddress("_pdtheta",&_pdtheta);
	fChain->SetBranchAddress("_rf_vtxchi2pdf",&_rf_vtxchi2pdf);
	fChain->SetBranchAddress("_rf_slope",&_rf_slope);
	//FVTX info
	fChain->SetBranchAddress("_chi2_fvtxmutr",&_chi2_fvtxmutr);
	fChain->SetBranchAddress("_nhits_fvtx",&_nhits_fvtx);
	fChain->SetBranchAddress("_nhits_fvtx_charge",&_nhits_fvtx_charge);
	fChain->SetBranchAddress("_hit_pattern",&_hit_pattern);
	fChain->SetBranchAddress("_chi2_fvtx",&_chi2_fvtx);
	fChain->SetBranchAddress("_px_fvtxmutr",&_px_fvtxmutr);
	fChain->SetBranchAddress("_py_fvtxmutr",&_py_fvtxmutr);
	fChain->SetBranchAddress("_pz_fvtxmutr",&_pz_fvtxmutr);
	fChain->SetBranchAddress("_dr_fvtx",&_dr_fvtx);
	fChain->SetBranchAddress("_dphi_fvtx",&_dphi_fvtx);
	fChain->SetBranchAddress("_dtheta_fvtx",&_dtheta_fvtx);
	fChain->SetBranchAddress("_mult_fvtxN",&_mult_fvtxN);
	fChain->SetBranchAddress("_mult_fvtxS",&_mult_fvtxS);

	fChain->SetBranchAddress("_projx_fvtxmutr",&_projx_fvtxmutr);
	fChain->SetBranchAddress("_projy_fvtxmutr",&_projy_fvtxmutr);

	fChain->SetBranchAddress("_xst1",&_xst1);
	fChain->SetBranchAddress("_yst1",&_yst1);
	fChain->SetBranchAddress("_xst2",&_xst2);
	fChain->SetBranchAddress("_yst2",&_yst2);
	fChain->SetBranchAddress("_xst3",&_xst3);
	fChain->SetBranchAddress("_yst3",&_yst3);

	int _prev_run = -1;
	bool bGOOD[narm] = {false, false};

	for (int ien=0; ien<nentries; ien++){

		if ( (ien%(nentries/100))==0 ){
			cout << setw(10) << int(ien*100.0/nentries + 0.1) << "\tpercent completed..." << endl;
		}

		fChain->GetEntry(ien);

		if ( _prev_run!=_runnumber ){
			bGOOD[0] = bGOOD[1] = true;
			if ( hbadrun[0]->GetBinContent(hbadrun[0]->FindBin(_runnumber))>0.5 ) bGOOD[0] = false;
			if ( hbadrun[1]->GetBinContent(hbadrun[1]->FindBin(_runnumber))>0.5 ) bGOOD[1] = false;

			if ( !bGOOD[0] ){
				cout << "ARM:0, BADRUN: " << _runnumber << endl;
			}
			if ( !bGOOD[1] ){
				cout << "ARM:1, BADRUN: " << _runnumber << endl;
			}

			_prev_run = _runnumber;
		}

		//Check good run
		if ( _arm==0 && !bGOOD[0] ) continue;
		if ( _arm==1 && !bGOOD[1] ) continue;

		if ( fabs(_bbcZ)>20.0 ) continue;
		if ( _arm==0 && _trchi2>15.0 ) continue;
		if ( _arm==1 && _trchi2>20.0 ) continue;
		if ( _lastgap>3.5 || _lastgap<1.5 ) continue;
		if ( _ntrhits<10.5 ) continue;
		if ( _DG0>cut_dg0[_arm] ) continue;
		if ( _DDG0>cut_ddg0[_arm] ) continue;
		if ( _rf_slope<0.10 ) continue;

		//MuTr acceptance cut
		float _phist1 = atan2(_yst1,_xst1);
		float _phist2 = atan2(_yst2,_xst2);
		float _phist3 = atan2(_yst3,_xst3);
		float _radst1 = sqrt(_xst1*_xst1 + _yst1*_yst1);
		float _radst2 = sqrt(_xst2*_xst2 + _yst2*_yst2);
		float _radst3 = sqrt(_xst3*_xst3 + _yst3*_yst3);

		if ( _arm==1 ){
			//if ( _phist2>0.35 && _phist2<0.75 ) continue;
			//if ( _phist2>0.77 && _phist2<1.20 && _radst2>100 ) continue;
			//if ( _phist3>+0.35 && _phist3<+0.77 ) continue;
			//if ( _phist3>+0.77 && _phist3<+1.15 && _radst3>160 && _radst3<275 ) continue;
			//if ( _phist3>-2.75 && _phist3<-1.92 && _radst3>237 && _radst3<275 ) continue;
			//if ( _phist2>-2.75 && _phist2<-1.92 && _radst2>140 && _radst2<163 ) continue;
			//if ( _phist1>-2.75 && _phist1<-2.00 && _radst1>76 && _radst1<80 ) continue;
			//if ( _phist1>+0.52 && _phist1<+1.11 && _radst1>60 && _radst1<75 ) continue;
			//if ( _phist2>-1.20 && _phist2<-0.90 ) continue;

			if ( _phist1>0.40 && _phist1<1.20 ) continue;
			if ( _phist2>0.40 && _phist2<1.20 ) continue;
			if ( _phist3>0.40 && _phist3<1.20 ) continue;
			if ( _phist1>-1.60 && _phist1<-0.88 ) continue;
			if ( _phist2>-1.60 && _phist2<-0.88 ) continue;
			if ( _phist3>-1.60 && _phist3<-0.88 ) continue;
			if ( _phist1>3.02 ) continue;
			if ( _phist2>3.02 ) continue;
			if ( _phist3>3.02 ) continue;
			if ( _phist1<-2.75 ) continue;
			if ( _phist2<-2.75 ) continue;
			if ( _phist3<-2.75 ) continue;
		}else{
			//if ( _phist2>0.00 && _phist2<0.40 ) continue;
			//if ( _phist1>1.55 && _phist1<1.90 ) continue;
			//if ( _phist2>1.55 && _phist2<1.90 ) continue;
			//if ( _phist3>1.55 && _phist3<1.90 ) continue;
			//if ( _phist3>-0.60 && _phist3<-0.40 ) continue;

			if ( _phist1>0 && _phist1<0.4 ) continue;
			if ( _phist2>0 && _phist2<0.4 ) continue;
			if ( _phist3>0 && _phist3<0.4 ) continue;
			if ( _phist1>1.55 && _phist1<1.9 ) continue;
			if ( _phist2>1.55 && _phist2<1.9 ) continue;
			if ( _phist3>1.55 && _phist3<1.9 ) continue;
			if ( _phist1>-0.60 && _phist1<-0.40 ) continue;
			if ( _phist2>-0.60 && _phist2<-0.40 ) continue;
			if ( _phist3>-0.60 && _phist3<-0.40 ) continue;
		}


		/*
		if ( _phist2<0 ) _phist2 += 2*TMath::Pi();
		if ( _phist3<0 ) _phist3 += 2*TMath::Pi();

		int _halfoctst2 = int(_phist2/(TMath::Pi()/8.));
		int _halfoctst3 = int(_phist3/(TMath::Pi()/8.));

		if ( _arm==0 && _halfoctst2==0 ) continue;
		if ( _arm==1 && _halfoctst2==1 ) continue;
		if ( _arm==1 && _halfoctst2==2 ) continue;
		if ( _arm==1 && _halfoctst2==15 ) continue;
		if ( _arm==1 && _halfoctst2==7 && _radst3<200 ) continue;
		if ( _arm==1 && _halfoctst2==8 && _radst3<200 ) continue;
		*/


		//Trigger cut - Run15pAu200 | Run15pAl200
		if ( trigger=="MU" ){
			if ( _arm==1 && !(_scaled_trigbit&(0x01000000)) ) continue;
			if ( _arm==0 && !(_scaled_trigbit&(0x02000000)) ) continue;
		}else if ( trigger=="MB" ){
			if ( !(_scaled_trigbit&(0x00000001)) && !(_scaled_trigbit&(0x00000002)) && !(_scaled_trigbit&(0x00000010)) ){
				//cout << "NON MB Triggered events" << endl;
				continue;
			}
		}

		if ( fabs(_eta)>1.4 && fabs(_eta)<(2.2+_arm*0.2) && fabs(_pz)>(3.5+0.5*(_lastgap-2)) && _pT>1.25 ){
			hRPHI_MUTR[_arm]->Fill(_phist2, _radst2);
		}

		if ( _nhits_fvtx<2.5 ) continue;
		if ( isnan(_px_fvtxmutr) || isnan(_py_fvtxmutr) || isnan(_pz_fvtxmutr) ) continue;
		if ( fabs(_px_fvtxmutr)<1e-10 || fabs(_py_fvtxmutr)<1e-10 || fabs(_pz_fvtxmutr)<1e-10 ) continue; 

		if ( _hit_pattern<256 && _nhits_fvtx_charge<2.5 ) continue;
		if ( _hit_pattern>256 && _nhits_fvtx_charge<2.5 ) continue;
		if ( _chi2_fvtxmutr>10.0 || _chi2_fvtxmutr<1e-10 ) continue;

		float pvalue = TMath::Prob(_chi2_fvtx*(2*_nhits_fvtx-5), 2*_nhits_fvtx-5);
		if ( pvalue<0.02 ) continue;

		TVector3 vec(_px_fvtxmutr,_py_fvtxmutr,_pz_fvtxmutr);
		_pT = vec.Pt();
		_eta = vec.Eta();
		float ptot = vec.Mag();

		if ( _arm==0 && (fabs(_eta)<1.4 || fabs(_eta)>2.2) ) continue;
		if ( bwide_eta ){
			if ( _arm==1 && (fabs(_eta)<1.4 || fabs(_eta)>2.4) ) continue;
		}else{
			if ( _arm==1 && (fabs(_eta)<1.4 || fabs(_eta)>2.2) ) continue;
		}

		int etabin = int((fabs(_eta)-1.2)/0.2);
		float pzcut = pz_cut[_arm][etabin]; 
		if ( _lastgap==2 ) pzcut -= 0.5;
		if ( _lastgap<3.5 && fabs(_pz)<pzcut ) continue;

		if ( _pT<pT_array[0] || _pT>pT_array[nptbin] ) continue;

		int ptbin = -999;
		for (int ipt=0; ipt<nptbin; ipt++){
			if ( _pT>pT_array[ipt] && _pT<=pT_array[ipt+1] ){
				ptbin = ipt;
				break;
			}
		}

		if ( ptbin<0 ){
			cout << ptbin << "\t" << _pT << endl;
			continue;
		}

		//if ( _rf_vtxchi2pdf>cut_vtxchi2[ptbin] ) continue;
		//if ( _pdtheta>cut_pdtheta[ptbin] ) continue;

		//if ( _rf_vtxchi2pdf>(cut_vtxchi2[_arm]-2*ptbin) ) continue;
		//if ( _rf_vtxchi2pdf>(cut_vtxchi2[_arm]-1*ptbin) ) continue;
		//if ( _rf_vtxchi2pdf>(cut_vtxchi2[_arm]-0.2*ptbin) ) continue;

		/*
		int centbin = -999;
		if ( (_mult_fvtxN+_mult_fvtxS)>30 ) continue;

		if ( (_mult_fvtxN+_mult_fvtxS)<=2 ) centbin = 4;
		else if ( (_mult_fvtxN+_mult_fvtxS)<=4 ) centbin = 3;
		else if ( (_mult_fvtxN+_mult_fvtxS)<=6 ) centbin = 2;
		else if ( (_mult_fvtxN+_mult_fvtxS)<=8 ) centbin = 1;
		else centbin = 0;

		if ( centbin<0 ){
			cout << centbin << "\t" << _mult_fvtxN+_mult_fvtxS << endl;
			continue;
		}
		*/

		float cut_dphi_fvtx = 0.025 + 0.20/ptot;
		float cut_dtheta_fvtx = 0.012 + 0.10/ptot;

		float cut_dphi_fvtx_loose = 1.5*cut_dphi_fvtx;
		float cut_dtheta_fvtx_loose = 1.5*cut_dtheta_fvtx;

		float cut_dphi_fvtx_tight = 0.5*cut_dphi_fvtx;
		float cut_dtheta_fvtx_tight = 0.5*cut_dtheta_fvtx;

		float w_trigeff_pT = ftrig_pT[_arm][_lastgap-2]->Eval(_pT); 
		float w_trigeff_eta = htrig_eta[_arm][_lastgap-2]->GetBinContent(htrig_eta[_arm][_lastgap-2]->FindBin(fabs(_eta)));

		if ( fabs(_dphi_fvtx)>cut_dphi_fvtx_loose ) continue;
		if ( fabs(_dtheta_fvtx)>cut_dtheta_fvtx_loose ) continue;

		hPT_ETA_loose[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta));
		hPT_ETA_loose[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta));

		if ( fabs(_dphi_fvtx)>cut_dphi_fvtx ) continue;
		if ( fabs(_dtheta_fvtx)>cut_dtheta_fvtx ) continue;

		float w_trigeff_all = w_trigeff_eta * w_trigeff_pT / ftrig_pT[_arm][_lastgap-2]->Eval(2.5);

		hPT_ETA[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta));
		hPT_ETA_cor_pT[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta), 1./w_trigeff_pT);
		hPT_ETA_cor_eta[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta), 1./w_trigeff_eta);
		hPT_ETA_cor_all[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta), 1./w_trigeff_all);
		hPT_ETA_fine[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta));
		hPT_ETA_fine_cor_pT[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta),1./w_trigeff_pT);
		hPT_ETA_fine_cor_eta[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta),1./w_trigeff_eta);
		hPT_ETA_fine_cor_all[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta),1./w_trigeff_all);
		hPT_mean[_arm][_lastgap-2][_charge]->Fill(_pT, _pT);

		hPT_ETA[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta));
		hPT_ETA_cor_pT[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta), 1./w_trigeff_pT);
		hPT_ETA_cor_eta[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta), 1./w_trigeff_eta);
		hPT_ETA_cor_all[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta), 1./w_trigeff_all);
		hPT_ETA_fine[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta));
		hPT_ETA_fine_cor_pT[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta),1./w_trigeff_pT);
		hPT_ETA_fine_cor_eta[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta),1./w_trigeff_eta);
		hPT_ETA_fine_cor_all[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta),1./w_trigeff_all);
		hPT_mean[_arm][_lastgap-2][2]->Fill(_pT, _pT);

		if ( _pT>1.5 ){
			hRUN[_arm][_lastgap-2]->Fill(_runnumber);
		}

		if ( _pT>1.25 ){
			hRPHI_FVTX[_arm]->Fill(atan2(_projy_fvtxmutr,_projx_fvtxmutr), sqrt(_projx_fvtxmutr*_projx_fvtxmutr + _projy_fvtxmutr*_projy_fvtxmutr));
		}


		if ( fabs(_dphi_fvtx)>cut_dphi_fvtx_tight ) continue;
		if ( fabs(_dtheta_fvtx)>cut_dtheta_fvtx_tight ) continue;

		hPT_ETA_tight[_arm][_lastgap-2][_charge]->Fill(_pT, fabs(_eta));
		hPT_ETA_tight[_arm][_lastgap-2][2]->Fill(_pT, fabs(_eta));

		/*
		//FVTX Cut
		if ( fabs(_fvtxZ)>15.0 ) continue;
		if ( fabs(_chi2_fvtxmutr)>8.0 ) continue;

		hPT_ETA_FVTX[_arm][_lastgap-2][_charge][centbin]->Fill(_pT, fabs(_eta));
		hPT_mean_FVTX[_arm][_lastgap-2][_charge][centbin]->Fill(_pT, _pT);

		hPT_ETA_FVTX[_arm][_lastgap-2][2][centbin]->Fill(_pT, fabs(_eta));
		hPT_mean_FVTX[_arm][_lastgap-2][2][centbin]->Fill(_pT, _pT);
		*/

	}//ien

	TFile *outfile;

	if ( bwide_eta ){
		outfile = new TFile(Form("%s%s_FVTX_hadron_histo_eta1224.root",dataset.c_str(),trigger.c_str()),"RECREATE");
	}else{
		outfile = new TFile(Form("%s%s_FVTX_hadron_histo_eta1222.root",dataset.c_str(),trigger.c_str()),"RECREATE");
	}
	for (int iarm=0; iarm<narm; iarm++){
		for (int igap=0; igap<ngap; igap++){
			for (int ichg=0; ichg<nchg; ichg++){
				hPT_ETA[iarm][igap][ichg]->Write();
				hPT_ETA_cor_pT[iarm][igap][ichg]->Write();
				hPT_ETA_cor_eta[iarm][igap][ichg]->Write();
				hPT_ETA_cor_all[iarm][igap][ichg]->Write();
				hPT_ETA_fine[iarm][igap][ichg]->Write();
				hPT_ETA_fine_cor_pT[iarm][igap][ichg]->Write();
				hPT_ETA_fine_cor_eta[iarm][igap][ichg]->Write();
				hPT_ETA_fine_cor_all[iarm][igap][ichg]->Write();
				hPT_mean[iarm][igap][ichg]->Write();
				hPT_ETA_tight[iarm][igap][ichg]->Write();
				hPT_ETA_loose[iarm][igap][ichg]->Write();
			}
			hRUN[iarm][igap]->Write();
		}
		hRPHI_MUTR[iarm]->Write();
		hRPHI_FVTX[iarm]->Write();
	}
	outfile->Close();

}
