#include "/phenix/u/shlim/Style.h"

void Draw_spectra_pAu(const bool bMultCut = false){

	const int narm = 2;
	const int nchg = 3;
	const int ngap = 2;
	const int netabin = 6;
	const int ncent = 6;
	const int nskip = 2;
	//const int ncent = 1;
	const float deta[narm] = {0.8, 1.0};

	const bool bSAVE = false;
	const bool bWRITE = false;

	int cent_array[ncent+1] = {0, 5, 10, 20, 40, 60, 84}; 
	float Ncoll[ncent] = {9.693, 8.379, 7.394, 6.085, 4.409, 2.637};
	float Npart[ncent] = {10.692, 9.399, 8.394, 7.085, 5.409, 3.637};
	float BiasF[ncent] = {0.858, 0.902, 0.935, 0.982, 1.026, 1.000};
	float Ncoll_err[ncent] = {0.0578, 0.0655, 0.0652, 0.0662, 0.0658, 0.0592};
	float BiasF_err[ncent] = {0.0093, 0.0067, 0.0043, 0.0051, 0.0127, 0.0675};
	const int nColor[ncent] = {1, 4, kGreen+2, 2, 6, 9};

	float BiasF_mb = 0.858;
	float BiasF_mb_err = 0.016;
	float Ncoll_mb = 4.7;
	float Ncoll_mb_err = 0.064;

	const float sys_frac_kpi = 0.07;
	const float sys_frac_shape[2] = {0.05, 0.03};
	const float sys_frac_muon[2] = {0.072, 0.027};
	const float sys_frac_fvtx[2] = {0.027, 0.023};
	const float sys_frac_proton = 0.05;
	
	float Nevt[ncent][narm] = {0.0};
	float Nevt_mb[narm] = {0.0};

	TH1F *hCentrality[narm];
	TFile *infile0 = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/10.runQA/00.eventhist/Run15pAu200_merged_eventhist_hadron_20180327.root","read");
	hCentrality[0] = (TH1F*)infile0->Get("hCentrality_S");
	hCentrality[1] = (TH1F*)infile0->Get("hCentrality_N");

	for (int icent=0; icent<ncent; icent++){
		cout << "CENTBIN: " << icent;
		for (int iarm=0; iarm<narm; iarm++){
			Nevt[icent][iarm] = hCentrality[iarm]->Integral(hCentrality[iarm]->FindBin(cent_array[icent]+1+0.001),hCentrality[iarm]->FindBin(cent_array[icent+1]+0.001));
			cout << ", ARM: " << iarm << ", N EVT: " << Nevt[icent][iarm];

			Nevt_mb[iarm] += Nevt[icent][iarm];
		}
		cout << endl;
	}
	//return;

	//TFile *infile_eff = new TFile("embedding_study_for_pp/Run15pp200_hadron_eff_output_set1_embed_eta1224.root","read");
	//TFile *infile_eff = new TFile("embedding_study_for_pp/Run15pAu200_CENT6084_hadron_eff_output_set1_embed_eta1224.root","read");
	TFile *infile_eff = new TFile("embedding_study_for_pp/Run15pAu200_CENT6084_hadron_eff_output_set2_embed_eta1224.root","read");
	TH1F *heff_pT[narm];
	TH1F *heff_etabin_pT[narm][6];
	TF1 *feff_pT[narm];
	TH1F *heff_eta[narm];
	for (int iarm=0; iarm<narm; iarm++){
		heff_pT[iarm] = (TH1F*)infile_eff->Get(Form("heff_pT_set0_arm%d",iarm));
		feff_pT[iarm] = (TF1*)infile_eff->Get(Form("feff_pT_set0_arm%d",iarm));
		heff_eta[iarm] = (TH1F*)infile_eff->Get(Form("heff_eta_set0_arm%d",iarm));

		for (int ieta=1; ieta<6-1+iarm; ieta++){
			heff_etabin_pT[iarm][ieta] = (TH1F*)infile_eff->Get(Form("heff_etabin_pT_set0_arm%d_%d",iarm,ieta));
			if ( !heff_etabin_pT[iarm][ieta] ){
				cout << "CAN NOT LOAD ACCxEFF HISTOGRAM ETABIN" << endl;
				exit(1);
			}
		}

		if ( !heff_pT[iarm] || !heff_eta[iarm] || !feff_pT[iarm] ){
			cout << "CAN NOT LOAD ACCxEFF HISTOGRAM" << endl;
			exit(1);
		}
	}
	//return;

	//Centrality dep.
	TFile *infile_eff2 = new TFile("Run15pAu200_hadron_eff_func.root","READ");
	TFile *infile_mis = new TFile("Run15pAu200_hadron_mismatch_frac.root","READ");
	TF1 *feff_pT_cent[ncent][narm];
	TH1F *heff_eta_cent[ncent][narm];
	TH1F *hmismatch_pT_cent[ncent][narm];
	for (int icent=0; icent<ncent; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			feff_pT_cent[icent][iarm] = (TF1*)infile_eff2->Get(Form("ffit_cent%d_arm%d",icent,iarm));
			heff_eta_cent[icent][iarm] = (TH1F*)infile_eff2->Get(Form("heff_rap_cent%d_arm%d",icent,iarm));
			hmismatch_pT_cent[icent][iarm] = (TH1F*)infile_mis->Get(Form("hmismatch_arm%d_centbin%d",iarm,icent));
		}
	}

	TFile *infile = new TFile("Run15pAu200MU_FVTX_hadron_histo.root","read");



	TH1F *hPT[ncent][narm][ngap][nchg];
	TH1F *hETA[ncent][narm][ngap][nchg];
	TH1F *hPT_etabin[ncent][narm][ngap][nchg][6];
	TH2F *hPT_ETA_cor_pT[ncent][narm][ngap][nchg];
	TH2F *hPT_ETA_cor_eta[ncent][narm][ngap][nchg];

	for (int icent=0; icent<ncent; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			for (int igap=0; igap<ngap; igap++){
				for (int ichg=0; ichg<nchg; ichg++){

					hPT_ETA_cor_pT[icent][iarm][igap][ichg] = (TH2F*)infile->Get(Form("hPT_ETA_cor_pT_arm%d_gap%d_chg%d_cent%02d%02d",iarm,igap+2,ichg,cent_array[icent],cent_array[icent+1]));
					hPT_ETA_cor_eta[icent][iarm][igap][ichg] = (TH2F*)infile->Get(Form("hPT_ETA_cor_eta_arm%d_gap%d_chg%d_cent%02d%02d",iarm,igap+2,ichg,cent_array[icent],cent_array[icent+1]));

					if ( !hPT_ETA_cor_pT[icent][iarm][igap][ichg] || !hPT_ETA_cor_eta[icent][iarm][igap][ichg] ){
						cout << "CAN NOT LOAD HISTOGRAMS!" << endl;
						exit(1);
					}//

					hPT[icent][iarm][igap][ichg] = (TH1F*)hPT_ETA_cor_pT[icent][iarm][igap][ichg]->ProjectionX(Form("hPT_cent%d_arm%d_gap%d_chg%d",icent,iarm,igap,ichg)); //1.4-2.2(2.4)

					int ptbin_min = hPT[icent][iarm][igap][ichg]->FindBin(2.5+0.01);
					int ptbin_max = hPT[icent][iarm][igap][ichg]->GetNbinsX();

					hETA[icent][iarm][igap][ichg] = (TH1F*)hPT_ETA_cor_eta[icent][iarm][igap][ichg]->ProjectionY(Form("hETA_cent%d_arm%d_gap%d_chg%d",icent,iarm,igap,ichg),ptbin_min,ptbin_max);
					hETA[icent][iarm][igap][ichg]->Rebin();

					for (int ieta=0; ieta<6; ieta++){
						hPT_etabin[icent][iarm][igap][ichg][ieta] = (TH1F*)hPT_ETA_cor_pT[icent][iarm][igap][ichg]->ProjectionX(Form("hPT_etabin_cent%d_arm%d_gap%d_chg%d_%d",icent,iarm,igap,ichg,ieta),2*ieta+1,2*ieta+2);
					}

				}//ichg
			}//igap
		}//iarm
	}//icent

	//return;

	//Get p+p
	TH1F *hyield_pp[narm];
	TH1F *hyield_sys_pp[narm];
	TH1F *hyield_pp_eta;
	TH1F *hyield_sys_pp_eta;
	TH1F *hyield_pp_eta_int;
	TH1F *hyield_sys_pp_eta_int;
	TFile *infile_pp = new TFile("Run15pp200_hadron_spectra_eta1224.root","read");
	hyield_pp[0] = (TH1F*)infile_pp->Get("hyield_pp_arm0");
	hyield_pp[1] = (TH1F*)infile_pp->Get("hyield_pp_arm1");
	hyield_sys_pp[0] = (TH1F*)infile_pp->Get("hyield_sys_pp_arm0");
	hyield_sys_pp[1] = (TH1F*)infile_pp->Get("hyield_sys_pp_arm1");
	hyield_pp_eta = (TH1F*)infile_pp->Get("hyield_pp_eta");
	hyield_sys_pp_eta = (TH1F*)infile_pp->Get("hyield_sys_pp_eta");
	hyield_pp_eta_int = (TH1F*)infile_pp->Get("hyield_pp_eta_int");
	hyield_sys_pp_eta_int = (TH1F*)infile_pp->Get("hyield_sys_pp_eta_int");

	int nptbin = 0;
	float pT_mean[20], pT_err[20], pT_syserr[20], pT_varbin[20];
	float yield[ncent][narm][20], yield_err[ncent][narm][20], yield_syserr[ncent][narm][20];
	float yield_etabin[ncent][narm][6][20], yield_etabin_combined[ncent][narm][20] = {0.};
	float RpA[ncent][narm][20], RpA_err[ncent][narm][20], RpA_syserr[ncent][narm][20];
	float yield_mb[narm][20] = {0.}, yield_mb_err[narm][20] = {0.}, yield_mb_syserr[narm][20] = {0.};
	float RpA_mb[narm][20] = {0.}, RpA_mb_err[narm][20] = {0.}, RpA_mb_syserr[narm][20] = {0.};

	TGraphErrors *gyield[ncent][narm];
	TGraphErrors *gyield_sys[ncent][narm];
	TGraphErrors *gyield_mb[narm];
	TGraphErrors *gyield_mb_sys[narm];
	TGraphErrors *gRpA[ncent][narm];
	TGraphErrors *gRpA_sys[ncent][narm];
	TGraphErrors *gRpA_mb[narm];
	TGraphErrors *gRpA_mb_sys[narm];

	for (int icent=0; icent<ncent; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			nptbin = hPT[icent][iarm][1][2]->GetNbinsX();

			for (int ipt=0; ipt<nptbin; ipt++){

				pT_mean[ipt] = hPT[icent][iarm][1][2]->GetBinCenter(ipt+1);
				pT_err[ipt] = 0.5*hPT[icent][iarm][1][2]->GetBinWidth(ipt+1);
				pT_varbin[ipt] = hPT[icent][iarm][1][2]->GetBinLowEdge(ipt+1);
				pT_syserr[ipt] = 0.125;

				yield[icent][iarm][ipt] = (hPT[icent][iarm][0][2]->GetBinContent(ipt+1) + hPT[icent][iarm][1][2]->GetBinContent(ipt+1));

				float err1 = hPT[icent][iarm][0][2]->GetBinContent(ipt+1);
				float err2 = hPT[icent][iarm][1][2]->GetBinContent(ipt+1);
				yield_err[icent][iarm][ipt] = sqrt(err1+err2);

				float eff_pT = 0.0;

				if ( pT_mean[ipt]<5.0 ){
					eff_pT = heff_pT[iarm]->GetBinContent(heff_pT[iarm]->FindBin(pT_mean[ipt]));
				}else{
					//eff_pT = feff_pT[iarm]->Eval(pT_mean[ipt]);
					eff_pT = heff_pT[iarm]->GetBinContent(heff_pT[iarm]->FindBin(pT_mean[ipt]));
				}

				eff_pT *= (feff_pT_cent[icent][iarm]->Eval(pT_mean[ipt])/feff_pT_cent[ncent-1][iarm]->Eval(pT_mean[ipt]));
				eff_pT *= (1-hmismatch_pT_cent[icent][iarm]->GetBinContent(hmismatch_pT_cent[icent][iarm]->FindBin(pT_mean[ipt])));
				//cout << "EFF SCALE CENT: " << icent << ", ARM: " << iarm << ", PT: " << ipt << ", A: " << (feff_pT_cent[icent][iarm]->Eval(pT_mean[ipt])/feff_pT_cent[ncent-1][iarm]->Eval(pT_mean[ipt])) << endl;
				//
				yield_mb[iarm][ipt] += yield[icent][iarm][ipt]/(2*pT_err[ipt]*deta[iarm]*eff_pT);
				yield_mb_err[iarm][ipt] += yield[icent][iarm][ipt];

				yield[icent][iarm][ipt] /= (2*pT_err[ipt]*deta[iarm]*eff_pT*Nevt[icent][iarm]);
				yield_err[icent][iarm][ipt] /= (2*pT_err[ipt]*deta[iarm]*eff_pT*Nevt[icent][iarm]);

				yield[icent][iarm][ipt] *= BiasF[icent];
				yield_err[icent][iarm][ipt] *= BiasF[icent];

				float syserr = sys_frac_kpi*sys_frac_kpi;
				syserr += sys_frac_shape[iarm]*sys_frac_shape[iarm];
				syserr += sys_frac_proton*sys_frac_proton;
				syserr += sys_frac_muon[iarm]*sys_frac_muon[iarm];
				syserr += sys_frac_fvtx[iarm]*sys_frac_fvtx[iarm];
			
				yield_syserr[icent][iarm][ipt] = yield[icent][iarm][ipt] * sqrt(syserr);

				yield_etabin_combined[icent][iarm][ipt] = 0.0;

				for (int ieta=1; ieta<6-1+iarm; ieta++){
					if ( ipt>1 ){
						yield_etabin[icent][iarm][ieta][ipt] = (hPT_etabin[icent][iarm][0][2][ieta]->GetBinContent(ipt+1) + hPT_etabin[icent][iarm][1][2][ieta]->GetBinContent(ipt+1));
						eff_pT = heff_etabin_pT[iarm][ieta]->GetBinContent(heff_etabin_pT[iarm][ieta]->FindBin(pT_mean[ipt]));
						eff_pT *= (feff_pT_cent[icent][iarm]->Eval(pT_mean[ipt])/feff_pT_cent[ncent-1][iarm]->Eval(pT_mean[ipt]));
						eff_pT *= (1-hmismatch_pT_cent[icent][iarm]->GetBinContent(hmismatch_pT_cent[icent][iarm]->FindBin(pT_mean[ipt])));

						yield_etabin[icent][iarm][ieta][ipt] /= (2*pT_err[ipt]*eff_pT*Nevt[icent][iarm]);
						yield_etabin[icent][iarm][ieta][ipt] *= BiasF[icent];
					}else{
						yield_etabin[icent][iarm][ieta][ipt] = 0.0;
					}

					yield_etabin_combined[icent][iarm][ipt] += yield_etabin[icent][iarm][ieta][ipt];
				}

				yield_etabin_combined[icent][iarm][ipt] /= deta[iarm];

				if ( pT_mean[ipt]<3.0 ){
					err1 = yield_err[icent][iarm][ipt] / yield[icent][iarm][ipt]; 
					err2 = yield_syserr[icent][iarm][ipt] / yield[icent][iarm][ipt]; 

					yield[icent][iarm][ipt] = yield_etabin_combined[icent][iarm][ipt];
					yield_err[icent][iarm][ipt] = yield[icent][iarm][ipt] * err1;
					yield_syserr[icent][iarm][ipt] = yield[icent][iarm][ipt] * err2;
				}

				//cout << "CENT: " << icent << ", ARM: " << iarm << ", PT: " << ipt << ", R: " << yield_etabin_combined[icent][iarm][ipt]/yield[icent][iarm][ipt] << endl;

				RpA[icent][iarm][ipt] = yield[icent][iarm][ipt] / hyield_pp[iarm]->GetBinContent(ipt+1) / Ncoll[icent];

				err1 = hyield_pp[iarm]->GetBinError(ipt+1)/hyield_pp[iarm]->GetBinContent(ipt+1);
				err2 = yield_err[icent][iarm][ipt] / yield[icent][iarm][ipt];
				RpA_err[icent][iarm][ipt] = RpA[icent][iarm][ipt] * sqrt(err1*err1 + err2*err2);

				err1 = hyield_sys_pp[iarm]->GetBinError(ipt+1)/hyield_sys_pp[iarm]->GetBinContent(ipt+1);
				err2 = yield_syserr[icent][iarm][ipt] / yield[icent][iarm][ipt];
				RpA_syserr[icent][iarm][ipt] = RpA[icent][iarm][ipt] * sqrt(err1*err1 + err2*err2);

			}//ipt

			gyield[icent][iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &yield[icent][iarm][nskip], 0, &yield_err[icent][iarm][nskip]);
			gyield[icent][iarm]->SetMarkerStyle(20);
			gyield[icent][iarm]->SetLineWidth(2);
			gyield[icent][iarm]->SetMarkerColor(nColor[icent]);
			gyield[icent][iarm]->SetLineColor(nColor[icent]);

			gyield_sys[icent][iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &yield[icent][iarm][nskip], &pT_syserr[nskip], &yield_syserr[icent][iarm][nskip]);
			gyield_sys[icent][iarm]->SetMarkerStyle(20);
			gyield_sys[icent][iarm]->SetLineWidth(2);
			gyield_sys[icent][iarm]->SetMarkerColor(nColor[icent]);
			gyield_sys[icent][iarm]->SetLineColor(nColor[icent]);
			gyield_sys[icent][iarm]->SetFillStyle(0);

			gRpA[icent][iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &RpA[icent][iarm][nskip], 0, &RpA_err[icent][iarm][nskip]);
			gRpA[icent][iarm]->SetMarkerStyle(20+4*iarm);
			gRpA[icent][iarm]->SetLineWidth(1);
			gRpA[icent][iarm]->SetMarkerColor(nColor[icent]);
			gRpA[icent][iarm]->SetLineColor(nColor[icent]);

			gRpA_sys[icent][iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &RpA[icent][iarm][nskip], &pT_syserr[nskip], &RpA_syserr[icent][iarm][nskip]);
			gRpA_sys[icent][iarm]->SetMarkerStyle(20+4*iarm);
			gRpA_sys[icent][iarm]->SetLineWidth(1);
			gRpA_sys[icent][iarm]->SetMarkerColor(nColor[icent]);
			gRpA_sys[icent][iarm]->SetLineColor(nColor[icent]);
			gRpA_sys[icent][iarm]->SetFillStyle(0);
		}//iarm
	}

	//return;

	for (int iarm=0; iarm<narm; iarm++){
		for (int ipt=0; ipt<nptbin; ipt++){
			float err = 1./sqrt(yield_mb_err[iarm][ipt]);
			yield_mb[iarm][ipt] /= Nevt_mb[iarm];
			yield_mb_err[iarm][ipt] = yield_mb[iarm][ipt] * err; 

			float syserr = sys_frac_kpi*sys_frac_kpi;
			syserr += sys_frac_shape[iarm]*sys_frac_shape[iarm];
			syserr += sys_frac_proton*sys_frac_proton;
			syserr += sys_frac_muon[iarm]*sys_frac_muon[iarm];
			syserr += sys_frac_fvtx[iarm]*sys_frac_fvtx[iarm];

			yield_mb_syserr[iarm][ipt] = yield_mb[iarm][ipt] * sqrt(syserr); 

			yield_mb[iarm][ipt] *= BiasF_mb;
			yield_mb_err[iarm][ipt] *= BiasF_mb;
			yield_mb_syserr[iarm][ipt] *= BiasF_mb;

			RpA_mb[iarm][ipt] = yield_mb[iarm][ipt] / hyield_pp[iarm]->GetBinContent(ipt+1) / Ncoll_mb;

			float err1 = hyield_pp[iarm]->GetBinError(ipt+1)/hyield_pp[iarm]->GetBinContent(ipt+1);
			float err2 = yield_mb_err[iarm][ipt] / yield_mb[iarm][ipt];
			RpA_mb_err[iarm][ipt] = RpA_mb[iarm][ipt] * sqrt(err1*err1 + err2*err2);

			err1 = hyield_sys_pp[iarm]->GetBinError(ipt+1)/hyield_sys_pp[iarm]->GetBinContent(ipt+1);
			err2 = yield_mb_syserr[iarm][ipt] / yield_mb[iarm][ipt];
			RpA_mb_syserr[iarm][ipt] = RpA_mb[iarm][ipt] * sqrt(err1*err1 + err2*err2);
		}

		gyield_mb[iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &yield_mb[iarm][nskip], 0, &yield_mb_err[iarm][nskip]);
		gyield_mb[iarm]->SetMarkerStyle(20);
		gyield_mb[iarm]->SetLineWidth(2);
		gyield_mb[iarm]->SetMarkerColor(1);
		gyield_mb[iarm]->SetLineColor(1);

		gyield_mb_sys[iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &yield_mb[iarm][nskip], &pT_syserr[nskip], &yield_mb_syserr[iarm][nskip]);
		gyield_mb_sys[iarm]->SetMarkerStyle(20);
		gyield_mb_sys[iarm]->SetLineWidth(2);
		gyield_mb_sys[iarm]->SetMarkerColor(1);
		gyield_mb_sys[iarm]->SetLineColor(1);
		gyield_mb_sys[iarm]->SetFillStyle(0);

		gRpA_mb[iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &RpA_mb[iarm][nskip], 0, &RpA_mb_err[iarm][nskip]);
		gRpA_mb[iarm]->SetMarkerStyle(20+iarm);
		gRpA_mb[iarm]->SetLineWidth(2);
		gRpA_mb[iarm]->SetMarkerColor(1);
		gRpA_mb[iarm]->SetLineColor(1);

		gRpA_mb_sys[iarm] = new TGraphErrors(nptbin-nskip, &pT_mean[nskip], &RpA_mb[iarm][nskip], &pT_syserr[nskip], &RpA_mb_syserr[iarm][nskip]);
		gRpA_mb_sys[iarm]->SetMarkerStyle(20+iarm);
		gRpA_mb_sys[iarm]->SetLineWidth(2);
		gRpA_mb_sys[iarm]->SetMarkerColor(1);
		gRpA_mb_sys[iarm]->SetLineColor(1);
		gRpA_mb_sys[iarm]->SetFillStyle(0);
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*800,400);
	c1->Divide(2,1);

	for (int iarm=0; iarm<narm; iarm++){
		c1->cd(iarm+1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1F*)gPad->DrawFrame(1.0,2e-8,10,2);
		SetHistoStyle("p_{T} (GeV/c)","d^{2}N / d#etadp_{T}");
		for (int icent=0; icent<ncent; icent++){
			gyield_sys[icent][iarm]->Draw("2");
			gyield[icent][iarm]->Draw("p");
		}
		gyield_mb_sys[iarm]->Draw("2");
		gyield_mb[iarm]->Draw("p");
	}

	//return;
	//
	TLine *line = new TLine(0.5,1,10.5,1);
	line->SetLineStyle(2);
	line->SetLineWidth(2);

	TLine *line_eta = new TLine(-3.2,1,3.2,1);
	line_eta->SetLineStyle(2);
	line_eta->SetLineWidth(2);

	//return;

	TCanvas *c2 = new TCanvas("c2","c2",1.1*800,400);
	c2->Divide(2,1);

	for (int iarm=0; iarm<narm; iarm++){
		c2->cd(iarm+1);
		SetPadStyle();

		htmp = (TH1F*)gPad->DrawFrame(1,0,10,3);
		SetHistoStyle("p_{T} (GeV/c)","R_{pA}");
		if ( iarm==1 ){
			TLegend *leg = new TLegend(0.25,0.65,0.9,0.9);
			leg->SetHeader("p+Au #sqrt{s_{NN}}=200 GeV, h^{#pm}, 1.2<#eta<2.4");
			leg->SetFillStyle(0);
			leg->SetNColumns(2);
		}
		for (int icent=0; icent<ncent; icent++){
			gRpA[icent][iarm]->Draw("p");
			if ( iarm==1 )
				leg->AddEntry(gRpA[icent][iarm],Form("%d-%d%c",cent_array[icent],cent_array[icent+1],'%'),"P");
		}
		if ( iarm==1 )
			leg->Draw();
	}

	TCanvas *c2_ = new TCanvas("c2_","c2_",400*3,400*2);
	TPad *pp2_[ncent];

	pp2_[0] = new TPad(Form("c0_0"),Form("c0_0"),0.08,0.52,0.38,1);
	pp2_[0]->Draw();
	pp2_[0]->cd();
	gPad->SetBottomMargin(0.000);
	gPad->SetRightMargin(0.000);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.001);
	htmp = (TH1F*)gPad->DrawFrame(0.5,0.005,10.5,3.5);
	SetHistoStyle("","R_{CP}");
	htmp->GetYaxis()->SetLabelSize(0.06);
	htmp->GetYaxis()->SetTitleSize(0.06);
	htmp->GetYaxis()->SetTitleOffset(1.15);

	line->Draw();

	c2_->cd();

	pp2_[1] = new TPad(Form("c0_1"),Form("c0_1"),0.3795,0.52,0.685,1);
	pp2_[1]->Draw();
	pp2_[1]->cd();
	gPad->SetBottomMargin(0.000);
	//gPad->SetRightMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.000);
	htmp = (TH1F*)gPad->DrawFrame(0.5,0.005,10.5,3.5);
	SetHistoStyle("","");

	line->Draw();

	c2_->cd();

	pp2_[2] = new TPad(Form("c0_2"),Form("c0_2"),0.6845,0.52,0.98,1);
	pp2_[2]->Draw();
	pp2_[2]->cd();
	gPad->SetBottomMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.001);
	htmp = (TH1F*)gPad->DrawFrame(0.5,0.005,10.5,3.5);
	SetHistoStyle("","");

	line->Draw();

	c2_->cd();

	pp2_[3] = new TPad(Form("c0_1"),Form("c0_1"),0.08,0.0,0.38,0.52);
	pp2_[3]->Draw();
	pp2_[3]->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.0);
	gPad->SetTopMargin(0.0);
	gPad->SetLeftMargin(0.00);
	htmp = (TH1F*)gPad->DrawFrame(0.5,0,10.5,3.495);
	SetHistoStyle("p_{T} (GeV/c)","R_{pA}");
	htmp->GetYaxis()->SetLabelSize(0.06);
	htmp->GetYaxis()->SetTitleSize(0.06);
	htmp->GetYaxis()->SetTitleOffset(1.15);
	htmp->GetXaxis()->SetLabelSize(0.07);
	htmp->GetXaxis()->SetTitleSize(0.07);

	line->Draw();

	c2_->cd();
	pp2_[4] = new TPad(Form("c0_1"),Form("c0_1"),0.3795,0.0,0.685,0.52);
	pp2_[4]->Draw();
	pp2_[4]->cd();
	gPad->SetBottomMargin(0.15);
	//gPad->SetRightMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.0);
	gPad->SetLeftMargin(0.000);
	htmp = (TH1F*)gPad->DrawFrame(0.5,0,10.5,3.495);
	SetHistoStyle("p_{T} (GeV/c)","");
	htmp->GetXaxis()->SetLabelSize(0.07);
	htmp->GetXaxis()->SetTitleSize(0.07);

	line->Draw();

	c2_->cd();
	pp2_[5] = new TPad(Form("c0_1"),Form("c0_1"),0.6845,0.0,0.98,0.52);
	pp2_[5]->Draw();
	pp2_[5]->cd();
	gPad->SetBottomMargin(0.15);
	//gPad->SetRightMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.0);
	gPad->SetLeftMargin(0.000);
	htmp = (TH1F*)gPad->DrawFrame(0.5,0,10.5,3.495);
	SetHistoStyle("p_{T} (GeV/c)","");
	htmp->GetXaxis()->SetLabelSize(0.07);
	htmp->GetXaxis()->SetTitleSize(0.07);

	line->Draw();

	c2_->cd();
	TGaxis *axis0 = new TGaxis(0.08,0.52,0.08,0.977,0.005,3.5,510);
	axis0->SetLabelSize(0.03);
	axis0->SetTitle("R_{pA}");
	axis0->SetTitleOffset(0.8);
	axis0->Draw();

	TGaxis *axis1 = new TGaxis(0.08,0.08,0.08,0.52,0.000,3.495,510);
	axis1->SetLabelSize(0.03);
	axis1->SetTitle("R_{pA}");
	axis1->SetTitleOffset(0.8);
	axis1->Draw();

	for (int icent=0; icent<ncent; icent++){
		pp2_[icent]->cd();

		float global_sys = sqrt(Ncoll_err[icent]*Ncoll_err[icent] + BiasF_err[icent]*BiasF_err[icent] + 0.101*0.101);
		TBox *bsys = new TBox(0.5,1-global_sys,0.8,1+global_sys);
		bsys->SetFillColor(1);
		bsys->SetLineColor(1);
		bsys->Draw();

		TLatex *tex = new TLatex(1.4,3.25,"p+Au#rightarrowh^{#pm}+X #sqrt{s_{NN}}=200 GeV");
		tex->SetTextSize(0.055);
		tex->Draw();

		TLatex *tex1 = new TLatex(1.4,3.00,Form("%d-%d%c centrality",cent_array[icent],cent_array[icent+1],'%'));
		tex1->SetTextSize(0.055);
		tex1->Draw();

		gRpA_sys[icent][0]->Draw("2");
		gRpA_sys[icent][1]->Draw("2");
		gRpA[icent][0]->Draw("p");
		gRpA[icent][1]->Draw("p");

		if ( icent<3 )
			TLegend *leg = new TLegend(0.30,0.65,0.90,0.78);
		else
			TLegend *leg = new TLegend(0.30,0.70,0.90,0.83);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		le = leg->AddEntry(gRpA[icent][0],"-2.2<#eta<-1.2 (Au-going)","P");
		le->SetTextSize(0.055);
		le = leg->AddEntry(gRpA[icent][1],"1.2<#eta<2.4 (p-going)","P");
		le->SetTextSize(0.055);
		leg->Draw();
	}

	//return;

	//int netabin = 0;
	float eta_mean[narm][20], eta_err[narm][20], eta_syserr[narm][20];
	float yield_eta[ncent][narm][20], yield_eta_err[ncent][narm][20], yield_eta_syserr[ncent][narm][20];
	float yield_mb_eta[narm][20] = {0.}, yield_mb_eta_err[narm][20] = {0.}, yield_mb_eta_syserr[narm][20] = {0.};
	float yield_int_eta[ncent][narm] = {0.}, yield_int_eta_err[ncent][narm] = {0.}, yield_int_eta_syserr[ncent][narm] = {0.};
	float eta_int_mean[narm] = {-1.8, 1.9}, eta_int_err[narm] = {0.1, 0.1};

	float RpA_eta[ncent][narm][20], RpA_eta_err[ncent][narm][20], RpA_eta_syserr[ncent][narm][20];
	float RpA_mb_eta[narm][20] = {0.}, RpA_mb_eta_err[narm][20] = {0.}, RpA_mb_eta_syserr[narm][20] = {0.};
	float RpA_int_eta[ncent][narm], RpA_int_eta_err[ncent][narm], RpA_int_eta_syserr[ncent][narm];

	TGraphErrors *gyield_eta[ncent][narm];
	TGraphErrors *gyield_mb_eta[narm];
	TGraphErrors *gRpA_eta[ncent][narm];
	TGraphErrors *gRpA_mb_eta[narm];

	TGraphErrors *gyield_sys_eta[ncent][narm];
	TGraphErrors *gyield_mb_sys_eta[narm];
	TGraphErrors *gRpA_sys_eta[ncent][narm];
	TGraphErrors *gRpA_mb_sys_eta[narm];

	TGraphErrors *gyield_eta_int[ncent];
	TGraphErrors *gyield_sys_eta_int[ncent];

	for (int icent=0; icent<ncent; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			//netabin = hETA[icent][iarm][1][2]->GetNbinsX();
			//cout << netabin << endl;

			for (int ieta=0; ieta<netabin; ieta++){

				float xx = hETA[icent][iarm][1][2]->GetBinCenter(ieta+1);

				eta_mean[iarm][ieta] = xx;
				eta_err[iarm][ieta] = 0.5*(hETA[icent][iarm][1][2]->GetBinWidth(ieta+1));
				eta_syserr[iarm][ieta] = 0.1;

				float eff_eta = heff_eta[iarm]->GetBinContent(heff_eta[iarm]->FindBin(xx));

				if ( eff_eta<1e-10 ){
					yield_eta[icent][iarm][ieta] = 0.0;
				}else{
					eff_eta *= (heff_eta_cent[icent][iarm]->GetBinContent(heff_eta_cent[icent][iarm]->FindBin(xx))/heff_eta_cent[ncent-1][iarm]->GetBinContent(heff_eta_cent[ncent-1][iarm]->FindBin(xx)));
					/*
					cout << "CENT: " << icent << ", ARM: " << iarm << ", ETA: " << ieta 
						<< ", A: " << (heff_eta_cent[icent][iarm]->GetBinContent(heff_eta_cent[icent][iarm]->FindBin(xx))/heff_eta_cent[ncent-1][iarm]->GetBinContent(heff_eta_cent[ncent-1][iarm]->FindBin(xx)))
						<< endl;
					*/
					yield_eta[icent][iarm][ieta] = hETA[icent][iarm][0][2]->GetBinContent(ieta+1)/eff_eta;
					yield_eta[icent][iarm][ieta] += hETA[icent][iarm][1][2]->GetBinContent(ieta+1)/eff_eta;
					yield_int_eta[icent][iarm] += yield_eta[icent][iarm][ieta];
				}

				if ( iarm==0 ) eta_mean[iarm][ieta] *= -1;

				float err = hETA[icent][iarm][0][2]->GetBinContent(ieta+1) + hETA[icent][iarm][1][2]->GetBinContent(ieta+1);
				yield_eta_err[icent][iarm][ieta] = sqrt(err)/(err)*yield_eta[icent][iarm][ieta];
				yield_int_eta_err[icent][iarm] += err;

				yield_mb_eta[iarm][ieta] += yield_eta[icent][iarm][ieta]/(2*eta_err[iarm][ieta]);
				yield_mb_eta_err[iarm][ieta] += yield_eta[icent][iarm][ieta];

				yield_eta[icent][iarm][ieta] /= (2*eta_err[iarm][ieta]*Nevt[icent][iarm]);
				yield_eta_err[icent][iarm][ieta] /= (2*eta_err[iarm][ieta]*Nevt[icent][iarm]);

				yield_eta[icent][iarm][ieta] *= BiasF[icent];
				yield_eta_err[icent][iarm][ieta] *= BiasF[icent];

				float syserr = sys_frac_kpi*sys_frac_kpi;
				syserr += sys_frac_shape[iarm]*sys_frac_shape[iarm];
				syserr += sys_frac_proton*sys_frac_proton;
				syserr += sys_frac_muon[iarm]*sys_frac_muon[iarm];
				syserr += sys_frac_fvtx[iarm]*sys_frac_fvtx[iarm];
				yield_int_eta_syserr[icent][iarm] = sqrt(syserr);

				yield_eta_syserr[icent][iarm][ieta] = yield_eta[icent][iarm][ieta] * sqrt(syserr);

				if ( eff_eta<1e-10 ){
					RpA_eta[icent][iarm][ieta] = 0;
					RpA_eta_err[icent][iarm][ieta] = 0;
					RpA_eta_syserr[icent][iarm][ieta] = 0;
				}else{
					float Y_pp = hyield_pp_eta->GetBinContent(hyield_pp_eta->FindBin(eta_mean[iarm][ieta]));
					float Y_pp_err = hyield_pp_eta->GetBinError(hyield_pp_eta->FindBin(eta_mean[iarm][ieta]));
					float Y_pp_syserr = hyield_sys_pp_eta->GetBinError(hyield_sys_pp_eta->FindBin(eta_mean[iarm][ieta]));

					RpA_eta[icent][iarm][ieta] = yield_eta[icent][iarm][ieta] / Y_pp / Ncoll[icent];

					float err1 = yield_eta_err[icent][iarm][ieta]/yield_eta[icent][iarm][ieta];
					float err2 = Y_pp_err / Y_pp;
					RpA_eta_err[icent][iarm][ieta] = RpA_eta[icent][iarm][ieta] * sqrt(err1*err1 + err2*err2);

					err1 = yield_eta_syserr[icent][iarm][ieta]/yield_eta[icent][iarm][ieta];
					err2 = Y_pp_syserr / Y_pp;
					RpA_eta_syserr[icent][iarm][ieta] = RpA_eta[icent][iarm][ieta] * sqrt(err1*err1 + err2*err2);
				}

			}//ieta

			yield_int_eta[icent][iarm] /= (deta[iarm]*Nevt[icent][iarm]);
			yield_int_eta_err[icent][iarm] = yield_int_eta[icent][iarm]*1./sqrt(yield_int_eta_err[icent][iarm]);
			yield_int_eta_syserr[icent][iarm] = yield_int_eta[icent][iarm]*yield_int_eta_syserr[icent][iarm];

			yield_int_eta[icent][iarm] *= BiasF[icent];
			yield_int_eta_err[icent][iarm] *= BiasF[icent];
			yield_int_eta_syserr[icent][iarm] *= BiasF[icent];

			float Y_pp = hyield_pp_eta_int->GetBinContent(hyield_pp_eta_int->FindBin(eta_int_mean[iarm]));
			float Y_pp_err = hyield_pp_eta_int->GetBinError(hyield_pp_eta_int->FindBin(eta_int_mean[iarm]));
			float Y_pp_syserr = hyield_sys_pp_eta_int->GetBinError(hyield_sys_pp_eta_int->FindBin(eta_int_mean[iarm]));

			RpA_int_eta[icent][iarm] = yield_int_eta[icent][iarm] / Y_pp / Ncoll[icent];

			float err1 = yield_int_eta_err[icent][iarm]/yield_int_eta[icent][iarm];
			float err2 = Y_pp_err / Y_pp;
			RpA_int_eta_err[icent][iarm] = RpA_int_eta[icent][iarm] * sqrt(err1*err1 + err2*err2);

			err1 = yield_int_eta_syserr[icent][iarm]/yield_int_eta[icent][iarm];
			err2 = Y_pp_syserr / Y_pp;
			RpA_int_eta_syserr[icent][iarm] = RpA_int_eta[icent][iarm] * sqrt(err1*err1 + err2*err2 + Ncoll_err[icent]*Ncoll_err[icent] + BiasF_err[icent]*BiasF_err[icent]);

			//cout << "CENT: " << icent << ", ARM: " << iarm << ", R: " << RpA_int_eta[icent][iarm] << ", ERR0: " << RpA_int_eta_err[icent][iarm] << ", ERR1: " << RpA_int_eta_syserr[icent][iarm] << endl;

			gyield_eta[icent][iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &yield_eta[icent][iarm][1], &eta_err[iarm][1], &yield_eta_err[icent][iarm][1]);
			gyield_eta[icent][iarm]->SetMarkerStyle(20);
			gyield_eta[icent][iarm]->SetLineWidth(2);
			gyield_eta[icent][iarm]->SetMarkerColor(nColor[icent]);
			gyield_eta[icent][iarm]->SetLineColor(nColor[icent]);

			gyield_sys_eta[icent][iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &yield_eta[icent][iarm][1], &eta_syserr[iarm][1], &yield_eta_syserr[icent][iarm][1]);
			gyield_sys_eta[icent][iarm]->SetMarkerStyle(20);
			gyield_sys_eta[icent][iarm]->SetLineWidth(2);
			gyield_sys_eta[icent][iarm]->SetMarkerColor(nColor[icent]);
			gyield_sys_eta[icent][iarm]->SetLineColor(nColor[icent]);
			gyield_sys_eta[icent][iarm]->SetFillStyle(0);

			gRpA_eta[icent][iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &RpA_eta[icent][iarm][1], 0, &RpA_eta_err[icent][iarm][1]);
			gRpA_eta[icent][iarm]->SetMarkerStyle(20);
			gRpA_eta[icent][iarm]->SetLineWidth(1);
			gRpA_eta[icent][iarm]->SetMarkerColor(nColor[icent]);
			gRpA_eta[icent][iarm]->SetLineColor(nColor[icent]);

			gRpA_sys_eta[icent][iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &RpA_eta[icent][iarm][1], &eta_syserr[iarm][1], &RpA_eta_syserr[icent][iarm][1]);
			gRpA_sys_eta[icent][iarm]->SetMarkerStyle(20);
			gRpA_sys_eta[icent][iarm]->SetLineWidth(1);
			gRpA_sys_eta[icent][iarm]->SetMarkerColor(nColor[icent]);
			gRpA_sys_eta[icent][iarm]->SetLineColor(nColor[icent]);
			gRpA_sys_eta[icent][iarm]->SetFillStyle(0);
		}//iarm

		gyield_eta_int[icent] = new TGraphErrors(narm, eta_int_mean, &yield_int_eta[icent][0], 0, &yield_int_eta_err[icent][0]);
		gyield_eta_int[icent]->SetMarkerStyle(25);
		gyield_eta_int[icent]->SetLineWidth(2);
		gyield_eta_int[icent]->SetMarkerColor(nColor[icent]);
		gyield_eta_int[icent]->SetLineColor(nColor[icent]);

		gyield_sys_eta_int[icent] = new TGraphErrors(narm, eta_int_mean, &yield_int_eta[icent][0], eta_int_err, &yield_int_eta_syserr[icent][0]);
		gyield_sys_eta_int[icent]->SetMarkerStyle(25);
		gyield_sys_eta_int[icent]->SetLineWidth(2);
		gyield_sys_eta_int[icent]->SetMarkerColor(nColor[icent]);
		gyield_sys_eta_int[icent]->SetLineColor(nColor[icent]);
		gyield_sys_eta_int[icent]->SetFillStyle(0);
	}//icent

	//return;

	for (int iarm=0; iarm<narm; iarm++){
		for (int ieta=0; ieta<netabin; ieta++){
			if ( yield_mb_eta_err[iarm][ieta]>0 ){
				float err = 1./sqrt(yield_mb_eta_err[iarm][ieta]);
				yield_mb_eta[iarm][ieta] /= Nevt_mb[iarm];
				yield_mb_eta_err[iarm][ieta] = yield_mb_eta[iarm][ieta] * err; 

				float syserr = sys_frac_kpi*sys_frac_kpi;
				syserr += sys_frac_shape[iarm]*sys_frac_shape[iarm];
				syserr += sys_frac_proton*sys_frac_proton;
				syserr += sys_frac_muon[iarm]*sys_frac_muon[iarm];
				syserr += sys_frac_fvtx[iarm]*sys_frac_fvtx[iarm];

				yield_mb_eta_syserr[iarm][ieta] = yield_mb_eta[iarm][ieta] * sqrt(syserr); 

				yield_mb_eta[iarm][ieta] *= BiasF_mb;
				yield_mb_eta_err[iarm][ieta] *= BiasF_mb;
				yield_mb_eta_syserr[iarm][ieta] *= BiasF_mb;

				float Y_pp = hyield_pp_eta->GetBinContent(hyield_pp_eta->FindBin(eta_mean[iarm][ieta]));
				float Y_pp_err = hyield_pp_eta->GetBinError(hyield_pp_eta->FindBin(eta_mean[iarm][ieta]));
				float Y_pp_syserr = hyield_sys_pp_eta->GetBinError(hyield_sys_pp_eta->FindBin(eta_mean[iarm][ieta]));

				RpA_mb_eta[iarm][ieta] = yield_mb_eta[iarm][ieta] / Y_pp / Ncoll_mb;

				float err1 = yield_mb_eta_err[iarm][ieta]/yield_mb_eta[iarm][ieta];
				float err2 = Y_pp_err / Y_pp;
				RpA_mb_eta_err[iarm][ieta] = RpA_mb_eta[iarm][ieta] * sqrt(err1*err1 + err2*err2);

				err1 = yield_mb_eta_syserr[iarm][ieta]/yield_mb_eta[iarm][ieta];
				err2 = Y_pp_syserr / Y_pp;
				RpA_mb_eta_syserr[iarm][ieta] = RpA_mb_eta[iarm][ieta] * sqrt(err1*err1 + err2*err2);

			}else{

				yield_mb_eta[iarm][ieta] = 0.;
				yield_mb_eta_err[iarm][ieta] = 0.;
				yield_mb_eta_syserr[iarm][ieta] = 0.;

				RpA_mb_eta[iarm][ieta] = 0.;
				RpA_mb_eta_err[iarm][ieta] = 0.;
				RpA_mb_eta_syserr[iarm][ieta] = 0.;

			}
		}

		gyield_mb_eta[iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &yield_mb_eta[iarm][1], 0, &yield_mb_eta_err[iarm][1]);
		gyield_mb_eta[iarm]->SetMarkerStyle(20);
		gyield_mb_eta[iarm]->SetLineWidth(2);
		gyield_mb_eta[iarm]->SetMarkerColor(1);
		gyield_mb_eta[iarm]->SetLineColor(1);

		gyield_mb_sys_eta[iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &yield_mb_eta[iarm][1], &eta_syserr[iarm][1], &yield_mb_eta_syserr[iarm][1]);
		gyield_mb_sys_eta[iarm]->SetMarkerStyle(20);
		gyield_mb_sys_eta[iarm]->SetLineWidth(2);
		gyield_mb_sys_eta[iarm]->SetMarkerColor(1);
		gyield_mb_sys_eta[iarm]->SetLineColor(1);
		gyield_mb_sys_eta[iarm]->SetFillStyle(0);

		gRpA_mb_eta[iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &RpA_mb_eta[iarm][1], 0, &RpA_mb_eta_err[iarm][1]);
		gRpA_mb_eta[iarm]->SetMarkerStyle(20);
		gRpA_mb_eta[iarm]->SetLineWidth(2);
		gRpA_mb_eta[iarm]->SetMarkerColor(1);
		gRpA_mb_eta[iarm]->SetLineColor(1);

		gRpA_mb_sys_eta[iarm] = new TGraphErrors(netabin-2+iarm, &eta_mean[iarm][1], &RpA_mb_eta[iarm][1], &eta_syserr[iarm][1], &RpA_mb_eta_syserr[iarm][1]);
		gRpA_mb_sys_eta[iarm]->SetMarkerStyle(20);
		gRpA_mb_sys_eta[iarm]->SetLineWidth(2);
		gRpA_mb_sys_eta[iarm]->SetMarkerColor(1);
		gRpA_mb_sys_eta[iarm]->SetLineColor(1);
		gRpA_mb_sys_eta[iarm]->SetFillStyle(0);

	}

	TCanvas *c3 = new TCanvas("c3","c3",1.1*400,400);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(-3.0,0,3.0,1.5e-2);
	SetHistoStyle();
	for (int icent=0; icent<ncent; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			gyield_sys_eta[icent][iarm]->Draw("2");
			gyield_eta[icent][iarm]->Draw("p same");
		}//iarm
		gyield_sys_eta_int[icent]->Draw("2");
		gyield_eta_int[icent]->Draw("p");
	}//icent
	gyield_mb_sys_eta[0]->Draw("2");
	gyield_mb_eta[0]->Draw("p same");
	gyield_mb_sys_eta[1]->Draw("2");
	gyield_mb_eta[1]->Draw("p same");

	//return;

	TCanvas *c4 = new TCanvas("c4","c4",1.1*400,400);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(-3.0,0,3.0,3.0);
	SetHistoStyle("#eta","R_{pA}");
	TLegend *leg = new TLegend(0.25,0.65,0.9,0.9);
	leg->SetHeader("p+Au #sqrt{s_{NN}}=200 GeV, h^{#pm}, 2.5<p_{T}<5");
	leg->SetFillStyle(0);
	leg->SetNColumns(2);
	for (int icent=0; icent<ncent; icent++){
		for (int iarm=0; iarm<narm; iarm++){
			gRpA_eta[icent][iarm]->Draw("p same");
			if ( iarm==1 )
				leg->AddEntry(gRpA_eta[icent][iarm],Form("%d-%d%c",cent_array[icent],cent_array[icent+1],'%'),"P");
		}//iarm
	}//icent

	leg->Draw();

	TCanvas *c4_ = new TCanvas("c4_","c4_",400*3,400*2);
	TPad *pp4_[ncent];

	pp4_[0] = new TPad(Form("c0_0"),Form("c0_0"),0.08,0.52,0.38,1);
	pp4_[0]->Draw();
	pp4_[0]->cd();
	gPad->SetBottomMargin(0.000);
	gPad->SetRightMargin(0.000);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.001);
	htmp = (TH1F*)gPad->DrawFrame(-3.2,0.005,3.2,3.5);
	SetHistoStyle("","R_{CP}");
	htmp->GetYaxis()->SetLabelSize(0.06);
	htmp->GetYaxis()->SetTitleSize(0.06);
	htmp->GetYaxis()->SetTitleOffset(1.15);

	line_eta->Draw();

	c4_->cd();

	pp4_[1] = new TPad(Form("c0_1"),Form("c0_1"),0.3795,0.52,0.685,1);
	pp4_[1]->Draw();
	pp4_[1]->cd();
	gPad->SetBottomMargin(0.000);
	//gPad->SetRightMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.000);
	htmp = (TH1F*)gPad->DrawFrame(-3.2,0.005,3.2,3.5);
	SetHistoStyle("","");

	line_eta->Draw();

	c4_->cd();

	pp4_[2] = new TPad(Form("c0_2"),Form("c0_2"),0.6845,0.52,0.98,1);
	pp4_[2]->Draw();
	pp4_[2]->cd();
	gPad->SetBottomMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.001);
	htmp = (TH1F*)gPad->DrawFrame(-3.2,0.005,3.2,3.5);
	SetHistoStyle("","");

	line_eta->Draw();

	c4_->cd();

	pp4_[3] = new TPad(Form("c0_1"),Form("c0_1"),0.08,0.0,0.38,0.52);
	pp4_[3]->Draw();
	pp4_[3]->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.0);
	gPad->SetTopMargin(0.0);
	gPad->SetLeftMargin(0.00);
	htmp = (TH1F*)gPad->DrawFrame(-3.2,0,3.2,3.495);
	SetHistoStyle("#eta","R_{pA}");
	htmp->GetYaxis()->SetLabelSize(0.07);
	htmp->GetYaxis()->SetTitleSize(0.06);
	htmp->GetYaxis()->SetTitleOffset(1.15);
	htmp->GetXaxis()->SetLabelSize(0.07);
	htmp->GetXaxis()->SetTitleSize(0.07);

	line_eta->Draw();

	c4_->cd();
	pp4_[4] = new TPad(Form("c0_1"),Form("c0_1"),0.3795,0.0,0.685,0.52);
	pp4_[4]->Draw();
	pp4_[4]->cd();
	gPad->SetBottomMargin(0.15);
	//gPad->SetRightMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.0);
	gPad->SetLeftMargin(0.000);
	htmp = (TH1F*)gPad->DrawFrame(-3.2,0,3.2,3.495);
	SetHistoStyle("#eta","");
	htmp->GetXaxis()->SetLabelSize(0.07);
	htmp->GetXaxis()->SetTitleSize(0.07);

	line_eta->Draw();

	c4_->cd();
	pp4_[5] = new TPad(Form("c0_1"),Form("c0_1"),0.6845,0.0,0.98,0.52);
	pp4_[5]->Draw();
	pp4_[5]->cd();
	gPad->SetBottomMargin(0.15);
	//gPad->SetRightMargin(0.000);
	gPad->SetRightMargin(0.003);
	gPad->SetTopMargin(0.0);
	gPad->SetLeftMargin(0.000);
	htmp = (TH1F*)gPad->DrawFrame(-3.2,0,3.2,3.495);
	SetHistoStyle("#eta","");
	htmp->GetXaxis()->SetLabelSize(0.07);
	htmp->GetXaxis()->SetTitleSize(0.07);

	line_eta->Draw();

	c4_->cd();
	TGaxis *axis0 = new TGaxis(0.08,0.52,0.08,0.977,0.005,3.5,510);
	axis0->SetLabelSize(0.03);
	axis0->SetTitle("R_{pA}");
	axis0->SetTitleOffset(0.8);
	axis0->Draw();

	TGaxis *axis1 = new TGaxis(0.08,0.08,0.08,0.52,0.000,3.495,510);
	axis1->SetLabelSize(0.03);
	axis1->SetTitle("R_{pA}");
	axis1->SetTitleOffset(0.8);
	axis1->Draw();

	for (int icent=0; icent<ncent; icent++){
		pp4_[icent]->cd();

		gRpA_sys_eta[icent][0]->Draw("2");
		gRpA_sys_eta[icent][1]->Draw("2");
		gRpA_eta[icent][0]->Draw("p");
		gRpA_eta[icent][1]->Draw("p");

		float global_sys = sqrt(Ncoll_err[icent]*Ncoll_err[icent] + BiasF_err[icent]*BiasF_err[icent] + 0.101*0.101);
		TBox *bsys = new TBox(3.0,1-global_sys,3.2,1+global_sys);
		bsys->SetFillColor(1);
		bsys->SetLineColor(1);
		bsys->Draw();

		TLatex *tex = new TLatex(-2.5,3.2,"p+Au#rightarrowh^{#pm}+X #sqrt{s_{NN}}=200 GeV");
		tex->SetTextSize(0.06);
		tex->Draw();

		TLatex *tex1 = new TLatex(-2.5,2.9,Form("2.5<p_{T}<5 GeV/c"));
		tex1->SetTextSize(0.06);
		tex1->Draw();

		TLatex *tex1 = new TLatex(-2.5,2.6,Form("%d-%d%c centrality",cent_array[icent],cent_array[icent+1],'%'));
		tex1->SetTextSize(0.06);
		tex1->Draw();

		if ( icent==0 ){
			TLatex *tex1 = new TLatex(-2.7,0.25,"Au-going");
			tex1->SetTextSize(0.06);
			tex1->Draw();

			TLatex *tex1 = new TLatex(1.3,0.25,"p-going");
			tex1->SetTextSize(0.06);
			tex1->Draw();
		}
	}

	//return;

	TCanvas *c5 = new TCanvas("c5","c5",1.1*800,400);
	c5->Divide(2,1);

	c5->cd(1);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(1.0,0,10.0,3.0);

	line->Draw();
	//bsys->Draw();
	SetHistoStyle("p_{T} (GeV/c)","R_{pA}");

	float global_sys = sqrt(Ncoll_mb_err*Ncoll_mb_err + BiasF_mb_err*BiasF_mb_err + 0.101*0.101);
	TBox *bsys = new TBox(9.7,1-global_sys,10.0,1+global_sys);
	bsys->SetFillColor(1);
	bsys->SetLineColor(1);
	bsys->Draw();

	gRpA_mb_sys[0]->Draw("2");
	gRpA_mb_sys[1]->Draw("2");
	gRpA_mb[0]->Draw("p");
	gRpA_mb[1]->Draw("p");

	c5->cd(2);
	SetPadStyle();

	htmp = (TH1F*)gPad->DrawFrame(-3.0,0,3.0,3.0);
	SetHistoStyle("#eta","R_{pA}");

	line_eta->Draw();

	float global_sys = sqrt(Ncoll_mb_err*Ncoll_mb_err + BiasF_mb_err*BiasF_mb_err + 0.101*0.101);
	TBox *bsys = new TBox(2.8,1-global_sys,3.0,1+global_sys);
	bsys->SetFillColor(1);
	bsys->SetLineColor(1);
	bsys->Draw();

	gRpA_mb_sys_eta[0]->Draw("2");
	gRpA_mb_sys_eta[1]->Draw("2");
	gRpA_mb_eta[0]->Draw("p");
	gRpA_mb_eta[1]->Draw("p");


	TGraphErrors *gRpA_cent[narm];
	TGraphErrors *gRpA_cent_sys[narm];

	for (int iarm=0; iarm<narm; iarm++){
		gRpA_cent[iarm] = new TGraphErrors(ncent);
		gRpA_cent[iarm]->SetMarkerStyle(20+iarm);
		gRpA_cent[iarm]->SetMarkerColor(1);
		gRpA_cent[iarm]->SetLineColor(1);
		gRpA_cent[iarm]->SetLineWidth(2);

		gRpA_cent_sys[iarm] = new TGraphErrors(ncent);
		gRpA_cent_sys[iarm]->SetMarkerStyle(20+iarm);
		gRpA_cent_sys[iarm]->SetMarkerColor(1);
		gRpA_cent_sys[iarm]->SetLineColor(1);
		gRpA_cent_sys[iarm]->SetLineWidth(2);
		gRpA_cent_sys[iarm]->SetFillStyle(0);
		for (int icent=0; icent<ncent; icent++){
			//gRpA_cent[iarm]->SetPoint(icent, Ncoll[icent], RpA_int_eta[icent][iarm]);
			gRpA_cent[iarm]->SetPoint(icent, Npart[icent], RpA_int_eta[icent][iarm]);
			gRpA_cent[iarm]->SetPointError(icent, 0, RpA_int_eta_err[icent][iarm]);

			//gRpA_cent_sys[iarm]->SetPoint(icent, Ncoll[icent], RpA_int_eta[icent][iarm]);
			gRpA_cent_sys[iarm]->SetPoint(icent, Npart[icent], RpA_int_eta[icent][iarm]);
			gRpA_cent_sys[iarm]->SetPointError(icent, 0.2, RpA_int_eta_syserr[icent][iarm]);
		}
	}

	TCanvas *c100 = new TCanvas("c100","c100",1.1*400,400);
	SetPadStyle();

	htmp = (TH1F*)gPad->DrawFrame(0,0,15,2.0);
	SetHistoStyle("<N_{coll}>","R_{pA}");

	gRpA_cent_sys[0]->Draw("2");
	gRpA_cent_sys[1]->Draw("2");

	gRpA_cent[0]->Draw("p");
	gRpA_cent[1]->Draw("p");

	if ( bSAVE ){
		c2_->cd();
		c2_->SaveAs("~/plots/Run15pp200_hadron/Run15pAu200_centbin_RpA_vs_pT.gif");
		c2_->SaveAs("~/plots/Run15pp200_hadron/Run15pAu200_centbin_RpA_vs_pT.pdf");

		c4_->cd();
		c4_->SaveAs("~/plots/Run15pp200_hadron/Run15pAu200_centbin_RpA_vs_eta.gif");
		c4_->SaveAs("~/plots/Run15pp200_hadron/Run15pAu200_centbin_RpA_vs_eta.pdf");
	}

	if ( bWRITE ){
		TFile *outfile = new TFile("outfile_pAu.root","recreate");
		for (int icent=0; icent<ncent; icent++){
			for (int iarm=0; iarm<narm; iarm++){
				gRpA_sys[icent][iarm]->Write(Form("pAu_gRpA_sys_cent%d_arm%d",icent,iarm));
				gRpA[icent][iarm]->Write(Form("pAu_gRpA_cent%d_arm%d",icent,iarm));
				gRpA_sys_eta[icent][iarm]->Write(Form("pAu_gRpA_sys_eta_cent%d_arm%d",icent,iarm));
				gRpA_eta[icent][iarm]->Write(Form("pAu_gRpA_eta_cent%d_arm%d",icent,iarm));
			}
		}

		for (int iarm=0; iarm<narm; iarm++){
			gRpA_mb_sys[iarm]->Write(Form("pAu_gRpA_mb_sys_arm%d",iarm));
			gRpA_mb[iarm]->Write(Form("pAu_gRpA_mb_arm%d",iarm));
			gRpA_mb_sys_eta[iarm]->Write(Form("pAu_gRpA_mb_sys_eta_arm%d",iarm));
			gRpA_mb_eta[iarm]->Write(Form("pAu_gRpA_mb_eta_arm%d",iarm));
			gRpA_cent_sys[iarm]->Write(Form("pAu_gRpA_sys_npart_arm%d",iarm));
			gRpA_cent[iarm]->Write(Form("pAu_gRpA_npart_arm%d",iarm));
		}
		outfile->Close();
	}

	return;



}
