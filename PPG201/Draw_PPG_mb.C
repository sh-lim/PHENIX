#include "/phenix/u/shlim/Style.h"

void Draw_PPG_mb(){

	gStyle->SetOptStat(0);
	gStyle->SetCanvasPreferGL(1);

	TFile *infile_epps16_pAu = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/16.BJpsi_ana/EPPS16_MB/charged_hadrons_EPPS16_Au_histos.root","read");
	TFile *infile_epps16_pAl = new TFile("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/16.BJpsi_ana/EPPS16_MB/charged_hadrons_EPPS16_Al_histos.root","read");
	TProfile *tmod_pAu_fwd[41];
	TProfile *tmod_pAu_bwd[41];
	TProfile *tmod_pAu_eta[41];
	TProfile *tmod_pAl_fwd[41];
	TProfile *tmod_pAl_bwd[41];
	TProfile *tmod_pAl_eta[41];

	for (int iset=0; iset<41; iset++){
		tmod_pAu_fwd[iset] = (TProfile*)infile_epps16_pAu->Get(Form("tmod_fwd_set%d",iset));
		tmod_pAu_bwd[iset] = (TProfile*)infile_epps16_pAu->Get(Form("tmod_bwd_set%d",iset));
		tmod_pAu_eta[iset] = (TProfile*)infile_epps16_pAu->Get(Form("tmod_eta_set%d",iset));

		tmod_pAl_fwd[iset] = (TProfile*)infile_epps16_pAl->Get(Form("tmod_fwd_set%d",iset));
		tmod_pAl_bwd[iset] = (TProfile*)infile_epps16_pAl->Get(Form("tmod_bwd_set%d",iset));
		tmod_pAl_eta[iset] = (TProfile*)infile_epps16_pAl->Get(Form("tmod_eta_set%d",iset));
	}

	float BiasF_pAu_mb_err = 0.016;
	float Ncoll_pAu_mb_err = 0.064;
	float global_err_pAu = sqrt(BiasF_pAu_mb_err*BiasF_pAu_mb_err + Ncoll_pAu_mb_err*Ncoll_pAu_mb_err + 0.1*0.1);

	float BiasF_pAl_mb_err = 0.025;
	float Ncoll_pAl_mb_err = 0.047;
	float global_err_pAl = sqrt(BiasF_pAl_mb_err*BiasF_pAl_mb_err + Ncoll_pAl_mb_err*Ncoll_pAl_mb_err + 0.1*0.1);

	TBox *sysbox_pAu_pT = new TBox(0,1-global_err_pAu,0.3,1+global_err_pAu);
	sysbox_pAu_pT->SetFillColor(1);
	sysbox_pAu_pT->SetLineColor(1);

	TBox *sysbox_pAl_pT = new TBox(0,1-global_err_pAl,0.3,1+global_err_pAl);
	sysbox_pAl_pT->SetFillColor(4);
	sysbox_pAl_pT->SetLineColor(4);

	TBox *sysbox_pAu_eta = new TBox(2.8,1-global_err_pAu,3.0,1+global_err_pAu);
	sysbox_pAu_eta->SetFillColor(1);
	sysbox_pAu_eta->SetLineColor(1);

	TBox *sysbox_pAl_eta = new TBox(2.8,1-global_err_pAl,3.0,1+global_err_pAl);
	sysbox_pAl_eta->SetFillColor(4);
	sysbox_pAl_eta->SetLineColor(4);

	TCanvas *c0_pAu = new TCanvas("c0_pAu","c0_pAu",1.1*800,400);
	c0_pAu->Divide(2,1);
	for (int iarm=0; iarm<2; iarm++){
		c0_pAu->cd(iarm+1);
		SetPadStyle();

		htmp = (TH1F*)gPad->DrawFrame(0,0,10,2);
		SetHistoStyle();

		for (int iset=0; iset<41; iset++){
			if ( iarm==0 ) tmod_pAu_bwd[iset]->Draw("same");
			else tmod_pAu_fwd[iset]->Draw("same");
		}
	}

	TCanvas *c0_pAl = new TCanvas("c0_pAl","c0_pAl",1.1*800,400);
	c0_pAl->Divide(2,1);
	for (int iarm=0; iarm<2; iarm++){
		c0_pAl->cd(iarm+1);
		SetPadStyle();

		htmp = (TH1F*)gPad->DrawFrame(0,0,10,2);
		SetHistoStyle();

		for (int iset=0; iset<41; iset++){
			if ( iarm==0 ) tmod_pAl_bwd[iset]->Draw("same");
			else tmod_pAl_fwd[iset]->Draw("same");
		}
	}
	//return;

	int epps16_eta_count = 0;
	double epps16_xx_eta[200];
	double epps16_RpAu_eta[200], epps16_RpAu_eta_up[200], epps16_RpAu_eta_dn[200];
	double epps16_RpAl_eta[200], epps16_RpAl_eta_up[200], epps16_RpAl_eta_dn[200];

	for (int ix=0; ix<tmod_pAu_eta[0]->GetNbinsX(); ix++){
		if ( tmod_pAu_eta[0]->GetBinContent(ix+1)<1e-5 ) continue;

		epps16_xx_eta[epps16_eta_count] = tmod_pAu_eta[0]->GetBinCenter(ix+1);
		epps16_RpAu_eta[epps16_eta_count] = tmod_pAu_eta[0]->GetBinContent(ix+1);
		epps16_RpAl_eta[epps16_eta_count] = tmod_pAl_eta[0]->GetBinContent(ix+1);

		epps16_RpAu_eta_up[epps16_eta_count] = epps16_RpAu_eta_dn[epps16_eta_count] = 0.0;
		epps16_RpAl_eta_up[epps16_eta_count] = epps16_RpAl_eta_dn[epps16_eta_count] = 0.0;

		for (int iset=0; iset<20; iset++){
			double deltadn = tmod_pAu_eta[2*iset+1]->GetBinContent(ix+1) - epps16_RpAu_eta[epps16_eta_count];
			double deltaup = tmod_pAu_eta[2*iset+2]->GetBinContent(ix+1) - epps16_RpAu_eta[epps16_eta_count];

			double dsigcemup = TMath::Max(deltaup, deltadn);
			if ( dsigcemup<0 ) dsigcemup = 0.0; 
			double dsigcemdn = TMath::Max(-deltaup, -deltadn);
			if ( dsigcemdn<0 ) dsigcemdn = 0.0;

			epps16_RpAu_eta_up[epps16_eta_count] += TMath::Power(dsigcemup, 2);
			epps16_RpAu_eta_dn[epps16_eta_count] += TMath::Power(dsigcemdn, 2);

			deltadn = tmod_pAl_eta[2*iset+1]->GetBinContent(ix+1) - epps16_RpAl_eta[epps16_eta_count];
			deltaup = tmod_pAl_eta[2*iset+2]->GetBinContent(ix+1) - epps16_RpAl_eta[epps16_eta_count];

			dsigcemup = TMath::Max(deltaup, deltadn);
			if ( dsigcemup<0 ) dsigcemup = 0.0; 
			dsigcemdn = TMath::Max(-deltaup, -deltadn);
			if ( dsigcemdn<0 ) dsigcemdn = 0.0;

			epps16_RpAl_eta_up[epps16_eta_count] += TMath::Power(dsigcemup, 2);
			epps16_RpAl_eta_dn[epps16_eta_count] += TMath::Power(dsigcemdn, 2);

		}//iset

		epps16_RpAu_eta_up[epps16_eta_count] = sqrt(epps16_RpAu_eta_up[epps16_eta_count]);
		epps16_RpAu_eta_dn[epps16_eta_count] = sqrt(epps16_RpAu_eta_dn[epps16_eta_count]);

		epps16_RpAl_eta_up[epps16_eta_count] = sqrt(epps16_RpAl_eta_up[epps16_eta_count]);
		epps16_RpAl_eta_dn[epps16_eta_count] = sqrt(epps16_RpAl_eta_dn[epps16_eta_count]);

		//cout << epps16_RpAu_eta[epps16_eta_count] << ", " << epps16_RpAu_eta_up[epps16_eta_count] << ", " << epps16_RpAu_eta_dn[epps16_eta_count] << endl;

		epps16_eta_count++;
	}//ix

	int epps16_pT_count[2] = {0};
	double epps16_xx_pT[2][100];
	double epps16_RpAu_pT[2][100], epps16_RpAu_pT_up[2][100], epps16_RpAu_pT_dn[2][100];
	double epps16_RpAl_pT[2][100], epps16_RpAl_pT_up[2][100], epps16_RpAl_pT_dn[2][100];

	for (int ix=0; ix<tmod_pAu_bwd[0]->GetNbinsX(); ix++){
		if ( tmod_pAu_bwd[0]->GetBinContent(ix+1)<1e-5 ) continue;

		epps16_xx_pT[0][epps16_pT_count[0]] = tmod_pAu_bwd[0]->GetBinCenter(ix+1);
		epps16_RpAu_pT[0][epps16_pT_count[0]] = tmod_pAu_bwd[0]->GetBinContent(ix+1);
		epps16_RpAl_pT[0][epps16_pT_count[0]] = tmod_pAl_bwd[0]->GetBinContent(ix+1);

		epps16_RpAu_pT_up[0][epps16_pT_count[0]] = epps16_RpAu_pT_dn[0][epps16_pT_count[0]] = 0.0;
		epps16_RpAl_pT_up[0][epps16_pT_count[0]] = epps16_RpAl_pT_dn[0][epps16_pT_count[0]] = 0.0;

		for (int iset=0; iset<20; iset++){
			double deltadn = tmod_pAu_bwd[2*iset+1]->GetBinContent(ix+1) - epps16_RpAu_pT[0][epps16_pT_count[0]];
			double deltaup = tmod_pAu_bwd[2*iset+2]->GetBinContent(ix+1) - epps16_RpAu_pT[0][epps16_pT_count[0]];

			double dsigcemup = TMath::Max(deltaup, deltadn);
			if ( dsigcemup<0 ) dsigcemup = 0.0; 
			double dsigcemdn = TMath::Max(-deltaup, -deltadn);
			if ( dsigcemdn<0 ) dsigcemdn = 0.0;

			epps16_RpAu_pT_up[0][epps16_pT_count[0]] += TMath::Power(dsigcemup, 2);
			epps16_RpAu_pT_dn[0][epps16_pT_count[0]] += TMath::Power(dsigcemdn, 2);

			deltadn = tmod_pAl_bwd[2*iset+1]->GetBinContent(ix+1) - epps16_RpAl_pT[0][epps16_pT_count[0]];
			deltaup = tmod_pAl_bwd[2*iset+2]->GetBinContent(ix+1) - epps16_RpAl_pT[0][epps16_pT_count[0]];

			dsigcemup = TMath::Max(deltaup, deltadn);
			if ( dsigcemup<0 ) dsigcemup = 0.0; 
			dsigcemdn = TMath::Max(-deltaup, -deltadn);
			if ( dsigcemdn<0 ) dsigcemdn = 0.0;

			epps16_RpAl_pT_up[0][epps16_pT_count[0]] += TMath::Power(dsigcemup, 2);
			epps16_RpAl_pT_dn[0][epps16_pT_count[0]] += TMath::Power(dsigcemdn, 2);

		}//iset

		epps16_RpAu_pT_up[0][epps16_pT_count[0]] = sqrt(epps16_RpAu_pT_up[0][epps16_pT_count[0]]);
		epps16_RpAu_pT_dn[0][epps16_pT_count[0]] = sqrt(epps16_RpAu_pT_dn[0][epps16_pT_count[0]]);

		epps16_RpAl_pT_up[0][epps16_pT_count[0]] = sqrt(epps16_RpAl_pT_up[0][epps16_pT_count[0]]);
		epps16_RpAl_pT_dn[0][epps16_pT_count[0]] = sqrt(epps16_RpAl_pT_dn[0][epps16_pT_count[0]]);

		epps16_pT_count[0]++;
	}//ix

	for (int ix=0; ix<tmod_pAu_fwd[0]->GetNbinsX(); ix++){
		if ( tmod_pAu_fwd[0]->GetBinContent(ix+1)<1e-5 ) continue;

		epps16_xx_pT[1][epps16_pT_count[1]] = tmod_pAu_fwd[0]->GetBinCenter(ix+1);
		epps16_RpAu_pT[1][epps16_pT_count[1]] = tmod_pAu_fwd[0]->GetBinContent(ix+1);
		epps16_RpAl_pT[1][epps16_pT_count[1]] = tmod_pAl_fwd[0]->GetBinContent(ix+1);

		epps16_RpAu_pT_up[1][epps16_pT_count[1]] = epps16_RpAu_pT_dn[1][epps16_pT_count[1]] = 0.0;
		epps16_RpAl_pT_up[1][epps16_pT_count[1]] = epps16_RpAl_pT_dn[1][epps16_pT_count[1]] = 0.0;

		for (int iset=0; iset<20; iset++){
			double deltadn = tmod_pAu_fwd[2*iset+1]->GetBinContent(ix+1) - epps16_RpAu_pT[1][epps16_pT_count[1]];
			double deltaup = tmod_pAu_fwd[2*iset+2]->GetBinContent(ix+1) - epps16_RpAu_pT[1][epps16_pT_count[1]];

			double dsigcemup = TMath::Max(deltaup, deltadn);
			if ( dsigcemup<0 ) dsigcemup = 0.0; 
			double dsigcemdn = TMath::Max(-deltaup, -deltadn);
			if ( dsigcemdn<0 ) dsigcemdn = 0.0;

			epps16_RpAu_pT_up[1][epps16_pT_count[1]] += TMath::Power(dsigcemup, 2);
			epps16_RpAu_pT_dn[1][epps16_pT_count[1]] += TMath::Power(dsigcemdn, 2);

			deltadn = tmod_pAl_fwd[2*iset+1]->GetBinContent(ix+1) - epps16_RpAl_pT[1][epps16_pT_count[1]];
			deltaup = tmod_pAl_fwd[2*iset+2]->GetBinContent(ix+1) - epps16_RpAl_pT[1][epps16_pT_count[1]];

			dsigcemup = TMath::Max(deltaup, deltadn);
			if ( dsigcemup<0 ) dsigcemup = 0.0; 
			dsigcemdn = TMath::Max(-deltaup, -deltadn);
			if ( dsigcemdn<0 ) dsigcemdn = 0.0;

			epps16_RpAl_pT_up[1][epps16_pT_count[1]] += TMath::Power(dsigcemup, 2);
			epps16_RpAl_pT_dn[1][epps16_pT_count[1]] += TMath::Power(dsigcemdn, 2);

		}//iset

		epps16_RpAu_pT_up[1][epps16_pT_count[1]] = sqrt(epps16_RpAu_pT_up[1][epps16_pT_count[1]]);
		epps16_RpAu_pT_dn[1][epps16_pT_count[1]] = sqrt(epps16_RpAu_pT_dn[1][epps16_pT_count[1]]);

		epps16_RpAl_pT_up[1][epps16_pT_count[1]] = sqrt(epps16_RpAl_pT_up[1][epps16_pT_count[1]]);
		epps16_RpAl_pT_dn[1][epps16_pT_count[1]] = sqrt(epps16_RpAl_pT_dn[1][epps16_pT_count[1]]);

		epps16_pT_count[1]++;
	}//ix

	//return;

	TGraphAsymmErrors *gepps16_pAu_eta = new TGraphAsymmErrors(epps16_eta_count, epps16_xx_eta, epps16_RpAu_eta, 0, 0, epps16_RpAu_eta_dn, epps16_RpAu_eta_up);
	gepps16_pAu_eta->SetLineColor(0);
	//gepps16_pAu_eta->SetFillColor(kGray+3);
	gepps16_pAu_eta->SetFillColor(15);
	gepps16_pAu_eta->SetFillStyle(3001);

	TGraphAsymmErrors *gepps16_pAl_eta = new TGraphAsymmErrors(epps16_eta_count, epps16_xx_eta, epps16_RpAl_eta, 0, 0, epps16_RpAl_eta_dn, epps16_RpAl_eta_up);
	gepps16_pAl_eta->SetLineColor(0);
	gepps16_pAl_eta->SetFillColor(kBlue-7);
	gepps16_pAl_eta->SetFillStyle(3001);

	TGraphAsymmErrors *gepps16_pAu_pT_bwd = new TGraphAsymmErrors(epps16_pT_count[0], &epps16_xx_pT[0][0], &epps16_RpAu_pT[0][0], 0, 0, &epps16_RpAu_pT_dn[0][0], &epps16_RpAu_pT_up[0][0]);
	gepps16_pAu_pT_bwd->SetLineColor(0);
	gepps16_pAu_pT_bwd->SetFillColorAlpha(kGreen+2, 0.7);
	//gepps16_pAu_pT_bwd->SetFillStyle(3354);
	gepps16_pAu_pT_bwd->SetFillStyle(1001);

	TGraphAsymmErrors *gepps16_pAu_pT_fwd = new TGraphAsymmErrors(epps16_pT_count[1], &epps16_xx_pT[1][0], &epps16_RpAu_pT[1][0], 0, 0, &epps16_RpAu_pT_dn[1][0], &epps16_RpAu_pT_up[1][0]);
	gepps16_pAu_pT_fwd->SetLineColor(0);
	gepps16_pAu_pT_fwd->SetFillColorAlpha(kRed-4, 0.5);
	//gepps16_pAu_pT_fwd->SetFillStyle(3345);
	gepps16_pAu_pT_fwd->SetFillStyle(3001);

	TGraphAsymmErrors *gepps16_pAl_pT_bwd = new TGraphAsymmErrors(epps16_pT_count[0], &epps16_xx_pT[0][0], &epps16_RpAl_pT[0][0], 0, 0, &epps16_RpAl_pT_dn[0][0], &epps16_RpAl_pT_up[0][0]);
	gepps16_pAl_pT_bwd->SetLineColor(0);
	gepps16_pAl_pT_bwd->SetFillColor(kGreen+2);
	//gepps16_pAl_pT_bwd->SetFillStyle(3354);
	gepps16_pAl_pT_bwd->SetFillStyle(3001);

	TGraphAsymmErrors *gepps16_pAl_pT_fwd = new TGraphAsymmErrors(epps16_pT_count[1], &epps16_xx_pT[1][0], &epps16_RpAl_pT[1][0], 0, 0, &epps16_RpAl_pT_dn[1][0], &epps16_RpAl_pT_up[1][0]);
	gepps16_pAl_pT_fwd->SetLineColor(0);
	gepps16_pAl_pT_fwd->SetFillColor(kRed-7);
	//gepps16_pAl_pT_fwd->SetFillStyle(3345);
	gepps16_pAl_pT_fwd->SetFillStyle(3001);

	TLegendEntry *le;

	const int narm = 2;
	const bool bSAVE = false;

	TFile *infile_pAl = new TFile("outfile_pAl.root","read");
	TFile *infile_pAu = new TFile("outfile_pAu.root","read");

	TGraphErrors *gpAu_pT_sys[narm];
	TGraphErrors *gpAu_pT[narm];
	TGraphErrors *gpAl_pT_sys[narm];
	TGraphErrors *gpAl_pT[narm];

	TGraphErrors *gpAu_eta_sys[narm];
	TGraphErrors *gpAu_eta[narm];
	TGraphErrors *gpAl_eta_sys[narm];
	TGraphErrors *gpAl_eta[narm];

	TGraphErrors *gpAu_npart_sys[narm];
	TGraphErrors *gpAu_npart[narm];
	TGraphErrors *gpAl_npart_sys[narm];
	TGraphErrors *gpAl_npart[narm];

	for (int iarm=0; iarm<narm; iarm++){
		gpAu_pT_sys[iarm] = (TGraphErrors*)infile_pAu->Get(Form("pAu_gRpA_mb_sys_arm%d",iarm));
		gpAu_pT[iarm] = (TGraphErrors*)infile_pAu->Get(Form("pAu_gRpA_mb_arm%d",iarm));
		gpAu_eta_sys[iarm] = (TGraphErrors*)infile_pAu->Get(Form("pAu_gRpA_mb_sys_eta_arm%d",iarm));
		gpAu_eta[iarm] = (TGraphErrors*)infile_pAu->Get(Form("pAu_gRpA_mb_eta_arm%d",iarm));
		gpAu_npart_sys[iarm] = (TGraphErrors*)infile_pAu->Get(Form("pAu_gRpA_sys_npart_arm%d",iarm));
		gpAu_npart[iarm] = (TGraphErrors*)infile_pAu->Get(Form("pAu_gRpA_npart_arm%d",iarm));

		gpAl_pT_sys[iarm] = (TGraphErrors*)infile_pAl->Get(Form("pAl_gRpA_mb_sys_arm%d",iarm));
		gpAl_pT[iarm] = (TGraphErrors*)infile_pAl->Get(Form("pAl_gRpA_mb_arm%d",iarm));
		gpAl_eta_sys[iarm] = (TGraphErrors*)infile_pAl->Get(Form("pAl_gRpA_mb_sys_eta_arm%d",iarm));
		gpAl_eta[iarm] = (TGraphErrors*)infile_pAl->Get(Form("pAl_gRpA_mb_eta_arm%d",iarm));
		gpAl_npart_sys[iarm] = (TGraphErrors*)infile_pAl->Get(Form("pAl_gRpA_sys_npart_arm%d",iarm));
		gpAl_npart[iarm] = (TGraphErrors*)infile_pAl->Get(Form("pAl_gRpA_npart_arm%d",iarm));

		gpAl_pT_sys[iarm]->SetLineColor(4);
		gpAl_pT_sys[iarm]->SetMarkerColor(4);
		gpAl_pT_sys[iarm]->SetLineWidth(1);
		gpAl_pT[iarm]->SetLineColor(4);
		gpAl_pT[iarm]->SetMarkerColor(4);
		gpAl_pT[iarm]->SetLineWidth(1);
		gpAl_eta_sys[iarm]->SetLineColor(4);
		gpAl_eta_sys[iarm]->SetMarkerColor(4);
		gpAl_eta_sys[iarm]->SetLineWidth(1);
		gpAl_eta[iarm]->SetLineColor(4);
		gpAl_eta[iarm]->SetMarkerColor(4);
		gpAl_eta[iarm]->SetLineWidth(1);
		gpAl_npart_sys[iarm]->SetLineColor(4);
		gpAl_npart_sys[iarm]->SetMarkerColor(4);
		gpAl_npart_sys[iarm]->SetLineWidth(1);
		gpAl_npart[iarm]->SetLineColor(4);
		gpAl_npart[iarm]->SetLineWidth(1);
		gpAl_npart[iarm]->SetMarkerColor(4);
		gpAl_npart[iarm]->SetMarkerSize(1.2);

		gpAu_pT[iarm]->SetLineWidth(1);
		gpAu_pT_sys[iarm]->SetLineWidth(1);
		gpAu_eta[iarm]->SetLineWidth(1);
		gpAu_eta_sys[iarm]->SetLineWidth(1);
		gpAu_npart[iarm]->SetMarkerSize(1.2);
		gpAu_npart[iarm]->SetLineWidth(2);
		gpAu_npart_sys[iarm]->SetLineWidth(3);
		gpAu_npart_sys[iarm]->SetLineWidth(1);
	}

	TLine *line = new TLine(0,1,10,1);
	line->SetLineWidth(2);
	line->SetLineStyle(2);

	TCanvas *c1_0 = new TCanvas("c1_0","c1_0",1.2*500,500);

	SetPadStyle();
	gPad->SetRightMargin(0.03);
	gPad->SetLeftMargin(0.15);
	htmp = (TH1F*)gPad->DrawFrame(0,0,10,3.0);
	SetHistoStyle("p_{T} (GeV/c)","R_{pA}");
	htmp->GetYaxis()->SetTitleOffset(1.1);
	line->Draw();
	gepps16_pAu_pT_bwd->Draw("l 3");
	gepps16_pAu_pT_fwd->Draw("l 3");
	gpAu_pT_sys[1]->SetMarkerStyle(24);
	gpAu_pT[1]->SetMarkerStyle(24);
	gpAu_pT_sys[0]->Draw("2");
	gpAu_pT_sys[1]->Draw("2");
	gpAu_pT[0]->Draw("p");
	gpAu_pT[1]->Draw("p");

	sysbox_pAu_pT->Draw();

	TLatex *tex = new TLatex(1.4,2.75,"p+Au#rightarrowh^{#pm}+X #sqrt{s_{NN}}=200 GeV");
	tex->SetTextSize(0.055);
	tex->Draw();

	TLatex *tex1 = new TLatex(1.4,2.50,"0-100% centrality");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLegend *leg = new TLegend(0.30,0.65,0.95,0.78);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	le = leg->AddEntry(gpAu_pT[0],"-2.2<#eta<-1.2 (Au-going)","P");
	le->SetTextSize(0.05);
	le = leg->AddEntry(gpAu_pT[1],"1.2<#eta<2.4 (p-going)","P");
	le->SetTextSize(0.05);
	leg->Draw();

	TLegend *leg = new TLegend(0.32,0.17,0.95,0.30);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	le = leg->AddEntry(gepps16_pAu_pT_bwd,"EPPS16+PYTHIA, -2.2<#eta<-1.2","F");
	le->SetTextSize(0.04);
	le = leg->AddEntry(gepps16_pAu_pT_fwd,"EPPS16+PYTHIA, 1.2<#eta<2.4","F");
	le->SetTextSize(0.04);
	leg->Draw();

	return;

	c1->cd(2);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(0,0,10,3.0);
	SetHistoStyle("p_{T} (GeV/c)","R_{pA}");
	line->Draw();
	gepps16_pAl_pT_bwd->Draw("l 3");
	gepps16_pAl_pT_fwd->Draw("l 3");
	gpAl_pT[1]->SetMarkerStyle(24);
	gpAl_pT_sys[1]->SetMarkerStyle(24);
	gpAl_pT_sys[0]->Draw("2");
	gpAl_pT_sys[1]->Draw("2");
	gpAl_pT[0]->Draw("p");
	gpAl_pT[1]->Draw("p");

	sysbox_pAl_pT->Draw();

	TLatex *tex = new TLatex(1.4,2.75,"p+Al#rightarrowh^{#pm}+X #sqrt{s_{NN}}=200 GeV");
	tex->SetTextSize(0.055);
	tex->Draw();

	TLatex *tex1 = new TLatex(1.4,2.50,"0-100% centrality");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLegend *leg = new TLegend(0.30,0.65,0.95,0.78);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	le = leg->AddEntry(gpAl_pT[0],"-2.2<#eta<-1.2 (Al-going)","P");
	le->SetTextSize(0.055);
	le = leg->AddEntry(gpAl_pT[1],"1.2<#eta<2.4 (p-going)","P");
	le->SetTextSize(0.055);
	leg->Draw();

	TLegend *leg = new TLegend(0.32,0.17,0.95,0.30);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	le = leg->AddEntry(gepps16_pAu_pT_bwd,"EPPS16+PYTHIA, -2.2<#eta<-1.2","F");
	le->SetTextSize(0.04);
	le = leg->AddEntry(gepps16_pAu_pT_fwd,"EPPS16+PYTHIA, 1.2<#eta<2.4","F");
	le->SetTextSize(0.04);
	leg->Draw();

	TLine *line_eta = new TLine(-3,1,3,1);
	line_eta->SetLineWidth(2);
	line_eta->SetLineStyle(2);

	TCanvas *c2 = new TCanvas("c2","c2",1.1*800*1.5,400*1.5);
	c2->Divide(2,1);

	c2->cd(1);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(-3,0,3,2.5);
	SetHistoStyle("#eta","R_{pA}");
	line_eta->Draw();
	gepps16_pAu_eta->Draw("l 3");
	gpAu_eta_sys[0]->Draw("2");
	gpAu_eta_sys[1]->Draw("2");
	gpAu_eta[0]->Draw("p");
	gpAu_eta[1]->Draw("p");

	sysbox_pAu_eta->Draw();

	TLatex *tex = new TLatex(-2,2.25,"p+Au#rightarrowh^{#pm}+X #sqrt{s_{NN}}=200 GeV");
	tex->SetTextSize(0.055);
	tex->Draw();

	TLatex *tex1 = new TLatex(-0.5,2.0,"0-100% centrality");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLatex *tex1 = new TLatex(-0.1,1.75,"2.5<p_{T}<5 GeV/c");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLatex *tex1 = new TLatex(-2.7,0.25,"Au-going");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLatex *tex1 = new TLatex(1.3,0.25,"p-going");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLegend *leg = new TLegend(0.50,0.60,0.95,0.65);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	le = leg->AddEntry(gepps16_pAu_eta,"EPPS16+PYTHIA","F");
	le->SetTextSize(0.045);
	leg->Draw();

	c2->cd(2);
	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(-3,0,3,2.5);
	SetHistoStyle("#eta","R_{pA}");
	line_eta->Draw();
	gepps16_pAl_eta->Draw("l 3");
	gpAl_eta_sys[0]->Draw("2");
	gpAl_eta_sys[1]->Draw("2");
	gpAl_eta[0]->Draw("p");
	gpAl_eta[1]->Draw("p");

	sysbox_pAl_eta->Draw();

	TLatex *tex = new TLatex(-2,2.25,"p+Al#rightarrowh^{#pm}+X #sqrt{s_{NN}}=200 GeV");
	tex->SetTextSize(0.055);
	tex->Draw();

	TLatex *tex1 = new TLatex(-0.5,2.0,"0-100% centrality");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLatex *tex1 = new TLatex(-0.1,1.75,"2.5<p_{T}<5 GeV/c");
	tex1->SetTextSize(0.055);
	tex1->Draw();

	TLatex *tex1 = new TLatex(-2.7,0.25,"Al-going");
	tex1->SetTextSize(0.055);
	tex1->SetTextColor(4);
	tex1->Draw();

	TLatex *tex1 = new TLatex(1.3,0.25,"p-going");
	tex1->SetTextSize(0.055);
	tex1->SetTextColor(4);
	tex1->Draw();

	TLegend *leg = new TLegend(0.50,0.60,0.95,0.65);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	le = leg->AddEntry(gepps16_pAl_eta,"EPPS16+PYTHIA","F");
	le->SetTextSize(0.045);
	leg->Draw();

	TCanvas *c3 = new TCanvas("c3","c3",1.2*800,800);

	SetPadStyle();
	htmp = (TH1F*)gPad->DrawFrame(0,0,15,3.5);
	SetHistoStyle("<N_{part}>","R_{pA}");

	TLine *line_npart = new TLine(0,1,15,1);
	line_npart->SetLineWidth(2);
	line_npart->SetLineStyle(2);
	line_npart->Draw();

	gpAu_npart_sys[1]->SetMarkerStyle(24);
	gpAu_npart[1]->SetMarkerStyle(24);
	gpAl_npart_sys[1]->SetMarkerSize(1.7);
	gpAl_npart_sys[1]->SetMarkerStyle(27);
	gpAl_npart[1]->SetMarkerSize(1.7);
	gpAl_npart[1]->SetMarkerStyle(27);
	gpAl_npart_sys[0]->SetMarkerSize(1.7);
	gpAl_npart_sys[0]->SetMarkerStyle(33);
	gpAl_npart[0]->SetMarkerSize(1.7);
	gpAl_npart[0]->SetMarkerStyle(33);

	gpAl_npart_sys[0]->Draw("2");
	gpAu_npart_sys[0]->Draw("2");
	gpAl_npart_sys[1]->Draw("2");
	gpAu_npart_sys[1]->Draw("2");
	gpAl_npart[0]->Draw("p");
	gpAu_npart[0]->Draw("p");
	gpAl_npart[1]->Draw("p");
	gpAu_npart[1]->Draw("p");

	TLatex *tex1 = new TLatex(1.5,2.3,"h^{#pm}, 2.5<p_{T}<5 GeV/c");
	tex1->SetTextSize(0.045);
	tex1->Draw();

	TLatex *tex1 = new TLatex(1.5,0.25,"10% global sys. uncertainty");
	tex1->SetTextSize(0.045);
	tex1->Draw();
	
	TLegend *leg = new TLegend(0.20,0.72,0.95,0.9);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	le = leg->AddEntry(gpAu_npart[0],"-2.2<#eta<-1.2, p+Au (Au-going)","P");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gpAl_npart[0],"-2.2<#eta<-1.2, p+Al (Al-going)","P");
	le->SetTextSize(0.045);
	le->SetTextColor(kBlue);
	le = leg->AddEntry(gpAu_npart[1],"1.2<#eta<2.4, p+Au (p-going)","P");
	le->SetTextSize(0.045);
	le = leg->AddEntry(gpAl_npart[1],"1.2<#eta<2.4, p+Al (p-going)","P");
	le->SetTextSize(0.045);
	le->SetTextColor(kBlue);
	leg->Draw();

	if ( bSAVE ){
		c1->cd();
		c1->SaveAs("~/plots/Run15pp200_hadron/Run15pAupAl200_mb_hadron_RpA_pT.gif");
		c1->SaveAs("~/plots/Run15pp200_hadron/Run15pAupAl200_mb_hadron_RpA_pT.pdf");

		c2->cd();
		c2->SaveAs("~/plots/Run15pp200_hadron/Run15pAupAl200_mb_hadron_RpA_eta.gif");
		c2->SaveAs("~/plots/Run15pp200_hadron/Run15pAupAl200_mb_hadron_RpA_eta.pdf");

		c3->cd();
		c3->SaveAs("~/plots/Run15pp200_hadron/Run15pAupAl200_mb_hadron_RpA_Npart.gif");
		c3->SaveAs("~/plots/Run15pp200_hadron/Run15pAupAl200_mb_hadron_RpA_Npart.pdf");
	}

}
