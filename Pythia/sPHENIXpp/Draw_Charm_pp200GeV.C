#include "/phenix/u/shlim/Style.h"

void Draw_Charm_pp200GeV(){

	float etarange = 1.0;

	const int nset = 5;

	const int nMarker[5] = {24, 25, 26, 20, 21};
	const int nColor[5] = {1, 2, 4, 6, 8};

	TFile *infile[nset];

	//infile[0] = new TFile("outfile_hist_pp200GeV_softqcd_grp002.root","read");
	//infile[1] = new TFile("outfile_hist_pp200GeV_softqcd_grp003.root","read");
	//infile[2] = new TFile("outfile_hist_pp200GeV_softqcd_grp004.root","read");
	//infile[3] = new TFile("outfile_hist_pp200GeV_softqcd_grp005.root","read");
	//infile[4] = new TFile("outfile_hist_pp200GeV_softqcd_grp006.root","read");

	infile[0] = new TFile("outfile_hist_pp200GeV_inelastic_grp000.root","read");
	infile[1] = new TFile("outfile_hist_pp200GeV_inelastic_grp001.root","read");
	infile[2] = new TFile("outfile_hist_pp200GeV_inelastic_grp002.root","read");
	infile[3] = new TFile("outfile_hist_pp200GeV_inelastic_grp003.root","read");
	infile[4] = new TFile("outfile_hist_pp200GeV_inelastic_grp004.root","read");

	char *dataset[nset] = {
		"default",
		"mode0",
		"mode2"
	};

	TCanvas *c0 = new TCanvas("c0","c0",1.2*3*250,2*250); 
	c0->Divide(3,2);

	TH1F *hxection[nset];
	float Xection_hardqcd[nset], nevent[nset];
	for (int iset=0; iset<nset; iset++){

		infile[iset]->cd();
		hxection[iset] = (TH1F*)infile[iset]->Get("hxection");

		c0->cd(iset+1);
		SetPadStyle();
		htmp = (TH1F*)hxection[iset];
		SetHistoStyle();
		hxection[iset]->Draw();

		Xection_hardqcd[iset] = hxection[iset]->GetMean();
		nevent[iset] = hxection[iset]->Integral();
	}

	//return;
	TH3F *hpT_eta_pi[nset];
	TH3F *hpT_eta_ka[nset];
	TH3F *hpT_eta_pr[nset];

	TH1D *heta_pi[nset];

	for (int iset=0; iset<nset; iset++){

		hpT_eta_pi[iset] = (TH3F*)infile[iset]->Get("hpT_eta_pi");
		hpT_eta_ka[iset] = (TH3F*)infile[iset]->Get("hpT_eta_ka");
		hpT_eta_pr[iset] = (TH3F*)infile[iset]->Get("hpT_eta_pr");

		hpT_eta_pi[iset]->Add(hpT_eta_ka[iset]);
		hpT_eta_pi[iset]->Add(hpT_eta_pr[iset]);

		heta_pi[iset] = (TH1D*)hpT_eta_pi[iset]->ProjectionY(Form("heta_pi_set%d",iset));
		heta_pi[iset]->Scale(1./(nevent[iset]*heta_pi[iset]->GetBinWidth(1)));
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.2*2*400,400);
	c1->Divide(2,1);

	c1->cd(1);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-3,0,3,5);
	SetHistoStyle("#eta","dN/d#eta");

	for (int iset=0; iset<nset; iset++){
		heta_pi[iset]->SetLineWidth(2);
		heta_pi[iset]->SetLineColor(nColor[iset]);
		heta_pi[iset]->Draw("same");
	}

	ifstream fdata;
	fdata.open("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/15.run15/18.dNdch/phobos_pp.txt");

	int count = 0;
	float ph_eta[100], ph_Y[100], ph_err1[100], ph_err2[100];

	while ( fdata >> ph_eta[count] >> ph_Y[count] >> ph_err1[count] >> ph_err2[count] ){
		count++;
	}

	TGraphAsymmErrors *gpp = new TGraphAsymmErrors(count,ph_eta,ph_Y,0,0,ph_err1,ph_err2);
	gpp->SetMarkerStyle(24);
	gpp->SetMarkerColor(1);
	gpp->SetLineColor(1);
	gpp->SetLineWidth(2);
	gpp->Draw("p");

	TLegend *leg = new TLegend(0.3,0.25,0.8,0.35);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(gpp,"PHOBOS, PRC 83, 024913 (2013)","P");
	leg->Draw();

	c1->cd(2);
	TLegend *leg = new TLegend(0.0,0.2,0.9,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("","PYTHIA v8.235 pp 200 GeV","");
	leg->AddEntry("","SoftQCD:inelastic","");
	leg->AddEntry(heta_pi[0],"case0: Monash, PDF:pSet=7, MultipartonInteractions:Kfactor=0.5","L");
	leg->AddEntry(heta_pi[1],"case1: Mode0, PDF:pSet=7, MultipartonInteractions:Kfactor=0.5","L");
	leg->AddEntry(heta_pi[2],"case2: Mode2, PDF:pSet=7, MultipartonInteractions:Kfactor=0.5","L");
	leg->AddEntry(heta_pi[3],"case3: Mode2, PDF:pSet=7, MultipartonInteractions:Kfactor=1.0","L");
	leg->AddEntry(heta_pi[4],"case4: Mode2, PDF:pSet=13, MultipartonInteractions:Kfactor=1.0","L");
	leg->Draw();

	//return;

	TH2F *hpT_y_Lc[nset];
	TH2F *hpT_y_D0[nset];

	TH1F *hpT_Lc[nset];
	TH1F *hpT_D0[nset];

	TH1F *hpT_R_Lc_D0[nset];
	TH1F *hpT_R_D0[nset];

	for (int iset=0; iset<nset; iset++){
		hpT_y_Lc[iset] = (TH2F*)infile[iset]->Get("hpT_y_Lc");
		hpT_y_D0[iset] = (TH2F*)infile[iset]->Get("hpT_y_D0");

		int etamin = hpT_y_Lc[iset]->GetYaxis()->FindBin(-etarange+0.0001);
		int etamax = hpT_y_D0[iset]->GetYaxis()->FindBin(+etarange-0.0001);

		hpT_Lc[iset] = (TH1F*)hpT_y_Lc[iset]->ProjectionX(Form("hpT_Lc_set%d",iset),etamin,etamax);
		hpT_D0[iset] = (TH1F*)hpT_y_D0[iset]->ProjectionX(Form("hpT_D0_set%d",iset),etamin,etamax);

		hpT_Lc[iset]->Rebin(2);
		hpT_D0[iset]->Rebin(2);
		hpT_Lc[iset]->Sumw2();
		hpT_D0[iset]->Sumw2();

		hpT_R_Lc_D0[iset] = (TH1F*)hpT_Lc[iset]->Clone(Form("hpT_R_Lc_D0_set%d",iset));
		hpT_R_Lc_D0[iset]->Divide(hpT_D0[iset]);

		hpT_R_D0[iset] = (TH1F*)hpT_D0[iset]->Clone(Form("hpT_R_D0_set%d",iset));
		hpT_R_D0[iset]->Divide(hpT_D0[0]);


	}//iset

	for (int iset=0; iset<nset; iset++){

		for (int ipt=0; ipt<hpT_D0[iset]->GetNbinsX(); ipt++){
			float dpt = hpT_D0[iset]->GetBinWidth(ipt+1);
			float pt = hpT_D0[iset]->GetBinCenter(ipt+1);

			float aa = 1./(2*TMath::Pi()*pt*dpt*2*etarange);

			float yy = hpT_D0[iset]->GetBinContent(ipt+1);
			float yy_err = hpT_D0[iset]->GetBinError(ipt+1);

			hpT_D0[iset]->SetBinContent(ipt+1, yy*aa/nevent[iset]*Xection_hardqcd[iset]/2.0);
			hpT_D0[iset]->SetBinError(ipt+1, yy_err*aa/nevent[iset]*Xection_hardqcd[iset]/2.0);
		}

	}

	TCanvas *c00 = new TCanvas("c00","c00",1.2*2*400,2*400);
	c00->Divide(2,2);

	c00->cd(1);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,10,1.2);
	SetHistoStyle("p_{T} (GeV/c)","#Lambda_{c}/D^{0}");

	for (int iset=0; iset<nset; iset++){
		hpT_R_Lc_D0[iset]->SetMarkerStyle(0);
		hpT_R_Lc_D0[iset]->SetMarkerColor(nColor[iset]);
		hpT_R_Lc_D0[iset]->SetLineColor(nColor[iset]);
		hpT_R_Lc_D0[iset]->SetLineWidth(2);
		hpT_R_Lc_D0[iset]->Draw("same");
	}

	c00->cd(3);
	SetPadStyle();
	gPad->SetLogy();
	htmp = (TH1D*)gPad->DrawFrame(0,5e-9,10,5e-2);
	SetHistoStyle("p_{T} (GeV/c)","Inv. Xection (mb)");
	for (int iset=0; iset<nset; iset++){
		hpT_D0[iset]->SetMarkerStyle(0);
		hpT_D0[iset]->SetMarkerColor(nColor[iset]);
		hpT_D0[iset]->SetLineColor(nColor[iset]);
		hpT_D0[iset]->SetLineWidth(2);
		hpT_D0[iset]->Draw("same");
	}

	c00->cd(4);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,10,2);
	SetHistoStyle("p_{T} (GeV/c)","Ratio to case0");
	for (int iset=0; iset<nset; iset++){
		hpT_R_D0[iset]->SetMarkerStyle(0);
		hpT_R_D0[iset]->SetMarkerColor(nColor[iset]);
		hpT_R_D0[iset]->SetLineColor(nColor[iset]);
		hpT_R_D0[iset]->SetLineWidth(2);
		hpT_R_D0[iset]->Draw("same");
	}

	float star_d0_xx[6] = {0.908, 1.57, 2.45, 3.44, 4.45, 5.45};
	float star_d0_yy[6] = {1.52e-2, 5.38e-3, 1.01e-3, 1.75e-4, 2.92e-5, 9.51e-6};
	float star_d0_yy_err[6] = {5.11e-3, 1.68e-3, 2.77e-4, 6.19e-5, 1.16e-5, 3.14e-6};

	for (int ipt=0; ipt<6; ipt++){
		star_d0_yy[ipt] *= 0.565;
		star_d0_yy_err[ipt] *= 0.565;
	}

	c00->cd(3);
	TGraphErrors *gstar_d0 = new TGraphErrors(6, star_d0_xx, star_d0_yy, 0, star_d0_yy_err);
	gstar_d0->SetMarkerStyle(25);
	gstar_d0->SetMarkerColor(1);
	gstar_d0->SetLineColor(1);
	gstar_d0->SetLineWidth(2);
	gstar_d0->Draw("p");

	TLegend *leg = new TLegend(0.3,0.8,0.8,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(gstar_d0,"STAR, D^{0} and D*, PRD 86, 072013 (2012)","P");
	leg->Draw();

	c00->cd(2);
	TLegend *leg = new TLegend(0.0,0.2,0.9,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("","PYTHIA v8.235 pp 200 GeV","");
	leg->AddEntry("","SoftQCD:inelastic","");
	leg->AddEntry(hpT_R_D0[0],"case0: Monash, PDF:pSet=7, MultipartonInteractions:Kfactor=0.5","L");
	leg->AddEntry(hpT_R_D0[1],"case1: Mode0, PDF:pSet=7, MultipartonInteractions:Kfactor=0.5","L");
	leg->AddEntry(hpT_R_D0[2],"case2: Mode2, PDF:pSet=7, MultipartonInteractions:Kfactor=0.5","L");
	leg->AddEntry(hpT_R_D0[3],"case3: Mode2, PDF:pSet=7, MultipartonInteractions:Kfactor=1.0","L");
	leg->AddEntry(hpT_R_D0[4],"case4: Mode2, PDF:pSet=13, MultipartonInteractions:Kfactor=1.0","L");
	leg->Draw();


	//leg->Draw();

	TFile *infile0 = new TFile("ppg223_refoldelectrons.root","read");
	TGraphErrors *gdata_cc = (TGraphErrors*)infile0->Get("cross_charm");
	TGraphErrors *gdata_bb = (TGraphErrors*)infile0->Get("cross_bottom");
	TGraphErrors *gdata_cb = (TGraphErrors*)infile0->Get("cross_sum");

	gdata_cb->SetMarkerStyle(20);
	gdata_cb->SetMarkerColor(1);
	gdata_cb->SetLineColor(1);

	gdata_cc->SetMarkerStyle(24);
	gdata_cc->SetMarkerColor(1);
	gdata_cc->SetLineColor(1);
	gdata_cc->SetLineWidth(2);

	gdata_bb->SetMarkerStyle(25);
	gdata_bb->SetMarkerColor(1);
	gdata_bb->SetLineColor(1);
	gdata_bb->SetLineWidth(2);

	TH3D *hpT_eta_cc[nset];
	TH3D *hpT_eta_bb[nset];
	TH3D *hpT_eta_cb[nset];

	TH1D *hpT_cc[nset];
	TH1D *hpT_bb[nset];
	TH1D *hpT_cb[nset];

	TH1D *hpT_R_cc[nset];
	TH1D *hpT_R_bb[nset];

	for (int iset=0; iset<nset; iset++){
		hpT_eta_cc[iset] = (TH3D*)infile[iset]->Get("hpT_eta_cc");
		hpT_eta_bb[iset] = (TH3D*)infile[iset]->Get("hpT_eta_bb");

		int etamin = (TH1D*)hpT_eta_cc[iset]->GetYaxis()->FindBin(-0.5+0.0001);
		int etamax = (TH1D*)hpT_eta_cc[iset]->GetYaxis()->FindBin(+0.5-0.0001);

		hpT_cc[iset] = (TH1D*)hpT_eta_cc[iset]->ProjectionX(Form("hpT_cc_set%d",iset),etamin,etamax,1,2);
		hpT_bb[iset] = (TH1D*)hpT_eta_bb[iset]->ProjectionX(Form("hpT_bb_set%d",iset),etamin,etamax,1,2);

		hpT_cc[iset]->Rebin(2);
		hpT_cc[iset]->Sumw2();

		hpT_bb[iset]->Rebin(2);
		hpT_bb[iset]->Sumw2();

		for (int ipt=0; ipt<hpT_cc[iset]->GetNbinsX(); ipt++){
			float dpt = hpT_cc[iset]->GetBinWidth(ipt+1);
			float pt = hpT_cc[iset]->GetBinCenter(ipt+1);

			float aa = 1./(2*TMath::Pi()*pt*dpt*2*0.5);

			float yy = hpT_cc[iset]->GetBinContent(ipt+1);
			float yy_err = hpT_cc[iset]->GetBinError(ipt+1);

			hpT_cc[iset]->SetBinContent(ipt+1, yy*aa/nevent[iset]*Xection_hardqcd[iset]/4.0);
			hpT_cc[iset]->SetBinError(ipt+1, yy_err*aa/nevent[iset]*Xection_hardqcd[iset]/4.0);
		}

		for (int ipt=0; ipt<hpT_bb[iset]->GetNbinsX(); ipt++){
			float dpt = hpT_bb[iset]->GetBinWidth(ipt+1);
			float pt = hpT_bb[iset]->GetBinCenter(ipt+1);

			float aa = 1./(2*TMath::Pi()*pt*dpt*2*0.5);

			float yy = hpT_bb[iset]->GetBinContent(ipt+1);
			float yy_err = hpT_bb[iset]->GetBinError(ipt+1);

			hpT_bb[iset]->SetBinContent(ipt+1, yy*aa/nevent[iset]*Xection_hardqcd[iset]/4.0);
			hpT_bb[iset]->SetBinError(ipt+1, yy_err*aa/nevent[iset]*Xection_hardqcd[iset]/4.0);
		}

		hpT_R_cc[iset] = (TH1D*)hpT_cc[iset]->Clone(Form("hpT_R_cc_set%d",iset));
		hpT_R_bb[iset] = (TH1D*)hpT_bb[iset]->Clone(Form("hpT_R_bb_set%d",iset));

		hpT_R_cc[iset]->Divide(hpT_cc[0]);
		hpT_R_bb[iset]->Divide(hpT_bb[0]);
	}

	TCanvas *c11 = new TCanvas("c11","c11",1.2*2*400,2*400);
	c11->Divide(2,2);

	c11->cd(1);
	SetPadStyle();
	gPad->SetLogy();
	htmp = (TH1D*)gPad->DrawFrame(0,2e-10,8,1e-2);
	SetHistoStyle("p_{T} (GeV/c)","Inv. Xection (mb)");

	for (int iset=0; iset<nset; iset++){
		hpT_cc[iset]->SetMarkerStyle(0);
		hpT_cc[iset]->SetMarkerColor(nColor[iset]);
		hpT_cc[iset]->SetLineColor(nColor[iset]);
		hpT_cc[iset]->SetLineWidth(2);
		hpT_cc[iset]->Draw("same");
	}
	gdata_cc->Draw("p");

	TLegend *leg = new TLegend(0.3,0.8,0.8,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(gdata_cc,"PHENIX, c#rightarrowe, arXiv:1901.08405","P");
	leg->Draw();

	c11->cd(2);
	SetPadStyle();
	gPad->SetLogy();
	htmp = (TH1D*)gPad->DrawFrame(0,1e-10,10,1e-1);
	SetHistoStyle("p_{T} (GeV/c)","Inv. Xection (mb)");

	for (int iset=0; iset<nset; iset++){
		hpT_bb[iset]->SetMarkerStyle(0);
		hpT_bb[iset]->SetMarkerColor(nColor[iset]);
		hpT_bb[iset]->SetLineColor(nColor[iset]);
		hpT_bb[iset]->SetLineWidth(2);
		hpT_bb[iset]->Draw("same");
	}
	gdata_bb->Draw("p");

	TLegend *leg = new TLegend(0.3,0.8,0.8,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(gdata_bb,"PHENIX, b#rightarrowe, arXiv:1901.08405","P");
	leg->Draw();

	c11->cd(3);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,8,2);
	SetHistoStyle("p_{T} (GeV/c)","Ratio to case0");
	for (int iset=0; iset<nset; iset++){
		hpT_R_cc[iset]->SetMarkerStyle(0);
		hpT_R_cc[iset]->SetMarkerColor(nColor[iset]);
		hpT_R_cc[iset]->SetLineColor(nColor[iset]);
		hpT_R_cc[iset]->SetLineWidth(2);
		hpT_R_cc[iset]->Draw("same");
	}

	c11->cd(4);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,8,2);
	SetHistoStyle("p_{T} (GeV/c)","Ratio to case0");
	for (int iset=0; iset<nset; iset++){
		hpT_R_bb[iset]->SetMarkerStyle(0);
		hpT_R_bb[iset]->SetMarkerColor(nColor[iset]);
		hpT_R_bb[iset]->SetLineColor(nColor[iset]);
		hpT_R_bb[iset]->SetLineWidth(2);
		hpT_R_bb[iset]->Draw("same");
	}

	c11->cd(2);
	TLegend *leg = new TLegend(0.3,0.55,0.95,0.8);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(hpT_R_D0[0],"case0: Monash, PDF:pSet=7, MPI:Kfactor=0.5","L");
	leg->AddEntry(hpT_R_D0[1],"case1: Mode0, PDF:pSet=7, MPI:Kfactor=0.5","L");
	leg->AddEntry(hpT_R_D0[2],"case2: Mode2, PDF:pSet=7, MPI:Kfactor=0.5","L");
	leg->AddEntry(hpT_R_D0[3],"case3: Mode2, PDF:pSet=7, MPI:Kfactor=1.0","L");
	leg->AddEntry(hpT_R_D0[4],"case4: Mode2, PDF:pSet=13, MPI:Kfactor=1.0","L");
	leg->Draw();

}
