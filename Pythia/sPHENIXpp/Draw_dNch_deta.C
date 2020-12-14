#include "Style.h"

void Draw_dNch_deta(){

	const int nset = 4;

	const int nMarker[5] = {24, 25, 26, 20, 21};
	const int nColor[5] = {1, 2, 4, 6, 8};

	TFile *infile[nset];

	infile[0] = new TFile("CTEQ_PDF_Kfactor1.0.root","read");
	infile[1] = new TFile("CTEQ_PDF_Kfactor0.5.root","read");
	infile[2] = new TFile("DefaultPDF_Kfactor1.0.root","read");
	infile[3] = new TFile("DefaultPDF_Kfactor0.8.root","read");

	TH1D *hdNch[nset];
	TH1D *hpT_pi[nset][2];
	TH1D *hpT_ka[nset][2];

	for (int iset=0; iset<nset; iset++){
		hdNch[iset] = (TH1D*)infile[iset]->Get("h_inc_eta");
		hdNch[iset]->SetLineColor(nColor[iset+1]);
		hdNch[iset]->SetMarkerColor(nColor[iset+1]);
		hdNch[iset]->SetLineWidth(2);

		for (int ichg=0; ichg<2; ichg++){
			hpT_pi[iset][ichg] = (TH1D*)infile[iset]->Get(Form("h_pi_pt_chg%d",ichg));
			hpT_ka[iset][ichg] = (TH1D*)infile[iset]->Get(Form("h_ka_pt_chg%d",ichg));
		}

		hpT_pi[iset][0]->Add(hpT_pi[iset][1]);
		hpT_pi[iset][0]->Scale(0.5);
		hpT_pi[iset][0]->SetLineColor(nColor[iset+1]);
		hpT_pi[iset][0]->SetMarkerColor(nColor[iset+1]);
		hpT_pi[iset][0]->SetMarkerStyle(20);
		hpT_pi[iset][0]->SetMarkerSize(0.7);
		hpT_pi[iset][0]->SetLineWidth(2);

		hpT_ka[iset][0]->Add(hpT_ka[iset][1]);
		hpT_ka[iset][0]->Scale(0.5);
		hpT_ka[iset][0]->SetLineColor(nColor[iset+1]);
		hpT_ka[iset][0]->SetMarkerColor(nColor[iset+1]);
		hpT_ka[iset][0]->SetMarkerStyle(20);
		hpT_ka[iset][0]->SetMarkerSize(0.7);
		hpT_ka[iset][0]->SetLineWidth(2);
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*500,500);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(-3,0,3,5);
	SetHistoStyle("#eta","dN/d#eta","",24,20);

	ifstream fdata;
	fdata.open("phobos_pp.txt");

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

	hdNch[0]->Draw("same");
	hdNch[1]->Draw("same");
	hdNch[2]->Draw("same");
	hdNch[3]->Draw("same");

	TLegend *leg = new TLegend(0.15,0.2,0.65,0.4);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035);
	leg->AddEntry("","pp 200 GeV","");
	leg->AddEntry(gpp,"PHOBOS","P");
	leg->AddEntry(hdNch[0],"PYTHIA8, CTEQ6L NLO, MPI Kfactor=1.0","P");
	leg->AddEntry(hdNch[1],"PYTHIA8, CTEQ6L NLO, MPI Kfactor=0.5","P");
	leg->AddEntry(hdNch[2],"PYTHIA8, NNPDF2.3 QCD+QED LO, MPI Kfactor=1.0","P");
	leg->AddEntry(hdNch[3],"PYTHIA8, NNPDF2.3 QCD+QED LO, MPI Kfactor=0.8","P");
	leg->Draw();

	TFile *infile1 = new TFile("identified_chd_pT.root","read");
	TGraphErrors *gpion[2];
	TGraphErrors *gkaon[2];

	for (int ichg=0; ichg<2; ichg++){
		gpion[ichg] = (TGraphErrors*)infile1->Get(Form("gpion_y0_chg%d",ichg));
		gkaon[ichg] = (TGraphErrors*)infile1->Get(Form("gkaon_y0_chg%d",ichg));
	}

	TGraphErrors *gpion_combined = new TGraphErrors;
	TGraphErrors *gpion_ratio[nset];
	TGraphErrors *gkaon_combined = new TGraphErrors;
	TGraphErrors *gkaon_ratio[nset];
	for (int iset=0; iset<nset; iset++){
		gpion_ratio[iset] = new TGraphErrors;
		gkaon_ratio[iset] = new TGraphErrors;
	}

	for (int ii=0; ii<gpion[0]->GetN(); ii++){
		double xx, yy0, yy1;
		gpion[0]->GetPoint(ii, xx, yy0);
		gpion[1]->GetPoint(ii, xx, yy1);
		float yy0_err = gpion[0]->GetErrorY(ii)/yy0;
		float yy1_err = gpion[1]->GetErrorY(ii)/yy1;

		float yy = (yy0 + yy1)/2;
		float yy_err = yy*sqrt(0.5*0.5*yy0_err*yy0_err + 0.5*0.5*yy1_err*yy1_err);
		gpion_combined->SetPoint(ii, xx, yy);
		gpion_combined->SetPointError(ii, 0, yy_err);

		for (int iset=0; iset<nset; iset++){
			int bin = hpT_pi[iset][0]->FindBin(xx);
			float yy2 = hpT_pi[iset][0]->GetBinContent(bin); 
			float yy2_err = hpT_pi[iset][0]->GetBinError(bin); 

			float ratio = yy2 / yy;
			float ratio_err = ratio*sqrt(pow(yy2_err/yy,2) + pow(yy_err/yy,2));
			gpion_ratio[iset]->SetPoint(ii, xx, ratio);
			gpion_ratio[iset]->SetPointError(ii, 0, ratio_err);

		}//iset
	}

	for (int ii=0; ii<gkaon[0]->GetN(); ii++){
		double xx, yy0, yy1;
		gkaon[0]->GetPoint(ii, xx, yy0);
		gkaon[1]->GetPoint(ii, xx, yy1);
		float yy0_err = gkaon[0]->GetErrorY(ii)/yy0;
		float yy1_err = gkaon[1]->GetErrorY(ii)/yy1;

		float yy = (yy0 + yy1)/2;
		float yy_err = yy*sqrt(0.5*0.5*yy0_err*yy0_err + 0.5*0.5*yy1_err*yy1_err);
		gkaon_combined->SetPoint(ii, xx, yy);
		gkaon_combined->SetPointError(ii, 0, yy_err);

		for (int iset=0; iset<nset; iset++){
			int bin = hpT_ka[iset][0]->FindBin(xx);
			float yy2 = hpT_ka[iset][0]->GetBinContent(bin); 
			float yy2_err = hpT_ka[iset][0]->GetBinError(bin); 

			float ratio = yy2 / yy;
			float ratio_err = ratio*sqrt(pow(yy2_err/yy,2) + pow(yy_err/yy,2));
			gkaon_ratio[iset]->SetPoint(ii, xx, ratio);
			gkaon_ratio[iset]->SetPointError(ii, 0, ratio_err);

		}//iset
	}

	TCanvas *c2 = new TCanvas("c2","c2",1.1*2*500,500);
	c2->Divide(2,1);

	{
		c2->cd(1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1D*)gPad->DrawFrame(0,1e-8,5,10);
		SetHistoStyle("p_{T} (GeV/c)","#frac{1}{2#pi p_{T}} #frac{d^{2}N}{dp_{T}d#eta} (GeV/c)^{-2}","",22,20);

		gpion_combined->SetMarkerStyle(24);
		gpion_combined->SetMarkerColor(1);
		gpion_combined->SetLineWidth(2);
		gpion_combined->Draw("p");
		for (int iset=0; iset<nset; iset++){
			hpT_pi[iset][0]->Draw("p same");
		}

		TLegend *leg = new TLegend(0.15,0.18,0.55,0.38);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->AddEntry("","(#pi^{+}+#pi^{-})/2, pp 200 GeV","");
		leg->AddEntry(gpp,"PHENIX+STAR","P");
		leg->AddEntry(hdNch[0],"PYTHIA8, CTEQ6L NLO, MPI Kfactor=1.0","P");
		leg->AddEntry(hdNch[1],"PYTHIA8, CTEQ6L NLO, MPI Kfactor=0.5","P");
		leg->AddEntry(hdNch[2],"PYTHIA8, NNPDF2.3 QCD+QED LO, MPI Kfactor=1.0","P");
		leg->AddEntry(hdNch[3],"PYTHIA8, NNPDF2.3 QCD+QED LO, MPI Kfactor=0.8","P");
		leg->Draw();

		c2->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,5,5);
		SetHistoStyle("p_{T} (GeV/c)","PYTHIA / DATA","",22,20);

		for (int iset=0; iset<nset; iset++){
			gpion_ratio[iset]->SetMarkerStyle(20);
			gpion_ratio[iset]->SetMarkerColor(nColor[iset+1]);
			gpion_ratio[iset]->SetLineColor(nColor[iset+1]);
			gpion_ratio[iset]->SetLineWidth(1);
			gpion_ratio[iset]->Draw("p");
		}

	}

	TCanvas *c3 = new TCanvas("c3","c3",1.1*2*500,500);
	c3->Divide(2,1);

	{
		c3->cd(1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1D*)gPad->DrawFrame(0,1e-9,5,1);
		SetHistoStyle("p_{T} (GeV/c)","#frac{1}{2#pi p_{T}} #frac{d^{2}N}{dp_{T}d#eta} (GeV/c)^{-2}","",22,20);

		gkaon_combined->SetMarkerStyle(24);
		gkaon_combined->SetMarkerColor(1);
		gkaon_combined->SetLineWidth(2);
		gkaon_combined->Draw("p");
		for (int iset=0; iset<nset; iset++){
			hpT_ka[iset][0]->Draw("p same");
		}

		TLegend *leg = new TLegend(0.15,0.18,0.55,0.38);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->AddEntry("","(#K^{+}+#K^{-})/2, pp 200 GeV","");
		leg->AddEntry(gpp,"PHENIX+STAR","P");
		leg->AddEntry(hdNch[0],"PYTHIA8, CTEQ6L NLO, MPI Kfactor=1.0","P");
		leg->AddEntry(hdNch[1],"PYTHIA8, CTEQ6L NLO, MPI Kfactor=0.5","P");
		leg->AddEntry(hdNch[2],"PYTHIA8, NNPDF2.3 QCD+QED LO, MPI Kfactor=1.0","P");
		leg->AddEntry(hdNch[3],"PYTHIA8, NNPDF2.3 QCD+QED LO, MPI Kfactor=0.8","P");
		leg->Draw();

		c3->cd(2);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(0,0,5,5);
		SetHistoStyle("p_{T} (GeV/c)","PYTHIA / DATA","",22,20);

		for (int iset=0; iset<nset; iset++){
			gkaon_ratio[iset]->SetMarkerStyle(20);
			gkaon_ratio[iset]->SetMarkerColor(nColor[iset+1]);
			gkaon_ratio[iset]->SetLineColor(nColor[iset+1]);
			gkaon_ratio[iset]->SetLineWidth(1);
			gkaon_ratio[iset]->Draw("p");
		}

	}


}
