#include "/phenix/u/shlim/Style.h"

void Draw_PPG_pAu(){

	const int narm = 2;
	const int ncent = 6;

	int cent_array[ncent+1] = {0, 5, 10, 20, 40, 60, 84}; 
	float Ncoll_err[ncent] = {0.0578, 0.0655, 0.0652, 0.0662, 0.0658, 0.0592};
	float BiasF_err[ncent] = {0.0093, 0.0067, 0.0043, 0.0051, 0.0127, 0.0675};
	const int nColor[ncent] = {1, 4, kGreen+2, 2, 6, 9};

	TLine *line = new TLine(0.5,1,10.5,1);
	line->SetLineStyle(2);
	line->SetLineWidth(2);

	TLine *line_eta = new TLine(-3.2,1,3.2,1);
	line_eta->SetLineStyle(2);
	line_eta->SetLineWidth(2);

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

		/*
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
		*/
	}


}
