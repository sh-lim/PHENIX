void DrawPythia(){

	gStyle->SetOptStat(0);

	TFile *infile = new TFile("outfile_pythia_jpsi.root","read");

	const int narm = 2;

	TH2D *h2d_fvtx_bbc[narm];
	TH2D *h2d_fvtx_svx[narm];

	TProfile *hprof_fvtx_bbc[narm];
	TProfile *hprof_fvtx_svx[narm];

	for (int iarm=0; iarm<narm; iarm++){

		h2d_fvtx_bbc[iarm] = (TH2D*)infile->Get(Form("h2d_fvtx_bbc_arm%d",iarm));
		h2d_fvtx_svx[iarm] = (TH2D*)infile->Get(Form("h2d_fvtx_svx_arm%d",iarm));

		hprof_fvtx_bbc[iarm] = (TProfile*)h2d_fvtx_bbc[iarm]->ProfileX();
		hprof_fvtx_svx[iarm] = (TProfile*)h2d_fvtx_svx[iarm]->ProfileX();

		hprof_fvtx_bbc[iarm]->SetMarkerStyle(20);
		hprof_fvtx_bbc[iarm]->SetMarkerColor(1);
		hprof_fvtx_bbc[iarm]->SetLineColor(1);

		hprof_fvtx_svx[iarm]->SetMarkerStyle(20);
		hprof_fvtx_svx[iarm]->SetMarkerColor(1);
		hprof_fvtx_svx[iarm]->SetLineColor(1);

	}//iarm

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*400,400);
	c1->Divide(2,1);

	{
		c1->cd(1);
		gPad->SetMargin(0.14,0.15,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,15,30);
		htmp->GetYaxis()->SetTitle("N_{ch}^{|#eta|<1}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		htmp->GetXaxis()->SetTitle("N_{ch}^{-2.2<#eta<-1.2}");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		h2d_fvtx_svx[0]->Draw("colz same");
		hprof_fvtx_svx[0]->Draw("p same");
		
		TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.045);
		leg->SetBorderSize(0);
		leg->SetTextFont(62);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","-2.2<y^{J/#psi}<-1.2","h");
		leg->Draw();
	}

	{
		c1->cd(2);
		gPad->SetMargin(0.14,0.15,0.12,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,15,10);
		htmp->GetYaxis()->SetTitle("N_{ch}^{-3.9<#eta<-3.1}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);

		htmp->GetXaxis()->SetTitle("N_{ch}^{-2.2<#eta<-1.2}");
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		h2d_fvtx_bbc[0]->Draw("colz same");
		hprof_fvtx_bbc[0]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.045);
		leg->SetBorderSize(0);
		leg->SetTextFont(62);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","-2.2<y^{J/#psi}<-1.2","h");
		leg->Draw();
	}


}
