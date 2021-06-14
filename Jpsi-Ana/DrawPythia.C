void DrawPythia(){

	gStyle->SetOptStat(0);

	TFile *infile = new TFile("outfile_pythia_jpsi.root","read");
	TFile *infile_mb = new TFile("outfile_pythia_mb.root","read");

	const int narm = 2;

	TH2D *h2d_fvtx_bbc[narm];
	TH2D *h2d_fvtx_bbc2[narm];
	TH2D *h2d_fvtx_svx[narm];
	TH2D *h2d_fvtx_fvtx[narm];

	TProfile *hprof_fvtx_bbc[narm];
	TProfile *hprof_fvtx_svx[narm];
	TProfile *hprof_fvtx_fvtx[narm];

	TH2D *h2d_mb_fvtx_bbc[narm];
	TH2D *h2d_mb_fvtx_bbc2[narm];
	TH2D *h2d_mb_fvtx_svx[narm];
	TH2D *h2d_mb_fvtx_fvtx[narm];

	TProfile *hprof_mb_fvtx_bbc[narm];
	TProfile *hprof_mb_fvtx_svx[narm];

	for (int iarm=0; iarm<narm; iarm++){

		h2d_fvtx_bbc[iarm] = (TH2D*)infile->Get(Form("h2d_fvtx_bbc_arm%d",iarm));
		h2d_fvtx_bbc2[iarm] = (TH2D*)infile->Get(Form("h2d_fvtx_bbc2_arm%d",iarm));
		h2d_fvtx_svx[iarm] = (TH2D*)infile->Get(Form("h2d_fvtx_svx_arm%d",iarm));
		h2d_fvtx_fvtx[iarm] = (TH2D*)infile->Get(Form("h2d_fvtx_fvtx_arm%d",iarm));

		hprof_fvtx_bbc[iarm] = (TProfile*)h2d_fvtx_bbc[iarm]->ProfileX();
		hprof_fvtx_svx[iarm] = (TProfile*)h2d_fvtx_svx[iarm]->ProfileX();
		hprof_fvtx_fvtx[iarm] = (TProfile*)h2d_fvtx_fvtx[iarm]->ProfileX();

		hprof_fvtx_bbc[iarm]->SetMarkerStyle(20);
		hprof_fvtx_bbc[iarm]->SetMarkerColor(1);
		hprof_fvtx_bbc[iarm]->SetLineColor(1);

		hprof_fvtx_svx[iarm]->SetMarkerStyle(20);
		hprof_fvtx_svx[iarm]->SetMarkerColor(1);
		hprof_fvtx_svx[iarm]->SetLineColor(1);

		hprof_fvtx_fvtx[iarm]->SetMarkerStyle(20);
		hprof_fvtx_fvtx[iarm]->SetMarkerColor(1);
		hprof_fvtx_fvtx[iarm]->SetLineColor(1);

		h2d_mb_fvtx_bbc[iarm] = (TH2D*)infile_mb->Get(Form("h2d_fvtx_bbc_arm%d",iarm));
		h2d_mb_fvtx_bbc2[iarm] = (TH2D*)infile_mb->Get(Form("h2d_fvtx_bbc2_arm%d",iarm));
		h2d_mb_fvtx_svx[iarm] = (TH2D*)infile_mb->Get(Form("h2d_fvtx_svx_arm%d",iarm));
		h2d_mb_fvtx_fvtx[iarm] = (TH2D*)infile_mb->Get(Form("h2d_fvtx_fvtx_arm%d",iarm));

		hprof_mb_fvtx_bbc[iarm] = (TProfile*)h2d_mb_fvtx_bbc[iarm]->ProfileX();
		hprof_mb_fvtx_svx[iarm] = (TProfile*)h2d_mb_fvtx_svx[iarm]->ProfileX();

		hprof_mb_fvtx_bbc[iarm]->SetMarkerStyle(20);
		hprof_mb_fvtx_bbc[iarm]->SetMarkerColor(1);
		hprof_mb_fvtx_bbc[iarm]->SetLineColor(1);

		hprof_mb_fvtx_svx[iarm]->SetMarkerStyle(20);
		hprof_mb_fvtx_svx[iarm]->SetMarkerColor(1);
		hprof_mb_fvtx_svx[iarm]->SetLineColor(1);

	}//iarm

	TCanvas *c1[narm];
	TCanvas *c2[narm];

	for (int iarm=1; iarm<narm; iarm++){

		c1[iarm] = new TCanvas(Form("c1_arm%d",iarm),Form("c1_arm%d",iarm),1.1*2*400,2*400);
		c1[iarm]->Divide(2,2);

		{
			c1[iarm]->cd(1);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,30);
			htmp->GetYaxis()->SetTitle("N_{ch}^{|#eta|<1}");
			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);

			if ( iarm==0 ){
				htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}else{
				htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_fvtx_svx[iarm]->Draw("colz same");
			//hprof_fvtx_svx[iarm]->Draw("p same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			if ( iarm==0 ){
				leg->AddEntry("","-2.2<y^{J/#psi}<-1.2","h");
			}else{
				leg->AddEntry("","1.2<y^{J/#psi}<2.2","h");
			}
			leg->Draw();
		}

		{
			c1[iarm]->cd(2);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,20);
			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);

			if ( iarm==0 ){
				htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
				htmp->GetYaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}else{
				htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
				htmp->GetYaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_fvtx_fvtx[iarm]->Draw("colz same");
			//hprof_fvtx_fvtx[iarm]->Draw("p same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			if ( iarm==0 ){
				leg->AddEntry("","-2.2<y^{J/#psi}<-1.2","h");
			}else{
				leg->AddEntry("","1.2<y^{J/#psi}<2.2","h");
			}
			leg->Draw();
		}

		{
			c1[iarm]->cd(3);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,10);
			if ( iarm==0 ){
					htmp->GetYaxis()->SetTitle("N_{ch}^{-3.9<#eta<-3.1}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}else{
					htmp->GetYaxis()->SetTitle("N_{ch}^{3.1<#eta<3.9}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}

			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_fvtx_bbc[iarm]->Draw("colz same");
			//hprof_fvtx_bbc[iarm]->Draw("p same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			if ( iarm==0 ){
				leg->AddEntry("","-2.2<y^{J/#psi}<-1.2","h");
			}else{
				leg->AddEntry("","1.2<y^{J/#psi}<2.2","h");
			}
			leg->Draw();
		}

		{
			c1[iarm]->cd(4);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,10);
			if ( iarm==0 ){
					htmp->GetYaxis()->SetTitle("N_{ch}^{3.1<#eta<3.9}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}else{
					htmp->GetYaxis()->SetTitle("N_{ch}^{-3.9<#eta<-3.1}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}

			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_fvtx_bbc2[iarm]->Draw("colz same");
			//hprof_fvtx_bbc[iarm]->Draw("p same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			if ( iarm==0 ){
				leg->AddEntry("","-2.2<y^{J/#psi}<-1.2","h");
			}else{
				leg->AddEntry("","1.2<y^{J/#psi}<2.2","h");
			}
			leg->Draw();
		}
	}

	//return;

	/*
	for (int iarm=0; iarm<narm; iarm++){

		c2[iarm] = new TCanvas(Form("c2_arm%d",iarm),Form("c2_arm%d",iarm),1.1*2*400,2*400);
		c2[iarm]->Divide(2,2);

		{
			c2[iarm]->cd(1);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,30);
			htmp->GetYaxis()->SetTitle("N_{ch}^{|#eta|<1}");
			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);

			if ( iarm==0 ){
				htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}else{
				htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_mb_fvtx_svx[iarm]->Draw("colz same");
			//hprof_mb_fvtx_svx[iarm]->Draw("p same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			leg->AddEntry("","SoftQCD:nonDiffractive","h");
			leg->Draw();
		}

		{
			c2[iarm]->cd(2);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,20);
			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);

			if ( iarm==0 ){
				htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
				htmp->GetYaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}else{
				htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
				htmp->GetYaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_mb_fvtx_fvtx[iarm]->Draw("colz same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			leg->AddEntry("","SoftQCD:nonDiffractive","h");
			leg->Draw();
		}

		{
			c2[iarm]->cd(3);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,10);
			if ( iarm==0 ){
					htmp->GetYaxis()->SetTitle("N_{ch}^{-3.9<#eta<-3.1}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}else{
					htmp->GetYaxis()->SetTitle("N_{ch}^{3.1<#eta<3.9}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}

			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_mb_fvtx_bbc[iarm]->Draw("colz same");
			//hprof_mb_fvtx_bbc[iarm]->Draw("p same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			leg->AddEntry("","SoftQCD:nonDiffractive","h");
			leg->Draw();
		}

		{
			c2[iarm]->cd(4);
			gPad->SetMargin(0.14,0.15,0.12,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,10);
			if ( iarm==0 ){
					htmp->GetYaxis()->SetTitle("N_{ch}^{3.1<#eta<3.9}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{-2.6<#eta<-1.2}");
			}else{
					htmp->GetYaxis()->SetTitle("N_{ch}^{-3.9<#eta<-3.1}");
					htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<#eta<2.6}");
			}

			htmp->GetYaxis()->SetTitleSize(0.05);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitleSize(0.05);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h2d_mb_fvtx_bbc2[iarm]->Draw("colz same");
			//hprof_fvtx_bbc[iarm]->Draw("p same");

			TLegend *leg = new TLegend(0.2,0.8,0.5,0.95);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.045);
			leg->SetBorderSize(0);
			leg->SetTextFont(62);
			leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
			leg->AddEntry("","SoftQCD:nonDiffractive","h");
			leg->Draw();
		}
	}
	*/

	TH1D *h1d_fvtx_same[2];
	TH1D *h1d_fvtx_oppo[2];
	TH1D *h1d_svx[2];
	TH1D *h1d_bbc_same[2];
	TH1D *h1d_bbc_oppo[2];

	for (int ii=0; ii<2; ii++){

		h1d_fvtx_same[ii] = (TH1D*)infile->Get(Form("h1d_fvtx_same_%d",ii)); 
		h1d_fvtx_oppo[ii] = (TH1D*)infile->Get(Form("h1d_fvtx_oppo_%d",ii)); 

		h1d_svx[ii] = (TH1D*)infile->Get(Form("h1d_svx_%d",ii)); 

		h1d_bbc_same[ii] = (TH1D*)infile->Get(Form("h1d_bbc_same_%d",ii)); 
		h1d_bbc_oppo[ii] = (TH1D*)infile->Get(Form("h1d_bbc_oppo_%d",ii)); 

		//h1d_fvtx_same[ii]->Rebin(2);
		//h1d_fvtx_oppo[ii]->Rebin(2);
		//h1d_svx[ii]->Rebin(2);
		//h1d_bbc_same[ii]->Rebin(2);
		//h1d_bbc_oppo[ii]->Rebin(2);

		h1d_fvtx_same[ii]->Sumw2();
		h1d_fvtx_oppo[ii]->Sumw2();
		h1d_svx[ii]->Sumw2();
		h1d_bbc_same[ii]->Sumw2();
		h1d_bbc_oppo[ii]->Sumw2();

	}

	TCanvas *c3 = new TCanvas("c3","c3",1.1*500,500);

	{
		gPad->SetMargin(0.18,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,25,0.1);
		htmp->GetYaxis()->SetTitle("B_{ll}^{#psi(2S)}N^{#psi(2S)}/B_{ll}^{J/#psi}N^{J/#psi}");
		htmp->GetYaxis()->SetTitleSize(0.06);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("N_{ch}");
		htmp->GetXaxis()->SetTitleSize(0.06);
		htmp->GetXaxis()->SetLabelSize(0.04);

		h1d_fvtx_same[1]->Divide(h1d_fvtx_same[0]);
		h1d_fvtx_oppo[1]->Divide(h1d_fvtx_oppo[0]);
		h1d_svx[1]->Divide(h1d_svx[0]);

		h1d_fvtx_same[1]->SetMarkerStyle(24);
		h1d_fvtx_same[1]->SetMarkerColor(1);
		h1d_fvtx_same[1]->SetLineColor(1);
		h1d_fvtx_same[1]->SetLineWidth(2);
		h1d_fvtx_same[1]->Draw("p same");

		h1d_fvtx_oppo[1]->SetMarkerStyle(24);
		h1d_fvtx_oppo[1]->SetMarkerColor(2);
		h1d_fvtx_oppo[1]->SetLineColor(2);
		h1d_fvtx_oppo[1]->SetLineWidth(2);
		h1d_fvtx_oppo[1]->Draw("p same");

		h1d_svx[1]->SetMarkerStyle(25);
		h1d_svx[1]->SetMarkerColor(4);
		h1d_svx[1]->SetLineColor(4);
		h1d_svx[1]->SetLineWidth(2);
		h1d_svx[1]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.7,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","1.2<y^{J/#psi}<2.2","h");
		leg->AddEntry(h1d_fvtx_same[1],"N_{ch} 1.2<#eta<2.6","p");
		leg->AddEntry(h1d_fvtx_oppo[1],"N_{ch} -2.6<#eta<-1.2","p");
		leg->AddEntry(h1d_svx[1],"N_{ch} |#eta|<1.0","p");
		leg->Draw();
	}
		
	TCanvas *c4 = new TCanvas("c4","c4",1.1*500,500);

	{
		gPad->SetMargin(0.18,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,15,0.1);
		htmp->GetYaxis()->SetTitle("B_{ll}^{#psi(2S)}N^{#psi(2S)}/B_{ll}^{J/#psi}N^{J/#psi}");
		htmp->GetYaxis()->SetTitleSize(0.06);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("N_{ch}");
		htmp->GetXaxis()->SetTitleSize(0.06);
		htmp->GetXaxis()->SetLabelSize(0.04);

		h1d_bbc_same[1]->Divide(h1d_bbc_same[0]);
		h1d_bbc_oppo[1]->Divide(h1d_bbc_oppo[0]);

		h1d_bbc_same[1]->SetMarkerStyle(26);
		h1d_bbc_same[1]->SetMarkerColor(kGreen+2);
		h1d_bbc_same[1]->SetLineColor(kGreen+2);
		h1d_bbc_same[1]->SetLineWidth(2);
		h1d_bbc_same[1]->Draw("p same");

		h1d_bbc_oppo[1]->SetMarkerStyle(26);
		h1d_bbc_oppo[1]->SetMarkerColor(6);
		h1d_bbc_oppo[1]->SetLineColor(6);
		h1d_bbc_oppo[1]->SetLineWidth(2);
		h1d_bbc_oppo[1]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.7,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","1.2<y^{J/#psi}<2.2","h");
		leg->AddEntry(h1d_bbc_same[1],"N_{ch} 3.1<#eta<3.9","p");
		leg->AddEntry(h1d_bbc_oppo[1],"N_{ch} -3.9<#eta<-3.1","p");
		leg->AddEntry("","","");
		leg->Draw();
	}

}
