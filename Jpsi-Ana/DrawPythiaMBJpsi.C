void DrawPythiaMBJpsi(){

	gStyle->SetOptStat(0);

	const int nset = 3;
	const int narm = 2;

	TFile *infile[nset];

	infile[0] = new TFile("outfile_pythiaMB_pp13TeV.root","read");
	infile[1] = new TFile("outfile_pythiaMB_pp200GeV.root","read");
	infile[2] = new TFile("outfile_pythiaMB_pp200GeV_Detroit.root","read");

	const float xmax[nset] = {60, 25, 25};

	TH1D *h1d_svx[nset][2];
	TH1D *h1d_scaled_svx[nset][2];

	TH1D *h1d_fvtx[nset][2];
	TH1D *h1d_scaled_fvtx[nset][2];

	TH1D *hmult_svx[nset];
	TH1D *hmult_fvtx[nset];

	TH1D *h1d_svx_ratio[nset];
	TH1D *h1d_fvtx_ratio[nset];

	for (int iset=0; iset<nset; iset++){
		for (int ii=0; ii<2; ii++){

			h1d_svx[iset][ii] = (TH1D*)infile[iset]->Get(Form("h1d_svx_%d",ii));
			h1d_scaled_svx[iset][ii] = (TH1D*)infile[iset]->Get(Form("h1d_scaled_svx_%d",ii));

			h1d_fvtx[iset][ii] = (TH1D*)infile[iset]->Get(Form("h1d_fvtx_%d",ii));
			h1d_scaled_fvtx[iset][ii] = (TH1D*)infile[iset]->Get(Form("h1d_scaled_fvtx_%d",ii));

			h1d_svx[iset][ii]->Sumw2();
			h1d_scaled_svx[iset][ii]->Sumw2();

			h1d_fvtx[iset][ii]->Sumw2();
			h1d_scaled_fvtx[iset][ii]->Sumw2();
		}

		hmult_svx[iset] = (TH1D*)infile[iset]->Get("hmult_svx");
		hmult_fvtx[iset] = (TH1D*)infile[iset]->Get("hmult_fvtx");

		float mean_ratio = h1d_svx[iset][1]->Integral()/h1d_svx[iset][0]->Integral();

		h1d_scaled_svx[iset][1]->Divide(h1d_scaled_svx[iset][0]);
		h1d_scaled_svx[iset][1]->Scale(1./mean_ratio);

		h1d_svx_ratio[iset] = (TH1D*)h1d_svx[iset][1]->Clone(Form("h1d_svx_ratio_%d",iset));
		h1d_svx_ratio[iset]->Divide(h1d_svx[iset][0]);
		h1d_svx_ratio[iset]->Scale(1./mean_ratio);

		h1d_svx[iset][0]->Scale(1./h1d_svx[iset][0]->Integral());
		h1d_svx[iset][1]->Scale(1./h1d_svx[iset][1]->Integral());

		h1d_scaled_fvtx[iset][1]->Divide(h1d_scaled_fvtx[iset][0]);
		h1d_scaled_fvtx[iset][1]->Scale(1./mean_ratio);

		h1d_fvtx_ratio[iset] = (TH1D*)h1d_fvtx[iset][1]->Clone(Form("h1d_fvtx_ratio_%d",iset));
		h1d_fvtx_ratio[iset]->Divide(h1d_fvtx[iset][0]);
		h1d_fvtx_ratio[iset]->Scale(1./mean_ratio);

		h1d_fvtx[iset][0]->Scale(1./h1d_fvtx[iset][0]->Integral());
		h1d_fvtx[iset][1]->Scale(1./h1d_fvtx[iset][1]->Integral());
	}


	TCanvas *c0[3];
	TCanvas *c1[3];
	
	for (int iset=0; iset<nset; iset++){
		c0[iset] = new TCanvas(Form("c0_%d",iset),Form("c0_%d",iset),1.1*2*500,500);
		c0[iset]->Divide(2,1);

		{
			c0[iset]->cd(1);
			gPad->SetMargin(0.15,0.03,0.14,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1e-5,xmax[iset],1.2*h1d_svx[iset][0]->GetMaximum());
			htmp->GetYaxis()->SetTitle("");
			htmp->GetYaxis()->SetTitleSize(0.045);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitle("N_{ch}^{|#eta|<1}");
			htmp->GetXaxis()->SetTitleSize(0.045);
			htmp->GetXaxis()->SetLabelSize(0.04);

			for (int ii=0; ii<2; ii++){
				h1d_svx[iset][ii]->SetLineColor(1+ii);
				h1d_svx[iset][ii]->Draw("same");
			}
		}

		{
			c0[iset]->cd(2);
			gPad->SetMargin(0.15,0.03,0.14,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.7,xmax[iset],1.3);
			htmp->GetYaxis()->SetTitle("#psi(2S) events / J/#psi events");
			htmp->GetYaxis()->SetTitleSize(0.045);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitle("N_{ch}^{|#eta|<1}");
			htmp->GetXaxis()->SetTitleSize(0.045);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h1d_svx_ratio[iset]->Draw("same");
		}

	}

	for (int iset=0; iset<nset; iset++){
		c1[iset] = new TCanvas(Form("c1_%d",iset),Form("c1_%d",iset),1.1*2*500,500);
		c1[iset]->Divide(2,1);

		{
			c1[iset]->cd(1);
			gPad->SetMargin(0.15,0.03,0.14,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1e-5,xmax[iset],1.2*h1d_fvtx[iset][0]->GetMaximum());
			htmp->GetYaxis()->SetTitle("");
			htmp->GetYaxis()->SetTitleSize(0.045);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<|#eta|<2.6}");
			htmp->GetXaxis()->SetTitleSize(0.045);
			htmp->GetXaxis()->SetLabelSize(0.04);

			for (int ii=0; ii<2; ii++){
				h1d_fvtx[iset][ii]->SetLineColor(1+ii);
				h1d_fvtx[iset][ii]->Draw("same");
			}
		}

		{
			c1[iset]->cd(2);
			gPad->SetMargin(0.15,0.03,0.14,0.05);

			TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.7,xmax[iset],1.3);
			htmp->GetYaxis()->SetTitle("#psi(2S) events / J/#psi events");
			htmp->GetYaxis()->SetTitleSize(0.045);
			htmp->GetYaxis()->SetLabelSize(0.04);
			htmp->GetXaxis()->SetTitle("N_{ch}^{1.2<|#eta|<2.6}");
			htmp->GetXaxis()->SetTitleSize(0.045);
			htmp->GetXaxis()->SetLabelSize(0.04);

			h1d_fvtx_ratio[iset]->Draw("same");
		}

	}

	return;

	TCanvas *c10 = new TCanvas("c10","c10",1.1*2*500,500);
	c10->Divide(2,1);
	{
		c10->cd(1);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.6,50,1.6);
		htmp->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J/#psi}/#LTN_{#psi(2S)}/N_{J/#psi}#GT");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta_{|#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleOffset(1.1);

		h1d_svx_ratio[0]->SetMarkerStyle(24);
		h1d_svx_ratio[0]->SetLineColor(1);
		h1d_svx_ratio[0]->SetMarkerColor(1);
		h1d_svx_ratio[0]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.68,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","J/#psi, #psi(2S), 2.5<|y|<4","h");
		leg->AddEntry("",Form("J/#psi event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[0][0]->GetMean()),"h");
		leg->AddEntry("",Form("#psi(2S) event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[0][1]->GetMean()),"h");
		//leg->AddEntry("",Form("#MB event, #LTdN_{ch}/d#eta#GT=%4.2f",hmult_svx[0]->GetMean()/2),"h");
		leg->Draw();

	}

	//return;

	{
		c10->cd(2);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.6,25,1.6);
		htmp->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J/#psi}/#LTN_{#psi(2S)}/N_{J/#psi}#GT");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta_{|#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleOffset(1.1);

		h1d_svx_ratio[1]->SetMarkerStyle(24);
		h1d_svx_ratio[1]->SetLineColor(1);
		h1d_svx_ratio[1]->SetMarkerColor(1);
		h1d_svx_ratio[1]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.68,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","J/#psi, #psi(2S), 1.2<|y|<2.2","h");
		leg->AddEntry("",Form("J/#psi event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[1][0]->GetMean()),"h");
		leg->AddEntry("",Form("#psi(2S) event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[1][1]->GetMean()),"h");
		//leg->AddEntry("",Form("#MB event, #LTdN_{ch}/d#eta#GT=%4.2f",hmult_svx[1]->GetMean()/2),"h");
		leg->Draw();
	}

	/*
	*/


	return;

	TCanvas *c3 = new TCanvas("c3","c3",1.1*2*500,500);
	c3->Divide(2,1);
	{
		c3->cd(1);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.6,6,1.6);
		htmp->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J/#psi}/#LTN_{#psi(2S)}/N_{J/#psi}#GT");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta/#LTdN_{ch}/d#eta#GT_{|#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleOffset(1.1);

		h1d_scaled_svx[0][1]->SetMarkerStyle(24);
		h1d_scaled_svx[0][1]->SetLineColor(1);
		h1d_scaled_svx[0][1]->SetMarkerColor(1);
		h1d_scaled_svx[0][1]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.68,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","J/#psi, #psi(2S), 2.5<|y|<4","h");
		leg->AddEntry("",Form("J/#psi event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[0][0]->GetMean()),"h");
		leg->AddEntry("",Form("#psi(2S) event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[0][1]->GetMean()),"h");
		//leg->AddEntry("",Form("#MB event, #LTdN_{ch}/d#eta#GT=%4.2f",hmult_svx[0]->GetMean()/2),"h");
		leg->Draw();

	}

	{
		c3->cd(2);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.6,6,1.6);
		htmp->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J/#psi}/#LTN_{#psi(2S)}/N_{J/#psi}#GT");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta/#LTdN_{ch}/d#eta#GT_{|#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleOffset(1.1);

		h1d_scaled_svx[1][1]->SetMarkerStyle(24);
		h1d_scaled_svx[1][1]->SetLineColor(1);
		h1d_scaled_svx[1][1]->SetMarkerColor(1);
		h1d_scaled_svx[1][1]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.68,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","J/#psi, #psi(2S), 1.2<|y|<2.2","h");
		leg->AddEntry("",Form("J/#psi event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[1][0]->GetMean()),"h");
		leg->AddEntry("",Form("#psi(2S) event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[1][1]->GetMean()),"h");
		//leg->AddEntry("",Form("#MB event, #LTdN_{ch}/d#eta#GT=%4.2f",hmult_svx[1]->GetMean()/2),"h");
		leg->Draw();
	}


	TCanvas *c2 = new TCanvas("c2","c2",1.1*2*500,500);
	c2->Divide(2,1);
	{
		c2->cd(1);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.6,50,1.6);
		htmp->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J/#psi}/#LTN_{#psi(2S)}/N_{J/#psi}#GT");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta_{|#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleOffset(1.1);

		h1d_svx_ratio[0]->SetMarkerStyle(24);
		h1d_svx_ratio[0]->SetLineColor(1);
		h1d_svx_ratio[0]->SetMarkerColor(1);
		h1d_svx_ratio[0]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.68,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 13 TeV","h");
		leg->AddEntry("","J/#psi, #psi(2S), 2.5<|y|<4","h");
		leg->AddEntry("",Form("J/#psi event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[0][0]->GetMean()),"h");
		leg->AddEntry("",Form("#psi(2S) event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[0][1]->GetMean()),"h");
		//leg->AddEntry("",Form("#MB event, #LTdN_{ch}/d#eta#GT=%4.2f",hmult_svx[0]->GetMean()/2),"h");
		leg->Draw();

	}

	{
		c2->cd(2);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.6,25,1.6);
		htmp->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J/#psi}/#LTN_{#psi(2S)}/N_{J/#psi}#GT");
		htmp->GetYaxis()->SetTitleSize(0.045);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta_{|#eta|<1}");
		htmp->GetXaxis()->SetTitleSize(0.045);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleOffset(1.1);

		h1d_svx_ratio[1]->SetMarkerStyle(24);
		h1d_svx_ratio[1]->SetLineColor(1);
		h1d_svx_ratio[1]->SetMarkerColor(1);
		h1d_svx_ratio[1]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.68,0.5,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","J/#psi, #psi(2S), 1.2<|y|<2.2","h");
		leg->AddEntry("",Form("J/#psi event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[1][0]->GetMean()),"h");
		leg->AddEntry("",Form("#psi(2S) event, #LTdN_{ch}/d#eta#GT=%4.2f",h1d_svx[1][1]->GetMean()),"h");
		//leg->AddEntry("",Form("#MB event, #LTdN_{ch}/d#eta#GT=%4.2f",hmult_svx[1]->GetMean()/2),"h");
		leg->Draw();
	}

	return;

	/*

	TCanvas *c3 = new TCanvas("c3","c3",1.1*2*500,500);
	c3->Divide(2,1);
	{
		c3->cd(1);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1e-5,35,1.2*hmult3[0][0]->GetMaximum());
		htmp->GetYaxis()->SetTitle("");
		htmp->GetYaxis()->SetTitleSize(0.06);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("N_{FVTXS+FVTXN+SVX}");
		htmp->GetXaxis()->SetTitleSize(0.06);
		htmp->GetXaxis()->SetLabelSize(0.04);

		for (int iset=0; iset<nset; iset++){
			hmult3[iset][0]->SetLineColor(1+iset);
			hmult3[iset][0]->Draw("same");
		}

		TLegend *leg = new TLegend(0.5,0.63,0.9,0.93);
		leg->SetFillStyle(0);
		leg->SetTextSize(0.04);
		leg->SetBorderSize(0);
		leg->AddEntry("","PYTHIA8 pp 200 GeV","h");
		leg->AddEntry("","Reconstruction level","h");
		leg->AddEntry(hfvtxz_mult3[0][0],"J/#psi, 1.2<|y|<2.2","p");
		leg->AddEntry("",Form("mean=%4.2f",hmult3[0][0]->GetMean()),"");
		leg->AddEntry(hfvtxz_mult3[1][0],"#psi(2S), 1.2<|y|<2.2","p");
		leg->AddEntry("",Form("mean=%4.2f",hmult3[1][0]->GetMean()),"");
		leg->Draw();
	}

	{
		c3->cd(2);
		gPad->SetMargin(0.15,0.03,0.14,0.05);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0.7,35,1.3);
		htmp->GetYaxis()->SetTitle("#psi(2S) events / J/#psi events");
		htmp->GetYaxis()->SetTitleSize(0.06);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("N_{FVTXS+FVTXN+SVX}");
		htmp->GetXaxis()->SetTitleSize(0.06);
		htmp->GetXaxis()->SetLabelSize(0.04);

		TH1D *hmult3_ratio = (TH1D*)hmult3[1][0]->Clone("hmult3_ratio");
		hmult3_ratio->Divide(hmult3[0][0]);
		hmult3_ratio->Draw("same");

	}
	*/

}
