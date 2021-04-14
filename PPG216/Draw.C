void Draw(){

	const int nset = 3;
	string system[nset] = {"pAu200", "dAu200", "HeAu200"};
	string legname[nset] = {"0-5% p+Au 200 GeV", "0-5% d+Au 200 GeV", "0-5% ^{3}He+Au 200 GeV"};

	int v2_np[nset] = {0}, v3_np[nset] = {0};

	float v2_pt[nset][30], v3_pt[nset][30];
	float v2[nset][30], v3[nset][30];
	float v2_err[nset][30], v3_err[nset][30];
	float v2_syserr_pos[nset][30], v3_syserr_pos[nset][30];
	float v2_syserr_neg[nset][30], v3_syserr_neg[nset][30];
	float v2_syserr_pt[nset][30], v3_syserr_pt[nset][30];

	char fname[300], buf[300];

	ifstream fdata;

	TGraphErrors *gv2[nset];
	TGraphErrors *gv3[nset];

	TGraphAsymmErrors *gv2_sys[nset];
	TGraphAsymmErrors *gv3_sys[nset];

	for (int iset=0; iset<nset; iset++){

		sprintf(fname, "v2_pt_%s_c0-5_asym.dat.txt", system[iset].c_str());
		fdata.open(fname);
		fdata.getline(buf, 300);
		int index = v2_np[iset];
		while ( fdata >> v2_pt[iset][index] >> v2[iset][index] >> v2_err[iset][index] >> v2_syserr_pos[iset][index] >> v2_syserr_neg[iset][index] ){
			v2_syserr_pt[iset][index] = 0.05;
			v2_np[iset]++;
			index = v2_np[iset];
		}
		fdata.close();

		sprintf(fname, "v3_pt_%s_c0-5_asym.dat.txt", system[iset].c_str());
		fdata.open(fname);
		fdata.getline(buf, 300);
		index = v3_np[iset];
		while ( fdata >> v3_pt[iset][index] >> v3[iset][index] >> v3_err[iset][index] >> v3_syserr_pos[iset][index] >> v3_syserr_neg[iset][index] ){
			v3_syserr_pt[iset][index] = 0.05;
			v3_np[iset]++;
			index = v3_np[iset];
		}
		fdata.close();

		gv2[iset] = new TGraphErrors(v2_np[iset], &v2_pt[iset][0], &v2[iset][0], 0, &v2_err[iset][0]);
		gv2[iset]->SetLineColor(1);
		gv2[iset]->SetLineWidth(2);
		gv2[iset]->SetMarkerStyle(20);

		gv3[iset] = new TGraphErrors(v3_np[iset], &v3_pt[iset][0], &v3[iset][0], 0, &v3_err[iset][0]);
		gv3[iset]->SetLineColor(4);
		gv3[iset]->SetLineWidth(2);
		gv3[iset]->SetMarkerStyle(20);
		gv3[iset]->SetMarkerColor(4);

		gv2_sys[iset] = new TGraphAsymmErrors(v2_np[iset], &v2_pt[iset][0], &v2[iset][0], &v2_syserr_pt[iset][0], &v2_syserr_pt[iset][0], &v2_syserr_neg[iset][0], &v2_syserr_pos[iset][0]);
		gv2_sys[iset]->SetLineColor(1);
		gv2_sys[iset]->SetLineWidth(1);
		gv2_sys[iset]->SetMarkerStyle(20);
		gv2_sys[iset]->SetFillStyle(0);

		gv3_sys[iset] = new TGraphAsymmErrors(v3_np[iset], &v3_pt[iset][0], &v3[iset][0], &v3_syserr_pt[iset][0], &v3_syserr_pt[iset][0], &v3_syserr_neg[iset][0], &v3_syserr_pos[iset][0]);
		gv3_sys[iset]->SetLineColor(4);
		gv3_sys[iset]->SetLineWidth(1);
		gv3_sys[iset]->SetMarkerStyle(20);
		gv3_sys[iset]->SetMarkerColor(4);
		gv3_sys[iset]->SetFillStyle(0);

	}//iset

	TCanvas *c1 = new TCanvas("c1","c1",1.1*3*400,400);
	c1->Divide(3,1);
	for (int iset=0; iset<nset; iset++){
		c1->cd(iset+1);
		gPad->SetMargin(0.14,0.05,0.12,0.05);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0, 0, 3, 0.2);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetYaxis()->SetTitle("v_{n}");

		gv2_sys[iset]->Draw("2");
		gv3_sys[iset]->Draw("2");
		gv2[iset]->Draw("P");
		gv3[iset]->Draw("P");

		TLegend *leg = new TLegend(0.2,0.65,0.5,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry("",legname[iset].c_str(),"h");
		leg->AddEntry(gv2[iset],"v_{2}","P");
		leg->AddEntry(gv3[iset],"v_{3}","P");
		leg->Draw();
	}

}
