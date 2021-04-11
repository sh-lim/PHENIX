#include "/phenix/u/shlim/Style.h"

void DrawQAHisto(){

	gStyle->SetOptStat(0);

	const int nset = 2;
	
	TFile *infile[nset];
	infile[0] = new TFile("Run14AuAu200_pro109_QA.root","read");
	infile[1] = new TFile("Run16AuAu200MB_pro111_QA_histos.root","read");

	string setname[nset] = {"Run14AuAu200", "Run16AuAu200"};

	TH2F *hvtx_xz_svx_only[nset];
	TH2F *hvtx_yz_svx_only[nset];

	TH2F *hvtx_xz_fvtx[nset];
	TH2F *hvtx_yz_fvtx[nset];

	TH2F *hvtx_zdiff_bbc_fvtx[nset];
	TH2F *hvtx_zdiff_svx_fvtx[nset];

	for (int iset=0; iset<nset; iset++){
		hvtx_xz_svx_only[iset] = (TH2F*)infile[iset]->Get("hvtx_xz_svx_only");
		hvtx_yz_svx_only[iset] = (TH2F*)infile[iset]->Get("hvtx_yz_svx_only");

		hvtx_xz_fvtx[iset] = (TH2F*)infile[iset]->Get("hvtx_xz_fvtx");
		hvtx_yz_fvtx[iset] = (TH2F*)infile[iset]->Get("hvtx_yz_fvtx");

		hvtx_zdiff_bbc_fvtx[iset] = (TH2F*)infile[iset]->Get("hvtx_zdiff_bbc_fvtx");
		hvtx_zdiff_svx_fvtx[iset] = (TH2F*)infile[iset]->Get("hvtx_zdiff_svx_fvtx");
	}


	TCanvas *c1[nset];

	for (int iset=0; iset<nset; iset++){
		c1[iset] = new TCanvas(Form("c1_set%d",iset),Form("c1_set%d",iset),1.2*2*400,2*400);
		c1[iset]->Divide(2,2);

		c1[iset]->cd(1);
		SetPadStyle(1);
		gPad->SetLogz();
		hvtx_xz_fvtx[iset]->Draw("colz");

		c1[iset]->cd(2);
		SetPadStyle(1);
		gPad->SetLogz();
		hvtx_yz_fvtx[iset]->Draw("colz");

		c1[iset]->cd(3);
		SetPadStyle(1);
		gPad->SetLogz();
		hvtx_xz_svx_only[iset]->Draw("colz");

		c1[iset]->cd(4);
		SetPadStyle(1);
		gPad->SetLogz();
		hvtx_yz_svx_only[iset]->Draw("colz");
	}

	TCanvas *c2 = new TCanvas("c2","c2",1.2*2*400,1*400);
	c2->Divide(2,1);

	for (int iset=0; iset<nset; iset++){
		c2->cd(iset+1);
		SetPadStyle(1);
		gPad->SetLogz();
		hvtx_zdiff_svx_fvtx[iset]->Draw("colz");

		htmp = (TH1D*)hvtx_zdiff_svx_fvtx[iset];
		SetHistoStyle("FVTX-z [cm]","SVX_PRECISE-z - FVTX-z [cm]","",0.05,0.04);
		htmp->SetTitle(setname[iset].c_str());

	}


}
