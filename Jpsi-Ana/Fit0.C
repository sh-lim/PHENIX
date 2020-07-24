using namespace RooFit;

void Fit0(){

	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadTopMargin(0.05);

	const int nset = 2;
	const int narm = 2;

	const char *pname[nset] = {"jpsi", "psip"};

	TFile *infile = new TFile("outfile_mass_pp_sim.root","read");

	TH1D *h1d_mass[nset][narm];
	TH1D *h1d_mass_fvtx[nset][narm];

	for (int iset=0; iset<nset; iset++){
		for (int iarm=0; iarm<narm; iarm++){

			h1d_mass[iset][iarm] = (TH1D*)infile->Get(Form("h1d_mass_%s_arm%d",pname[iset],iarm));
			h1d_mass_fvtx[iset][iarm] = (TH1D*)infile->Get(Form("h1d_mass_fvtx_%s_arm%d",pname[iset],iarm));

		}
	}

	TCanvas *c1 = new TCanvas("c1","c1",1.2*2*400,2*400);
	c1->Divide(2,2);

	for (int iset=0; iset<nset; iset++){
		for (int iarm=0; iarm<narm; iarm++){

			c1->cd(narm*iset + iarm + 1);
			gPad->SetLogy();

			TH1D *hdata = (TH1D*)h1d_mass[iset][iarm]->Clone();

			float max = hdata->GetMaximum();
			TH1F *htmp = (TH1F*)gPad->DrawFrame(2.0,5e-1,6.0,3.0*max);
			htmp->GetXaxis()->SetTitle("mass (GeV/c^{2})");

			RooRealVar x("x","mass",2.0,6.0);
			RooDataHist rdh_all("rdh_all","hist",RooArgSet(x),hdata);

			float par_mean = (iset==0) ? 3.1 : 3.7;
			RooRealVar meanCB("meanCB","meanCB",par_mean,3.0,3.8);
			RooRealVar sigmaCB("sigmaCB","sigmaCB",0.15,0.05,0.50);
			RooRealVar alpha("alpha_jpsi","alpha_jpsi", 0.75, 0.0, 5.0);
			RooRealVar n("n_jpsi","n_jpsi", 5.0, 1.0, 10.0);
			RooCBShape CB("CB","CB_jpsi",x,meanCB,sigmaCB,alpha,n);

			float nNorm = hdata->Integral();
			RooRealVar nCB("nCB","nCB",0.95*nNorm,0,1.05*nNorm);

			RooAddPdf pdf_all("pdf_all","pdf_all",RooArgList(CB),RooArgList(nCB));

			RooFitResult *r_all = pdf_all.fitTo(rdh_all,Save(kTRUE),Verbose(kFALSE),Range(2.0,4.0),PrintLevel(-1));
			//RooFitResult *r_all = pdf_all.fitTo(rdh_all,Save(kTRUE),Verbose(kFALSE),Range(2.0,6.0));

			RooPlot* xframe = x.frame(Title("Data Fit"),Range(2.0,6.0));
			rdh_all.plotOn(xframe,MarkerStyle(20));
			CB.plotOn(xframe,LineColor(4),Normalization(nCB.getValV(),RooAbsReal::NumEvent));
			xframe->Draw("same");


		}
	}

	//return;

	TCanvas *c2 = new TCanvas("c2","c2",1.2*2*400,2*400);
	c2->Divide(2,2);

	for (int iset=0; iset<nset; iset++){
		for (int iarm=0; iarm<narm; iarm++){

			c2->cd(narm*iset + iarm + 1);
			gPad->SetLogy();

			TH1D *hdata = (TH1D*)h1d_mass_fvtx[iset][iarm]->Clone();

			float max = hdata->GetMaximum();
			TH1F *htmp = (TH1F*)gPad->DrawFrame(2.0,5e-1,6.0,3.0*max);
			htmp->GetXaxis()->SetTitle("mass (GeV/c^{2})");

			RooRealVar x("x","mass",2.0,6.0);
			RooDataHist rdh_all("rdh_all","hist",RooArgSet(x),hdata);

			float par_mean = (iset==0) ? 3.1 : 3.7;
			RooRealVar meanCB("meanCB","meanCB",par_mean,3.0,3.8);
			RooRealVar sigmaCB("sigmaCB","sigmaCB",0.15,0.05,0.50);
			RooRealVar alpha("alpha_jpsi","alpha_jpsi", 0.75, 0.0, 5.0);
			RooRealVar n("n_jpsi","n_jpsi", 5.0, 1.0, 10.0);
			RooCBShape CB("CB","CB_jpsi",x,meanCB,sigmaCB,alpha,n);

			float nNorm = hdata->Integral();
			RooRealVar nCB("nCB","nCB",0.95*nNorm,0,1.05*nNorm);

			RooAddPdf pdf_all("pdf_all","pdf_all",RooArgList(CB),RooArgList(nCB));

			RooFitResult *r_all = pdf_all.fitTo(rdh_all,Save(kTRUE),Verbose(kFALSE),Range(2.0,4.0),PrintLevel(-1));

			RooPlot* xframe = x.frame(Title("Data Fit"),Range(2.0,6.0));
			rdh_all.plotOn(xframe,MarkerStyle(20));
			CB.plotOn(xframe,LineColor(4),Normalization(nCB.getValV(),RooAbsReal::NumEvent));
			xframe->Draw("same");


		}
	}

}
