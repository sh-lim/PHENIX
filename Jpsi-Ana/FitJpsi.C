using namespace RooFit;

void FitJpsi(){

	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadTopMargin(0.05);

	const int nchg = 3;

	const float fit_min = 2.0;
	const float fit_max = 5.0;

	TFile *infile = new TFile("matched_event_masshisto.root","read");

	TH1D *h1d_mass[nchg];

	h1d_mass[0] = (TH1D*)infile->Get("m2");
	h1d_mass[1] = (TH1D*)infile->Get("m1");
	h1d_mass[2] = (TH1D*)infile->Get("m3");

	TCanvas *c1 = new TCanvas("c1","c1",1.1*500,1*500);
	gPad->SetLogy();

	TH1D *hdata = (TH1D*)h1d_mass[0]->Clone();

	float max = hdata->GetMaximum();
	float min = hdata->GetBinContent(hdata->FindBin(fit_max));

	TH1F *htmp = (TH1F*)gPad->DrawFrame(fit_min,0.01*min,fit_max,3.0*max);
	htmp->GetXaxis()->SetTitle("mass (GeV/c^{2})");

	RooRealVar x("x","mass",fit_min,fit_max);
	RooDataHist rdh_all("rdh_all","hist",RooArgSet(x),hdata);

	//combinatorial background
	TH1D *hcomb = (TH1D*)h1d_mass[1]->Clone();
	hcomb->Add(h1d_mass[2]);
	RooDataHist rdh_comb("rdh_comb","Comb hist",RooArgSet(x),hcomb);
	RooHistPdf pdf_comb("pdf_comb","Comb pdf",RooArgSet(x),rdh_comb);

	//Jpsi
	RooRealVar meanCB_jpsi("meanCB_jpsi","meanCB_jpsi",3.1,3.0,3.25);
	RooRealVar sigmaCB_jpsi("sigmaCB_jpsi","sigmaCB_jpsi",0.15,0.10,0.30);
	RooRealVar alpha_jpsi("alpha_jpsi","alpha_jpsi", 0.75, 0.0, 5.0);
	RooRealVar n_jpsi("n_jpsi","n_jpsi", 5, 1.0, 10.0);
	RooCBShape CB_jpsi("CB_jpsi","CB_jpsi",x,meanCB_jpsi,sigmaCB_jpsi,alpha_jpsi,n_jpsi);

	RooRealVar fracG_jpsi("fracG_jpsi","fracG_jpsi",0.03,0,0.5);
	RooRealVar sigmaG_jpsi("sigmaG_jpsi","sigmaG_jpsi",0.2,0.1,1.0);
	RooGaussian G_jpsi("G_jpsi","G_jpsi",x,meanCB_jpsi,sigmaG_jpsi);

	RooAddPdf pdf_jpsi("pdf_jpsi","pdf_jpsi",RooArgList(G_jpsi,CB_jpsi),fracG_jpsi);

	//Psip
	RooFormulaVar meanCB_psip("meanCB_psip","@0 + 0.589*@0/3.0969",RooArgSet(meanCB_jpsi));
	RooFormulaVar sigmaCB_psip("sigmaCB_psip","@0*1.15",RooArgSet(sigmaCB_jpsi));
	RooCBShape CB_psip("CB_psip","CB_psip",x,meanCB_psip,sigmaCB_psip,alpha_jpsi,n_jpsi);

	RooRealVar sigmaG_psip("sigmaG_psip","sigmaG_psip",0.3,0.1,1.0);
	RooGaussian G_psip("G_psip","G_psip",x,meanCB_psip,sigmaG_jpsi);

	RooAddPdf pdf_psip("pdf_psip","pdf_psip",RooArgList(G_psip,CB_psip),fracG_jpsi);

	//Correlated BKG
	RooRealVar exp_alpha("exp_alpha","exp_alpha",-1,-5,5);
	RooExponential pdf_exp("pdf_exp","Correlated Background",x,exp_alpha);

	RooRealVar par_eff("par_eff","par_eff",2.50);
	RooFormulaVar eff("eff","(TMath::Erf((@0-@1)/1.0)+1)",RooArgSet(x,par_eff));

	RooEffProd pdf_exp_eff("pdf_exp_eff","model with efficiency",pdf_exp,eff);

	//Add pdf
	float nNorm = hdata->Integral(hdata->FindBin(fit_min+0.001),hdata->FindBin(fit_max-0.001));
	float nComb = hcomb->Integral(hcomb->FindBin(fit_min+0.001),hcomb->FindBin(fit_max-0.001));
	RooRealVar njpsi("njpsi","njpsi",0.9*nNorm,0,nNorm);
	RooRealVar npsip("npsip","npsip",0.05*nNorm,0,0.1*nNorm);
	RooRealVar ncomb("ncomb","ncomb",nComb);
	RooRealVar nexp("nexp","nexp",0.1*nNorm,0,0.3*nNorm);
	//RooAddPdf pdf_all("pdf_all","pdf_all",RooArgList(pdf_jpsi,pdf_psip,pdf_exp,pdf_comb),RooArgList(njpsi,npsip,nexp,ncomb));
	RooAddPdf pdf_all("pdf_all","pdf_all",RooArgList(pdf_jpsi,pdf_psip,pdf_exp_eff,pdf_comb),RooArgList(njpsi,npsip,nexp,ncomb));

	//Fit
	RooFitResult *r_all = pdf_all.fitTo(rdh_all,Save(kTRUE),Verbose(kFALSE),Range(fit_min,fit_max));

	//Draw
	RooPlot* xframe = x.frame(Title("Data Fit"),Range(fit_min,fit_max));
	rdh_all.plotOn(xframe,MarkerStyle(20));
	pdf_all.plotOn(xframe,LineColor(2),Normalization(nNorm,RooAbsReal::NumEvent),Range(fit_min,fit_max));
	pdf_comb.plotOn(xframe,LineColor(2),Normalization(nComb,RooAbsReal::NumEvent),Range(fit_min,fit_max));
	//pdf_exp.plotOn(xframe,LineColor(8),Normalization(nexp.getValV(),RooAbsReal::NumEvent),Range(fit_min,fit_max));
	pdf_exp_eff.plotOn(xframe,LineColor(8),Normalization(nexp.getValV(),RooAbsReal::NumEvent),Range(fit_min,fit_max));
	pdf_jpsi.plotOn(xframe,LineColor(4),Normalization(njpsi.getValV(),RooAbsReal::NumEvent),Range(fit_min,fit_max));
	pdf_psip.plotOn(xframe,LineColor(4),Normalization(npsip.getValV(),RooAbsReal::NumEvent),Range(fit_min,fit_max));
	xframe->Draw("same");

}
