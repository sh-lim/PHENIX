#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRatioPlot.h"
#include "cmath"
#include "TString.h"
#include "TF1.h"
#include "TF2.h"
#include "iostream"
#include "list"
#include "fstream"
#include "TProfile.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "RooFit.h"
#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <TCanvas.h>
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include <fstream>
#include <iostream>
#include <string>
using namespace RooFit;

void jpsi_fitting_code()
{

	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadTopMargin(0.05);

	const int narm = 2;
	const int nchg = 3;
	const float fit_min = 2.0;
	const float fit_max = 5.0;

	float Njpsi[narm];
	float Npsip[narm];

	float Njpsi_err[narm];
	float Npsip_err[narm];

	float Fpsip[narm];
	float Fpsip_err[narm];

	float pairs;
	float par_sigma_sys[narm];
	float par_jpsi_sigma[narm];
	const float par_fracG_jpsi[narm] = {2.31037e-01, 1.97724e-01};
	const float par_sigmaG_jpsi[narm] = {2.13686e-01, 2.33722e-01};
	const float par_sigmaCB_jpsi[narm] = {1.01434e-01, 1.02440e-01};
	const float par_meanCB_jpsi[narm] = {3.13118e+00, 3.16754e+00};
	const float par_sigmaG_2nd_jpsi[narm] = {2.0910329, 2.2736887};
	// krista
	const float par_n_jpsi[narm] = {7.57221e+00, 8.21088e+00};
	const float par_alpha_jpsi[narm] = {1.00229e+00, 9.50002e-01};

	float svxpairss[8];
	float svxpairsn[8];
	float dimupairss[8];
	float dimupairsn[8];

	// for same, oppo
	TFile *infile_mixed = new TFile("out_fvtx_mixed_evt.root", "read");
	TFile *infile_same = new TFile("out_fvtx_evt.root", "read");

	// for svx

	// TFile *infile_mixed = new TFile("out_svx_mixed_evt.root", "read");
	// TFile *infile_same = new TFile("out_svx_evt.root", "read");

	TH1D *h1d_mass_same[narm][nchg];
	TH1D *h1d_mass_mixed[narm][nchg];

	TCanvas *c1 = new TCanvas("lower_south", "lower_south", 1.2 * 3 * 500, 2 * 500);
	c1->Divide(4, 2);
	TCanvas *c2 = new TCanvas("lower_north", "lower_north", 1.2 * 3 * 500, 2 * 500);
	c2->Divide(4, 2);
	TCanvas *c3 = new TCanvas("uppper_south", "upper_south", 1.2 * 2 * 500, 1 * 500);
	c3->Divide(2, 1);
	TCanvas *c4 = new TCanvas("upper_north", "upper_north", 1.2 * 2 * 500, 1 * 500);
	c4->Divide(2, 1);

	const float par_corr_a[narm] = {-1.96190e-01, 1.01851e-02};
	const float par_corr_b[narm] = {3.90264e-01, 3.15392e-01};
	const float par_corr_d[narm] = {2.40817e+00, 2.78569e+00};
	const float par_corr_e[narm] = {3.74928e+00, 3.65850e+00};
	for (int mult = 0; mult < 8; mult++)
	{

		// for (int iarm=0; iarm<narm; iarm++){
		for (int ichg = 0; ichg < nchg; ichg++)
		{
			//h1d_mass_same[0][ichg] = (TH1D *)infile_same->Get(Form("same%d0%d", mult, ichg + 1));
			//h1d_mass_same[1][ichg] = (TH1D *)infile_same->Get(Form("same%d1%d", mult, ichg + 1));
			//h1d_mass_mixed[0][ichg] = (TH1D *)infile_mixed->Get(Form("msx%d", ichg + 1));
			//h1d_mass_mixed[1][ichg] = (TH1D *)infile_mixed->Get(Form("mnx%d", ichg + 1));

			h1d_mass_same[0][ichg] = (TH1D *)infile_same->Get(Form("oppo%d0%d", mult, ichg + 1));
			h1d_mass_same[1][ichg] = (TH1D *)infile_same->Get(Form("oppo%d1%d", mult, ichg + 1));
			h1d_mass_mixed[0][ichg] = (TH1D *)infile_mixed->Get(Form("msx%d", ichg + 1));
			h1d_mass_mixed[1][ichg] = (TH1D *)infile_mixed->Get(Form("mnx%d", ichg + 1));
		}

		h1d_mass_mixed[0][1]->Add(h1d_mass_mixed[0][0]);
		h1d_mass_mixed[0][1]->Add(h1d_mass_mixed[0][2]);
		h1d_mass_mixed[1][1]->Add(h1d_mass_mixed[1][0]);
		h1d_mass_mixed[1][1]->Add(h1d_mass_mixed[1][2]);

		for (int iarm = 0; iarm < narm; iarm++)
		{
			if (iarm == 0)
			{
				c1->cd(mult + 1);
			}
			else
			{
				c2->cd(mult + 1);
			}
			gPad->SetLogy();
			gPad->SetTicks();
			TH1D *hdata = (TH1D *)h1d_mass_same[iarm][1]->Clone();

			float max = hdata->GetMaximum();
			float min = hdata->GetBinContent(hdata->FindBin(fit_max));

			TH1F *htmp = (TH1F *)gPad->DrawFrame(fit_min, 1e-4 * max, fit_max, 1.5 * max);

			htmp->GetXaxis()->SetTitle("mass (GeV/c^{2})");

			RooRealVar x("x", "mass", fit_min, fit_max);
			RooDataHist rdh_all("rdh_all", "hist", RooArgSet(x), hdata);

			// combinatorial background
			TH1D *hcomb = (TH1D *)h1d_mass_mixed[iarm][0]->Clone();
			// TH1D *hcomb = (TH1D*)h1d_mass_same[iarm][0]->Clone(Form("hLS_%d",iarm));
			// hcomb->Add(h1d_mass_mixed[iarm][2]);
			RooDataHist rdh_comb("rdh_comb", "Comb hist", RooArgSet(x), hcomb);
			RooHistPdf pdf_comb("pdf_comb", "Comb pdf", RooArgSet(x), rdh_comb);

			float nNorm = hdata->Integral(hdata->FindBin(fit_min + 0.001), hdata->FindBin(fit_max - 0.001));
			float nNorm_pp = h1d_mass_same[iarm][0]->Integral(hdata->FindBin(fit_min + 0.001), hdata->FindBin(fit_max - 0.001));
			float nNorm_mm = h1d_mass_same[iarm][2]->Integral(hdata->FindBin(fit_min + 0.001), hdata->FindBin(fit_max - 0.001));
			float nComb = 2 * sqrt(nNorm_pp * nNorm_mm);

			// Jpsi
			RooRealVar meanCB_jpsi("meanCB_jpsi", "meanCB_jpsi", 3.1, 3.0, 3.25);
			RooRealVar sigmaCB_jpsi("sigmaCB_jpsi", "sigmaCB_jpsi", 0.12, 0.08, 0.30);
			RooRealVar alpha_jpsi("alpha_jpsi", "alpha_jpsi", par_alpha_jpsi[iarm]);
			RooRealVar n_jpsi("n_jpsi", "n_jpsi", par_n_jpsi[iarm]);
			RooCBShape CB_jpsi("CB_jpsi", "CB_jpsi", x, meanCB_jpsi, sigmaCB_jpsi, alpha_jpsi, n_jpsi);

			RooRealVar fracG_jpsi("fracG_jpsi", "fracG_jpsi", par_fracG_jpsi[iarm]);
			RooRealVar scale_sigmaG_jpsi("scale_sigmaG_jpsi", "scale_sigmaG_jpsi", par_sigmaG_2nd_jpsi[iarm]);
			RooFormulaVar sigmaG_jpsi("sigmaG_jpsi", "@0*@1", RooArgSet(sigmaCB_jpsi, scale_sigmaG_jpsi));
			RooGaussian G_jpsi("G_jpsi", "G_jpsi", x, meanCB_jpsi, sigmaG_jpsi);

			RooAddPdf pdf_jpsi("pdf_jpsi", "pdf_jpsi", RooArgList(G_jpsi, CB_jpsi), fracG_jpsi);

			// Psip
			RooFormulaVar meanCB_psip("meanCB_psip", "@0 + 0.589*@0/3.0969", RooArgSet(meanCB_jpsi));
			RooRealVar sigma_sys("sigma_sys", "sigma_sys", par_sigma_sys[iarm]);
			RooFormulaVar sigmaCB_psip("sigmaCB_psip", "@0*1.15", RooArgSet(sigmaCB_jpsi));

			RooCBShape CB_psip("CB_psip", "CB_psip", x, meanCB_psip, sigmaCB_psip, alpha_jpsi, n_jpsi);

			RooFormulaVar sigmaG_psip("sigmaG_psip", "@0*1.15", RooArgSet(sigmaG_jpsi));
			RooGaussian G_psip("G_psip", "G_psip", x, meanCB_psip, sigmaG_psip);

			RooAddPdf pdf_psip("pdf_psip", "pdf_psip", RooArgList(G_psip, CB_psip), fracG_jpsi);
			// Signal pdf
			RooRealVar fpsip("fpsip", "fpsip", 0.05, 0, 1);
			RooAddPdf pdf_sig("pdf_sig", "pdf_sig", RooArgList(pdf_psip, pdf_jpsi), fpsip);

			// Correlated BKG
			RooRealVar exp_alpha("exp_alpha", "exp_alpha", -2, -5, 5);
			RooExponential pdf_exp("pdf_exp", "Correlated Background", x, exp_alpha);
			RooRealVar corr_a("corr_a", "corr_a", par_corr_a[iarm]);
			RooRealVar corr_b("corr_b", "corr_b", par_corr_b[iarm]);
			RooRealVar corr_c("corr_c", "corr_c", 0.001, 0, 1);
			RooRealVar corr_d("corr_d", "corr_d", par_corr_d[iarm]);
			RooRealVar corr_e("corr_e", "corr_e", par_corr_e[iarm] * 2, 0, 10); // ori
			RooGenericPdf pdf_corr("pdf_corr", "corr_c/(pow(exp(-corr_a*x - corr_b*x*x) + x/corr_d,corr_e))", RooArgSet(x, corr_a, corr_b, corr_c, corr_d, corr_e));
			//  orignal fitting
			// Add pdf
			RooRealVar nsig("nsig", "nsig", 0.85 * nNorm, 0, nNorm);
			RooRealVar ncomb("ncomb", "ncomb", nComb);
			RooRealVar nexp("nexp", "nexp", 0.1 * nNorm, 0, 0.3 * nNorm);
			RooAddPdf pdf_all("pdf_all", "pdf_all", RooArgList(pdf_sig, pdf_corr, pdf_comb), RooArgList(nsig, nexp, ncomb));

			// Fit
			RooFitResult *r_all = pdf_all.fitTo(rdh_all, Save(kTRUE), Verbose(kFALSE), Range(fit_min, fit_max));
			Njpsi[iarm] = nsig.getValV() * (1 - fpsip.getValV());
			Npsip[iarm] = nsig.getValV() * fpsip.getValV();
			Fpsip[iarm] = fpsip.getValV();
			Fpsip_err[iarm] = fpsip.getError();

			RooPlot *xframe = x.frame(Title("Data Fit"), Range(fit_min, fit_max));
			pdf_all.plotOn(xframe, Name("Total"), LineColor(45), LineWidth(2), Normalization(nNorm, RooAbsReal::NumEvent), Range(fit_min - 0.01, fit_max));
			pdf_comb.plotOn(xframe, Name("comb_BKG"), LineColor(2), LineWidth(2), Normalization(nComb, RooAbsReal::NumEvent), Range(fit_min - 0.01, fit_max));
			// pdf_exp.plotOn(xframe,LineColor(8),Normalization(nexp.getValV(),RooAbsReal::NumEvent),Range(fit_min,fit_max));
			pdf_corr.plotOn(xframe, Name("corr_BKG"), LineColor(8), LineWidth(2), Normalization(nexp.getValV(), RooAbsReal::NumEvent), Range(fit_min, fit_max));
			//pdf_sig.plotOn(xframe, Name("Signal"), LineColor(4), Normalization(nsig.getValV(), RooAbsReal::NumEvent), Range(fit_min, fit_max));
			pdf_jpsi.plotOn(xframe, Name("Jpsi(1S)"), LineColor(4), LineWidth(2), Normalization(Njpsi[iarm], RooAbsReal::NumEvent), Range(fit_min, fit_max));
			pdf_psip.plotOn(xframe, Name("Psi(2S)"), LineColor(6), LineWidth(2), Normalization(Npsip[iarm], RooAbsReal::NumEvent), Range(fit_min, fit_max));
			//G_jpsi.plotOn(xframe, Name("Sec_Gaussian(1S)"), LineColor(40), LineStyle(10), Normalization(Njpsi[iarm] * fracG_jpsi.getValV(), RooAbsReal::NumEvent), Range(fit_min, fit_max));
			//G_psip.plotOn(xframe, Name("Psi(122S)"), LineColor(40), LineStyle(10), Normalization(Npsip[iarm] * fracG_jpsi.getValV(), RooAbsReal::NumEvent), Range(fit_min, fit_max));
			rdh_all.plotOn(xframe, MarkerStyle(20));
			xframe->Draw("same");

			cout << fit_min << fit_max << endl;
			Fpsip[iarm] = fpsip.getValV();
			Fpsip_err[iarm] = fpsip.getError();

			TLegend *leg = new TLegend(0.7, 0.93, 0.55, 0.6);
			leg->SetFillStyle(0);
			leg->SetTextSize(0.04);
			leg->SetBorderSize(0);
			leg->SetTextFont(42);

			leg->AddEntry("", "PHENIX pp 200 GeV", "h");

			if (iarm == 0)
			{
				leg->AddEntry("", "MuTr+FVTX south arm", "h");
			}
			else if (iarm == 1)
			{
				leg->AddEntry("", "MuTr+FVTX north arm", "h");
			}

			leg->AddEntry("", "all fix parameter", "h");
			leg->AddEntry("Jpsi(1S)", "J/#psi", "l");
			leg->AddEntry("Psi(2S)", "#psi(2S)", "l");
			leg->AddEntry("corr_BKG", "corr bkg", "l");
			leg->AddEntry("comb_BKG", "comb bkg", "l");
			leg->AddEntry("Sec_Gaussian(1S)", "Second_Gaussian", "l");
			leg->AddEntry("jpsi_CB", "J/#psi_and_#psi_CB shape", "l");
			leg->AddEntry("", Form("#psi(2S)/J/#psi ratio : %4.4f", Fpsip[iarm]), "h");
			leg->AddEntry("", Form("ratio err : %4.4f", Fpsip_err[iarm]), "h");
			leg->Draw();

			if ( iarm==1 && mult==0 ){

				c4->cd(1);
				gPad->SetMargin(0.13, 0.01, 0.13, 0.03);
				gPad->SetTicks();
				gPad->SetLogy();

				TH1F *htmp = (TH1F *)gPad->DrawFrame(fit_min, 1e-4 * max, fit_max, 10 * max);
				htmp->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (GeV/c^{2})");
				htmp->GetXaxis()->SetLabelSize(0.05);
				htmp->GetXaxis()->SetTitleSize(0.06);
				htmp->GetXaxis()->SetNdivisions(9,4,0);
				htmp->GetXaxis()->SetTitleOffset(1.0);
				htmp->GetYaxis()->SetTitle("Raw counts / (100 MeV/c^{2})");
				htmp->GetYaxis()->SetLabelSize(0.05);
				htmp->GetYaxis()->SetTitleSize(0.06);
				htmp->GetYaxis()->SetTitleOffset(1.1);
				htmp->GetYaxis()->SetNdivisions(9,5,0);
				xframe->Draw("same");

				TLegend *leg = new TLegend(0.16, 0.94-0.06*4, 0.4, 0.94);
				leg->SetFillStyle(0);
				leg->SetTextSize(0.05);
				leg->SetBorderSize(0);
				leg->AddEntry("", "p+p #sqrt{s} = 200 GeV", "h");
				leg->AddEntry("", "1.2 < y < 2.2", "h");
				//leg->AddEntry("", "N_{ch}/#LTN_{ch}#GT = 0.5", "h");
				leg->AddEntry("", "0 < N_{ch}/#LTN_{ch}#GT < 0.6", "h");
				leg->AddEntry("", "N_{ch}: #LT|#Delta#eta|#GT = 3.4", "h");
				leg->Draw();

				leg = new TLegend(0.58, 0.94-0.055*5, 0.93, 0.94);
				leg->SetFillStyle(0);
				leg->SetTextSize(0.05);
				leg->SetBorderSize(0);
				leg->AddEntry("Total", "Total fit", "l");
				leg->AddEntry("Jpsi(1S)", "J/#psi", "l");
				leg->AddEntry("Psi(2S)", "#psi(2S)", "l");
				leg->AddEntry("corr_BKG", "Correlated bkg", "l");
				leg->AddEntry("comb_BKG", "Mixed events bkg", "l");
				leg->Draw();

			}else if ( iarm==1 && mult==6 ){

				c4->cd(2);
				gPad->SetMargin(0.13, 0.01, 0.13, 0.03);
				gPad->SetTicks();
				gPad->SetLogy();

				TH1F *htmp = (TH1F *)gPad->DrawFrame(fit_min, 1e-4 * max, fit_max, 10.0 * max);
				htmp->GetXaxis()->SetTitle("#mu^{+}#mu^{-} mass (GeV/c^{2})");
				htmp->GetXaxis()->SetLabelSize(0.05);
				htmp->GetXaxis()->SetTitleSize(0.06);
				htmp->GetXaxis()->SetNdivisions(9,4,0);
				htmp->GetXaxis()->SetTitleOffset(1.0);
				htmp->GetYaxis()->SetTitle("Raw counts / (100 MeV/c^{2})");
				htmp->GetYaxis()->SetLabelSize(0.05);
				htmp->GetYaxis()->SetTitleSize(0.06);
				htmp->GetYaxis()->SetTitleOffset(1.1);
				htmp->GetYaxis()->SetNdivisions(9,5,0);
				xframe->Draw("same");

				TLegend *leg = new TLegend(0.16, 0.94-0.06*4, 0.4, 0.94);
				leg->SetFillStyle(0);
				leg->SetTextSize(0.05);
				leg->SetBorderSize(0);
				leg->AddEntry("", "p+p #sqrt{s} = 200 GeV", "h");
				leg->AddEntry("", "1.2 < y < 2.2", "h");
				//leg->AddEntry("", "N_{ch}/#LTN_{ch}#GT = 4.2", "h");
				leg->AddEntry("", "2.2 < N_{ch}/#LTN_{ch}#GT < 3", "h");
				leg->AddEntry("", "N_{ch}: #LT|#Delta#eta|#GT = 3.4", "h");
				leg->Draw();

				leg = new TLegend(0.58, 0.94-0.055*5, 0.93, 0.94);
				leg->SetFillStyle(0);
				leg->SetTextSize(0.05);
				leg->SetBorderSize(0);
				leg->AddEntry("Total", "Total fit", "l");
				leg->AddEntry("Jpsi(1S)", "J/#psi", "l");
				leg->AddEntry("Psi(2S)", "#psi(2S)", "l");
				leg->AddEntry("corr_BKG", "Correlated bkg", "l");
				leg->AddEntry("comb_BKG", "Mixed events bkg", "l");
				leg->Draw();

			}
		}//iarm
		for (int iarm = 0; iarm < narm; iarm++)
		{
			cout
				<< "ARM: " << iarm
				<< ", F: " << Fpsip[iarm] << " +/-" << Fpsip_err[iarm]
				<< endl;
		}
	}//mult

	return;
	c1->SaveAs("svxafixjpsisouth.png");
	c2->SaveAs("svxafixjpsinorth.png");
	c3->SaveAs("svxdfixjpsisouth.png");
	c4->SaveAs("svxdfixjpsinorth.png");
	}
