#include "Style.h"

void Draw_Kpi_pAu(){

	TFile *infile = new TFile("dAu200_kpi_ratio_chg_func.root","read");
	TFile *infile1 = new TFile("runcuts_Run15pAu200_COMBINED_QGSP_BERT_TIGHT.root","read");
	TFile *infile2 = new TFile("runcuts_Run15pAu200_COMBINED_FTFP_BERT_TIGHT.root","read");
	TFile *infile3 = new TFile("runcuts_Run15pp200_COMBINED_QGSP_BIC_TIGHT.root","read");

	TH1D *hPT_KPI[2][2]; //[arm][chg]
	hPT_KPI[0][0] = (TH1D*)infile->Get("hkpi_r_dAu_chg0_etabin06");
	hPT_KPI[0][1] = (TH1D*)infile->Get("hkpi_r_dAu_chg1_etabin06");
	hPT_KPI[1][0] = (TH1D*)infile->Get("hkpi_r_dAu_chg0_etabin23");
	hPT_KPI[1][1] = (TH1D*)infile->Get("hkpi_r_dAu_chg1_etabin23");

	TH1D *hreco1_PION[2][2][2]; //[arm][gap][chg]
	TH1D *hreco1_KAON[2][2][2]; //[arm][gap][chg]
	TH1D *hreco2_PION[2][2][2]; //[arm][gap][chg]
	TH1D *hreco2_KAON[2][2][2]; //[arm][gap][chg]
	TH1D *hreco3_PION[2][2][2]; //[arm][gap][chg]
	TH1D *hreco3_KAON[2][2][2]; //[arm][gap][chg]

	for (int iarm=0; iarm<2; iarm++){
		for (int igap=0; igap<2; igap++){
			for (int ichg=0; ichg<2; ichg++){
				hreco1_PION[iarm][igap][ichg] = (TH1D*)infile1->Get(Form("n_pion_arm%d_gap%d_chg%d",iarm,igap+2,ichg));
				hreco1_KAON[iarm][igap][ichg] = (TH1D*)infile1->Get(Form("n_kaon_arm%d_gap%d_chg%d",iarm,igap+2,ichg));
				hreco1_PION[iarm][igap][ichg]->Rebin(2);
				hreco1_KAON[iarm][igap][ichg]->Rebin(2);

				hreco2_PION[iarm][igap][ichg] = (TH1D*)infile2->Get(Form("n_pion_arm%d_gap%d_chg%d",iarm,igap+2,ichg));
				hreco2_KAON[iarm][igap][ichg] = (TH1D*)infile2->Get(Form("n_kaon_arm%d_gap%d_chg%d",iarm,igap+2,ichg));
				hreco2_PION[iarm][igap][ichg]->Rebin(2);
				hreco2_KAON[iarm][igap][ichg]->Rebin(2);

				hreco3_PION[iarm][igap][ichg] = (TH1D*)infile3->Get(Form("n_pion_arm%d_gap%d_chg%d",iarm,igap+2,ichg));
				hreco3_KAON[iarm][igap][ichg] = (TH1D*)infile3->Get(Form("n_kaon_arm%d_gap%d_chg%d",iarm,igap+2,ichg));
				hreco3_PION[iarm][igap][ichg]->Rebin(2);
				hreco3_KAON[iarm][igap][ichg]->Rebin(2);
			}
		}
	}

	for (int ichg=0; ichg<2; ichg++){
		hreco1_PION[0][0][ichg]->Add(hreco1_PION[0][1][ichg]);
		hreco1_PION[0][0][ichg]->Add(hreco1_PION[1][0][ichg]);
		hreco1_PION[0][0][ichg]->Add(hreco1_PION[1][1][ichg]);

		hreco1_KAON[0][0][ichg]->Add(hreco1_KAON[0][1][ichg]);
		hreco1_KAON[0][0][ichg]->Add(hreco1_KAON[1][0][ichg]);
		hreco1_KAON[0][0][ichg]->Add(hreco1_KAON[1][1][ichg]);

		for (int ipt=0; ipt<hreco1_PION[0][0][ichg]->GetNbinsX(); ipt++){
			float dpt = hreco1_PION[0][0][ichg]->GetBinWidth(ipt+1);
			float yy = hreco1_PION[0][0][ichg]->GetBinContent(ipt+1);
			float yy_err = hreco1_PION[0][0][ichg]->GetBinError(ipt+1);

			hreco1_PION[0][0][ichg]->SetBinContent(ipt+1, yy/dpt);
			if ( yy_err/yy<0.1 ){
				hreco1_PION[0][0][ichg]->SetBinError(ipt+1, yy_err*2/dpt);
			}else{
				hreco1_PION[0][0][ichg]->SetBinError(ipt+1, yy_err/dpt);
			}

			dpt = hreco1_KAON[0][0][ichg]->GetBinWidth(ipt+1);
			yy = hreco1_KAON[0][0][ichg]->GetBinContent(ipt+1);
			yy_err = hreco1_KAON[0][0][ichg]->GetBinError(ipt+1);

			hreco1_KAON[0][0][ichg]->SetBinContent(ipt+1, yy/dpt);
			if ( yy_err/yy<0.1 ){
				hreco1_KAON[0][0][ichg]->SetBinError(ipt+1, yy_err*2/dpt);
			}else{
				hreco1_KAON[0][0][ichg]->SetBinError(ipt+1, yy_err/dpt);
			}
		}//ipt
	}//ichg

	for (int ichg=0; ichg<2; ichg++){
		hreco2_PION[0][0][ichg]->Add(hreco2_PION[0][1][ichg]);
		hreco2_PION[0][0][ichg]->Add(hreco2_PION[1][0][ichg]);
		hreco2_PION[0][0][ichg]->Add(hreco2_PION[1][1][ichg]);

		hreco2_KAON[0][0][ichg]->Add(hreco2_KAON[0][1][ichg]);
		hreco2_KAON[0][0][ichg]->Add(hreco2_KAON[1][0][ichg]);
		hreco2_KAON[0][0][ichg]->Add(hreco2_KAON[1][1][ichg]);

		for (int ipt=0; ipt<hreco2_PION[0][0][ichg]->GetNbinsX(); ipt++){
			float dpt = hreco2_PION[0][0][ichg]->GetBinWidth(ipt+1);
			float yy = hreco2_PION[0][0][ichg]->GetBinContent(ipt+1);
			float yy_err = hreco2_PION[0][0][ichg]->GetBinError(ipt+1);

			hreco2_PION[0][0][ichg]->SetBinContent(ipt+1, yy/dpt);
			if ( yy_err/yy<0.1 ){
				hreco2_PION[0][0][ichg]->SetBinError(ipt+1, yy_err*2/dpt);
			}else{
				hreco2_PION[0][0][ichg]->SetBinError(ipt+1, yy_err/dpt);
			}

			dpt = hreco2_KAON[0][0][ichg]->GetBinWidth(ipt+1);
			yy = hreco2_KAON[0][0][ichg]->GetBinContent(ipt+1);
			yy_err = hreco2_KAON[0][0][ichg]->GetBinError(ipt+1);

			hreco2_KAON[0][0][ichg]->SetBinContent(ipt+1, yy/dpt);
			if ( yy_err/yy<0.1 ){
				hreco2_KAON[0][0][ichg]->SetBinError(ipt+1, yy_err*2/dpt);
			}else{
				hreco2_KAON[0][0][ichg]->SetBinError(ipt+1, yy_err/dpt);
			}
		}//ipt
	}//ichg

	for (int ichg=0; ichg<2; ichg++){
		hreco3_PION[0][0][ichg]->Add(hreco3_PION[0][1][ichg]);
		hreco3_PION[0][0][ichg]->Add(hreco3_PION[1][0][ichg]);
		hreco3_PION[0][0][ichg]->Add(hreco3_PION[1][1][ichg]);

		hreco3_KAON[0][0][ichg]->Add(hreco3_KAON[0][1][ichg]);
		hreco3_KAON[0][0][ichg]->Add(hreco3_KAON[1][0][ichg]);
		hreco3_KAON[0][0][ichg]->Add(hreco3_KAON[1][1][ichg]);

		for (int ipt=0; ipt<hreco3_PION[0][0][ichg]->GetNbinsX(); ipt++){
			float dpt = hreco3_PION[0][0][ichg]->GetBinWidth(ipt+1);
			float yy = hreco3_PION[0][0][ichg]->GetBinContent(ipt+1);
			float yy_err = hreco3_PION[0][0][ichg]->GetBinError(ipt+1);

			hreco3_PION[0][0][ichg]->SetBinContent(ipt+1, yy/dpt);
			if ( yy_err/yy<0.1 ){
				hreco3_PION[0][0][ichg]->SetBinError(ipt+1, yy_err*2/dpt);
			}else{
				hreco3_PION[0][0][ichg]->SetBinError(ipt+1, yy_err/dpt);
			}

			dpt = hreco3_KAON[0][0][ichg]->GetBinWidth(ipt+1);
			yy = hreco3_KAON[0][0][ichg]->GetBinContent(ipt+1);
			yy_err = hreco3_KAON[0][0][ichg]->GetBinError(ipt+1);

			hreco3_KAON[0][0][ichg]->SetBinContent(ipt+1, yy/dpt);
			if ( yy_err/yy<0.1 ){
				hreco3_KAON[0][0][ichg]->SetBinError(ipt+1, yy_err*2/dpt);
			}else{
				hreco3_KAON[0][0][ichg]->SetBinError(ipt+1, yy_err/dpt);
			}
		}//ipt
	}//ichg

	//return;

	TH1D *hreco1_ratio[2];
	TH1D *hreco2_ratio[2];
	TH1D *hreco3_ratio[2];

	TCanvas *c0_1 = new TCanvas("c0_1","c0_1",1.2*2*400,400);
	c0_1->Divide(2,1);

	TCanvas *c0_2 = new TCanvas("c0_2","c0_2",1.2*2*400,400);
	c0_2->Divide(2,1);

	TCanvas *c0_3 = new TCanvas("c0_3","c0_3",1.2*2*400,400);
	c0_3->Divide(2,1);

	for (int ichg=0; ichg<2; ichg++){

		c0_1->cd(ichg+1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1D*)gPad->DrawFrame(1.0, 0.1, 10.0, 2.0*hreco1_PION[0][0][0]->GetMaximum());
		SetHistoStyle();

		hreco1_PION[0][0][ichg]->SetLineColor(1);
		hreco1_PION[0][0][ichg]->SetLineWidth(2);
		hreco1_PION[0][0][ichg]->Draw("same");

		hreco1_KAON[0][0][ichg]->SetLineColor(2);
		hreco1_KAON[0][0][ichg]->SetLineWidth(2);
		hreco1_KAON[0][0][ichg]->Draw("same");

		//TF1 *fpion = new TF1(Form("fpion_chg%d",ichg),"[0]*(exp([1]*x) + (x/[2]))^[3]",1.4,10.0);
		//fpion->SetParameters(hreco1_PION[0][0][ichg]->GetMaximum(), 0.1, 1.1, -10.0);
		TF1 *fpion = new TF1(Form("fpion_chg%d",ichg),"[0]*(exp([1]*x+[2]*x*x) + (x/[3]))^[4]",1.4,10.0);
		fpion->SetParameters(hreco1_PION[0][0][ichg]->GetMaximum(), -0.1, +0.1, 2.0, -8.0);
		hreco1_PION[0][0][ichg]->Fit(fpion,"R0Q");
		hreco1_PION[0][0][ichg]->Fit(fpion,"R0Q");
		hreco1_PION[0][0][ichg]->Fit(fpion,"R0Q");
		fpion->SetLineWidth(2);
		fpion->SetLineColor(1);
		fpion->Draw("same");

		//TF1 *fkaon = new TF1(Form("fkaon_chg%d",ichg),"[0]*(exp([1]*x) + (x/[2]))^[3]",1.4,10.0);
		//fkaon->SetParameters(hreco1_KAON[0][0][ichg]->GetMaximum(), 0.1, 1.1, -10.0);
		TF1 *fkaon = new TF1(Form("fkaon_chg%d",ichg),"[0]*(exp([1]*x + [2]*x*x) + (x/[3]))^[4]",1.4,10.0);
		fkaon->SetParameters(0.1*hreco1_KAON[0][0][ichg]->GetMaximum(), -1.0, +0.1, 2.0, -8.0);
		hreco1_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		hreco1_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		hreco1_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		fkaon->SetLineWidth(2);
		fkaon->SetLineColor(2);
		fkaon->Draw("same");

		hreco1_ratio[ichg] = (TH1D*)hreco1_KAON[0][0][ichg]->Clone(Form("hreco1_ratio_chg%d",ichg));
		hreco1_ratio[ichg]->Divide(hreco1_PION[0][0][ichg]);
		hreco1_ratio[ichg]->SetLineWidth(4);
		hreco1_ratio[ichg]->SetLineColor(ichg+1);
		hreco1_ratio[ichg]->SetLineStyle(2);

		for (int ipt=0; ipt<hreco1_ratio[ichg]->GetNbinsX(); ipt++){

			hreco1_ratio[ichg]->SetBinError(ipt+1, 0);

			float xx = hreco1_ratio[ichg]->GetBinCenter(ipt+1);
			if ( xx<1.4 ) continue;
			float yy = fkaon->Eval(xx) / fpion->Eval(xx);

			hreco1_ratio[ichg]->SetBinContent(ipt+1, yy);
		}

	}//

	//return;

	for (int ichg=0; ichg<2; ichg++){

		c0_2->cd(ichg+1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1D*)gPad->DrawFrame(1.0, 0.1, 10.0, 2.0*hreco2_PION[0][0][0]->GetMaximum());
		SetHistoStyle();

		hreco2_PION[0][0][ichg]->SetLineColor(1);
		hreco2_PION[0][0][ichg]->SetLineWidth(2);
		hreco2_PION[0][0][ichg]->Draw("same");

		hreco2_KAON[0][0][ichg]->SetLineColor(2);
		hreco2_KAON[0][0][ichg]->SetLineWidth(2);
		hreco2_KAON[0][0][ichg]->Draw("same");

		//TF1 *fpion = new TF1(Form("fpion_chg%d",ichg),"[0]*(exp([1]*x) + (x/[2]))^[3]",1.4,10.0);
		//fpion->SetParameters(hreco2_PION[0][0][ichg]->GetMaximum(), 0.1, 1.1, -10.0);
		TF1 *fpion = new TF1(Form("fpion_chg%d",ichg),"[0]*(exp([1]*x + [2]*x*x) + (x/[3]))^[4]",1.4,10.0);
		fpion->SetParameters(hreco2_PION[0][0][ichg]->GetMaximum(), -0.1, +0.1, 1.1, -8.0);
		hreco2_PION[0][0][ichg]->Fit(fpion,"R0Q");
		hreco2_PION[0][0][ichg]->Fit(fpion,"R0Q");
		hreco2_PION[0][0][ichg]->Fit(fpion,"R0Q");
		fpion->SetLineWidth(2);
		fpion->SetLineColor(1);
		fpion->Draw("same");

		//TF1 *fkaon = new TF1(Form("fkaon_chg%d",ichg),"[0]*(exp([1]*x) + (x/[2]))^[3]",1.4,10.0);
		//fkaon->SetParameters(hreco2_KAON[0][0][ichg]->GetMaximum(), 0.1, 1.1, -10.0);
		TF1 *fkaon = new TF1(Form("fkaon_chg%d",ichg),"[0]*(exp([1]*x + [2]*x*x) + (x/[3]))^[4]",1.4,10.0);
		fkaon->SetParameters(0.01*hreco2_KAON[0][0][ichg]->GetMaximum(), -1.0, +0.1, 1.1, -10.0);
		hreco2_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		hreco2_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		hreco2_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		fkaon->SetLineWidth(2);
		fkaon->SetLineColor(2);
		fkaon->Draw("same");

		hreco2_ratio[ichg] = (TH1D*)hreco2_KAON[0][0][ichg]->Clone(Form("hreco2_ratio_chg%d",ichg));
		hreco2_ratio[ichg]->Divide(hreco2_PION[0][0][ichg]);
		hreco2_ratio[ichg]->SetLineWidth(4);
		hreco2_ratio[ichg]->SetLineColor(ichg+1);
		hreco2_ratio[ichg]->SetLineStyle(2);

		for (int ipt=0; ipt<hreco2_ratio[ichg]->GetNbinsX(); ipt++){

			hreco2_ratio[ichg]->SetBinError(ipt+1, 0);

			float xx = hreco2_ratio[ichg]->GetBinCenter(ipt+1);
			if ( xx<1.4 ) continue;
			float yy = fkaon->Eval(xx) / fpion->Eval(xx);

			hreco2_ratio[ichg]->SetBinContent(ipt+1, yy);
		}

	}//

	//return;

	for (int ichg=0; ichg<2; ichg++){

		c0_3->cd(ichg+1);
		SetPadStyle();
		gPad->SetLogy();

		htmp = (TH1D*)gPad->DrawFrame(1.0, 0.1, 10.0, 2.0*hreco3_PION[0][0][0]->GetMaximum());
		SetHistoStyle();

		hreco3_PION[0][0][ichg]->SetLineColor(1);
		hreco3_PION[0][0][ichg]->SetLineWidth(2);
		hreco3_PION[0][0][ichg]->Draw("same");

		hreco3_KAON[0][0][ichg]->SetLineColor(2);
		hreco3_KAON[0][0][ichg]->SetLineWidth(2);
		hreco3_KAON[0][0][ichg]->Draw("same");

		//TF1 *fpion = new TF1(Form("fpion_chg%d",ichg),"[0]*(exp([1]*x) + (x/[2]))^[3]",1.4,10.0);
		//fpion->SetParameters(hreco3_PION[0][0][ichg]->GetMaximum(), 0.1, 1.1, -10.0);
		TF1 *fpion = new TF1(Form("fpion_chg%d",ichg),"[0]*(exp([1]*x + [2]*x*x) + (x/[3]))^[4]",1.4,10.0);
		fpion->SetParameters(hreco3_PION[0][0][ichg]->GetMaximum(), -0.1, +0.1, 1.1, -8.0);
		hreco3_PION[0][0][ichg]->Fit(fpion,"R0Q");
		hreco3_PION[0][0][ichg]->Fit(fpion,"R0Q");
		hreco3_PION[0][0][ichg]->Fit(fpion,"R0Q");
		fpion->SetLineWidth(2);
		fpion->SetLineColor(1);
		fpion->Draw("same");

		//TF1 *fkaon = new TF1(Form("fkaon_chg%d",ichg),"[0]*(exp([1]*x) + (x/[2]))^[3]",1.4,10.0);
		//fkaon->SetParameters(hreco3_KAON[0][0][ichg]->GetMaximum(), 0.1, 1.1, -10.0);
		TF1 *fkaon = new TF1(Form("fkaon_chg%d",ichg),"[0]*(exp([1]*x + [2]*x*x) + (x/[3]))^[4]",1.4,10.0);
		fkaon->SetParameters(0.01*hreco3_KAON[0][0][ichg]->GetMaximum(), -1.0, +0.1, 1.1, -8.0);
		hreco3_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		hreco3_KAON[0][0][ichg]->Fit(fkaon,"R0Q");
		hreco3_KAON[0][0][ichg]->Fit(fkaon,"R0");
		fkaon->SetLineWidth(2);
		fkaon->SetLineColor(2);
		fkaon->Draw("same");

		hreco3_ratio[ichg] = (TH1D*)hreco3_KAON[0][0][ichg]->Clone(Form("hreco3_ratio_chg%d",ichg));
		hreco3_ratio[ichg]->Divide(hreco3_PION[0][0][ichg]);
		hreco3_ratio[ichg]->SetLineWidth(4);
		hreco3_ratio[ichg]->SetLineColor(ichg+1);
		hreco3_ratio[ichg]->SetLineStyle(2);

		for (int ipt=0; ipt<hreco3_ratio[ichg]->GetNbinsX(); ipt++){

			hreco3_ratio[ichg]->SetBinError(ipt+1, 0);

			float xx = hreco3_ratio[ichg]->GetBinCenter(ipt+1);
			if ( xx<1.4 ) continue;
			float yy = fkaon->Eval(xx) / fpion->Eval(xx);

			hreco3_ratio[ichg]->SetBinContent(ipt+1, yy);
		}

	}//

	TGraphErrors *greco_ratio[2];

	for (int ichg=0; ichg<2; ichg++){
		greco_ratio[ichg] = new TGraphErrors;
		greco_ratio[ichg]->SetLineColor(ichg+1);
		greco_ratio[ichg]->SetLineStyle(2);
		greco_ratio[ichg]->SetLineWidth(4);
		greco_ratio[ichg]->SetFillColorAlpha(ichg+1,0.3);

		for (int ipt=0; ipt<hreco1_ratio[ichg]->GetNbinsX(); ipt++){
			float xx = hreco1_ratio[ichg]->GetBinCenter(ipt+1);
			float xx_err = hreco1_ratio[ichg]->GetBinWidth(ipt+1)/2;
			float yy1 = hreco1_ratio[ichg]->GetBinContent(ipt+1);
			float yy2 = hreco2_ratio[ichg]->GetBinContent(ipt+1);
			greco_ratio[ichg]->SetPoint(ipt, xx, (yy1+yy2)/2);
			greco_ratio[ichg]->SetPointError(ipt, xx_err, fabs(yy1-yy2)/2);
		}//ipt
	}//ichg

	//return;

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	SetPadStyle();
	gPad->SetRightMargin(0.03);
	gPad->SetLeftMargin(0.13);

	htmp = (TH1D*)gPad->DrawFrame(1.5,0,8.0,2.0);
	SetHistoStyle("p_{T} (GeV/c)","K/#pi","",24,24);
	htmp->GetYaxis()->SetTitleOffset(1.0);
	htmp->GetXaxis()->SetTitleOffset(1.2);
	htmp->GetXaxis()->SetLabelSize(24);
	htmp->GetXaxis()->SetTitleSize(28);
	htmp->GetYaxis()->SetLabelSize(24);
	htmp->GetYaxis()->SetTitleSize(32);

	for (int iarm=1; iarm<2; iarm++){
		for (int ichg=0; ichg<2; ichg++){
			hPT_KPI[iarm][ichg]->SetLineWidth(4);
			hPT_KPI[iarm][ichg]->SetLineColor(ichg+1);
			hPT_KPI[iarm][ichg]->Draw("same");
		}
	}

	greco_ratio[0]->Draw("c 3");
	greco_ratio[1]->Draw("c 3");

	TLegend *leg = new TLegend(0.5,0.6,0.95,0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry("","Simulation for p+Au","h");
	leg->AddEntry(hPT_KPI[1][0],"Generation, K^{-}/#pi^{-}","L");
	leg->AddEntry(hPT_KPI[1][1],"Generation, K^{+}/#pi^{+}","L");
	leg->AddEntry(greco_ratio[0],"Reconstruction, K^{-}/#pi^{-}","LF");
	leg->AddEntry(greco_ratio[1],"Reconstruction, K^{+}/#pi^{+}","LF");
	leg->Draw();



}
