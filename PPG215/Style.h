#include <TH1.h>
#include <TPad.h>
#include <TStyle.h>
#include <string>

TH1 *htmp;

void SetPadStyle(bool log = false){
	gPad->SetLeftMargin(0.20);
	if ( log ){
		gPad->SetRightMargin(0.12);
	}else{
		gPad->SetRightMargin(0.05);
	}
	gPad->SetTopMargin(0.07);
	gPad->SetBottomMargin(0.15);

	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);

	return;
}

void SetHistoStyle(std::string xtitle="", std::string ytitle="", std::string ztitle="",float size_title=24, float size_label=24){
	htmp->SetTitleFont(63);
	htmp->GetYaxis()->SetLabelFont(63);
	htmp->GetYaxis()->SetTitleFont(63);
	htmp->GetYaxis()->SetTitleOffset(1.3);
	//htmp->GetYaxis()->SetLabelFont(102);
	//htmp->GetYaxis()->SetTitleFont(102);
	htmp->GetYaxis()->SetLabelSize(0.9*size_label);
	htmp->GetYaxis()->SetTitleSize(size_title);
	htmp->GetYaxis()->SetTitle(ytitle.c_str());
	htmp->GetYaxis()->SetNdivisions(9,5,0);
	htmp->GetXaxis()->SetLabelFont(63);
	htmp->GetXaxis()->SetTitleFont(63);
	//htmp->GetXaxis()->SetLabelFont(102);
	//htmp->GetXaxis()->SetTitleFont(102);
	htmp->GetXaxis()->SetLabelSize(0.9*size_label);
	htmp->GetXaxis()->SetTitleSize(size_title);
	htmp->GetXaxis()->SetTitle(xtitle.c_str());
	//htmp->GetXaxis()->SetNdivisions(505);
	htmp->GetXaxis()->SetNdivisions(9,2,0);

	htmp->GetZaxis()->SetLabelFont(63);
	htmp->GetZaxis()->SetTitleFont(63);
	htmp->GetZaxis()->SetTitleOffset(1.3);
	//htmp->GetZaxis()->SetLabelSize(20);
	//htmp->GetZaxis()->SetTitleSize(22);
	htmp->GetZaxis()->SetLabelSize(size_label);
	htmp->GetZaxis()->SetTitleSize(size_title);
	htmp->GetZaxis()->SetTitle(ztitle.c_str());

	return;
}

