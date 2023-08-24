#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TProfile.h>
#include <TH2.h>

#include <fstream>
#include <iostream>

using namespace std;

void ScanEventTree01(){

	TFile *infileQA = new TFile("/home/shlim/Work/PHENIX/Code/dimuon/Run15pp200_badrun_for_hadron_20161231.root","read");
	TH1D *hbadrun[2];
	hbadrun[0] = (TH1D*)infileQA->Get("hBADRUN_arm0");
	hbadrun[1] = (TH1D*)infileQA->Get("hBADRUN_arm1");

	//TFile *infile = new TFile("/home/ojh7976/phenix_analysis/realevent_fvtx_mult/realoutput.root","read");
	//TProfile *hprofPre = (TProfile*)infile->Get("mult_allfvtx_both_1030_mean");

	int Run_Number;
	unsigned int trigbit_scaled;

	float bbcZ;
	float fvtxZ;

	int mult_fvtxN_eta1226, mult_fvtxN_eta1030, mult_fvtxN;
	int mult_fvtxS_eta1226, mult_fvtxS_eta1030, mult_fvtxS;
	int mult_svx, mult_svx_eta10;

	ifstream flist;
	flist.open("file.lst");

	char fname[500];

	TProfile *hprof_fvtxS = new TProfile("hprof_fvtxS","",30,-30,30);
	TProfile *hprof_fvtxS_eta1226 = new TProfile("hprof_fvtxS_eta1226","",30,-30,30);
	TProfile *hprof_fvtxS_eta1030 = new TProfile("hprof_fvtxS_eta1030","",30,-30,30);

	TProfile *hprof_fvtxN = new TProfile("hprof_fvtxN","",30,-30,30);
	TProfile *hprof_fvtxN_eta1226 = new TProfile("hprof_fvtxN_eta1226","",30,-30,30);
	TProfile *hprof_fvtxN_eta1030 = new TProfile("hprof_fvtxN_eta1030","",30,-30,30);

	TProfile *hprof_fvtxSN = new TProfile("hprof_fvtxSN","",30,-30,30);
	TProfile *hprof_fvtxSN_eta1226 = new TProfile("hprof_fvtxSN_eta1226","",30,-30,30);
	TProfile *hprof_fvtxSN_eta1030 = new TProfile("hprof_fvtxSN_eta1030","",30,-30,30);

	TProfile *hprof_svx = new TProfile("hprof_svx","",30,-30,30);
	TProfile *hprof_svx_eta10 = new TProfile("hprof_svx_eta10","",30,-30,30);

	/*
	TH2D *h2d_fvtxN = new TH2D("h2d_fvtxN","",20,-20,20,31,-0.5,30.5);
	TH2D *h2d_fvtxN_eta1226 = new TH2D("h2d_fvtxN_eta1226","",20,-20,20,31,-0.5,30.5);
	TH2D *h2d_fvtxN_eta1030 = new TH2D("h2d_fvtxN_eta1030","",20,-20,20,31,-0.5,30.5);

	TH2D *h2d_fvtxS = new TH2D("h2d_fvtxS","",20,-20,20,31,-0.5,30.5);
	TH2D *h2d_fvtxS_eta1226 = new TH2D("h2d_fvtxS_eta1226","",20,-20,20,31,-0.5,30.5);
	TH2D *h2d_fvtxS_eta1030 = new TH2D("h2d_fvtxS_eta1030","",20,-20,20,31,-0.5,30.5);

	TH2D *h2d_fvtxSN = new TH2D("h2d_fvtxSN","",20,-20,20,31,-0.5,30.5);
	TH2D *h2d_fvtxSN_eta1226 = new TH2D("h2d_fvtxSN_eta1226","",20,-20,20,31,-0.5,30.5);
	TH2D *h2d_fvtxSN_eta1030 = new TH2D("h2d_fvtxSN_eta1030","",20,-20,20,31,-0.5,30.5);

	TH2D *h2d_scaled_fvtxSN_eta1030 = new TH2D("h2d_scaled_fvtxSN_eta1030","",20,-20,20,100,0.0,10.0);
	*/

	while ( flist >> fname ){

		cout << "OPEN: " << fname << endl;

		TFile *infile = new TFile(fname,"read");

		TTree *T = (TTree*)infile->Get("evt");

		if ( !T ){
			infile->Close();
			delete infile;
			continue;
		}

		T->SetBranchAddress("Run_Number",&Run_Number);
		T->SetBranchAddress("trigbit_scaled",&trigbit_scaled);
		T->SetBranchAddress("bbcZ",&bbcZ);
		T->SetBranchAddress("fvtxZ",&fvtxZ);
		T->SetBranchAddress("mult_fvtxN_eta1226",&mult_fvtxN_eta1226);
		T->SetBranchAddress("mult_fvtxN_eta1030",&mult_fvtxN_eta1030);
		T->SetBranchAddress("mult_fvtxN",&mult_fvtxN);
		T->SetBranchAddress("mult_fvtxS_eta1226",&mult_fvtxS_eta1226);
		T->SetBranchAddress("mult_fvtxS_eta1030",&mult_fvtxS_eta1030);
		T->SetBranchAddress("mult_fvtxS",&mult_fvtxS);
		T->SetBranchAddress("mult_svx",&mult_svx);
		T->SetBranchAddress("mult_svx_eta10",&mult_svx_eta10);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);

			int bin = hbadrun[0]->FindBin(Run_Number);
			if ( hbadrun[0]->GetBinContent(bin) || hbadrun[1]->GetBinContent(bin) ) continue;

			if ( fabs(bbcZ)>=20 ) continue;

			hprof_fvtxS->Fill(bbcZ, mult_fvtxS);
			hprof_fvtxS_eta1226->Fill(bbcZ, mult_fvtxS_eta1226);
			hprof_fvtxS_eta1030->Fill(bbcZ, mult_fvtxS_eta1030);

			hprof_fvtxN->Fill(bbcZ, mult_fvtxN);
			hprof_fvtxN_eta1226->Fill(bbcZ, mult_fvtxN_eta1226);
			hprof_fvtxN_eta1030->Fill(bbcZ, mult_fvtxN_eta1030);

			hprof_fvtxSN->Fill(bbcZ, mult_fvtxS+mult_fvtxN);
			hprof_fvtxSN_eta1226->Fill(bbcZ, mult_fvtxS_eta1226+mult_fvtxN_eta1226);
			hprof_fvtxSN_eta1030->Fill(bbcZ, mult_fvtxS_eta1030+mult_fvtxN_eta1030);

			hprof_svx->Fill(bbcZ, mult_svx);
			hprof_svx_eta10->Fill(bbcZ, mult_svx_eta10);

			/*
			h2d_fvtxS->Fill(bbcZ, mult_fvtxS);
			h2d_fvtxS_eta1226->Fill(bbcZ, mult_fvtxS_eta1226);
			h2d_fvtxS_eta1030->Fill(bbcZ, mult_fvtxS_eta1030);

			h2d_fvtxN->Fill(bbcZ, mult_fvtxN);
			h2d_fvtxN_eta1226->Fill(bbcZ, mult_fvtxN_eta1226);
			h2d_fvtxN_eta1030->Fill(bbcZ, mult_fvtxN_eta1030);

			h2d_fvtxSN->Fill(bbcZ, mult_fvtxS+mult_fvtxN);
			h2d_fvtxSN_eta1226->Fill(bbcZ, mult_fvtxS_eta1226+mult_fvtxN_eta1226);
			h2d_fvtxSN_eta1030->Fill(bbcZ, mult_fvtxS_eta1030+mult_fvtxN_eta1030);

			float mean = hprofPre->GetBinContent(hprofPre->FindBin(bbcZ));

			h2d_scaled_fvtxSN_eta1030->Fill(bbcZ, (mult_fvtxS_eta1030+mult_fvtxN_eta1030)/mean);
			*/
		}

		//if ( fabs(bbcZ)>30.0 ) continue;


		infile->Close();
		delete infile;
		continue;

	}//

	flist.close();

	TFile *outfile = new TFile("outfileScanEventTree.root","recreate");
	hprof_fvtxS->Write();
	hprof_fvtxS_eta1226->Write();
	hprof_fvtxS_eta1030->Write();

	hprof_fvtxN->Write();
	hprof_fvtxN_eta1226->Write();
	hprof_fvtxN_eta1030->Write();

	hprof_fvtxSN->Write();
	hprof_fvtxSN_eta1226->Write();
	hprof_fvtxSN_eta1030->Write();

	hprof_svx->Write();
	hprof_svx_eta10->Write();

	/*
	h2d_fvtxS->Write();
	h2d_fvtxS_eta1226->Write();
	h2d_fvtxS_eta1030->Write();

	h2d_fvtxN->Write();
	h2d_fvtxN_eta1226->Write();
	h2d_fvtxN_eta1030->Write();

	h2d_fvtxSN->Write();
	h2d_fvtxSN_eta1226->Write();
	h2d_fvtxSN_eta1030->Write();

	h2d_scaled_fvtxSN_eta1030->Write();
	*/

	outfile->Close();


}
