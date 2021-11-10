void ScanPythiaMB(){

	ifstream flist;
	flist.open("file_pythiaMB.lst");

	char fname[500];

	int np;
	int p_id[1000];
	float p_pt[1000], p_eta[1000], p_phi[1000];

	TH2D *h2d_fvtx_bbc[2];
	TH2D *h2d_fvtx_bbc2[2];
	TH2D *h2d_fvtx_svx[2];
	TH2D *h2d_fvtx_fvtx[2];

	for (int iarm=0; iarm<2; iarm++){
		h2d_fvtx_bbc[iarm] = new TH2D(Form("h2d_fvtx_bbc_arm%d",iarm),"",50,0,50,50,0,50);
		h2d_fvtx_bbc2[iarm] = new TH2D(Form("h2d_fvtx_bbc2_arm%d",iarm),"",50,0,50,50,0,50);
		h2d_fvtx_svx[iarm] = new TH2D(Form("h2d_fvtx_svx_arm%d",iarm),"",50,0,50,50,0,50);
		h2d_fvtx_fvtx[iarm] = new TH2D(Form("h2d_fvtx_fvtx_arm%d",iarm),"",50,0,50,50,0,50);
	}
	
	while ( flist >> fname ){

		TFile *infile = new TFile(fname,"read");
		cout << "OPEN: " << infile->GetName() << endl;

		TTree *T = (TTree*)infile->Get("T");
		T->SetBranchAddress("np",&np);
		T->SetBranchAddress("p_id",p_id);
		T->SetBranchAddress("p_pt",p_pt);
		T->SetBranchAddress("p_eta",p_eta);
		T->SetBranchAddress("p_phi",p_phi);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);

			int mult_svx = 0;
			int mult_fvtxs = 0;
			int mult_fvtxn = 0;
			int mult_bbcs = 0;
			int mult_bbcn = 0;

			for (int ip=0; ip<np; ip++){

				int pid = p_id[ip];
				float eta = p_eta[ip];

				if ( !(abs(pid)==211 || abs(pid)==321 || abs(pid)==2212) ) continue;

				if ( fabs(eta)<1.0 ){
					mult_svx++;
				}else if ( eta>1.2 && eta<2.6 ){
					mult_fvtxn++;
				}else if ( eta>3.1 && eta<3.9 ){
					mult_bbcn++;
				}else if ( eta>-2.6 && eta<-1.2 ){
					mult_fvtxs++;
				}else if ( eta>-3.9 && eta<-3.1 ){
					mult_bbcs++;
				}
			}//ip

			if ( mult_bbcs<1 || mult_bbcn<1 ) continue;

			h2d_fvtx_bbc[0]->Fill(mult_fvtxs, mult_bbcs);
			h2d_fvtx_bbc2[0]->Fill(mult_fvtxs, mult_bbcn);
			h2d_fvtx_svx[0]->Fill(mult_fvtxs, mult_svx);
			h2d_fvtx_fvtx[0]->Fill(mult_fvtxs, mult_fvtxn);

			h2d_fvtx_bbc[1]->Fill(mult_fvtxn, mult_bbcn);
			h2d_fvtx_bbc2[1]->Fill(mult_fvtxn, mult_bbcs);
			h2d_fvtx_svx[1]->Fill(mult_fvtxn, mult_svx);
			h2d_fvtx_fvtx[1]->Fill(mult_fvtxn, mult_fvtxs);

		}//

		infile->Close();
		delete infile;

	}//

	TFile *outfile = new TFile("outfile_pythia_mb.root","recreate");

	for (int iarm=0; iarm<2; iarm++){
		h2d_fvtx_bbc[iarm]->Write();
		h2d_fvtx_bbc2[iarm]->Write();
		h2d_fvtx_svx[iarm]->Write();
		h2d_fvtx_fvtx[iarm]->Write();
	}

	outfile->Close();


}
