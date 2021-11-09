void ScanPythiaJpsi(int opt=0){

	ifstream flist;
	if ( opt==0 ){
		flist.open("file_pythiaMB_13TeV.lst");
	}else if ( opt==1 ){
		flist.open("file_pythiaCharm_13TeV.lst");
	}else if ( opt==2 ){
		flist.open("file_pythiaMB_200GeV.lst");
	}

	char fname[500];

	int np;
	int p_id[1000];
	float p_pt[1000], p_eta[1000], p_phi[1000];

	TH1D *h1d_svx[2];
	TH1D *h1d_scaled_svx[2];

	for (int ii=0; ii<2; ii++){
		h1d_svx[ii] = new TH1D(Form("h1d_svx_%d",ii),"",200,0,200);
		h1d_scaled_svx[ii] = new TH1D(Form("h1d_scaled_svx_%d",ii),"",20,0,10);
	}
	
	TH1D *hjpsi_rap = new TH1D("hjpsi_rap","",100,-5,5);

	TH1D *hmult_svx = new TH1D("hmult_svx","",1000,0,100);

	while ( flist >> fname ){

		TFile *infile = new TFile(fname,"read");
		cout << "OPEN: " << infile->GetName() << endl;

		TH1D *hevent = (TH1D*)infile->Get("hevent");
		TH1D *heta = (TH1D*)infile->Get("heta");

		float mb_mult_svx = heta->Integral(heta->FindBin(-1+0.001), heta->FindBin(+1-0.001));
		mb_mult_svx /= hevent->Integral();

		hmult_svx->Fill(mb_mult_svx);

		//cout << mb_mult_svx << endl;

		TTree *T = (TTree*)infile->Get("T");
		T->SetBranchAddress("np",&np);
		T->SetBranchAddress("p_id",p_id);
		T->SetBranchAddress("p_pt",p_pt);
		T->SetBranchAddress("p_eta",p_eta);
		T->SetBranchAddress("p_phi",p_phi);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);

			float jpsi_rap = p_eta[0];

			int iarm = -1;
			if ( jpsi_rap>1.2 && jpsi_rap<2.2 ) iarm = 1;
			else if ( jpsi_rap>-2.2 && jpsi_rap<-1.2 ) iarm = 0;
			else if ( jpsi_rap>2.5 && jpsi_rap<4.0 ) iarm = 2;
			else if ( jpsi_rap>-4.0 && jpsi_rap<-2.5 ) iarm = 2;

			if ( iarm<0 ) continue;

			int jpsi_index = -1;
			if ( p_id[0]==443 ) jpsi_index = 0;
			else if ( p_id[0]==100443 ) jpsi_index = 1;

			if ( jpsi_index<0 ) continue;

			hjpsi_rap->Fill(jpsi_rap);

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

			if ( opt==0 ){

				if ( iarm==2 ){
					h1d_svx[jpsi_index]->Fill(mult_svx/2.);
					h1d_scaled_svx[jpsi_index]->Fill(mult_svx/mb_mult_svx);
				}

			}else if ( opt==1 ){

				if ( iarm==0 || iarm==1 ){
					h1d_svx[jpsi_index]->Fill(mult_svx/2.);
					h1d_scaled_svx[jpsi_index]->Fill(mult_svx/mb_mult_svx);
				}

			}

		}//

		infile->Close();
		delete infile;

	}//

	TFile *outfile;

	if ( opt==0 ){
		outfile = new TFile("outfile_pythiaMB_pp13TeV.root","recreate");
	}else if ( opt==1 ){
		outfile = new TFile("outfile_pythiaCharm_pp13TeV.root","recreate");
	}else if ( opt==2 ){
		outfile = new TFile("outfile_pythiaMB_pp200GeV.root","recreate");
	}

	hjpsi_rap->Write();
	hmult_svx->Write();

	for (int ii=0; ii<2; ii++){
		h1d_svx[ii]->Write();
		h1d_scaled_svx[ii]->Write();
	}

	outfile->Close();


}
