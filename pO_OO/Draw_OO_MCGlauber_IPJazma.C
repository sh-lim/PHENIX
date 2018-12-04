#include "/phenix/u/shlim/Style.h"

void draw_hist_OO_IP(){

	gStyle->SetOptStat(0);

	const float cut_mult[2] = {68.25, 67.5};
	int N_event_mb[2] = {0.0};
	float N_eccG2[2] = {0.0}, N_eccG3[2] = {0.0}, N_eccG4[2] = {0.0};
	int N_event[2] = {0}, N_event2[2] = {0};
	float N_eccG2_0005[2] = {0.0}, N_eccG3_0005[2] = {0.0}, N_eccG4_0005[2] = {0.0};
	float N_eccJ2_0005[2] = {0.0}, N_eccJ3_0005[2] = {0.0}, N_eccJ4_0005[2] = {0.0};
	float N_mult[2] = {0.0}, N_mult2[2] = {0.0};

	const int nset = 2;

	TFile *infile1[2];
	infile1[0] = new TFile("outfile_OO_6370GeV.root","read");
	infile1[1] = new TFile("outfile_OO_7000GeV_IPJazma.root","read");

	TTree *Tout[2];

	TProfile *hprof_v2_b[nset];
	TProfile *hprof_v3_b[nset];
	TProfile *hprof_v4_b[nset];
	TProfile *hprof_v2_b2[nset];
	TProfile *hprof_v3_b2[nset];
	TProfile *hprof_v4_b2[nset];
	TProfile *hprof_v2_pT[nset];
	TProfile *hprof_v3_pT[nset];
	TProfile *hprof_v4_pT[nset];
	TProfile *hprof_v2_pT_0005[nset];
	TProfile *hprof_v3_pT_0005[nset];
	TProfile *hprof_v4_pT_0005[nset];
	TProfile *hprof_v2_pT_sp[nset];
	TProfile *hprof_v3_pT_sp[nset];
	TProfile *hprof_v4_pT_sp[nset];
	TProfile *hprof_v2_ecc[nset];
	TProfile *hprof_v3_ecc[nset];
	TProfile *hprof_v4_ecc[nset];
	TProfile *hprof_v2_ecc2[nset];
	TProfile *hprof_v3_ecc2[nset];
	TProfile *hprof_v4_ecc2[nset];
	TProfile *hprof_v2_mult[nset];
	TProfile *hprof_v3_mult[nset];
	TProfile *hprof_v4_mult[nset];
	TProfile *hprof_v2_mult2[nset];
	TProfile *hprof_v3_mult2[nset];
	TProfile *hprof_v4_mult2[nset];
	TProfile *hprof_ecc2_b[nset];
	TProfile *hprof_ecc3_b[nset];
	TProfile *hprof_ecc4_b[nset];

	TProfile *hprof_eccj2_b[nset];
	TProfile *hprof_eccj3_b[nset];
	TProfile *hprof_eccj4_b[nset];

	TProfile *hprof_ecc2_mult[nset];
	TProfile *hprof_ecc3_mult[nset];
	TProfile *hprof_ecc4_mult[nset];

	TProfile *hprof_v2_b_ecc[nset];
	TProfile *hprof_v3_b_ecc[nset];
	TProfile *hprof_v4_b_ecc[nset];

	TProfile *hprof_v2_pT_ecc[nset];
	TProfile *hprof_v3_pT_ecc[nset];
	TProfile *hprof_v4_pT_ecc[nset];
	TProfile *hprof_v2_pT_ecc_0005[nset];
	TProfile *hprof_v3_pT_ecc_0005[nset];
	TProfile *hprof_v4_pT_ecc_0005[nset];

	TProfile *hprof_mult_b[nset];
	TProfile *hprof_npart_b[nset];

	TProfile *hprof_mult_npart[nset];

	TH1D *hmult[nset];

	float b_B, b_Mult;
	float b_Npart;
	float b_Ecc2G, b_Ecc3G, b_Ecc4G;
	float b_Ecc2J, b_Ecc3J, b_Ecc4J;
	float b_pt[25], b_v2[25], b_v3[25], b_v4[25];

	for (int iset=0; iset<nset; iset++){
		Tout[iset] = (TTree*)infile1[iset]->Get("Tout");
		Tout[iset]->SetBranchAddress("B",&b_B);
		Tout[iset]->SetBranchAddress("Mult",&b_Mult);
		Tout[iset]->SetBranchAddress("Npart",&b_Npart);
		Tout[iset]->SetBranchAddress("Ecc2G",&b_Ecc2G);
		Tout[iset]->SetBranchAddress("Ecc3G",&b_Ecc3G);
		Tout[iset]->SetBranchAddress("Ecc4G",&b_Ecc4G);
		if ( iset==1 ){
			Tout[iset]->SetBranchAddress("Ecc2J",&b_Ecc2J);
			Tout[iset]->SetBranchAddress("Ecc3J",&b_Ecc3J);
			Tout[iset]->SetBranchAddress("Ecc4J",&b_Ecc4J);
		}
		Tout[iset]->SetBranchAddress("pT",b_pt);
		Tout[iset]->SetBranchAddress("v2",b_v2);
		Tout[iset]->SetBranchAddress("v3",b_v3);
		Tout[iset]->SetBranchAddress("v4",b_v4);

		hmult[iset] = new TH1D(Form("hmult_set%d",iset),"",100,0,100);

		hprof_v2_b[iset] = new TProfile(Form("hprof_v2_b_set%d",iset),"",12,0,12);
		hprof_v3_b[iset] = new TProfile(Form("hprof_v3_b_set%d",iset),"",12,0,12);
		hprof_v4_b[iset] = new TProfile(Form("hprof_v4_b_set%d",iset),"",12,0,12);
		hprof_v2_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_v2_b[iset]->SetMarkerColor(iset+1);
		hprof_v2_b[iset]->SetLineColor(iset+1);
		hprof_v3_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_v3_b[iset]->SetMarkerColor(iset+1);
		hprof_v3_b[iset]->SetLineColor(iset+1);
		hprof_v4_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_v4_b[iset]->SetMarkerColor(iset+1);
		hprof_v4_b[iset]->SetLineColor(iset+1);

		hprof_v2_b2[iset] = new TProfile(Form("hprof_v2_b2_set%d",iset),"",12,0,12);
		hprof_v3_b2[iset] = new TProfile(Form("hprof_v3_b2_set%d",iset),"",12,0,12);
		hprof_v4_b2[iset] = new TProfile(Form("hprof_v4_b2_set%d",iset),"",12,0,12);
		hprof_v2_b2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v2_b2[iset]->SetMarkerColor(iset+1);
		hprof_v2_b2[iset]->SetLineColor(iset+1);
		hprof_v3_b2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v3_b2[iset]->SetMarkerColor(iset+1);
		hprof_v3_b2[iset]->SetLineColor(iset+1);
		hprof_v4_b2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v4_b2[iset]->SetMarkerColor(iset+1);
		hprof_v4_b2[iset]->SetLineColor(iset+1);

		hprof_v2_pT[iset] = new TProfile(Form("hprof_v2_pT_set%d",iset),"",25,0,0.16*25);
		hprof_v3_pT[iset] = new TProfile(Form("hprof_v3_pT_set%d",iset),"",25,0,0.16*25);
		hprof_v4_pT[iset] = new TProfile(Form("hprof_v4_pT_set%d",iset),"",25,0,0.16*25);
		hprof_v2_pT[iset]->SetMarkerStyle(20+4*iset);
		hprof_v2_pT[iset]->SetMarkerColor(iset+1);
		hprof_v2_pT[iset]->SetLineColor(iset+1);
		hprof_v3_pT[iset]->SetMarkerStyle(20+4*iset);
		hprof_v3_pT[iset]->SetMarkerColor(iset+1);
		hprof_v3_pT[iset]->SetLineColor(iset+1);
		hprof_v4_pT[iset]->SetMarkerStyle(20+4*iset);
		hprof_v4_pT[iset]->SetMarkerColor(iset+1);
		hprof_v4_pT[iset]->SetLineColor(iset+1);

		hprof_v2_pT_0005[iset] = new TProfile(Form("hprof_v2_pT_0005_set%d",iset),"",25,0,0.16*25);
		hprof_v3_pT_0005[iset] = new TProfile(Form("hprof_v3_pT_0005_set%d",iset),"",25,0,0.16*25);
		hprof_v4_pT_0005[iset] = new TProfile(Form("hprof_v4_pT_0005_set%d",iset),"",25,0,0.16*25);
		hprof_v2_pT_0005[iset]->SetMarkerStyle(20+4*iset);
		hprof_v2_pT_0005[iset]->SetMarkerColor(iset+1);
		hprof_v2_pT_0005[iset]->SetLineColor(iset+1);
		hprof_v3_pT_0005[iset]->SetMarkerStyle(20+4*iset);
		hprof_v3_pT_0005[iset]->SetMarkerColor(iset+1);
		hprof_v3_pT_0005[iset]->SetLineColor(iset+1);
		hprof_v4_pT_0005[iset]->SetMarkerStyle(20+4*iset);
		hprof_v4_pT_0005[iset]->SetMarkerColor(iset+1);
		hprof_v4_pT_0005[iset]->SetLineColor(iset+1);

		hprof_v2_pT_sp[iset] = new TProfile(Form("hprof_v2_pT_sp_set%d",iset),"",25,0,0.16*25);
		hprof_v3_pT_sp[iset] = new TProfile(Form("hprof_v3_pT_sp_set%d",iset),"",25,0,0.16*25);
		hprof_v4_pT_sp[iset] = new TProfile(Form("hprof_v4_pT_sp_set%d",iset),"",25,0,0.16*25);
		hprof_v2_pT_sp[iset]->SetMarkerStyle(21+4*iset);
		hprof_v2_pT_sp[iset]->SetMarkerColor(4);
		hprof_v2_pT_sp[iset]->SetLineColor(4);
		hprof_v3_pT_sp[iset]->SetMarkerStyle(21+4*iset);
		hprof_v3_pT_sp[iset]->SetMarkerColor(4);
		hprof_v3_pT_sp[iset]->SetLineColor(4);
		hprof_v4_pT_sp[iset]->SetMarkerStyle(21+4*iset);
		hprof_v4_pT_sp[iset]->SetMarkerColor(4);
		hprof_v4_pT_sp[iset]->SetLineColor(4);

		hprof_v2_mult[iset] = new TProfile(Form("hprof_v2_mult_set%d",iset),"",10,0,60);
		hprof_v3_mult[iset] = new TProfile(Form("hprof_v3_mult_set%d",iset),"",10,0,60);
		hprof_v4_mult[iset] = new TProfile(Form("hprof_v4_mult_set%d",iset),"",10,0,60);
		hprof_v2_mult[iset]->SetMarkerStyle(20+4*iset);
		hprof_v2_mult[iset]->SetMarkerColor(iset+1);
		hprof_v2_mult[iset]->SetLineColor(iset+1);
		hprof_v3_mult[iset]->SetMarkerStyle(20+4*iset);
		hprof_v3_mult[iset]->SetMarkerColor(iset+1);
		hprof_v3_mult[iset]->SetLineColor(iset+1);
		hprof_v4_mult[iset]->SetMarkerStyle(20+4*iset);
		hprof_v4_mult[iset]->SetMarkerColor(iset+1);
		hprof_v4_mult[iset]->SetLineColor(iset+1);

		hprof_v2_mult2[iset] = new TProfile(Form("hprof_v2_mult2_set%d",iset),"",10,0,60);
		hprof_v3_mult2[iset] = new TProfile(Form("hprof_v3_mult2_set%d",iset),"",10,0,60);
		hprof_v4_mult2[iset] = new TProfile(Form("hprof_v4_mult2_set%d",iset),"",10,0,60);
		hprof_v2_mult2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v2_mult2[iset]->SetMarkerColor(iset+1);
		hprof_v2_mult2[iset]->SetLineColor(iset+1);
		hprof_v3_mult2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v3_mult2[iset]->SetMarkerColor(iset+1);
		hprof_v3_mult2[iset]->SetLineColor(iset+1);
		hprof_v4_mult2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v4_mult2[iset]->SetMarkerColor(iset+1);
		hprof_v4_mult2[iset]->SetLineColor(iset+1);

		hprof_v2_ecc[iset] = new TProfile(Form("hprof_v2_ecc_set%d",iset),"",10,0,1);
		hprof_v3_ecc[iset] = new TProfile(Form("hprof_v3_ecc_set%d",iset),"",10,0,1);
		hprof_v4_ecc[iset] = new TProfile(Form("hprof_v4_ecc_set%d",iset),"",10,0,1);
		hprof_v2_ecc[iset]->SetMarkerStyle(20+4*iset);
		hprof_v2_ecc[iset]->SetMarkerColor(iset+1);
		hprof_v2_ecc[iset]->SetLineColor(iset+1);
		hprof_v3_ecc[iset]->SetMarkerStyle(20+4*iset);
		hprof_v3_ecc[iset]->SetMarkerColor(iset+1);
		hprof_v3_ecc[iset]->SetLineColor(iset+1);
		hprof_v4_ecc[iset]->SetMarkerStyle(20+4*iset);
		hprof_v4_ecc[iset]->SetMarkerColor(iset+1);
		hprof_v4_ecc[iset]->SetLineColor(iset+1);

		hprof_v2_ecc2[iset] = new TProfile(Form("hprof_v2_ecc2_set%d",iset),"",10,0,1);
		hprof_v3_ecc2[iset] = new TProfile(Form("hprof_v3_ecc2_set%d",iset),"",10,0,1);
		hprof_v4_ecc2[iset] = new TProfile(Form("hprof_v4_ecc2_set%d",iset),"",10,0,1);
		hprof_v2_ecc2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v2_ecc2[iset]->SetMarkerColor(iset+1);
		hprof_v2_ecc2[iset]->SetLineColor(iset+1);
		hprof_v3_ecc2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v3_ecc2[iset]->SetMarkerColor(iset+1);
		hprof_v3_ecc2[iset]->SetLineColor(iset+1);
		hprof_v4_ecc2[iset]->SetMarkerStyle(21+4*iset);
		hprof_v4_ecc2[iset]->SetMarkerColor(iset+1);
		hprof_v4_ecc2[iset]->SetLineColor(iset+1);

		hprof_ecc2_b[iset] = new TProfile(Form("hprof_ecc2_b_set%d",iset),"",12,0,12);
		hprof_ecc3_b[iset] = new TProfile(Form("hprof_ecc3_b_set%d",iset),"",12,0,12);
		hprof_ecc4_b[iset] = new TProfile(Form("hprof_ecc4_b_set%d",iset),"",12,0,12);
		hprof_ecc2_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_ecc2_b[iset]->SetMarkerColor(iset+1);
		hprof_ecc2_b[iset]->SetLineColor(iset+1);
		hprof_ecc3_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_ecc3_b[iset]->SetMarkerColor(iset+1);
		hprof_ecc3_b[iset]->SetLineColor(iset+1);
		hprof_ecc4_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_ecc4_b[iset]->SetMarkerColor(iset+1);
		hprof_ecc4_b[iset]->SetLineColor(iset+1);

		hprof_eccj2_b[iset] = new TProfile(Form("hprof_eccj2_b_set%d",iset),"",12,0,12);
		hprof_eccj3_b[iset] = new TProfile(Form("hprof_eccj3_b_set%d",iset),"",12,0,12);
		hprof_eccj4_b[iset] = new TProfile(Form("hprof_eccj4_b_set%d",iset),"",12,0,12);
		hprof_eccj2_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_eccj2_b[iset]->SetMarkerColor(iset+1);
		hprof_eccj2_b[iset]->SetLineColor(iset+1);
		hprof_eccj3_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_eccj3_b[iset]->SetMarkerColor(iset+1);
		hprof_eccj3_b[iset]->SetLineColor(iset+1);
		hprof_eccj4_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_eccj4_b[iset]->SetMarkerColor(iset+1);
		hprof_eccj4_b[iset]->SetLineColor(iset+1);

		hprof_ecc2_mult[iset] = new TProfile(Form("hprof_ecc2_mult_set%d",iset),"",10,0,60);
		hprof_ecc3_mult[iset] = new TProfile(Form("hprof_ecc3_mult_set%d",iset),"",10,0,60);
		hprof_ecc4_mult[iset] = new TProfile(Form("hprof_ecc4_mult_set%d",iset),"",10,0,60);
		hprof_ecc2_mult[iset]->SetMarkerStyle(20+4*iset);
		hprof_ecc2_mult[iset]->SetMarkerColor(iset+1);
		hprof_ecc2_mult[iset]->SetLineColor(iset+1);
		hprof_ecc3_mult[iset]->SetMarkerStyle(20+4*iset);
		hprof_ecc3_mult[iset]->SetMarkerColor(iset+1);
		hprof_ecc3_mult[iset]->SetLineColor(iset+1);
		hprof_ecc4_mult[iset]->SetMarkerStyle(20+4*iset);
		hprof_ecc4_mult[iset]->SetMarkerColor(iset+1);
		hprof_ecc4_mult[iset]->SetLineColor(iset+1);

		hprof_mult_b[iset] = new TProfile(Form("hprof_mult_b_set%d",iset),"",12,0,12);
		hprof_mult_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_mult_b[iset]->SetMarkerColor(iset+1);
		hprof_mult_b[iset]->SetLineColor(iset+1);
		hprof_mult_b[iset]->SetErrorOption("i");

		hprof_npart_b[iset] = new TProfile(Form("hprof_npart_b_set%d",iset),"",12,0,12);
		hprof_npart_b[iset]->SetMarkerStyle(20+4*iset);
		hprof_npart_b[iset]->SetMarkerColor(iset+1);
		hprof_npart_b[iset]->SetLineColor(iset+1);
		hprof_npart_b[iset]->SetErrorOption("i");

		hprof_mult_npart[iset] = new TProfile(Form("hprof_mult_npart_set%d",iset),"",50,0,50);
		hprof_mult_npart[iset]->SetMarkerStyle(20+4*iset);
		hprof_mult_npart[iset]->SetMarkerColor(iset+1);
		hprof_mult_npart[iset]->SetLineColor(iset+1);
		hprof_mult_npart[iset]->SetErrorOption("i");

		int nentries = Tout[iset]->GetEntries();
		for (int ien=0; ien<nentries; ien++){

			Tout[iset]->GetEntry(ien);

			hmult[iset]->Fill(b_Mult);

			hprof_v2_b[iset]->Fill(b_B, b_v2[12], b_Mult);
			hprof_v3_b[iset]->Fill(b_B, b_v3[12], b_Mult);
			hprof_v4_b[iset]->Fill(b_B, b_v4[12], b_Mult);

			hprof_v2_b2[iset]->Fill(b_B, b_v2[6], b_Mult);
			hprof_v3_b2[iset]->Fill(b_B, b_v3[6], b_Mult);
			hprof_v4_b2[iset]->Fill(b_B, b_v4[6], b_Mult);

			hprof_v2_mult[iset]->Fill(b_Mult, b_v2[12], b_Mult);
			hprof_v3_mult[iset]->Fill(b_Mult, b_v3[12], b_Mult);
			hprof_v4_mult[iset]->Fill(b_Mult, b_v4[12], b_Mult);

			hprof_v2_mult2[iset]->Fill(b_Mult, b_v2[6], b_Mult);
			hprof_v3_mult2[iset]->Fill(b_Mult, b_v3[6], b_Mult);
			hprof_v4_mult2[iset]->Fill(b_Mult, b_v4[6], b_Mult);

			for (int ipt=0; ipt<25; ipt++){
				hprof_v2_pT[iset]->Fill(b_pt[ipt], b_v2[ipt], b_Mult);
				hprof_v3_pT[iset]->Fill(b_pt[ipt], b_v3[ipt], b_Mult);
				hprof_v4_pT[iset]->Fill(b_pt[ipt], b_v4[ipt], b_Mult);

				if ( b_Mult>cut_mult[iset] ){
					hprof_v2_pT_0005[iset]->Fill(b_pt[ipt], b_v2[ipt], b_Mult);
					hprof_v3_pT_0005[iset]->Fill(b_pt[ipt], b_v3[ipt], b_Mult);
					hprof_v4_pT_0005[iset]->Fill(b_pt[ipt], b_v4[ipt], b_Mult);

					if ( ipt==0 ){
						N_event[iset]++;
						N_eccG2_0005[iset] += b_Ecc2G;
						N_eccG3_0005[iset] += b_Ecc3G;
						N_eccG4_0005[iset] += b_Ecc4G;

						N_eccJ2_0005[iset] += b_Ecc2J;
						N_eccJ3_0005[iset] += b_Ecc3J;
						N_eccJ4_0005[iset] += b_Ecc4J;
						N_mult[iset] += b_Mult;
					}
				}

				if ( b_Mult>42.0 && b_Mult<53.5 ){
					hprof_v2_pT_sp[iset]->Fill(b_pt[ipt], b_v2[ipt], b_Mult);
					hprof_v3_pT_sp[iset]->Fill(b_pt[ipt], b_v3[ipt], b_Mult);
					hprof_v4_pT_sp[iset]->Fill(b_pt[ipt], b_v4[ipt], b_Mult);
					if ( ipt==0 ){
						N_event2[iset]++;
						N_mult2[iset] += b_Mult;
					}
				}
			}//ipt

			hprof_v2_ecc[iset]->Fill(b_Ecc2G, b_v2[12], b_Mult);
			hprof_v3_ecc[iset]->Fill(b_Ecc3G, b_v3[12], b_Mult);
			hprof_v4_ecc[iset]->Fill(b_Ecc4G, b_v4[12], b_Mult);

			hprof_v2_ecc2[iset]->Fill(b_Ecc2G, b_v2[6], b_Mult);
			hprof_v3_ecc2[iset]->Fill(b_Ecc3G, b_v3[6], b_Mult);
			hprof_v4_ecc2[iset]->Fill(b_Ecc4G, b_v4[6], b_Mult);

			hprof_ecc2_b[iset]->Fill(b_B, b_Ecc2G);
			hprof_ecc3_b[iset]->Fill(b_B, b_Ecc3G);
			hprof_ecc4_b[iset]->Fill(b_B, b_Ecc4G);

			hprof_eccj2_b[iset]->Fill(b_B, b_Ecc2J);
			hprof_eccj3_b[iset]->Fill(b_B, b_Ecc3J);
			hprof_eccj4_b[iset]->Fill(b_B, b_Ecc4J);

			hprof_ecc2_mult[iset]->Fill(b_Mult, b_Ecc2G);
			hprof_ecc3_mult[iset]->Fill(b_Mult, b_Ecc3G);
			hprof_ecc4_mult[iset]->Fill(b_Mult, b_Ecc4G);

			hprof_mult_b[iset]->Fill(b_B, b_Mult);
			hprof_npart_b[iset]->Fill(b_B, b_Npart);

			hprof_mult_npart[iset]->Fill(b_Npart, b_Mult);

			N_event_mb[iset]++;
			N_eccG2[iset] += b_Ecc2G;
			N_eccG3[iset] += b_Ecc3G;
			N_eccG4[iset] += b_Ecc4G;

		}//ien

		hprof_v2_b_ecc[iset] = (TProfile*)hprof_v2_b[iset]->Clone(Form("hprof_v2_b_ecc_set%d",iset));
		hprof_v3_b_ecc[iset] = (TProfile*)hprof_v3_b[iset]->Clone(Form("hprof_v3_b_ecc_set%d",iset));
		hprof_v4_b_ecc[iset] = (TProfile*)hprof_v4_b[iset]->Clone(Form("hprof_v4_b_ecc_set%d",iset));

		hprof_v2_b_ecc[iset]->Divide(hprof_ecc2_b[iset]);
		hprof_v3_b_ecc[iset]->Divide(hprof_ecc3_b[iset]);
		hprof_v4_b_ecc[iset]->Divide(hprof_ecc4_b[iset]);

		hprof_v2_pT_ecc[iset] = (TProfile*)hprof_v2_pT[iset]->Clone(Form("hprof_v2_pT_ecc_set%d",iset));
		hprof_v3_pT_ecc[iset] = (TProfile*)hprof_v3_pT[iset]->Clone(Form("hprof_v3_pT_ecc_set%d",iset));
		hprof_v4_pT_ecc[iset] = (TProfile*)hprof_v4_pT[iset]->Clone(Form("hprof_v4_pT_ecc_set%d",iset));
		hprof_v2_pT_ecc[iset]->Scale(1./(N_eccG2[iset]/N_event_mb[iset]));
		hprof_v3_pT_ecc[iset]->Scale(1./(N_eccG3[iset]/N_event_mb[iset]));
		hprof_v4_pT_ecc[iset]->Scale(1./(N_eccG4[iset]/N_event_mb[iset]));

		hprof_v2_pT_ecc_0005[iset] = (TProfile*)hprof_v2_pT_0005[iset]->Clone(Form("hprof_v2_pT_ecc_0005_set%d",iset));
		hprof_v3_pT_ecc_0005[iset] = (TProfile*)hprof_v3_pT_0005[iset]->Clone(Form("hprof_v3_pT_ecc_0005_set%d",iset));
		hprof_v4_pT_ecc_0005[iset] = (TProfile*)hprof_v4_pT_0005[iset]->Clone(Form("hprof_v4_pT_ecc_0005_set%d",iset));
		if ( iset==0 ){
			hprof_v2_pT_ecc_0005[iset]->Scale(1./(N_eccG2_0005[iset]/N_event[iset]));
			hprof_v3_pT_ecc_0005[iset]->Scale(1./(N_eccG3_0005[iset]/N_event[iset]));
			hprof_v4_pT_ecc_0005[iset]->Scale(1./(N_eccG4_0005[iset]/N_event[iset]));
		}else{
			hprof_v2_pT_ecc_0005[iset]->Scale(1./(N_eccJ2_0005[iset]/N_event[iset]));
			hprof_v3_pT_ecc_0005[iset]->Scale(1./(N_eccJ3_0005[iset]/N_event[iset]));
			hprof_v4_pT_ecc_0005[iset]->Scale(1./(N_eccJ4_0005[iset]/N_event[iset]));
		}


		cout << N_mult[iset]/N_event[iset] << " " << N_event[iset] << " " << N_mult2[iset]/N_event2[iset] << " " << N_event2[iset] << endl;
		cout << "mb e2 : " << N_eccG2[iset]/N_event_mb[iset] << ", e3: " << N_eccG3[iset]/N_event_mb[iset] << ", e4: " << N_eccG4[iset]/N_event_mb[iset] << endl;
		cout << "0-5% e2 : " << N_eccG2_0005[iset]/N_event[iset] << ", e3: " << N_eccG3_0005[iset]/N_event[iset] << ", e4: " << N_eccG4_0005[iset]/N_event[iset] << endl;
		cout << "0-5% je2 : " << N_eccJ2_0005[iset]/N_event[iset] << ", je3: " << N_eccJ3_0005[iset]/N_event[iset] << ", je4: " << N_eccJ4_0005[iset]/N_event[iset] << endl;
	}//

	//return;


	/*
	TCanvas *c2 = new TCanvas("c2","c2",1.1*3*350,2*350);
	c2->Divide(3,2);
	c2->cd(1);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.25);
	SetHistoStyle("b (fm)", "v_{2}");
	hprof_v2_b[0]->Draw("p same");
	hprof_v2_b[1]->Draw("p same");

	TLegend *leg = new TLegend(0.4,0.75,0.9,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry("","O+O 7 TeV","");
	leg->AddEntry(hprof_v2_b[0],"Full, p_{T}=2 GeV/c","P");
	leg->AddEntry(hprof_v2_b[1],"WS, p_{T}=2 GeV/c","P");
	leg->Draw();

	c2->cd(2);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.1);
	SetHistoStyle("b (fm)", "v_{3}");
	hprof_v3_b[0]->Draw("p same");
	hprof_v3_b[1]->Draw("p same");
	leg->Draw();

	c2->cd(3);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.05);
	SetHistoStyle("b (fm)", "v_{4}");
	hprof_v4_b[0]->Draw("p same");
	hprof_v4_b[1]->Draw("p same");
	leg->Draw();

	c2->cd(4);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.70);
	SetHistoStyle("b (fm)", "v_{2}/<#varepsilon_{2}>");
	hprof_v2_b_ecc[0]->Draw("p same");
	hprof_v2_b_ecc[1]->Draw("p same");

	c2->cd(5);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.50);
	SetHistoStyle("b (fm)", "v_{3}/<#varepsilon_{3}>");
	hprof_v3_b_ecc[0]->Draw("p same");
	hprof_v3_b_ecc[1]->Draw("p same");

	c2->cd(6);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.20);
	SetHistoStyle("b (fm)", "v_{4}/<#varepsilon_{4}>");
	hprof_v4_b_ecc[0]->Draw("p same");
	hprof_v4_b_ecc[1]->Draw("p same");
	*/

	/*
	TCanvas *c7 = new TCanvas("c7","c7",1.1*2*350,1*350);
	c7->Divide(2,1);

	c7->cd(1);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,100);
	SetHistoStyle("b (fm)","<dN_{ch}/dy>");
	hprof_mult_b[0]->Draw("p same");
	hprof_mult_b[1]->Draw("p same");

	TLegend *leg = new TLegend(0.3,0.7,0.9,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(hprof_mult_b[0],"Full","P");
	leg->AddEntry(hprof_mult_b[1],"WS","P");
	leg->Draw();

	c7->cd(2);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,50);
	SetHistoStyle("b (fm)","<N_{part}>");
	hprof_npart_b[0]->Draw("p same");
	hprof_npart_b[1]->Draw("p same");
	leg->Draw();

	//return;
	
	for (int ibin=0; ibin<hmult[1]->GetNbinsX(); ibin++){
		cout << ibin+1 << ", " << hmult[1]->Integral(1,ibin+1)/hmult[1]->Integral() << endl;
	}
	*/

	TCanvas *c3 = new TCanvas("c3","c3",1.1*3*350,350);
	c3->Divide(3,1);
	c3->cd(1);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.7);
	SetHistoStyle("b (fm)", "#LT#varepsilon_{2}#GT");
	hprof_ecc2_b[0]->Draw("p same");
	hprof_ecc2_b[1]->Draw("p same");
	hprof_eccj2_b[1]->SetMarkerColor(4);
	hprof_eccj2_b[1]->SetLineColor(4);
	hprof_eccj2_b[1]->Draw("p same");

	c3->cd(2);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.7);
	SetHistoStyle("b (fm)", "#LT#varepsilon_{3}#GT");
	hprof_ecc3_b[0]->Draw("p same");
	hprof_ecc3_b[1]->Draw("p same");
	hprof_eccj3_b[1]->SetMarkerColor(4);
	hprof_eccj3_b[1]->SetLineColor(4);
	hprof_eccj3_b[1]->Draw("p same");

	TLegend *leg = new TLegend(0.2,0.7,0.9,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(hprof_ecc3_b[0],"eccgaus, Sanghoon","p");
	leg->AddEntry(hprof_ecc3_b[1],"eccgaus, Jamie","p");
	leg->AddEntry(hprof_eccj3_b[1],"eccjazma, Jamie","p");
	leg->Draw();

	c3->cd(3);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,12,0.7);
	SetHistoStyle("b (fm)", "#LT#varepsilon_{4}#GT");
	hprof_ecc4_b[0]->Draw("p same");
	hprof_ecc4_b[1]->Draw("p same");
	hprof_eccj4_b[1]->SetMarkerColor(4);
	hprof_eccj4_b[1]->SetLineColor(4);
	hprof_eccj4_b[1]->Draw("p same");

	//return;
	
	TCanvas *c1 = new TCanvas("c1","c1",1.1*3*350,1*350);
	c1->Divide(3,1);

	hmult[0]->Rebin(2);
	hmult[1]->Rebin(2);

	c1->cd(1);
	SetPadStyle();
	gPad->SetLogy();
	htmp = (TH1D*)gPad->DrawFrame(0,1,100,2*hmult[0]->GetMaximum());
	SetHistoStyle("dN_{ch}/dy","");
	hmult[0]->SetLineColor(1);
	hmult[0]->Draw("same");
	hmult[1]->SetLineColor(2);
	hmult[1]->Draw("same");

	TLegend *leg = new TLegend(0.1,0.75,0.9,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry("","O+O 7 TeV","");
	leg->AddEntry(hmult[0],Form("MC Glauber, <dN_{ch}/dy> MB=%4.1f, 0-5%c=%4.1f",hmult[0]->GetMean(),'%',N_mult[0]/N_event[0]),"");
	leg->AddEntry(hmult[1],Form("IP-JAZMA, <dN_{ch}/dy> MB=%4.1f, 0-5%c=%4.1f",hmult[1]->GetMean(),'%',N_mult[1]/N_event[1]),"");
	leg->Draw();

	c1->cd(2);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.25);
	SetHistoStyle("p_{T} (GeV/c)","v_{n}");
	//htmp->GetXaxis()->SetTitleOffset(1.8);
	//htmp->GetYaxis()->SetTitleOffset(1.8);

	for (int iset=0; iset<nset; iset++){
		hprof_v2_pT_0005[iset]->SetMarkerColor(1);
		hprof_v2_pT_0005[iset]->SetLineColor(1);
		hprof_v3_pT_0005[iset]->SetMarkerColor(2);
		hprof_v3_pT_0005[iset]->SetLineColor(2);
		hprof_v4_pT_0005[iset]->SetMarkerColor(4);
		hprof_v4_pT_0005[iset]->SetLineColor(4);
		hprof_v2_pT_0005[iset]->Draw("p same");
		hprof_v3_pT_0005[iset]->Draw("p same");
		hprof_v4_pT_0005[iset]->Draw("p same");
	}

	TLegend *leg = new TLegend(0.20,0.7,0.95,0.9);
	leg->SetFillStyle(0);
	leg->SetNColumns(3);
	leg->SetHeader("O+O 7 TeV, 0-5%");
	leg->AddEntry(hprof_v2_pT_0005[0],"MC Glauber, v_{2}","P");
	leg->AddEntry(hprof_v3_pT_0005[0],"MC Glauber, v_{3}","P");
	leg->AddEntry(hprof_v4_pT_0005[0],"MC Glauber, v_{4}","P");
	leg->AddEntry(hprof_v2_pT_0005[1],"IP-JAZMA, v_{2}","P");
	leg->AddEntry(hprof_v3_pT_0005[1],"IP-JAZMA, v_{3}","P");
	leg->AddEntry(hprof_v4_pT_0005[1],"IP-JAZMA, v_{4}","P");
	leg->Draw();

	c1->cd(3);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.7);
	SetHistoStyle("p_{T} (GeV/c)","v_{n}/#LT#varepsilon_{n}#GT");
	//htmp->GetXaxis()->SetTitleOffset(1.8);
	//htmp->GetYaxis()->SetTitleOffset(1.8);

	for (int iset=0; iset<nset; iset++){
		hprof_v2_pT_ecc_0005[iset]->SetMarkerColor(1);
		hprof_v2_pT_ecc_0005[iset]->SetLineColor(1);
		hprof_v3_pT_ecc_0005[iset]->SetMarkerColor(2);
		hprof_v3_pT_ecc_0005[iset]->SetLineColor(2);
		hprof_v4_pT_ecc_0005[iset]->SetMarkerColor(4);
		hprof_v4_pT_ecc_0005[iset]->SetLineColor(4);
		hprof_v2_pT_ecc_0005[iset]->Draw("p same");
		hprof_v3_pT_ecc_0005[iset]->Draw("p same");
		hprof_v4_pT_ecc_0005[iset]->Draw("p same");
	}


	return;

	c1->cd(5);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,3,1.0);
	SetHistoStyle("p_{T} (GeV/c)","v_{n}/<#varepsilon_{n}>");

	for (int iset=0; iset<nset; iset++){
		hprof_v2_pT_ecc[iset]->SetMarkerColor(1);
		hprof_v2_pT_ecc[iset]->SetLineColor(1);
		hprof_v3_pT_ecc[iset]->SetMarkerColor(2);
		hprof_v3_pT_ecc[iset]->SetLineColor(2);
		hprof_v4_pT_ecc[iset]->SetMarkerColor(4);
		hprof_v4_pT_ecc[iset]->SetLineColor(4);
		hprof_v2_pT_ecc[iset]->Draw("p same");
		hprof_v3_pT_ecc[iset]->Draw("p same");
		hprof_v4_pT_ecc[iset]->Draw("p same");
	}

	TLegend *leg = new TLegend(0.25,0.65,0.93,0.9);
	leg->SetFillStyle(0);
	leg->SetNColumns(3);
	leg->SetHeader("O+O 7 TeV, 0-100%");
	leg->AddEntry(hprof_v2_pT[0],"Full, v_{2}/<#varepsilon_{2}>","P");
	leg->AddEntry(hprof_v3_pT[0],"Full, v_{3}/<#varepsilon_{3}>","P");
	leg->AddEntry(hprof_v4_pT[0],"Full, v_{4}/<#varepsilon_{4}>","P");
	leg->AddEntry(hprof_v2_pT[1],"WS, v_{2}/<#varepsilon_{2}>","P");
	leg->AddEntry(hprof_v3_pT[1],"WS, v_{3}/<#varepsilon_{3}>","P");
	leg->AddEntry(hprof_v4_pT[1],"WS, v_{4}/<#varepsilon_{4}>","P");
	leg->Draw();

	c1->cd(6);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,3,1.0);
	SetHistoStyle("p_{T} (GeV/c)","v_{n}/<#varepsilon_{n}>");

	for (int iset=0; iset<nset; iset++){
		hprof_v2_pT_ecc_0005[iset]->SetMarkerColor(1);
		hprof_v2_pT_ecc_0005[iset]->SetLineColor(1);
		hprof_v3_pT_ecc_0005[iset]->SetMarkerColor(2);
		hprof_v3_pT_ecc_0005[iset]->SetLineColor(2);
		hprof_v4_pT_ecc_0005[iset]->SetMarkerColor(4);
		hprof_v4_pT_ecc_0005[iset]->SetLineColor(4);
		hprof_v2_pT_ecc_0005[iset]->Draw("p same");
		hprof_v3_pT_ecc_0005[iset]->Draw("p same");
		hprof_v4_pT_ecc_0005[iset]->Draw("p same");
	}

	TLegend *leg = new TLegend(0.25,0.65,0.93,0.9);
	leg->SetFillStyle(0);
	leg->SetNColumns(3);
	leg->SetHeader("O+O 7 TeV, 0-5%");
	leg->AddEntry(hprof_v2_pT[0],"Full, v_{2}/<#varepsilon_{2}>","P");
	leg->AddEntry(hprof_v3_pT[0],"Full, v_{3}/<#varepsilon_{3}>","P");
	leg->AddEntry(hprof_v4_pT[0],"Full, v_{4}/<#varepsilon_{4}>","P");
	leg->AddEntry(hprof_v2_pT[1],"WS, v_{2}/<#varepsilon_{2}>","P");
	leg->AddEntry(hprof_v3_pT[1],"WS, v_{3}/<#varepsilon_{3}>","P");
	leg->AddEntry(hprof_v4_pT[1],"WS, v_{4}/<#varepsilon_{4}>","P");
	leg->Draw();

	TCanvas *c8 = new TCanvas("c8","c8",1.1*500,500);
	SetPadStyle();
	gPad->SetLogy();
	gPad->SetLeftMargin(0.15);

	TH1D *hmult_sonic = (TH1D*)hmult[0]->Clone("hmult_sonic");
	hmult_sonic->SetLineWidth(2);
	hmult_sonic->Scale(1./hmult_sonic->Integral());

	TH1D *hmult_sonic_0005 = (TH1D*)hmult_sonic->Clone("hmult_sonic_0005");
	for (int ibin=0; ibin<hmult_sonic_0005->GetNbinsX(); ibin++){
		float val = hmult_sonic->Integral(1,ibin+1);
		if ( val<0.95 )
			hmult_sonic_0005->SetBinContent(ibin+1, 0);
	}

	htmp = (TH1D*)gPad->DrawFrame(0,1,300,2*hmult_sonic->GetMaximum());
	SetHistoStyle("<dN_{ch}/dy>");
	hmult_sonic_0005->SetFillStyle(3001);
	hmult_sonic_0005->SetFillColor(15);
	hmult_sonic_0005->SetLineWidth(0);
	hmult_sonic_0005->Draw("same");
	hmult_sonic->Draw("same");

	TFile *infile00 = new TFile("oo_7tev_1.0.root","read");
	TH1D *hmult_ampt = (TH1D*)infile00->Get("h_mult");
	hmult_ampt->SetLineWidth(2);
	hmult_ampt->SetLineColor(2);
	hmult_ampt->Scale(1./hmult_ampt->Integral());

	TH1D *hmult_ampt_0005 = (TH1D*)hmult_ampt->Clone("hmult_ampt_0005");
	for (int ibin=0; ibin<hmult_ampt_0005->GetNbinsX(); ibin++){
		float val = hmult_ampt->Integral(1,ibin+1);
		if ( val<0.95 )
			hmult_ampt_0005->SetBinContent(ibin+1, 0);
	}

	hmult_ampt_0005->SetFillStyle(3001);
	hmult_ampt_0005->SetFillColor(kRed-7);
	hmult_ampt_0005->SetLineColor(kRed-7);
	hmult_ampt_0005->SetLineWidth(0);
	hmult_ampt_0005->Draw("same");
	hmult_ampt->Draw("same");

	TLegend *leg = new TLegend(0.35,0.6,0.9,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry("","O+O 7 TeV","");
	leg->AddEntry(hmult_sonic,Form("SONIC, 0-100%c, mean=%4.1f",'%',hmult_sonic->GetMean()),"L");
	leg->AddEntry(hmult_sonic_0005,Form("SONIC, 0-5%c, mean=%4.1f",'%',hmult_sonic_0005->GetMean()),"F");
	leg->AddEntry(hmult_ampt,Form("AMPT, 0-100%c, mean=%4.1f",'%',hmult_ampt->GetMean()),"L");
	leg->AddEntry(hmult_ampt_0005,Form("AMPT, 0-5%c, mean=%4.1f",'%',hmult_ampt_0005->GetMean()),"F");
	leg->Draw();

	TCanvas *c9 = new TCanvas("c9","c9",1.1*400,400);
	SetPadStyle();

	htmp = (TH1D*)gPad->DrawFrame(0,0,40,200);
	SetHistoStyle("<N_{part}>","<dN_{ch}/dy>");

	hprof_mult_npart[0]->Draw("p same");

	TF1 *fmult = new TF1("fmult","5.8*(x-2)+4.7",0,50);
	fmult->Draw("same");

	TGraphErrors *gpPb = new TGraphErrors();
	gpPb->SetPoint(0, 8.02, 20.10);
	gpPb->SetMarkerStyle(25);
	gpPb->SetMarkerSize(1.5);
	gpPb->Draw("p same");

	TGraphErrors *gXeXe = new TGraphErrors();
	gXeXe->SetPoint(0, 5.14, 13.3);
	gXeXe->SetPoint(1, 10.5, 32.0);
	gXeXe->SetPoint(3, 19.7, 64.7);
	gXeXe->SetPoint(4, 34.1, 118.0);
	gXeXe->SetMarkerStyle(27);
	gXeXe->SetMarkerSize(2.0);
	gXeXe->Draw("p same");

	TLegend *leg = new TLegend(0.25,0.7,0.7,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(hprof_mult_npart[0],"O+O 7 TeV, SONIC","P");
	leg->AddEntry(gXeXe,"Xe+Xe 5.44 TeV, ALICE","P");
	leg->AddEntry(gpPb,"p+Pb 8.16 TeV, CMS","P");
	leg->AddEntry(fmult,"5.8#times(N_{part}-2)+4.7","L");
	leg->Draw();

	return;


	{
		c1->cd(3*iset+2);


		c1->cd(3*iset+3);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,0,3.5,0.25);
		SetHistoStyle("p_{T} (GeV/c)","v_{n}");
		hprof_v2_pT_0005[iset]->SetMarkerColor(1);
		hprof_v2_pT_0005[iset]->SetLineColor(1);
		hprof_v3_pT_0005[iset]->SetMarkerColor(2);
		hprof_v3_pT_0005[iset]->SetLineColor(2);
		hprof_v4_pT_0005[iset]->SetMarkerColor(4);
		hprof_v4_pT_0005[iset]->SetLineColor(4);
		hprof_v2_pT_0005[iset]->Draw("p same");
		hprof_v3_pT_0005[iset]->Draw("p same");
		hprof_v4_pT_0005[iset]->Draw("p same");

		TLegend *leg = new TLegend(0.2,0.65,0.8,0.9);
		leg->SetFillStyle(0);
		leg->AddEntry("","O+O 7 TeV, 0-5%","");
		leg->AddEntry(hprof_v2_pT[0],"v_{2}","P");
		leg->AddEntry(hprof_v3_pT[0],"v_{3}","P");
		leg->AddEntry(hprof_v4_pT[0],"v_{4}","P");
		leg->Draw();

	}

	//return;

	//return;

	//return;

	/*
	TCanvas *c4 = new TCanvas("c4","c4",1.1*3*350,350);
	c4->Divide(3,1);
	c4->cd(1);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,1,0.25);
	SetHistoStyle("<#varepsilon_{2}>","v_{2}, p_{T}=2 GeV/c");
	hprof_v2_ecc[0]->Draw("p same");
	hprof_v2_ecc[1]->Draw("p same");
	hprof_v2_ecc2[0]->Draw("p same");
	hprof_v2_ecc2[1]->Draw("p same");

	c4->cd(2);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,1,0.1);
	SetHistoStyle("<#varepsilon_{3}>","v_{3}, p_{T}=2 GeV/c");
	hprof_v3_ecc[0]->Draw("p same");
	hprof_v3_ecc[1]->Draw("p same");
	hprof_v3_ecc2[0]->Draw("p same");
	hprof_v3_ecc2[1]->Draw("p same");

	c4->cd(3);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,1,0.05);
	SetHistoStyle("<#varepsilon_{4}>","v_{4}, p_{T}=2 GeV/c");
	hprof_v4_ecc[0]->Draw("p same");
	hprof_v4_ecc[1]->Draw("p same");
	hprof_v4_ecc2[0]->Draw("p same");
	hprof_v4_ecc2[1]->Draw("p same");
	*/

	/*
	TCanvas *c5 = new TCanvas("c5","c5",1.1*3*350,350);
	c5->Divide(3,1);
	c5->cd(1);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,60,0.25);
	SetHistoStyle("dN_{ch}/dy","v_{2}");
	hprof_v2_mult[0]->Draw("p same");
	hprof_v2_mult2[0]->Draw("p same");

	TLegend *leg = new TLegend(0.2,0.75,0.7,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry("","O+O 7 TeV","");
	leg->AddEntry(hprof_v2_mult2[0],"p_{T}=1 GeV/c","P");
	leg->AddEntry(hprof_v2_mult[0],"p_{T}=2 GeV/c","P");
	leg->Draw();

	c5->cd(2);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,60,0.1);
	SetHistoStyle("dN_{ch}/dy","v_{3}");
	hprof_v3_mult[0]->Draw("p same");
	hprof_v3_mult2[0]->Draw("p same");
	leg->Draw();

	c5->cd(3);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,60,0.05);
	SetHistoStyle("dN_{ch}/dy","v_{4}");
	hprof_v4_mult[0]->Draw("p same");
	hprof_v4_mult2[0]->Draw("p same");
	leg->Draw();
	*/
	
	//return;
	//
	TCanvas *c7 = new TCanvas("c7","c7",1.1*3*350,350);
	c7->Divide(3,1);
	c7->cd(1);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,60,0.7);
	SetHistoStyle("dN_{ch}/dy", "<#varepsilon_{2}>");
	hprof_ecc2_mult[0]->Draw("p same");
	hprof_ecc2_mult[1]->Draw("p same");

	c7->cd(2);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,60,0.7);
	SetHistoStyle("dN_{ch}/dy", "<#varepsilon_{3}>");
	hprof_ecc3_mult[0]->Draw("p same");
	hprof_ecc3_mult[1]->Draw("p same");

	c7->cd(3);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,60,0.7);
	SetHistoStyle("dN_{ch}/dy", "<#varepsilon_{4}>");
	hprof_ecc4_mult[0]->Draw("p same");
	hprof_ecc4_mult[1]->Draw("p same");

}
