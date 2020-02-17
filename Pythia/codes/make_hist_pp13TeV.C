#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>

using namespace std;

void make_hist_pp13TeV(const char *fname="file.lst"){

	//CMS pt range
	//const float pt_min = 0.3;
	//const float pt_max = 3.0;
	//ATLAS pt range
	const float pt_min = 0.5;
	const float pt_max = 5.0;
	const float eta_max = 2.5;
	const float deta_max = 5.0;
	//const float eta_max = 0.9;
	//const float deta_max = 1.8;
	
	const float deta_cut = 2.0;
	const float const_pi = TMath::Pi();

	const int nmult = 20;

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int npart;
	int part_pid[1000];
	float part_eta[1000], part_phi[1000], part_pt[1000];

	TH2D *h2d_same_dphi_deta[nmult];
	TH2D *h2d_mixed_dphi_deta[nmult];

	TProfile *hv22pc_mult = new TProfile("hv22pc_mult","",nmult,0,10*nmult,-1,1);

	for (int im=0; im<nmult; im++){
		h2d_same_dphi_deta[im] = new TH2D(Form("h2d_same_dphi_deta_m%02d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		h2d_mixed_dphi_deta[im] = new TH2D(Form("h2d_mixed_dphi_deta_m%02d",im),"",96,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		//h2d_same_dphi_deta[im] = new TH2D(Form("h2d_same_dphi_deta_m%02d",im),"",30*4,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
		//h2d_mixed_dphi_deta[im] = new TH2D(Form("h2d_mixed_dphi_deta_m%02d",im),"",30*4,-const_pi/2,3*const_pi/2,100,-5.0,5.0);
	}

	TH1D *hmult_all = new TH1D("hmult_all","",500,0,500);
	TH1D *hmult_mid = new TH1D("hmult_mid","",500,0,500);
	TH1D *hmult_fwd = new TH1D("hmult_fwd","",500,0,500);
	TH1D *hntrig_mid = new TH1D("hntrig_mid","",200,0,200);
	TH1D *hntrig_mixed_mid = new TH1D("hntrig_mixed_mid","",200,0,200);

	TH2D *heta_pt_all = new TH2D("heta_pt_all","",100,-5,5,100,0,10);
	TH2D *heta_pt_ana = new TH2D("heta_pt_ana","",100,-5,5,100,0,10);

	vector<float> vec_pt[nmult];
	vector<float> vec_phi[nmult];
	vector<float> vec_eta[nmult];

	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("T");
		T->SetBranchAddress("np",&npart);
		T->SetBranchAddress("part_pid",part_pid);
		T->SetBranchAddress("part_eta",part_eta);
		T->SetBranchAddress("part_phi",part_phi);
		T->SetBranchAddress("part_pt",part_pt);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			hmult_all->Fill(npart);

			int nmult_mid = 0, nmult_fwd = 0;

			//count multiplicity for bin
			for (int ip=0; ip<npart; ip++){
				if ( fabs(part_eta[ip])<2.5 && part_pt[ip]>0.4 ){
					nmult_mid++;
				}
				heta_pt_all->Fill(part_eta[ip], part_pt[ip]);
			}

			hmult_mid->Fill(nmult_mid);

			int imult = int(nmult_mid/10);
			if ( imult>=nmult ) continue;

			for (int ip=0; ip<npart; ip++){
				//int ip_ptbin = int(part_pt[ip]/0.5);
				//if ( ip_ptbin>=20 ) continue;

				if ( part_pt[ip]<pt_min || part_pt[ip]>pt_max ) continue;
				if ( fabs(part_eta[ip])>eta_max ) continue;

				hntrig_mid->Fill(nmult_mid);
				heta_pt_ana->Fill(part_eta[ip], part_pt[ip]);

				for (int jp=0; jp<npart; jp++){
					if ( ip==jp ) continue;
					if ( part_pt[jp]<pt_min || part_pt[jp]>pt_max ) continue;
					if ( fabs(part_eta[jp])>eta_max ) continue;

					float dphi = part_phi[jp] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = part_eta[jp] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_same_dphi_deta[imult]->Fill(dphi, deta);

					if ( fabs(deta)>deta_cut ){
						hv22pc_mult->Fill(nmult_mid, cos(2*(part_phi[jp] - part_phi[ip])));
					}

				}//jp

				if ( vec_pt[imult].size()>0 ){
					hntrig_mixed_mid->Fill(nmult_mid);
				}

				//mixed events
				for (unsigned int jj=0; jj<vec_pt[imult].size(); jj++){
					if ( vec_pt[imult][jj]<pt_min || vec_pt[imult][jj]>pt_max ) continue;
					if ( fabs(vec_eta[imult][jj])>eta_max ) continue;

					float dphi = vec_phi[imult][jj] - part_phi[ip];
					if ( dphi<-const_pi/2 ) dphi += 2*const_pi;
					else if ( dphi>3*const_pi/2 ) dphi -= 2*const_pi;

					float deta = vec_eta[imult][jj] - part_eta[ip];
					if ( fabs(deta)>deta_max ) continue;

					h2d_mixed_dphi_deta[imult]->Fill(dphi, deta);
				}//jj
			}//ip


			//fill mixed event pool
			vec_pt[imult].clear();
			vec_phi[imult].clear();
			vec_eta[imult].clear();

			for (int ip=0; ip<npart; ip++){
				if ( part_pt[ip]<pt_min || part_pt[ip]>pt_max ) continue;
				if ( fabs(part_eta[ip])>eta_max ) continue;
				vec_pt[imult].push_back(part_pt[ip]);
				vec_phi[imult].push_back(part_phi[ip]);
				vec_eta[imult].push_back(part_eta[ip]);
			}

		}//ien


		delete infile;

	}//

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	for (int im=0; im<nmult; im++){
		h2d_same_dphi_deta[im]->Write();
		h2d_mixed_dphi_deta[im]->Write();
	}

	heta_pt_all->Write();
	heta_pt_ana->Write();
	hmult_all->Write();
	hmult_mid->Write();
	hntrig_mid->Write();
	hntrig_mixed_mid->Write();
	hv22pc_mult->Write();
	outfile->Close();

}
