// main03.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how different processes can be selected and studied.
// All input is specified in the main03.cmnd file.
// Also illustrated output to be plotted by Python/Matplotlib/pyplot.

#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

//#include "Pythia8Plugins/FastJet3.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
	pythia.readFile("mainEx04.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	//float f_scale;
	//float f_bMPI;
	//int i_nMPI;

	int i_np;
	int i_mult_svx;
	int i_mult_fvtxs1030;
	int i_mult_fvtxs1226;
	int i_mult_fvtxn1030;
	int i_mult_fvtxn1226;
	int i_mult_bbcs;
	int i_mult_bbcn;

	int i_p_id[1000];
	float f_p_pt[1000];
	float f_p_eta[1000];
	float f_p_phi[1000];

	/*
	int i_p_id[10];
	float f_p_pt[10];
	float f_p_eta[10];
	float f_p_phi[10];
	*/

	// Tree output
	auto T = new TTree("T","Pythia event");
	//T->Branch("scale",&f_scale,"scale/F");
	//T->Branch("bMPI",&f_bMPI,"bMPI/F");
	//T->Branch("nMPI",&i_nMPI,"nMPI/I");
	//T->Branch("multV0M",&i_mult_v0m,"multV0M/I");
	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");

	TH1D *hnsvx_mb = new TH1D("hnsvx_mb","",200,0,200);
	TH1D *hnfvtxs1030_mb = new TH1D("hnfvtxs1030_mb","",200,0,200);
	TH1D *hnfvtxs1226_mb = new TH1D("hnfvtxs1226_mb","",200,0,200);
	TH1D *hnfvtxn1030_mb = new TH1D("hnfvtxn1030_mb","",200,0,200);
	TH1D *hnfvtxn1226_mb = new TH1D("hnfvtxn1226_mb","",200,0,200);

	TH1D *hnjpsi = new TH1D("hnjpsi","",10,0,10);
	TH2D *heta_pt = new TH2D("heta_pt","",100,-5,5,100,0,10);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

		//f_scale = pythia.event.scale();
		//f_bMPI = pythia.info.bMPI();
		//i_nMPI = pythia.info.nMPI();

		//multiplicity calculation
		/*
		i_mult_v0m = 0;

		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) ) continue;

			float tmp_eta = event[i].eta();

			if ( (tmp_eta>2.8 && tmp_eta<5.1) || (tmp_eta>-3.7 && tmp_eta<-1.7) ){
				i_mult_v0m++;
			}

		}

		hmult_v0m->Fill(i_mult_v0m);
		*/

		//BBC trig.
		//multiplicity distribution
		i_mult_svx = 0;
		i_mult_fvtxs1030 = i_mult_fvtxn1030 = 0;
		i_mult_fvtxs1226 = i_mult_fvtxn1226 = 0;
		i_mult_bbcs = i_mult_bbcn = 0;
		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal() && event[i].isCharged()) ) continue;

			float tmp_eta = event[i].eta();

			if ( fabs(tmp_eta)<1.0 ){
				i_mult_svx++;
			}else if ( tmp_eta>1.0 && tmp_eta<3.0 ){
				i_mult_fvtxn1030++;
				if ( tmp_eta>1.2 && tmp_eta<2.6 ){
					i_mult_fvtxn1226++;
				}
			}else if ( tmp_eta>-3.0 && tmp_eta<-1.0 ){
				i_mult_fvtxs1030++;
				if ( tmp_eta>-2.6 && tmp_eta<-1.2 ){
					i_mult_fvtxs1226++;
				}
			}else if ( tmp_eta>3.0 && tmp_eta<3.9 ){
				i_mult_bbcn++;
			}else if ( tmp_eta>-3.9 && tmp_eta<-3.0 ){
				i_mult_bbcs++;
			}
		}//

		if ( i_mult_bbcs==0 || i_mult_bbcn==0 ) continue;

		hnsvx_mb->Fill(i_mult_svx);
		hnfvtxs1030_mb->Fill(i_mult_fvtxs1030);
		hnfvtxs1226_mb->Fill(i_mult_fvtxs1226);
		hnfvtxn1030_mb->Fill(i_mult_fvtxn1030);
		hnfvtxn1226_mb->Fill(i_mult_fvtxn1226);

		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal() && event[i].isCharged()) ) continue;

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();

			if ( fabs(tmp_eta)>5.0 ) continue;

			heta_pt->Fill(tmp_eta, tmp_pt);
		}


		i_np = 0;
		set<int> set_id;
		set<int> set_pid;
		set<int>::iterator iter;

		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal() && event[i].isCharged()) ) continue;

			int tmp_id = event[i].id();
			//if ( abs(tmp_id)!=13 ) continue;
			if ( !(abs(tmp_id)==13 || abs(tmp_id)==11) ) continue;

			int tmp_pid = event[event[i].mother1()].id();
			//if ( tmp_pid!=443 ) continue;
			if ( !(tmp_pid==443 || tmp_pid==100443) ) continue;

			set_id.insert(i);
			set_pid.insert(event[i].mother1());
		}

		hnjpsi->Fill(set_pid.size());
		if ( set_pid.size()!=1 ) continue;

		/*
		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) ) continue;

			int tmp_id = event[i].id();

			if ( abs(tmp_id)==211 || abs(tmp_id)==321 || abs(tmp_id)==2212 ){
				heta->Fill(event[i].eta());
			}

			if ( !(abs(tmp_id)==11 || abs(tmp_id)==13) ) continue;

			int tmp_pid = event[event[i].mother1()].id();

			//if ( tmp_pid==553 || tmp_pid==100553 || tmp_pid==200553 ){
			if ( tmp_pid==443 || tmp_pid==100443 ){
				set_pid.insert(event[i].mother1());
			}
		}
		*/

		/*
		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) ) continue;

			int tmp_id = event[i].id();

			if ( !(abs(tmp_id)==5132 || abs(tmp_id)==5232) ) continue;

			i_p_id[i_np] = tmp_id;
			f_p_pt[i_np] = event[i].pT();
			f_p_eta[i_np] = event[i].eta();
			f_p_phi[i_np] = event[i].phi();

			i_np++;
		}
		*/
		//if ( i_np<1 ) continue;

		for (iter=set_pid.begin(); iter!=set_pid.end(); iter++){
			int tmp_index = *iter;

			int tmp_id = event[tmp_index].id();
			float tmp_pt = event[tmp_index].pT();
			float tmp_eta = event[tmp_index].eta();
			float tmp_phi = event[tmp_index].phi();
			//bool tmp_charge = event[tmp_index].isCharged();

			i_p_id[i_np] = tmp_id;
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;

			i_np++;
		}

		for (iter=set_id.begin(); iter!=set_id.end(); iter++){
			int tmp_index = *iter;

			int tmp_id = event[tmp_index].id();
			float tmp_pt = event[tmp_index].pT();
			float tmp_eta = event[tmp_index].eta();
			float tmp_phi = event[tmp_index].phi();
			//bool tmp_charge = event[tmp_index].isCharged();

			i_p_id[i_np] = tmp_id;
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;

			i_np++;
		}

		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal() && event[i].isCharged()) ) continue;
			//if ( !(event[i].isFinal()) ) continue;

			int tmp_id = event[i].id();
			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();
			//bool tmp_charge = event[i].isCharged();

			if ( fabs(tmp_eta)>5.0 ) continue;

			i_p_id[i_np] = tmp_id;
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;

			i_np++;

			if ( i_np>=1000 ) break;
		}

		T->Fill();

		for (int ii=0; ii<1000; ii++){
			i_p_id[ii] = 0;
			f_p_pt[ii] = f_p_eta[ii] = f_p_phi[ii] = -999;
		}

  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");

	hnsvx_mb->Write();
	hnfvtxs1030_mb->Write();
	hnfvtxs1226_mb->Write();
	hnfvtxn1030_mb->Write();
	hnfvtxn1226_mb->Write();

	hnjpsi->Write();
	heta_pt->Write();
	T->Write();
	outfile->Close();

  // Done.
  return 0;
}
