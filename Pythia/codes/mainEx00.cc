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

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("mainEx00.cfg");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

	float f_scale;
	float f_bMPI;
	int i_nMPI;
	int i_np;

	int i_p_id[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];

	const float const_pt_cut = 0.2;

	// Tree output
	auto T = new TTree("T","Pythia event");
	T->Branch("scale",&f_scale,"scale/F");
	T->Branch("bMPI",&f_bMPI,"bMPI/F");
	T->Branch("nMPI",&i_nMPI,"nMPI/I");
	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

		f_scale = pythia.event.scale();
		f_bMPI = pythia.info.bMPI();
		i_nMPI = pythia.info.nMPI();

		/*
		cout 
			<< pythia.event.scale() << " " 
			<< pythia.info.bMPI() << " " 
			<< pythia.info.nMPI() << " " 
			//<< pythia.info.nISR() << " " 
			//<< pythia.info.nFSRinProc() << " " 
			//<< pythia.info.nFSRinRes() << " " 
			<< endl;
		*/

		i_np = 0;
    for (int i = 0; i < event.size(); ++i) {

			if ( !(event[i].isFinal()) || !(event[i].isCharged()) ) continue;

			i_p_id[i_np] = event[i].id();

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();

			if ( fabs(tmp_eta)>6.0 ) continue;
			if ( fabs(tmp_eta)<1.5 && tmp_pt<const_pt_cut ) continue;

			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;

			/*
			int status = event[i].status();
			int id = event[i].id();

			double xprod = event[i].xProd();
			double yprod = event[i].yProd();
			double zprod = event[i].zProd();
			double tprod = event[i].tProd();
			double rapidity_tau = (tprod-zprod)<1e-10 ? 0 : 0.5*log((tprod+zprod)/(tprod-zprod));

			double px = event[i].px();
			double py = event[i].py();
			double pz = event[i].pz();
			double ee = event[i].e();
			double rapidity = (ee-fabs(pz))<1e-10 ? 0 : 0.5*log((ee+pz)/(ee-pz));
			*/

			i_np++;

		}

		T->Fill();

  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");
	T->Write();
	outfile->Close();

  // Done.
  return 0;
}
