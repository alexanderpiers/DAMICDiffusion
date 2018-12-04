#include "muonFilter.h"
#include "analysis.h"
#include <string>

using namespace std;

int main(int argc, char **argv){

	// Get tree of muon data
	TFile* f = new TFile("~/large_muon_set.root");
	TTree* t = (TTree*)f->Get("clusters_tree");
	cout << "Number of tree entries: " << t->GetEntries() <<  endl;
/*	
	// Create TApplication object and run analysis
	const char *outfile = "~/DAMICDiffusion/rootfile/sigma_v_depth_10um.root";
	//const char *outfile = "~/DAMICDiffusion/rootfile/e_flucuations100.root";
	// Iterate over a range of energies and plot
	int estart=3; int eend=8; double dz=10;
	TGraph *tg;
	for(int e=estart; e<=eend; e++){
		tg = sigmaVDepth(t, dz, 0, 450, e, e+1);
		// Create the name for the graph
		char *histname = new char[200];
		sprintf(histname, "sd%i_%ikev", e, e+1);
		savegraph(tg, outfile, histname);
		delete histname;
	}
	
	for(int i=0; i<100; i++){
		
		TH1D *h1 = new TH1D();
		h1 = dedxFluctuation(t, i);
		char *histname = new char[200];
		sprintf(histname, "event_%i", i);
		savehist(h1, outfile, histname);
		delete histname;
		delete h1;
	}*/

	// Plot low depth distribution of charge. 
	/*int estart=2, eend=9;
	const char *outfile = "~/DAMICDiffusion/rootfile/low_depth_distI.root";
	for(int i=estart; i<eend; i++){
		TH1D *h1 = new TH1D();
		h1 = histDistance(t, 0, 5, true, i, i+1, true, false);
		char *histname = new char[100];
		sprintf(histname, "depth_0_5um_%i_%ikev", i, i+1);
		savehist(h1, outfile, histname);
		delete histname; delete h1;

	}*/

	// Filter out higher energy segments of track and save to file
	const char* outfilename = "~/DAMICDiffusion/rootfile/dedxfilter_8kev.root";
	double energyThreshold = 8;
	dedxFilterTree(t, outfilename, energyThreshold);
	
	cout << "Analysis successfully run." << endl;

	return 0;
}

void newMuonRootFile(){

	TChain* c = new TChain("clusters_tree");
	c->Add("/gscratch/damic/data/uchicago/processed/D3500/SbBe_2015-04-10_FullBe/root/Image1*.root");
	muonFilter(c, "~/large_muon_set.root", 500, 0.995);
	return;
}





