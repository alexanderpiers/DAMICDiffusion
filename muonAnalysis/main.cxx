#include "muonFilter.h"
#include "analysis.h"
#include <string>

using namespace std;

void newMuonRootFile();

int main(int argc, char **argv){

	// Get tree of muon data
	TFile * f = new TFile("~/DAMICDiffusion/rootfiles/dedxfilter_8kev.root");
	TTree * t = (TTree*)f->Get("dedxFilterTree");
	// TFile * f = new TFile("~/large_muon_set.root");
	// TTree * t = (TTree*)f->Get("clusters_tree");
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
	resampleMuonTrack(t, "~/DAMICDiffusion/rootfiles/testresample.root", 5, true);
	// Plot low distribution of charge for filtered tree. 
	/*int estart=3, eend=7;
  	int zstart=0, zincrement=25, zend=500;
	const char *outfile = "~/DAMICDiffusion/rootfiles/charge_dist_dedxfilt8kev_newEFilter.root";
	for(int i=estart; i<eend; i++){
		cout << i << endl;
		for(int j=zstart; j<zend; j+=zincrement){
			TH1D *h1 = new TH1D();
			h1 = histDistance(t, j, j+zincrement,true, true, i, i+1, true, false);
			char *histname = new char[100];
			sprintf(histname, "depth_%i_%ium_%i_%ikev", j, (j+zincrement), i, (i+1));
			cout << histname << endl;
			savehist(h1, outfile, histname);
			delete histname; delete h1;
		}
	}*/
	
	// Filter out higher energy segments of track and save to file
	// const char* outfilename = "~/DAMICDiffusion/rootfiles/dedxfilter_6kev.root";
	// double energyThreshold = 6;
	// dedxFilterTree(t, outfilename, energyThreshold);

	// Create set of muon tracks from all images
	// newMuonRootFile();
	
	cout << "Analysis successfully run." << endl;

	return 0;
}

void newMuonRootFile(){

	TChain* c = new TChain("clusters_tree");
	c->Add("/gscratch/damic/data/uchicago/processed/D3500/SbBe_2015-04-10_FullBe/root/Image1*.root");
	muonFilter(c, "~/large_muon_set.root", 500, 0.995);
	return;
}





