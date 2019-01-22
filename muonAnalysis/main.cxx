#include "muonFilter.h"
#include "analysis.h"
#include "newdedxFilterTree.h"
#include "muonFit.h"
#include <string>

using namespace std;

void newMuonRootFile();
void saveAllMuonTracks();

int main(int argc, char **argv){

	// Get tree of muon data
	TFile * f = new TFile("~/DAMICDiffusion/rootfiles/dedxfilter_8kev_fine.root");
	TTree * t = (TTree*)f->Get("dedxFilterTreeFine");

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

	// resampleMuonTrack(t, "~/DAMICDiffusion/rootfiles/dedxfilter_6kev_fine.root", 5, true);
	// Plot low distribution of charge for filtered tree. 
	// int estart=3, eend=7;
 //  	int zstart=0, zincrement=20, zend=500;
	// const char *outfile = "~/DAMICDiffusion/rootfiles/charge_dist_dedxfilt8kev_trackThresh_fine.root";
	// for(int i=estart; i<eend; i++){
	// 	cout << i << endl;
	// 	for(int j=zstart; j<zend; j+=zincrement){
	// 		TH1D *h1 = new TH1D();
	// 		h1 = histDistance(t, j, j+zincrement,true, true, i, i+1, true, false, 1);
	// 		char *histname = new char[100];
	// 		sprintf(histname, "depth_%i_%ium_%i_%ikev", j, (j+zincrement), i, (i+1));
	// 		cout << histname << endl;
	// 		savehist(h1, outfile, histname);
	// 		delete histname; delete h1;
	// 	}
	// }

	// Plot sigma v depth for multiple energies
	// int estart=3; int eend=7;
	// const char *outfileGraphError = "~/DAMICDiffusion/rootfiles/sigma_depth_dedx8kev.root";
	// for(int i=estart; i<eend; i++){
	// 	char *graphname = new char[100];
	// 	sprintf(graphname, "sigma_depth_%i_%ikeV", i, i+1);
	// 	TGraphErrors *outgraph = sigmaVDepthFile("~/DAMICDiffusion/rootfiles/charge_dist_dedxfilt8kev_newEFilter.root", i, i+1);
	// 	savegraph(outgraph, outfileGraphError, graphname);
	// 	delete graphname; delete outgraph;
	// }
	
	// Filter out higher energy segments of track and save to file
	// const char* outfilename = "~/DAMICDiffusion/rootfiles/dedxfilter_6kev.root";
	// double energyThreshold = 6;
	// dedxFilterTree(t, outfilename, energyThreshold);

	// Create set of muon tracks from all images
	// newMuonRootFile();

	// Create Set of "good" muon clusterss
	// TFile *muonClusterFile = new TFile("~/large_muon_set.root");
	// TTree *muonClusterTree = (TTree*)muonClusterFile->Get("clusters_tree");
	// muonFilterNoDelta(muonClusterTree, "~/test_muon_no_delta.root", 500, 0.995, 225, 20);
	// newMuonRootFile();
	// createDedxTree();
	// createGoodMuonTree();
	// plotAllMuonTracks("~/goodMuons.root");
	// saveAllMuonTracks();

	saveMuonFitToFile("~/DAMICDiffusion/rootfiles/likelihoodFit/muon15kev_0999ccf_fits.root", "~/goodMuons15keV_0999ccf.root");

	
	cout << "Analysis successfully run." << endl;

	return 0;
}

void newMuonRootFile(){

	TChain* c = new TChain("clusters_tree");
	c->Add("/gscratch/damic/data/uchicago/processed/D3500/SbBe_2015-04-10_FullBe/root/Image1*.root");
	muonFilter(c, "~/test_muon_no_delta_0999ccf.root", 500, 0.999);
	// muonFilterNoDelta(c, "~/test_muon_no_delta.root", 500, 0.995, 225, 20);

	return;
}

void saveAllMuonTracks(){

	const char *outfile = "~/DAMICDiffusion/rootfiles/all_muon_tracks_20kev_0999ccf.root";
	const char *infile = "~/goodMuons20keV_0999ccf.root";

	TFile *f = new TFile(infile);
	TTree *t = (TTree*)f->Get("clusters_tree");

	TArrayD *x;
	TArrayD *y;
	TArrayD *q;
	TArrayD xcopy, ycopy, qcopy;
	double *xArr, *yArr;
	int n;
	TParameter<double> *ccf;
	t->SetBranchAddress("curve_correlation_factor", &ccf);
	t->SetBranchAddress("pixel_x", &x);
	t->SetBranchAddress("pixel_y", &y);
	t->SetBranchAddress("pixel_val", &q);

	char *histname = new char[100];

	for (int i = 0; i < t->GetEntries(); i++)
	{
		cout << i << endl;
		t->GetEntry(i);
		xcopy = *x;
		ycopy = *y;
		qcopy = *q;

		TH2D* h2 = new TH2D();
		h2 = plot2DTrack(&xcopy, &ycopy, &qcopy, false);
		h2->SetOption("colz");
		sprintf(histname, "muon%i", i);
		savehist(h2, outfile, histname);
		delete h2;
	}



	return;
}





