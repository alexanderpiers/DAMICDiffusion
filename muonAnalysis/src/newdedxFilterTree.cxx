#include "newdedxFilterTree.h"

using namespace std;

#ifndef analysis_H
#define analysis_H

double getArrayMax(double* arr, int n){
	// Gets and returns the maximum value from the array
	double max = 0;

	for(int i=0; i<n; i++){
		if(arr[i] > max) max=arr[i];
	}
	
	return max;
}


double getArrayMin(double* arr, int n){
	// Gets and returns the maximum value from the array
	double min;

	for(int i=0; i<n; i++){
		if(i == 0) min = arr[i];
		if(arr[i] < min) min = arr[i];
	}
	
	return min;
}

#endif

// Create a second tree to track index and dedx values because for some reason 
// thingss are not working as I think they should
void createDedxTree(){

	TFile *muonClusterFile = new TFile("~/test_muon_no_delta_0999ccf.root", "UPDATE");	
	TTree *muonClusterTree = (TTree*)muonClusterFile->Get("clusters_tree");

	// Define input variables
	double minEnergy = 500, minccf = 0.995, maxdEdx = 30;
	int minTrackLength = 225;
	const char *outfile = "~/muon_no_delta_0999ccf_dedx.root";
	cout << "Creating new root file: " << outfile << endl;
	TFile *goodMuonFile = new TFile(outfile, "RECREATE");
	cout << "Creating new tree..." << endl;
	cout << "Parameters. Minimum cluster energy: " << minEnergy << " Minimum curve correlation factor: " << minccf << " Minimum track length: " << minTrackLength << " Max dEdx: " << maxdEdx << endl;

	TTree *goodMuonTree = new TTree("dedx", "dedx");

	// Create and set variables of the new tree
	double eFluctuationMax;
	int idx;
	goodMuonTree->Branch("energyFluctuationMaximum", &eFluctuationMax);
	goodMuonTree->Branch("idx", &idx);



	for(int i=0; i<muonClusterTree->GetEntries(); i++){
		TH1D * h1 = new TH1D();	
		h1 = dedxFluctuation(muonClusterTree, i);
		eFluctuationMax = h1->GetMaximum();
		idx = i;
		cout << eFluctuationMax << endl;
		goodMuonTree->Fill();
		delete h1;

	}

	goodMuonTree->Write();
	delete goodMuonFile;


	return;
}

// Uses the tree of muon tracks
void createGoodMuonTree(){

	double minccf = 0.995;
	double minTrackLength = 175;
	double minEnergy = 500;
	double maxdedx = 12;

	// Read in the required trees
	TFile *muonClusterFile = new TFile("~/test_muon_no_delta_0999ccf.root");
	TTree *muonClusterTree = (TTree*)muonClusterFile->Get("clusters_tree");
	TFile *dedxFile = new TFile("~/muon_no_delta_0999ccf_dedx.root");
	TTree *dedxTree = (TTree*)dedxFile->Get("dedx");

	// Clone old tree
	TFile *newMuonClusterFile = new TFile("~/goodMuons12keV_0999ccf.root", "RECREATE");
	TTree *newMuonClusterTree = muonClusterTree->CloneTree(0);

	// Set tree variables
	TParameter<double> *ccf;
	TParameter<double> *chargetot;
	TParameter<double> *trackLength;
	double efMax;
	muonClusterTree->SetBranchAddress("curve_correlation_factor", &ccf);
	muonClusterTree->SetBranchAddress("charge_total", &chargetot);
	muonClusterTree->SetBranchAddress("track_length", &trackLength);
	dedxTree->SetBranchAddress("energyFluctuationMaximum", &efMax);

	int nentries = muonClusterTree->GetEntries();

	for(int i=0; i<nentries; i++)
	{
		muonClusterTree->GetEntry(i);
		dedxTree->GetEntry(i);
		// cout << trackLength->GetVal() << endl;

		if(trackLength->GetVal() > minTrackLength && efMax < maxdedx){
			cout << efMax << endl;
			newMuonClusterTree->Fill();
		}


	}

	// Write the data to file
	newMuonClusterTree->Write();
	cout << "Total Good Muon entries: " << newMuonClusterTree->GetEntries() << endl;
	// delete 


	return;
}

void plotAllMuonTracks(const char* infile, int i){

	// TApplication myapp("myapp", 0, 0);
	// TCanvas *c = new TCanvas("c", "c", 800, 600);
	TFile *f = new TFile(infile);
	TTree *t = (TTree*)f->Get("clusters_tree");

	// Set address of parameters we need
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


	for(int i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		xcopy = *x;
		ycopy = *y;
		qcopy = *q;
		

		xArr = xcopy.GetArray(); yArr = ycopy.GetArray();	
		n = x->GetSize();

		// Get the limits of the track
		double xmin = getArrayMin(xArr, n);
		double xmax = getArrayMax(xArr, n);
		double ymin = getArrayMin(yArr, n);
		double ymax = getArrayMax(yArr, n);
		// Create the 2D histogram
		TH2D* h2 = new TH2D("tr", "track", xmax-xmin, xmin, xmax, ymax-ymin, ymin, ymax);

		// Fill the histogram
		h2->FillN(n, xArr, yArr, qcopy.GetArray());
		h2->Draw("colz");

		cout << "Enter for next histogram: " << ccf->GetVal();
		cin.get();
		delete h2;
	}
	return;
}