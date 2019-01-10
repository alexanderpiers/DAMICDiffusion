#include <string>

// Plotting function for comparing the difference between the histogram charge distribution when filtering out the dedx in two methods:
// 	--Total track dedx - so summing all the energies of a segment and dividing by the segment length
//  --Individual pixel dedx - using the pixel dedx coming from the filtering by dedx length
TCanvas * compareTrackDedx2PixelDedx(int dlow=200, int dhigh=225, int elow=5, int ehigh=6){
	
	// Read files that contain the charge distribution histograms.
	TFile *fDelta = new TFile("~/DAMICDiffusion/rootfile/charge_dist_deltarays_56kev.root");
	TFile *fNoDelta = new TFile("~/DAMICDiffusion/rootfile/charge_dist_nodeltarays_56kev.root");
	TFile *fDedx = new TFile("~/DAMICDiffusion/rootfile/charge_dist_dedx_6kev_56kev.root");
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	// Create TCanvas and add histograms to it.
	TCanvas *c = new TCanvas("c","c", 800, 600);
	c->SetLogy();
	double frange = 1.5;

	// Histogram name
	char histname[100];
	sprintf(histname, "depth_%i_%ium_%i_%ikev", dlow, dhigh, elow, ehigh);

	// Plot Histograms
	// Total Track Dedx

	return 
}