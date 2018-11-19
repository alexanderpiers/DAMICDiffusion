#include "analysis.h"

using namespace std;


// Plots the 2D histogam with 
TH2D* histEnergyvDistance(TTree *tree, double zmin, double zmax, bool resolveDeltaRay){
	// Make definitions
	TH2D *h2 = new TH2D();
	TArrayD *x, *y, *q;
	double *xx, *yy, *qq;
	double conversionFactor = 10300/6.4;
	int n, zcount;
	TF1 *tf;
	position ip;
	double dedx = 0;
	double vx, vy, sx, sy, inter, slope, projection, xmin, xmax, length, dedxMax=1., projMax=1., projArrayMax;
	// Set the appropriate branch address
	tree->SetBranchAddress("pixel_x", &x);
	tree->SetBranchAddress("pixel_y", &y);
	tree->SetBranchAddress("pixel_val", &q);
	int nTracks = tree->GetEntries();

	// Iterate over all track entries
	for(int i=0; i<nTracks; i++)
	{
		tree->GetEntry(i);
		xx = x->GetArray(); yy = y->GetArray(); qq = q->GetArray();
		n = x->GetSize();


		// Create array to store values in acceptable z range
		double *projArray = new double[10000];
		double *qArray = new double[10000];
		zcount = 0;	
		dedx = 0;

		// Conver pixel value to energy
		convertVal2Energy(qq, n, conversionFactor);

		getDistanceFromTrack(xx, yy, qq, n, zmin, zmax, resolveDeltaRay, projArray, qArray, dedx, zcount);
		
		// Set axis of the histogram
		projArrayMax = getArrayMax(projArray, zcount);
		if(projArrayMax > projMax){
			projMax = int(projArrayMax);
			h2->SetBins(2*projMax, 0, projMax, dedxMax, 0, dedxMax);
		}

		if(dedx > dedxMax){
			dedxMax = int(dedx);
			h2->SetBins(2*projMax, 0, projMax, dedxMax, 0, dedxMax);
		}
		
		// Create an array to hold the dedx values for the histogram
		double * dedxArray = new double[zcount];
		std::fill_n(dedxArray, zcount, dedx);
		h2->FillN(zcount, projArray, dedxArray, qArray);

		// Delete unecessary arrays
		delete projArray; delete qArray; delete dedxArray;

	}

	// Set histogram properties
	h2->SetStats(false);
	char titleString[100];
	sprintf(titleString, "Energy and spread hist for depth range %.1f to %0.1f um", zmin,  zmax);
	h2->SetTitle(titleString);
	h2->GetYaxis()->SetTitle("#frac{#DeltaE}{#Deltas}  (keV/Pixel)");
	h2->GetXaxis()->SetTitle("Distance from Track (Pixels)");
	h2->Draw("colz");
	return h2;


}

TH1D* histDistance(TTree *tree, double zmin, double zmax, bool energyFilt, double emin, double emax, bool deltaRayRejection, bool draw){

	// Define parameters
	TH1D *h1 = new TH1D("t", "Spread of Muon Tracks", 25, 0, 5);
	TStyle *gStyle = new TStyle();
	TArrayD *x = new TArrayD(); TArrayD *y = new TArrayD(); TArrayD *q = new TArrayD();
	double *xx, *yy, *qq;
	double conversionFactor = 10300/6.4;
	int n, zcount;
	TF1 *tf;
	position ip;
	bool deltaray;
	double vx, vy, sx, sy, inter, slope, projection, xmin, xmax, length, dedx;
	double projArray[100000];
	double qEnergyArray[100000];
	// Set the appropriate branch address
	tree->SetBranchAddress("pixel_x", &x);
	tree->SetBranchAddress("pixel_y", &y);
	tree->SetBranchAddress("pixel_val", &q);
	int nTracks = tree->GetEntries();
	// Iterate over all track entries
	for(int i=0; i<nTracks; i++){
	
		tree->GetEntry(i);
		xx = x->GetArray(); yy = y->GetArray(); qq = q->GetArray();
		n = x->GetSize();

		// Conver q to energy via conversion factor
		for(int k=0; k<n; k++){
			qq[k] /= conversionFactor;
		}	
			
		zcount = 0;
		dedx = 0;
		
		// Finds the distance of each point in the track. Returns projArry and qEnergyArray with the correct depth filter		
		getDistanceFromTrack(xx, yy, qq, n, zmin, zmax, deltaRayRejection, projArray, qEnergyArray, dedx, zcount);
		// Fill the histogram
		if(energyFilt){
			if(dedx > emin && dedx < emax) h1->FillN(zcount, projArray, qEnergyArray);
		} else{
			h1->FillN(zcount, projArray, qEnergyArray);
		}

	}

	// Create the histogram title
	char titlestring[250];
	if(energyFilt){
		sprintf(titlestring, "Diffusion of charge between z=%.1f and z=%.1f um and e=%.1f and e=%.1f keV", zmin, zmax, emin, emax);
	}else{
		sprintf(titlestring, "Diffusion of charge between z=%.1f and z=%.1f", zmin, zmax);
	}

	// Set histogram parameters
	h1->SetTitle(titlestring);
	h1->GetXaxis()->SetTitle("Distance (pixels)");
	h1->GetYaxis()->SetTitle("Amount of charge (keV)");
	h1->SetLineWidth(3);
	// Gaussian fit over the low diffusion regime
	TF1* f = new TF1("f1", "gaus", 0, 4);
	f->SetParameter(1,0);
	h1->Fit("f1", "QR0", 0, 1.5);
	gStyle->SetOptFit();
	if(draw){
		TCanvas* c = new TCanvas("c", "charge_diffusion", 800, 600);
		c->SetLogy();
		h1->Draw("HIST");
		h1->GetFunction("f1")->Draw("same");
	}
	delete x; delete y; delete q;
	return h1;
}

void convertVal2Energy(double *q, int n,  double conversionFactor){

	for(int i=0; i<n; i++){
		q[i] /= conversionFactor;
	}
	return;
}


// Plot the diffusion sigma as a function of depth
// Inputs:
// TTree *tree - tree containing muon tracks
// double deltaZ - zmax-zmin of each track, so the extent of the z direction of each point
// double emin - minimum dE/dx
// double emax - maximum dE/dx
TGraph* sigmaVDepth(TTree *tree, double deltaZ, double zstart,  double emin, double emax){
	
	// Define parameters
	double z[500];
	double sigma[500];
	int tGraphCnt = 0;
	int zmin, zmax;
	int nSlices = int((CCDWidth-zstart)/deltaZ);

	// Iterate over all slices, getting the diffusion spread
	for(int i=0; i<nSlices; i++){
		zmin = i*deltaZ + zstart; zmax = (i+1)*deltaZ + zstart;
		cout << zmin << "-" << zmax << endl;
		TH1D* h = new TH1D();
		h = histDistance(tree, zmin, zmax, true, emin, emax, true, true);
		z[i] = (zmin + zmax)/2;
		cout << zmin << "-" << zmax << ":" << endl;
		cout << h->GetFunction("f1")->GetParameter(2) << endl;
		sigma[i] =  h->GetFunction("f1")->GetParameter(2);  
		delete h;	
	}
	TGraph *tg = new TGraph(nSlices, z, sigma);
	return tg;

}

TGraph* plotSigmaVDepth(){
	int n = 20;
	double z[] = {12.5, 37.5, 62.5, 87.5, 112.5, 137.5, 162.5, 187.5, 212.5, 237.5, 262.5, 287.5, 312.5, 337.5, 362.5, 387.5, 412.5, 437.5, 462.5, 487.5};
	double sig45[] = {0.377, 0.43223, 0.53102, 0.4719, 0.4251, 0.5233, 0.5016, 0.5030, 0.5453, 0.5595, 0.5304, 0.6107, 0.6571, 0.7569, 0.6596, 0.6359, 0.6345, 0.6230, 0.7075, 0.8592};
	//double sig23[] = {};
	for(int i=0; i<n; i++){
		sig45[i] *= 15;
	}

	TGraph *tg = new TGraph(n, z, sig45);
	tg->SetMarkerSize(1);
	tg->SetMarkerStyle(5);
	tg->SetTitle("Diffusion #\sigma_{x} vs. Depth of Interaction");
	tg->GetXaxis()->SetTitle("Depth (#\mum)");
	tg->GetYaxis()->SetTitle("#\sigma_{x} (#\mum)");
	//tg->GetYaxis()->SetRangeUser(0, 1);
	tg->Draw("AP");

	// Fit expected function to data
	double q = 1.6e-19; // C
	double T = 130.; // K 
	double kb = 1.38e-11; // um^2 kg / K s^2 
	double epsSi = 11.68 * 8.85e-18; // C / Vum 
	double Vapp = 128.;
	double rho = 3.9e-1*q; // um^-3 
	double conversion =1.e-6;
	
	// TF1 *f = new TF1("sig_analytic", "TMath::Sqrt( [0]/[5]*TMath::Log([1]-x/([2]/[5] + [3])))*[4]", 1, CCDWidth);
	TF1 *f = new TF1("sig_analytic", "TMath::Sqrt( [0]*TMath::Log([2]-x*[1]))", 0, CCDWidth);	
	
	f->SetParameter(0, -2*kb*T*epsSi/(rho*q));
	f->SetParameter(1, 1/(epsSi*Vapp/(rho*CCDWidth)+0.5*CCDWidth));
	f->FixParameter(2, 1);
//	f->FixParameter(3, conversion);
	// Set Parameters
	//f->SetParameter(0, -2*kb*T*epsSi/q);
	//f->FixParameter(1, 1);
	//f->FixParameter(2, epsSi*Vapp/CCDWidth);
	//f->FixParameter(3, 0.5*CCDWidth);
	//f->FixParameter(4, conversion);
	//f->SetParameter(5, rho);
//	f->Draw("A");
	//return f;	
	tg->Fit(f, "R");
	cout << "Rho: " << f->GetParameter(5) << endl;
	cout << f->Eval(100) << endl;
	cout << f->Eval(200) << endl;

	return tg;
}
void saveAllHist(TTree *tree){
	int n = tree->GetEntries();
	TArrayD *x, *y, *q;
	double *xx, *yy, *qq;
	tree->SetBranchAddress("pixel_x", &x);
	tree->SetBranchAddress("pixel_y", &y);
	tree->SetBranchAddress("pixel_val", &q);
	int m = 0;
	for(int i=0; i<n; i++){
		// Get the current entry from the tree
		tree->GetEntry(i);
		xx = x->GetArray(); yy = y->GetArray(); m = x->GetSize();

		TH2D *h2 = new TH2D();	
		h2 = plot2DTrackDepth(x, y, false);
		TF1 * fit = new TF1();
		fit = fitMuonLine(xx, yy, m);
		// Create canvase and image
		TCanvas *c = new TCanvas;
		h2->Draw("colz");
		fit->Draw("same");
		TImage *img = TImage::Create();
		img->FromPad(c);

		// Create the filename
		char *titlestr = new char[100];
		sprintf(titlestr, "/usr/lusers/apiers/images/depth_track%i.png", i);

		// Write to file
		img->WriteImage(titlestr);
		
		// Cleanup pointers
		delete h2;
		delete fit;
		delete c;
		delete titlestr;
		delete img;
	
	}
	

}

// Saves a 1D histogram to filename
void savehist(TH1D *h, const char *filename, char *histname){
	
	TFile *f = new TFile(filename, "UPDATE");
	h->Write(histname);
	delete f;

	return;


}

void savehist(TH2D *h, const char *filename, char *histname){

	TFile *f = new TFile(filename, "UPDATE");
	h->Write(histname);
	delete f;

	return;
}

void savegraph(TGraph *graph, const char *filename, char* objname){

	TFile *f = new TFile(filename, "UPDATE");
	graph->Write(objname);
	delete f;

	return;
}