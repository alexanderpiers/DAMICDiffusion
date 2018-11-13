#include "muonFilter.c"
#include <cmath>
#include <TImage>
#include <algorithm>


// Plots the 2D histogam with 
TH2D* histEnergyvDistance(TTree *tree, double zmin=400, double zmax=425, bool resolveDeltaRay=true){
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

TH1D* histDistance(TTree *tree, double zmin=400, double zmax=450, bool energyFilt=false, double emin=2., double emax=4., bool deltaRayRejection=true){

	// Define parameters
	TCanvas* c = new TCanvas("c", "charge_diffusion", 800, 600);
	c->SetLogy();
	TH1D *h1 = new TH1D("t", "Spread of Muon Tracks", 25, 0, 5);
	TArrayD *x, *y, *q;
	double *xx, *yy, *qq;
	double conversionFactor = 10300/6.4;
	int n, zcount;
	TF1 *tf;
	position ip;
	double dedx = 0;
	bool deltaray;
	double vx, vy, sx, sy, inter, slope, projection, xmin, xmax, length, dedx;
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
		double *projArray = new double[100000];
		double *qEnergyArray = new double[100000];
		
		// Finds the distance of each point in the track. Returns projArry and qEnergyArray with the correct depth filter		
		getDistanceFromTrack(xx, yy, qq, n, zmin, zmax, deltaRayRejection, projArray, qEnergyArray, dedx, zcount);

		// Fill the histogram
		if(energyFilt){
			if(dedx > emin && dedx < emax) h1->FillN(zcount, projArray, qEnergyArray);
		} else{
			h1->FillN(zcount, projArray, qEnergyArray);
		}

		// Delete arrays
		delete projArray; delete qEnergyArray;


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
	h1->Fit("f1", "QR", 0, 1.5);
	gStyle->SetOptFit();
	h1->Draw("HIST");
	h1->GetFunction("f1")->Draw("same");
	return h1;
}

// Helper function to do the computation calculating the distance from the best fit line. Boolean option to exclude delta rays with threshold cutoff of distance
// Inputs:
// double *x - array of x positions
// double *y - array of y positions
// double *q - array of the value of the CCD
// int n - number of entries in the x and y arrays
// double zmin - minimum z threshold
// double zmax - maximum z threshold
// bool resolveDeltaRay - if true, reject events where delta rays are detected
// 
// Outputs:
// double *proj - array of distance from muon track line
// double *q - array of the amount of charge in each 

void  getDistanceFromTrack(double *x, double *y, double *q, int n, double zmin, double zmax,  bool resolveDeltaRay=true, double *proj, double *qEnergy, double &dedx, int &zcount){
	
	bool deltaray = false;
	double sx, sy, projection;
	// Get information about the line
	TF1* tf = fitMuonLine(x, y, n);
	double slope = -1/tf->GetParameter(1);

	double *z = getZ(x, y, n);
	double vx = 1/sqrt(1 + pow(slope,2));
	double vy = slope/sqrt(1 + pow(slope,2));
	position ip = getInitialPosition(x, y, n);
	
	// Calculate the length of the line segment
	double xmin = getArrayMin(x, n); double xmax = getArrayMax(x, n);
	double length = sqrt(pow(xmax-xmin,2) + pow(tf->Eval(xmax)-tf->Eval(xmin),2))*(zmax-zmin)/CCDWidth;

	// Iterate over all x,y pairs in the

	for(int j=0; j<n; j++){
		sx = x[j] - ip.x + 0.5;
		sy = y[j] - ip.y + 0.5;
		projection = vx*sx + vy*sy;
		
		// if in the depth range, add to histogram
		if(z[j] > zmin && z[j] < zmax && abs(projection) < 4){
			proj[zcount] = abs(projection);
			qEnergy[zcount] = q[j];
			dedx += q[j];
			zcount++;
		}else if(z[j] > zmin && z[j] < zmax && abs(projection) >= 4){
			deltaray = true;
			if(!resolveDeltaRay){
				proj[zcount] = abs(projection);
				qEnergy[zcount] = q[j];
				dedx += q[j];
				zcount++;
			}
		}

	}

	// Find the energy per  unit length
	dedx /= length;

	// If delta ray detected set the zcount to 0. Then when we fill the histogram
	if(deltaray && resolveDeltaRay) zcount=0;

	return;

}

void convertVal2Energy(double *q, int n,  double conversionFactor=10300/6.4){

	for(int i=0; i<n; i++){
		q[i] /= conversionFactor;
	}
	return;
}

void saveAllHist(TTree *tree){
	int n = tree->GetEntries();
	TArrayD *x, *y, *z;
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
