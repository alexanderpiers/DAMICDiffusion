#include "muonFilter.c"
#include <cmath>
#include <TImage>
#include <algorithm>

// Plots the 2D histogam with 
TH2D* histEnergyvDistance(TTree *tree, double zmin=400, double zmax=425){
	// Make definitions
	TH2D *h2 = new TH2D();
	TArrayD *x, *y, *q, *z;
	double conversionFactor = 10300/6.4;
	int n, zcount;
	TF1 *tf;
	position ip;
	double dedx = 0;
	double vx, vy, sx, sy, inter, slope, projection, xmin, xmax, length, dedxMax=1., projMax=1.;
	// Set the appropriate branch address
	tree->SetBranchAddress("pixel_x", &x);
	tree->SetBranchAddress("pixel_y", &y);
	tree->SetBranchAddress("pixel_val", &q);
	int nTracks = tree->GetEntries();
	// Iterate over all track entries
	for(int i=0; i<nTracks; i++)
	{
		tree->GetEntry(i);
		tf = fitMuonLine(x, y);
		inter = tf->GetParameter(0);
		slope = -1/tf->GetParameter(1);
		z = getZ(x, y);
		vx = 1/sqrt(1 + pow(slope,2));
		vy = slope/sqrt(1 + pow(slope,2));
		ip = getInitialPosition(x, y);
		// Calculate the length of the line segment
		xmin = getArrayMin(x); xmax = getArrayMax(x);
		length = sqrt(pow(xmax-xmin,2) + pow(tf->Eval(xmax)-tf->Eval(xmin),2))*(zmax-zmin)/CCDWidth;
		n = x->GetSize();

		// Iterate over all pixel values to get the total energy of the slice
		dedx = 0;
	//	for(int j=0; j<n; j++){
	//		dedx += (*q)[j];
	//	}
	//	dedx /= (conversionFactor*length);

		// Create array to store values in acceptable z range
		double *projArray = new double[10000];
		double *qArray = new double[10000];
		zcount = 0;

		for(int j=0; j<n; j++){
			sx = (*x)[j] - ip.x + 0.5;
			sy = (*y)[j] - ip.y + 0.5;
			projection = vx*sx + vy*sy;
		//	projection = abs(inter + slope * ((*x)[j] + 0.5) - (*y)[j] - 0.5)/sqrt(1 + pow(slope,2));		
			// if in the depth range, add to histogram
			if((*z)[j] > zmin && (*z)[j] < zmax){
				// cout << "dist away: " << projection<< endl;
				// Set hist limits appropriately
				if(abs(projection) > projMax){
				projMax = int(abs(projection));
				h2->SetBins(2*projMax, 0, projMax, dedxMax, 0, dedxMax);
				}
				dedx += (*q)[j]/conversionFactor;
				projArray[zcount] = abs(projection);
				qArray[zcount] = (*q)[j]/conversionFactor;
 				zcount++;
				//h2->Fill(abs(projection), dedx, (*q)[j]/conversionFactor);
			}
		}
		// Create dedx array and then fill all values in hist
		double * dedxArray = new double[zcount];
		dedx /= length;
		if(dedx > dedxMax){
			dedxMax = int(dedx);
			cout << dedxMax << endl;
			h2->SetBins(2*projMax, 0, projMax, dedxMax, 0, dedxMax);
		}
		std::fill_n(dedxArray, zcount, dedx);
		h2->FillN(zcount, projArray, dedxArray, qArray);

		// Delete unecessary arrays
		delete projArray; delete qArray; delete dedxArray;

	}
	h2->SetStats(false);
	char titleString[100];
	sprintf(titleString, "Energy and spread hist for depth range %.1f to %0.1f um", zmin,  zmax);
	h2->SetTitle(titleString);
	h2->Draw("colz");

	return h2;


}

TH1D* histDistance(TTree *tree, double zmin=400, double zmax=450, bool energyFilt=false, double emin=2., double emax=4.){

	// Define parameters
	TCanvas* c = new TCanvas("c", "charge_diffusion", 800, 600);
	c->SetLogy();
	TH1D *h1 = new TH1D("t", "Spread of Muon Tracks", 25, 0, 5);
	TArrayD *x, *y, *q, *z;
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
		tf = fitMuonLine(x, y);
		slope = -1/tf->GetParameter(1);
		z = getZ(x, y);
		vx = 1/sqrt(1 + pow(slope,2));
		vy = slope/sqrt(1 + pow(slope,2));
		ip = getInitialPosition(x, y);
		
		// Calculate the length of the line segment
		xmin = getArrayMin(x); xmax = getArrayMax(x);
		length = sqrt(pow(xmax-xmin,2) + pow(tf->Eval(xmax)-tf->Eval(xmin),2))*(zmax-zmin)/CCDWidth;
		n = x->GetSize();
		
		deltaray = false;
		zcount = 0;
		dedx = 0;
		double *zAcceptableArray = new double[100000];
		double *qAcceptableArray = new double[100000];
		for(int j=0; j<n; j++){
			sx = (*x)[j] - ip.x + 0.5;
			sy = (*y)[j] - ip.y + 0.5;
			projection = vx*sx + vy*sy;
		
			// if in the depth range, add to histogram
			if((*z)[j] > zmin && (*z)[j] < zmax && abs(projection) < 4){
				zAcceptableArray[zcount] = abs(projection);
				qAcceptableArray[zcount] = (*q)[j]/conversionFactor;
				dedx += (*q)[j]/conversionFactor;
				zcount++;
			}else if((*z)[j] > zmin && (*z)[j] < zmax && abs(projection) >= 4){
				deltaray = true;
			}

		}
		dedx /= length;
		// Fill the histogram if no delta rays detected
		if(!deltaray){
			if(energyFilt){
				if(dedx > emin && dedx < emax) h1->FillN(zcount, zAcceptableArray, qAcceptableArray);
			}else {
				h1->FillN(zcount, zAcceptableArray, qAcceptableArray);
			}	
		}


		// Delete arrays
		delete zAcceptableArray; delete qAcceptableArray;


	}
	char titlestring[250];
	if(energyFilt){
		sprintf(titlestring, "Diffusion of charge between z=%.1f and z=%.1f um and e=%.1f and e=%.1f keV", zmin, zmax, emin, emax);
	}else{
		sprintf(titlestring, "Diffusion of charge between z=%.1f and z=%.1f", zmin, zmax);
	}
	h1->SetTitle(titlestring);
	h1->SetXTitle("Distance (pixels)");
	h1->SetYTitle("Amount of charge (keV)");
	h1->SetLineWidth(3);
	//h1->SetStats(false);
	TF1* f = new TF1("f1", "gaus", 0, 4);
	f->SetParameter(1,0);
	h1->Fit("f1", "QR", 0, 1.5);
	gStyle->SetOptFit();
	h1->Draw("HIST");
	h1->GetFunction("f1")->Draw("same");
	return h1;
}

void saveAllHist(TTree *tree){
	int n = tree->GetEntries();
	TArrayD *x, *y, *z;
	tree->SetBranchAddress("pixel_x", &x);
	tree->SetBranchAddress("pixel_y", &y);
	tree->SetBranchAddress("pixel_val", &q);
	 
	for(int i=0; i<n; i++){
		// Get the current entry from the tree
		tree->GetEntry(i);

		TH2D *h2 = new TH2D();
		TF1 *fit = new TF1();	
		h2 = plot2DTrack(x, y, q, fit, false);

		// Create canvase and image
		TCanvas *c = new TCanvas;
		h2->Draw("colz");
		fit->Draw("same");
		TImage *img = TImage::Create();
		img->FromPad(c);

		// Create the filename
		char *titlestr = new char[100];
		sprintf(titlestr, "/usr/lusers/apiers/images/track%i.png", i);

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
