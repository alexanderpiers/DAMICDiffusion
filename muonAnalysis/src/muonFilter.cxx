#include "muonFilter.h"

using namespace std;

void muonFilter(TChain* chain, char const * outfile, double minEnergy, double minccf)
{
	cout << "Creating new root file: " << outfile << endl;
	TFile* muonData = new TFile(outfile, "recreate");
	cout << "Createing new tree ... " << endl;
	
	// Need to create a copy of the chain and then convert to tree
	TChain* newchain = (TChain*)chain->CloneTree(0);
	TTree* newtree = newchain->GetTree();

	int nentries = int(chain->GetEntries());
	cout << "number of entries: " << nentries << endl;
	int ientry = 0, nbytes = 0, nb = 0;
	
	// Create variables used for the filtering. Set Branch address accordingly
	TParameter<double> * ccf;
	TParameter<double> * chargetot;
	chain->SetBranchAddress("curve_correlation_factor", &ccf);
	chain->SetBranchAddress("charge_total", &chargetot);

	// Loop over all the entries in the TChain and see which ones satisfy conditions for muon
	for(int i=1; i<nentries; i++)
	{
		ientry = chain->LoadTree(i);
		if(ientry < 0) break;
		nb = chain->GetEntry(i);
		nbytes += nb;
		
		// If we satisfy the rejection criteria for a muon, fill the tree
		if( abs(ccf->GetVal()) >  minccf && chargetot->GetVal()/10300*6.4 > minEnergy)
		{
		newtree->Fill();
		
		}
	}
	
	// Write data to file	
	newtree->Write();
	// muonData->Close();
		
	// Print how many entries in the new tree there are
	cout << "Total muon entries: " << newtree->GetEntries() << endl; 

}
TH2D* plot2DTrackDepth(TArrayD *x, TArrayD *y, bool plot){
	
	// Convert TArrayD to regular array
	double *xArr = x->GetArray(); double *yArr = y->GetArray();	
	int n = x->GetSize();
	double *z;

	// Get the limits of the track
	double xmin = getArrayMin(xArr, n);
	double xmax = getArrayMax(xArr, n);
	double ymin = getArrayMin(yArr, n);
	double ymax = getArrayMax(yArr, n);

	// Create the 2D histogram
	TH2D* h2 = new TH2D("tr", "track", ymax-ymin, xmin, xmin+ymax-ymin, ymax-ymin, ymin, ymax);
	
	// Get Depth of the track
	TF1 * tf =  fitMuonLine(xArr, yArr, n);
	getZ(xArr, yArr, z, n);
	h2->FillN(n, xArr, yArr, z);
	
	if(plot){	
		h2->Draw("colz");
		tf->Draw("same");
	}
	
	return h2;

}
// Plot 2D histogram of track
TH2D* plot2DTrack(TArrayD *x, TArrayD *y, TArrayD *q, bool plot){
	
	// Convert TArrayD to regular array
	double *xArr = x->GetArray(); double *yArr = y->GetArray();	
	int n = x->GetSize();

	// Get the limits of the track
	double xmin = getArrayMin(xArr, n);
	double xmax = getArrayMax(xArr, n);
	double ymin = getArrayMin(yArr, n);
	double ymax = getArrayMax(yArr, n);
	// Create the 2D histogram
	TH2D* h2 = new TH2D("tr", "track", xmax-xmin, xmin, xmax, ymax-ymin, ymin, ymax);

	// Fill the histogram
	h2->FillN(n, xArr, yArr, q->GetArray());

	// Do the fit
	TF1 *fit = fitMuonLine(xArr, yArr, n);
	if(plot){
		h2->Draw("colz");
		fit->Draw("same");
	}
	return h2; 

}
void  pixel2pos(double* pixel, int n){
	
	for(int i=0; i<n; i++){
		pixel[i] +=  0.5;
	}
}

void pos2pixel(double* pixel, int n){
	for(int i=0; i<n; i++){
		pixel[i] -= 0.5;
	}
}

TF1 * fitMuonLine(double *x, double *y, int n){
	
	double xmin = getArrayMin(x, n);
	double xmax = getArrayMax(x, n);
	
	// Create TGraph object
	pixel2pos(x, n); pixel2pos(y, n);
	TGraph* gr = new TGraph(n, x, y);
	pos2pixel(x, n); pos2pixel(y, n);
	// gr->Draw("apsame");
	// Create TF1 Object
	TF1* tf = new TF1("linear_fit","[0] + [1]*x", xmin, xmax);
	
	// Estimate slope and intercept
	double slopeEstimate = (y[int(3*n/4)]-y[int(n/4)])/(x[int(3*n/4)]-x[int(n/4)]);
	double interceptEstimate = y[0] - slopeEstimate*x[0];
	tf->SetParameter(interceptEstimate, slopeEstimate);
	
	// Fit the function tf to the data
	gr->Fit(tf, "SQR"); 	
	tf->SetLineColor(1);
	tf->SetLineWidth(4);	
	delete gr;
	return tf;
}

void  getZ(double *x, double *y, double *z, int n){

	double xmin = getArrayMin(x, n);
	double xmax = getArrayMax(x, n);
	position ip = getInitialPosition(x, y, n);
	double xi = ip.x; double yi = ip.y;

	// Get the best fit of the data
	TF1* fit;
	fit = fitMuonLine(x, y, n);

	// Calculate the slope in the z direction
	double slope = fit->GetParameter(1);
	double xrange = xmax - xmin;
	double yrange = xrange*slope;
	double length = sqrt(pow(xrange,2) + pow(yrange,2));
	double zslope = CCDWidth / length;

	// Project pixel values onto the line. All vectors are computed relative to the inital point of the 
	// s is the vector we are projecting
	double vx, vy, sx, sy, r;
	vx = 1/sqrt(1 + pow(slope,2));
	vy = slope/sqrt(1 + pow(slope,2));
	for(int i=0; i<n; i++){
		sx = x[i] + 0.5 - xi;
		sy = y[i] + 0.5 - yi;
		r = vx*sx + vy*sy;
		z[i] = abs(r)*zslope;
	}
	delete fit;
	return;

}
position::position(){
	this->x = 0.;
	this->y = 0.;
}

position::position(double x, double y){
	this->x = x;
	this->y = y;
}

position getInitialPosition(double *x, double *y, int n){

	double xmin = getArrayMin(x, n);
	double xmax = getArrayMax(x, n);
	double xmean = xmin + (xmax - xmin)/2;
	
	// Compute the weighted x position
	double xw = 0;
	for(int i=0; i<n; i++){
		xw += x[i];
	}	
	xw /= n;

	// Gets the slope track because it is necessary for distiguishing to inital values
	TF1 * tf = fitMuonLine(x, y, n);
	double slope = tf->GetParameter(1);

	position ip;
	if(xw > xmean){
		ip.x = getArrayMin(x, n);
		ip.y = tf->Eval(ip.x);
	}else if (xw < xmean){
		ip.x = getArrayMax(x, n);
		ip.y = tf->Eval(ip.x);
	}else{
		ip.x = x[0];
		ip.y = y[0];
	}  
	delete tf;	
	return ip;
}

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

// Helper function to do the computation calculating the distance from the best fit line. Boolean option to exclude delta rays with threshold cutoff of distance
// Inputs:
// double *x - array of x positions
// double *y - array of y positions
// double *q - array of the value of the CCD
// int n - number of entries in the x and y arrays
// double zmin - minimum z threshold
// double zmax - maximum z threshold
// double dx - spacing between each grid point of data measured in number of pixels 
// bool resolveDeltaRay - if true, reject events where delta rays are detected
// 
// Outputs:
// double *proj - array of distance from muon track line
// double *q - array of the amount of charge in each pixel
// double dedx - pass by reference value of the mean dedx of the line segment
// int zcount - pass by reference value of the number of pixels that are in the correct zrange

void getDistanceFromTrack(double *x, double *y, double *q, int n, double zmin, double zmax, bool resolveDeltaRay, double *proj, double *qEnergy, double &dedx, int &zcount, double dx){
	
	bool deltaray = false;
	double sx, sy, projection;
	// Get information about the line
	TF1* tf = fitMuonLine(x, y, n);
	double slope = -1/tf->GetParameter(1);
	double *z = new double[n];
	getZ(x, y, z, n);
	double vx = 1/sqrt(1 + pow(slope,2));
	double vy = slope/sqrt(1 + pow(slope,2));
	position ip = getInitialPosition(x, y, n);
	
	// Calculate the length of the line segment
	double xmin = getArrayMin(x, n); double xmax = getArrayMax(x, n);
	double length = sqrt(pow(xmax-xmin,2) + pow(tf->Eval(xmax)-tf->Eval(xmin),2))*(zmax-zmin)/CCDWidth;
	// Iterate over all x,y pairs in the cluster
	for(int j=0; j<n; j++){
		sx = x[j] - ip.x + dx/2;
		sy = y[j] - ip.y + dx/2;
		projection = vx*sx + vy*sy;
		// if in the depth range, add to histogram
		if(z[j] > zmin && z[j] < zmax){
			proj[zcount] = projection;
			qEnergy[zcount] = q[j];
			dedx += q[j];
			zcount++;
		}
		if(abs(projection) >= 4){
			deltaray = true;
		}

	}
	delete z; delete tf;	
	// Find the energy per  unit length
	dedx /= length;
	// If delta ray detected set the zcount to 0. Then when we fill the histogram those segments are ignored
	if(deltaray && resolveDeltaRay) zcount=0;

	return;
}

// Overloaded function from above with minor changes. Everything not listed below has the same value/purpose
// Inputs:
// double *dedxIn - array of the dedx values associated with each pixel. Comes from the dedx filtering
// Outputs:
// double *dedxOut - output array of dedx values associated with each pixel
void getDistanceFromTrack(double *x, double *y, double *q, double *dedxIn, int n, double zmin, double zmax, bool resolveDeltaRay, double *proj, double *qEnergy, double *dedxOut, int &zcount, double dx){

	bool deltaray = false;
	double sx, sy, projection;
	// Get information about the line
	TF1* tf = fitMuonLine(x, y, n);
	double slope = -1/tf->GetParameter(1);
	double *z = new double[n];
	getZ(x, y, z, n);
	double vx = 1/sqrt(1 + pow(slope,2));
	double vy = slope/sqrt(1 + pow(slope,2));
	position ip = getInitialPosition(x, y, n);
	
	// Calculate the length of the line segment
	double xmin = getArrayMin(x, n); double xmax = getArrayMax(x, n);
	double length = sqrt(pow(xmax-xmin,2) + pow(tf->Eval(xmax)-tf->Eval(xmin),2))*(zmax-zmin)/CCDWidth;
	// Iterate over all x,y pairs in the cluster
	for(int j=0; j<n; j++){
		sx = x[j] - ip.x + dx/2;
		sy = y[j] - ip.y + dx/2;
		projection = vx*sx + vy*sy;
		// if in the depth range, add to histogram
		if(z[j] > zmin && z[j] < zmax){
			proj[zcount] = projection;
			qEnergy[zcount] = q[j];
			dedxOut[zcount] = dedxIn[j];
			zcount++;
		}
		if(abs(projection) >= 4){
			deltaray = true;
		}

	}
	delete z; delete tf;	
	
	// If delta ray detected set the zcount to 0. Then when we fill the histogram those segments are ignored
	if(deltaray && resolveDeltaRay) zcount=0;

	return;

}

// Wrapper function to extract x, y, and q arrays from a tree entry
void getXYE(TTree *tree, int i, double *xx, double *yy, double *qq, int &n){

	// Define parameters
	TArrayD *x = new TArrayD(); TArrayD *y = new TArrayD(); TArrayD *q = new TArrayD();
	double conversionFactor = 10300/6.4;
	
	// Set the appropriate branch address
	tree->SetBranchAddress("pixel_x", &x);
	tree->SetBranchAddress("pixel_y", &y);
	tree->SetBranchAddress("pixel_val", &q);


	tree->GetEntry(i);

	// Conver TArrayD to double arrays
	n = x->GetSize();
	for(int j=0; j<n; j++){
		xx[j] = x->GetAt(j); 
		yy[j] = y->GetAt(j); 
		qq[j] = q->GetAt(j) / conversionFactor; 

	}
	
	// Remove TArrayD
	delete x; delete y; delete q;
	return;
}
