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
	int projMax=1, projMin=-1, dedxMax=1;
	double vx, vy, sx, sy, inter, slope, projection, xmin, xmax, length, projArrayMax, projArrayMin;
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
		projArrayMin = getArrayMin(projArray, zcount);
		if(projArrayMax > projMax){
			projMax = int(projArrayMax);
			h2->SetBins(2*(projMax-projMin), -projMin, projMax, dedxMax, 0, dedxMax);
		}
		if(projArrayMin < projMin){
			projMin = int(projArrayMin);
			h2->SetBins(2*(projMax-projMin), -projMin, projMax, dedxMax, 0, dedxMax);
		}
		if(dedx > dedxMax){
			dedxMax = int(dedx);
			h2->SetBins(2*(projMax-projMin), -projMin, projMax, dedxMax, 0, dedxMax);
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

TH1D* histDistance(TTree *tree, double zmin, double zmax, bool dedxFilt, bool energyFilt, double emin, double emax, bool deltaRayRejection, bool draw){

	// Define parameters
	TH1D *h1 = new TH1D("t", "Spread of Muon Tracks", 80, -4, 4);
	TStyle *gStyle = new TStyle();

	// Defining variables
	
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
	TF1* f = new TF1("f1", "gaus", -1.5, 1.5);
	f->SetParameter(1,0);
	h1->Fit("f1", "QI0");
	if(draw){
		gStyle->SetOptFit();
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
TGraph* sigmaVDepth(TTree *tree, double deltaZ, double zstart, double zend, double emin, double emax){
	
	// Define parameters
	double z[500];
	double sigma[500];
	double sigmaerr[500];
	int tGraphCnt = 0;
	int zmin, zmax;
	int nSlices = int((zend-zstart)/deltaZ);

	// Iterate over all slices, getting the diffusion spread
	for(int i=0; i<nSlices; i++){
		zmin = i*deltaZ + zstart; zmax = (i+1)*deltaZ + zstart;
		TH1D* h = new TH1D();
		h = histDistance(tree, zmin, zmax, true, emin, emax, true, false);
		z[i] = (zmin + zmax)/2;
		cout << zmin << "-" << zmax << ":" << endl;
		sigma[i] =  h->GetFunction("f1")->GetParameter(2);
		sigmaerr[i] = h->GetFunction("f1")->GetParError(2);  
		delete h;	
	}
	TGraphErrors *tg = new TGraphErrors(nSlices, z, sigma, NULL, sigmaerr);
	
	// Get fit the graph	
	// Fit expected function to data
	double q = 1.6e-19; // C
	double T = 130.; // K 
	double kb = 1.38e-11; // um^2 kg / K s^2 
	double epsSi = 11.68 * 8.85e-18; // C / Vum 
	double Vapp = 128.;
	double rho = 3.9e-1*q; // um^-3 
	double conversion =1.e-6;
	
	// TF1 *f = new TF1("sig_analytic", "TMath::Sqrt( [0]/[5]*TMath::Log([1]-x/([2]/[5] + [3])))*[4]", 1, CCDWidth);
	TF1 *f = new TF1("sig_analytic", "[3]*TMath::Sqrt( [0]*TMath::Log([2]-x*[1]))", 0, CCDWidth);	
	
	f->SetParameter(0, -2*kb*T*epsSi/(rho*q));
	f->SetParameter(1, 1/(epsSi*Vapp/(rho*CCDWidth)+0.5*CCDWidth));
	f->FixParameter(2, 1);
	f->FixParameter(3, conversion);
	//return fit pointer;	
	char titlestr[200];
	sprintf(titlestr, "#sigma_{x} vs Depth. E_{min}=%.0f, E_{max}=%.0f.", emin, emax);
	tg->SetTitle(titlestr);
	tg->GetXaxis()->SetTitle("Depth (#mum)");
	tg->GetYaxis()->SetTitle("#sigma_{x} (#mum)");
	tg->Fit(f, "SR0");
	
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
	tg->SetTitle("Diffusion #sigma_{x} vs. Depth of Interaction");
	tg->GetXaxis()->SetTitle("Depth (#mum)");
	tg->GetYaxis()->SetTitle("#sigma_{x} (#mum)");
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
	
	TF1 *f = new TF1("sig_analytic", "[3]*TMath::Sqrt( [0]*TMath::Log([2]-x*[1]))", 0, CCDWidth);	
	
	f->SetParameter(0, -2*kb*T*epsSi/(rho*q));
	f->SetParameter(1, 1/(epsSi*Vapp/(rho*CCDWidth)+0.5*CCDWidth));
	f->FixParameter(2, 1);
	f->FixParameter(3, conversion);
	tg->Fit(f, "R");
	cout << "Rho: " << f->GetParameter(5) << endl;
	cout << f->Eval(100) << endl;
	cout << f->Eval(200) << endl;

	return tg;
}

// Function to look at the variation of de/dx as a function of track length
TH1D* dedxFluctuation(TTree *tree, int i){
	
	// Define parameters
	int n = 0;
	double *x = new double[30000];
	double *y = new double[30000];
	double *energy = new double[30000];
	
	// Get array of x, y pixels
	getXYE(tree, i, x, y, energy, n);
	
	// Get slope of line and find track length
	TF1 *tf = fitMuonLine(x, y, n);
	double slope = tf->GetParameter(1);
	double xrange = getArrayMax(x, n) - getArrayMin(x, n);
	double yrange = xrange*slope;
	double length = sqrt(pow(xrange, 2) + pow(yrange, 2));

	position ip = getInitialPosition(x, y, n);
	double xi = ip.x; double yi = ip.y;

	// Create histogram
	TH1D *h = new TH1D("h", "dE/dx v Position", int(length), 0, length);
	TStyle *gStyle = new TStyle();

	// Iterate over all points, find the position and energy, and add it to the histogram
	double sx, sy, r;
	double vx = 1 / sqrt(1 + pow(slope, 2));
	double vy = slope / sqrt(1 + pow(slope, 2));
	for(int i=0; i<n; i++){
		sx = x[i] + 0.5 - xi;
		sy = y[i] + 0.5 - yi;
		r = vx*sx + vy*sy;

		h->Fill(abs(r), energy[i]);
	}

	// Setting histogram parameters
	h->SetTitle("Approx. #frac{dE}{dx} vs Position along the track");
	h->GetXaxis()->SetTitle("Position Along Track (pixels)");
	h->GetYaxis()->SetTitle("#frac{dE}{dx} (keV/pixel)");
	gStyle->SetOptStat(0);
	return h; 
}

void dedxFilterTree(TTree *tree, const char *outfile, double dedxThresh, int startIdx){

	TFile *outfileRT = new TFile(outfile, "RECREATE");
	TTree *dedxTree = new TTree("dedxFilterTree","filter by dedx");

	// Define parameters to be used for the branches of the new tree
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> q;
	vector<double> dedx;
	int trackID, trackLength;
	double xi, yi, slope, slopeNew, xiNew, yiNew;
	position ip, ipNew;


	// Add branches to the new tree
	dedxTree->Branch("trackID", &trackID);
	dedxTree->Branch("x", &x);
	dedxTree->Branch("y", &y);
	dedxTree->Branch("z", &z);
	dedxTree->Branch("q", &q);
	dedxTree->Branch("dedx", &dedx);
	dedxTree->Branch("xi", &xi);
	dedxTree->Branch("yi", &yi);
	dedxTree->Branch("slope", &slope);
	dedxTree->Branch("trackLength", &trackLength);
	dedxTree->Branch("slopeNew", &slopeNew);
	dedxTree->Branch("xiNew", &xiNew);
	dedxTree->Branch("yiNew", &yiNew);
	
	// Set Parameters for the current treei
	double *xxold = new double[20000];
	double *yyold = new double[20000];
	double *qqold = new double[20000];
	int nold = 0;
	int nTracks = tree->GetEntries();
	
	// Creating other variables needed for storing "good" events
	TF1 *tf, *tfNew;
	int nGoodSlices;
	double delta;
	vector<double> goodDepth;
	vector<double> sliceDedx;

	// Iterate over all entries, adding values to the new trees
	for(int i=startIdx; i<nTracks; i++){
		getXYE(tree, i, xxold, yyold, qqold, nold);
		// Get a histogram of the energy fluctuation
		TH1D * energyfluc = new TH1D();
		energyfluc = dedxFluctuation(tree, i);
		trackLength = energyfluc->GetNbinsX();
		// Iterate over the histogram and find list of good depths
		// Definition of good is current value and ajacent values are less than threshold
		for(int j=0; j<trackLength; j++){
		
			if(j == 0){
				if(energyfluc->GetBinContent(j) < dedxThresh && energyfluc->GetBinContent(j+1) < dedxThresh){
				goodDepth.push_back((j+0.5)*CCDWidth/trackLength);
				sliceDedx.push_back(energyfluc->GetBinContent(j));
				}
			}else if(j == (trackLength -1)){
				if(energyfluc->GetBinContent(j) < dedxThresh && energyfluc->GetBinContent(j-1) < dedxThresh){
				goodDepth.push_back((j+0.5)*CCDWidth/trackLength);
				sliceDedx.push_back(energyfluc->GetBinContent(j));
				}
			}else{
				if(energyfluc->GetBinContent(j) < dedxThresh && energyfluc->GetBinContent(j+1) < dedxThresh &&energyfluc->GetBinContent(j-1) < dedxThresh){
				goodDepth.push_back((j+0.5)*CCDWidth/trackLength);
				sliceDedx.push_back(energyfluc->GetBinContent(j));
				}
			}
		}
		// Calculate the +/- range in depth. Delta D/2
		delta = 0.5*CCDWidth/trackLength;
		// Compare the depths of the actual pixels to the acceptable ones
		// Put good depths in the new tree
	
		double *zzold = new double[nold];
		getZ(xxold, yyold, zzold, nold);
		for(int j=0; j<nold; j++){
			
			// Compare each pixel to the acceptable ranges
			for(int k=0; k<goodDepth.size(); k++){
				if(zzold[j] < (goodDepth[k] + delta) && zzold[j] > (goodDepth[k] - delta)){
					x.push_back(xxold[j]);
					y.push_back(yyold[j]);
					z.push_back(zzold[j]);
					q.push_back(qqold[j]);
					dedx.push_back(sliceDedx[k]);
					break;
				}

			}
		}

		// Set remaining branch parameters
		// Simple parameters
		trackID = i;
		ip = getInitialPosition(xxold, yyold, nold);
		xi = ip.x; yi = ip.y;
		tf = fitMuonLine(xxold, yyold, nold);
		slope = tf->GetParameter(1);

		// Find the new best fit line with the remaining 
		tfNew = fitMuonLine(&x[0], &y[0], x.size());
		slopeNew = tfNew->GetParameter(1);
		ipNew = getInitialPosition(&x[0], &y[0], x.size());
		xiNew = ipNew.x; yiNew = ipNew.y;

		// Fill the tree 
		dedxTree->Fill();

		// Clear all the necessary vectors
		delete zzold;
		delete energyfluc;
		x.clear();
		y.clear();
		z.clear();
		q.clear();
		dedx.clear();
		goodDepth.clear(); sliceDedx.clear();
	}
	dedxTree->Write();
	outfileRT->Close();
	return;
}

void testdedxFilter(TTree *tree, int i=0){

	gStyle->SetPalette(55);
	TCanvas *c = new TCanvas("c","c", 700, 500);
	c->Divide(1,3);
	gStyle->SetOptStat(0);
	c->cd(1);
	TArrayD *x; TArrayD *y; TArrayD *q;
	tree->SetBranchAddress("pixel_x", &x);
	tree->SetBranchAddress("pixel_y", &y);
	tree->SetBranchAddress("pixel_val", &q);

	tree->GetEntry(i);
	TH2D *rawTrack = plot2DTrack(x, y, q, false);
	rawTrack->Draw("colz");

	c->cd(2);
	TH1D *fluc = dedxFluctuation(tree, i);
	fluc->Draw("hist");

	c->cd(3);
	const char *filename = "./test.root";
	dedxFilterTree(tree, filename, 8, i);
	TFile *fdedx = new TFile(filename);
	TTree *tdedx = (TTree*)fdedx->Get("dedxFilterTree");
	tdedx->Draw("y:x","q","colz");
	TH1D *filtTrack = (TH1D*)gROOT->FindObject("htemp");
	double xmin=rawTrack->GetXaxis()->GetXmin();
	double xmax=rawTrack->GetXaxis()->GetXmax();
	double ymin=rawTrack->GetYaxis()->GetXmin();
	double ymax=rawTrack->GetYaxis()->GetXmax();
	filtTrack->SetBins(rawTrack->GetNbinsX(), xmin, xmax, rawTrack->GetNbinsY(), ymin, ymax);
	filtTrack->Draw("colz");

	return;
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
