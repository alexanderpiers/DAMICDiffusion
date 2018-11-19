#ifndef muonFilter_H
#define muonFilter_H

#include <iostream>
#include <cmath>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TImage.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TChain.h"
#include "TParameter.h"
#include "TArrayD.h"


const double CCDWidth = 500; // CCD width in um

void muonFilter(TChain* chain, char const * outfile, double minEnergy=500., double minccf=0.99);

TH2D* plot2DTrackDepth(TArrayD *x, TArrayD *y, bool plot=true);
	
TH2D* plot2DTrack(TArrayD *x, TArrayD *y, TArrayD *q, bool plot=true);
	
void  pixel2pos(double* pixel, int n);

void pos2pixel(double* pixel, int n);

TF1 * fitMuonLine(double *x, double *y, int n);

void  getZ(double *x, double *y, double *z, int n);

struct position{
	double x;
	double y;
};

position getInitialPosition(double *x, double *y, int n);

double getArrayMax(double* arr, int n);

double getArrayMin(double* arr, int n);

void  getDistanceFromTrack(double *x, double *y, double *q, int n, double zmin, double zmax,  bool resolveDeltaRay, double *proj, double *qEnergy, double &dedx, int &zcount);

#endif
