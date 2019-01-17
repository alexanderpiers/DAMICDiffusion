#ifndef muonFilter_H
#define muonFilter_H

#include <iostream>
#include <cmath>
#include <string>
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

#include "settings.h"

void muonFilter(TChain* chain, char const * outfile, double minEnergy=500., double minccf=0.99);

void muonFilterNoDelta(TTree *chain, char const *outfile, double minEnergy, double minccf, double minTrackLength, double maxdEdx);

TH2D* plot2DTrackDepth(TArrayD *x, TArrayD *y, bool plot=true);
	
TH2D* plot2DTrack(TArrayD *x, TArrayD *y, TArrayD *q, bool plot=true);
	
void pixel2pos(double* pixel, int n, double dx=1);

void pos2pixel(double* pixel, int n, double dx=1);

TF1 * fitMuonLine(double *x, double *y, int n, double dx=1);

void getXYE(TTree *tree, int i, double *xx, double *yy, double *qq, int &n);

void getZ(double *x, double *y, double *z, int n, double dx=1);

struct position{
	double x;
	double y;

	position();
	position(double x, double y);
};

position getInitialPosition(double *x, double *y, int n, double dx=1);

double getArrayMax(double* arr, int n);

double getArrayMin(double* arr, int n);

void getDistanceFromTrack(double *x, double *y, double *q, int n, double zmin, double zmax, bool resolveDeltaRay, double *proj, double *qEnergy, double &dedx, int &zcount, double dx=1);

void getDistanceFromTrack(double *x, double *y, double *q, double *dedxIn, int n, double zmin, double zmax, bool resolveDeltaRay, double *proj, double *qEnergy, double *dedxOut, int &zcount, double dx=1);

#endif
