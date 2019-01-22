#ifndef analysis_H
#define analysis_H

#include "muonFilter.h"
#include "settings.h"
#include "TStyle.h"
#include "TObject.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TKey.h"
#include "TCollection.h"
#include <vector>

void convertVal2Energy(double *q, int n,  double conversionFactor=10300/6.4);

TH2D* histEnergyvDistance(TTree *tree, double zmin=400, double zmax=425, bool resolveDeltaRay=true);

TH1D* histDistance(TTree *tree, double zmin=400, double zmax=450, bool dedxFilt=true, bool energyFilt=false, double emin=2., double emax=4., bool deltaRayRejection=true, bool draw=true, double dx=1, int trackThreshold=250);

TGraph* sigmaVDepth(TTree *tree, double deltaZ, double zstart, double zend,  double emin=3., double emax=4., bool dedxFilt=true, bool deltarayFilt=true);	

TGraphErrors* sigmaVDepthFile(const char *filename, double emin, double emax);

TGraph* plotSigmaVDepth();

TH1D* dedxFluctuation(TTree *tree, int i);

TH1D* dedxFluctuation(double *x, double *y, double *e, int n);

void dedxFilterTree(TTree *tree, const char *outfile, double dedxThresh, int startIdx=0);

void resampleMuonTrack(TTree *tree, const char *outfile, const int sampleRatio=5, bool dedxFilter=true);

void saveAllHist(TTree *tree);

void savehist(TH1D *h, const char *filename, char *histname);

void savehist(TH2D *h, const char *filename, char *histname);

void savegraph(TGraph *graph, const char *filename, char* objname);

#endif
