#include "muonFilter.h"
#include "TStyle.h"
#include "TObject.h"
#include "TGraphErrors.h"
#include <vector>

void convertVal2Energy(double *q, int n,  double conversionFactor=10300/6.4);

TH2D* histEnergyvDistance(TTree *tree, double zmin=400, double zmax=425, bool resolveDeltaRay=true);

TH1D* histDistance(TTree *tree, double zmin=400, double zmax=450, bool dedxFilt=true, bool energyFilt=false, double emin=2., double emax=4., bool deltaRayRejection=true, bool draw=true);

TGraph* sigmaVDepth(TTree *tree, double deltaZ, double zstart, double zend,  double emin=3., double emax=4.);	

TGraph* plotSigmaVDepth();

TH1D* dedxFluctuation(TTree *tree, int i);

void dedxFilterTree(TTree *tree, const char *outfile, double dedxThresh, int startIdx=0);

void saveAllHist(TTree *tree);

void savehist(TH1D *h, const char *filename, char *histname);

void savehist(TH2D *h, const char *filename, char *histname);

void savegraph(TGraph *graph, const char *filename, char* objname);
