#ifndef muonFit_H
#define muonFit_H

#include "TH2F.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include "TGraph.h"
#include "TF1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TArray.h"
#include "TChain.h"
#include "TParameter.h"
#include "TFile.h"
#include <fstream>

const Double_t zd = 675; //thickness of the CCD in microns
const Double_t ps = 15; //pixel size in um
const Double_t cyield = 70; //yield (in electrons per um) of region with ~constant deposition


TH2F* ArraysToTH2F(TArrayD* x, TArrayD* y, TArrayD* val);

Double_t maxk(Double_t obs, Double_t ne, Double_t c);

Double_t fitf(Double_t z, Double_t* par);

double TrackLogLikelihood(const double *xx);

void muonFit(TH2F* cl, TTree *tMuon = nullptr, Double_t *b = nullptr, Double_t c = 5.5E-4, Double_t sadc_ = 28, Double_t sthr_ = 2);

void muonFitMany(TChain* t, TTree *tMuon = nullptr, Double_t *b = nullptr, Double_t c = 5.5E-4, Double_t sadc_ = 18, Double_t sthr_ = 3);

void saveMuonFitToFile(const char* outfile, const char *infile, Double_t c = 5.5E-4, Double_t sadc_ = 18, Double_t sthr_ = 3);

#endif