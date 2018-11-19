#include "muonFilter.h"
#include "analysis.h"

using namespace std;

int main(int argc, char **argv){

	// Get tree of muon data
	TFile* f = new TFile("~/muontracks2.root");
	TTree* t = (TTree*)f->Get("clusters_tree");
	cout << "Number of tree entries: " << t->GetEntries() <<  endl;
	
	// Create TApplication object and run analysis
	const char *outfile = "~/DAMICDiffusion/rootfile/sigma_v_depth_25um.root";

	// Iterate over a range of energies and plot
	int estart=4; int eend=5; double dz=10;
	TGraph *tg;
	for(int e=estart; e<=eend; e++){
		tg = sigmaVDepth(t, dz, 0, 500, e, e+1);
		// Create the name for the graph
		char *histname = new char[200];
		sprintf(histname, "sd_%i_%ikev", e, e+1);
		savegraph(tg, outfile, histname);
		delete histname;
	}

	cout << "Analysis successfully run." << endl;

	return 0;
}







