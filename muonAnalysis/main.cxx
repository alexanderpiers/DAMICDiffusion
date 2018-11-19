#include "muonFilter.h"
#include "analysis.h"

using namespace std;

int main(int argc, char **argv){

	// Get tree of muon data
	TFile* f = new TFile("~/muontracks.root");
	TTree* t = (TTree*)f->Get("clusters_tree");
	cout << "Number of tree entries: " << t->GetEntries() <<  endl;
	
	// Create TApplication object and run analysis
	const char *outfile = "~/DAMICDiffusion/rootfile/sigma_v_depth_25um.root";

	// Iterate over a range of energies and plot
	int estart=2; int eend=10; double dz=25;
	TGraph *tg;
	cout << "1" <<endl;
	for(int e=estart; e<=eend; e++){
		tg = sigmaVDepth(t, dz, 0, e, e+1);
		cout << "2" << endl;
		// Create the name for the graph
		char *histname = new char[200];
		sprintf(histname, "sd_%i_%ikev", e, e+1);
		cout << "3" << endl;
		savegraph(tg, outfile, histname);
		cout << "4" << endl;
		delete histname;
	}

	cout << "Analysis successfully run." << endl;

	return 0;
}







