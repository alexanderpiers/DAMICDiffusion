

void plotManyTracks(int start=0){

	const char *infile = "/home/apiers/mnt/mox/DAMICDiffusion/rootfiles/snolab/hist/muon_tracks_35kev_0997ccf.root";
	TFile *f = new TFile(infile, "READ");
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	char histogramName[100];
	int n=6;
	TCanvas *c = new TCanvas("c", "c", 800, 600);
	c->Divide(n,n);

	for(int i=0; i<n*n; i++){
		cout << i << endl;
		sprintf(histogramName, "muon%i", i+start);
		TH2D *h2 = new TH2D();
		h2 = (TH2D*)f->Get(histogramName);
		h2->GetXaxis()->SetLabelSize(0);
		h2->GetYaxis()->SetLabelSize(0);
		c->cd(i+1);
		h2->Draw();

	}

	return;
}