{
	
	TFile *f = TFile::Open("../../rootfiles/dedxfilter_8kev.root", "READ");
	if ((!f) || (f->IsZombie())) {delete f; return;} //
	TTree * t; f->GetObject("dedxFilterTree", t);
	TH2D * h2 = new TH2D("h2", "h2", 4000, 0, 4000, 2000, 0, 2000);
	
	vector<double> *x;
	vector<double> *y;
	t->SetBranchAddress("x", &x);
	t->SetBranchAddress("y", &y);


}