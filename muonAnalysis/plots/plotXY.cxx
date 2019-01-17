void plotXY(vector<double> *x, vector<double> *y, TH2D *h2){

	// TFile *f = TFile::Open("../../rootfiles/dedxfilter_8kev.root", "READ");
	// if ((!f) || (f->IsZombie())) {delete f; return;} //
	// TTree * t; f->GetObject("dedxFilterTree", t);
	// TH2D * h2 = new TH2D("h2", "h2", 4000, 0, 4000, 2000, 0, 2000);
	
	// vector<double> *x;
	// vector<double> *y;
	// cout << -1 << endl;
	// t->SetBranchAddress("x", &x);
	// t->SetBranchAddress("y", &y);
	// cout << "test" << endl;
	// t->GetEntry(35);

	double xmin, ymin, xmax, ymax;

	xmin = *min_element(x->begin(), x->end());
	xmax = *max_element(x->begin(), x->end());
	ymin = *min_element(y->begin(), y->end());
	ymax = *max_element(y->begin(), y->end());

	h2->GetXaxis()->SetRangeUser(xmin, xmax);
	h2->GetYaxis()->SetRangeUser(ymin, ymax);
	cout << 1 << endl;

	// char *logicString = new char[100];
	// // sprintf(logicString, "z && (trackID == %i)", i);
	// cout << 2 << endl;
	// t->Draw("y:x>>h2", "z*(trackID==35)", "colz");

	// delete f;
	// delete h2;
	// delete x;
	// delete y;
	return;

}