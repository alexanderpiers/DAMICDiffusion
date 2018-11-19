{
	gROOT->ProcessLine(".L ~/DAMICDiffusion/muonAnalysis/src/muonFilter.cxx");
	gROOT->ProcessLine(".L ~/DAMICDiffusion/muonAnalysis/src/analysis.cxx");

	TFile* f = new TFile("~/muontracks2.root");
	TTree* t = f->Get("clusters_tree");

	//TArrayD *x, *y, *q;
	//t->SetBranchAddress("pixel_x", &x);
	//t->SetBranchAddress("pixel_y", &y);
	//t->SetBranchAddress("pixel_val", &q);
	
	// Finding lines with good slopes
//	TF1* tf;
//	double s;
//	for(int i=0; i<t->GetEntries(); i++){
//		t->GetEntry(i);
//		tf = fitMuonLine(x,y);
//		s = abs(tf->GetParameter(1));
//		if(s > 0.99 &&  s < 1.01){
//		cout << "slope: " << s << " index: " << i << endl;
//
//		}
//	}	


}
