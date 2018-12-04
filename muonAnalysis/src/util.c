
{
	gInterpreter->AddIncludePath("~/DAMICDiffusion/muonAnalysis/include/");
	gROOT->ProcessLine(".L ~/DAMICDiffusion/muonAnalysis/src/muonFilter.cxx");
	gROOT->ProcessLine(".L ~/DAMICDiffusion/muonAnalysis/src/analysis.cxx");
	
	TFile* f = new TFile("~/muontracks2.root");
	TTree* t = f->Get("clusters_tree");

	TArrayD *xTA, *yTA, *qTA;
	t->SetBranchAddress("pixel_x", &xTA);
	t->SetBranchAddress("pixel_y", &yTA);
	t->SetBranchAddress("pixel_val", &qTA);
	
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
