void plotDiffusionParameters(double ANeutron=220, double bNeutron=5.6E-4){

	const double zd = 500;
	gStyle->SetOptFit(1);
	string infile = "../../rootfiles/likelihoodFit/muonTree15kev_0999ccf_fits.root";
	TFile *f = new TFile(infile.c_str());
	TTree *t = (TTree*)f->Get("muonFit");


	// Create and fill histograms from tree data
	TH1D *sigmamax = new TH1D("sigmamax","#sigma_{max}", 200, 0, 10);
	TH1D *a = new TH1D("a", "a", 100, -10, 10);

	double smaxVal, aVal;
	t->SetBranchAddress("sigmamax", &smaxVal);
	t->SetBranchAddress("a", &aVal);

	for(int i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		sigmamax->Fill(smaxVal);
		a->Fill(aVal);
	}

	// Plot sigma max and fit
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	sigmamax->SetLineWidth(2);
	TF1 *sigmamaxFit = new TF1("sigmamaxFit", "gaus", 7.5, 10);
	// TF1 *sigmamaxFit = new TF1("sigmamaxFit", "[0]*TMath::Exp(-[1])*TMath::Power([1],x)/TMath::Gamma(x+1)", 0, 10);
	// sigmamaxFit->SetParameter(1, 8.5);
	// sigmamaxFit->SetParameter(0, 260);
	sigmamax->Fit("sigmamaxFit", "QIR");

	
	sigmamax->Draw();
	sigmamaxFit->Draw("same");

	// Plot a and fit
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	a->SetLineWidth(2);
	TF1 *fA = new TF1("aFit", "gaus", -3.5, 1.5);
	// TF1 *fA = new TF1("aFit", "[0]*TMath::Exp(-[1])*TMath::Power([1],-x)/TMath::Gamma(TMath::Abs(x)+1)", -3, 0);
	// fA->SetParameter(1, 0.9);
	// fA->SetParameter(0, 200);
	a->Fit("aFit", "QIR");
	
	c1->cd(2);
	a->Draw();
	fA->Draw("same");
	cout << sigmamaxFit->GetParameter(1) << endl;
	cout << fA->GetParameter(1) << endl;
	
	// Compute diffusion parameters
	double b = TMath::Exp(-fA->GetParameter(1)) / zd;
	double A = -TMath::Power(sigmamaxFit->GetParameter(1), 2) / TMath::Log(1 - b * zd);
	cout << "A: " << A << endl;
	cout << "b: " << b << endl;

	// Plot comparison of diffusion parameters
	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);

	// neutron data diffusion
	TF1 *fn = new TF1("fNeutron", "TMath::Sqrt( [0]*TMath::Log([2]-x*[1]))", 0, zd);		
	fn->SetParameter(0, -ANeutron);
	fn->SetParameter(1, bNeutron);
	fn->SetParameter(2, 1);
	fn->SetLineWidth(2);
	fn->SetLineColor(kRed);
	fn->Draw();
	fn->SetTitle("Diffusion Model");

	// muon data diffusion
	TF1 *fm = new TF1("fMuon", "TMath::Sqrt( [0]*TMath::Log([2]-x*[1]))", 0, zd);		
	fm->SetParameter(0, -A);
	fm->SetParameter(1, b);
	fm->SetParameter(2, 1);
	fm->SetLineWidth(2);
	fm->SetLineColor(kBlue);
	fm->Draw("same");

	// Make Legend
	TLegend *leg = new TLegend(0.42, 0.15, 0.88, 0.35);
	char fMuonLeg[200], fNeutronLeg[200];
	sprintf(fMuonLeg, "Muon Fit. A=%.1e #mum^{2}, b=%.1e #mum^{-1}", A, b);
	sprintf(fNeutronLeg, "Neutron Fit. A=%.1e #mum^{2}, b=%.1e #mum^{-1}", ANeutron, bNeutron);
	leg->AddEntry("fMuon", fMuonLeg);
	leg->AddEntry("fNeutron", fNeutronLeg);
	leg->Draw();

	return;
}