void plotDiffusionParameters(char *infile, double ANeutron=220, double bNeutron=5.6E-4){

	const double zd = 500;
	gStyle->SetOptFit(1);

	TFile *f = new TFile(infile);
	TH1D *sigmamax = (TH1D*)f->Get("sigmamax");
	TH1D *a = (TH1D*)f->Get("a");

	// Plot sigma max and fit
	sigmamax->SetLineWidth(2);
	TF1 *sigmamaxFit = new TF1("sigmamaxFit", "gaus", 7.5, 10);
	sigmamax->Fit("sigmamaxFit", "QIR");

	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	sigmamax->Draw();
	sigmamaxFit->Draw("same");

	// Plot a and fit
	// a->SetLineWidth(2);
	TF1 *fA = new TF1("aFit", "gaus", -3.5, 1.5);
	a->Fit("aFit", "QIR");
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	c1->cd(2);
	a->Draw();
	aFit->Draw("same");
	
	// Compute diffusion parameters
	double b = exp(fA->GetParameter(1)) / zd;
	double A = -pow(sigmamaxFit->GetParameter(1), 2) / log(1 - b * zd);
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