double skewedGaussian(double *x, double *par){
	// Defining the formula for a skewed gaussian. f(x) = 2*phi(x)*Phi(a*x).
	// https://en.wikipedia.org/wiki/Skew_normal_distribution
	double xprime = (x[0] - par[1]) / par[2];
	double fitVal = par[0]*TMath::Exp((-0.5 * TMath::Power(xprime, 2))) * (1 + TMath::Erf(par[3] * xprime / TMath::Sqrt(2)));
	return fitVal;
}

void plotSkewedGaussian(double N=1, double mean=0, double sigma=1, double alpha=0){

	TF1 *func = new TF1("testSkewed", skewedGaussian, -3, 3, 4);
	func->SetParameters(N, mean, sigma, alpha);
	func->Draw();
	double distMean = mean + sigma * alpha * TMath::Sqrt(2/(TMath::Pi() * (1 + alpha*alpha)));
	cout << "Mean: " <<  distMean << endl;;
	return;
}

double computeSkewedMean(double mean, double sigma, double alpha){

	return mean + sigma * alpha * TMath::Sqrt(2/(TMath::Pi() * (1 + alpha*alpha)));
}

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
	TF1 *sigmamaxFit = new TF1("sigmamaxFit", "gaus", 8, 10);
	// TF1 *sigmamaxFit = new TF1("sigmamaxFit", "[0]*TMath::Exp(-[1])*TMath::Power([1],x)/TMath::Gamma(x+1)", 0, 10);
	// sigmamaxFit->SetParameter(1, 8.5);
	// sigmamaxFit->SetParameter(0, 260);
	sigmamax->Fit("sigmamaxFit", "QIR0");

	// Create skewed gaussian function to fit
	TF1 *sigmamaxSkew = new TF1("sigmaSkew", skewedGaussian, 8, 10, 4);
	sigmamaxSkew->SetParameters(sigmamaxFit->GetParameter(0), sigmamaxFit->GetParameter(1), sigmamaxFit->GetParameter(2), 2);
	sigmamaxSkew->SetParNames("N", "mean", "sigma", "skew");
	sigmamax->Fit("sigmaSkew", "QIR");

	
	sigmamax->Draw();
	// sigmamaxFit->Draw("same");
	sigmamaxSkew->Draw("same");

	// Plot a and fit
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	a->SetLineWidth(2);
	TF1 *fA = new TF1("aFit", "gaus", -3.5, 1.5);
	a->Fit("aFit", "QIR");

	// Make the skewed guassian

	TF1 *fASkew = new TF1("aSkew", skewedGaussian, -4, 1, 4);
	fASkew->SetParameters(fA->GetParameter(0), fA->GetParameter(1), fA->GetParameter(2), 1);
	fASkew->SetParNames("N", "mean", "sigma", "skew");
	a->Fit("aSkew", "QIR");
	
	a->Draw();
	fASkew->Draw("same");
	cout << sigmamaxFit->GetParameter(1) << endl;	
	cout << computeSkewedMean(fASkew->GetParameter(1), fASkew->GetParameter(2), fASkew->GetParameter(3)) << endl;
	cout << computeSkewedMean(sigmamaxSkew->GetParameter(1), sigmamaxSkew->GetParameter(2), sigmamaxSkew->GetParameter(3)) << endl;
	// cout << fASkew->GetMaximumX() << endl;
	
	// Compute diffusion parameters
	double b = TMath::Exp(computeSkewedMean(fASkew->GetParameter(1), fASkew->GetParameter(2), fASkew->GetParameter(3))) / zd;
	// double A = -TMath::Power(sigmamaxFit->GetParameter(1), 2) / TMath::Log(1 - b * zd);
	double A = -TMath::Power(8.48, 2) / TMath::Log(1 - b * zd);
	cout << "A: " << A << endl;
	cout << "b: " << b << endl;

	// Plot comparison of diffusion parameters
	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);

	// neutron data diffusion
	TF1 *fn = new TF1("fNeutron", "TMath::Sqrt([0]*TMath::Log([2]-x*[1]))", 0, 2.5*zd);		
	fn->SetParameter(0, -ANeutron);
	fn->SetParameter(1, bNeutron);
	fn->SetParameter(2, 1);
	fn->SetLineWidth(2);
	fn->SetLineColor(kRed);
	fn->Draw();
	fn->SetTitle("Diffusion Model");

	// muon data diffusion
	TF1 *fm = new TF1("fMuon", "TMath::Sqrt([0]*TMath::Log([2]-x*[1]))", 0, 2.5*zd);		
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