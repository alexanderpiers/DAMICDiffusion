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

void plotDiffusionParameters(bool fixA=false, double energyCorrection=0.97, double expectedSigmaMax=1, double ANeutron=220, double bNeutron=5.6E-4){

	const double zd = 675;
	gStyle->SetOptFit(1);
	string infile = "../../rootfiles/snolab/likelihoodFit/muonTree40kev_0997ccf_fits_subext_exp.root";
	// string infile = "../../rootfiles/likelihoodFit/muonTree15kev_0999ccf_fits.root";
	TFile *f = new TFile(infile.c_str());
	TTree *t = (TTree*)f->Get("muonFit");

	// Define ranges for each fit
	double sigmaRangeMin = 10;
	double sigmaRangeMax = 30;
	double sigmaFitMin = 14.5;
	double sigmaFitMax = 17;
	double aFitMin = 0;
	double aFitMax = 1.5;

	// double sigmaRangeMin = -5;
	// double sigmaRangeMax = 15;
	// double sigmaFitMin = 7;
	// double sigmaFitMax = 10;
	// double aFitMin = -3;
	// double aFitMax = 1.5;


	// Create and fill histograms from tree data
	TH1D *sigmamax = new TH1D("sigmamax","#sigma_{max}", 200, sigmaRangeMin, sigmaRangeMax);
	TH1D *a = new TH1D("a", "a", 100, -5, 5);

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
	TF1 *sigmamaxFit = new TF1("sigmamaxFit", "landau", sigmaFitMin, sigmaFitMax);
	sigmamax->Draw();

	// TF1 *sigmamaxFit = new TF1("sigmamaxFit", "[0]*TMath::Exp(-[1])*TMath::Power([1],x)/TMath::Gamma(x+1)", 0, 10);
	// sigmamaxFit->SetParameter(1, 8.5);
	// sigmamaxFit->SetParameter(0, 260);
	sigmamax->Fit("sigmamaxFit", "QIR0");
	sigmamaxFit->Draw("same");
	cout << sigmamaxFit->GetParameter(0) << "   " << sigmamaxFit->GetParameter(1) << endl;

	// Create skewed gaussian function to fit
	TF1 *sigmamaxSkew = new TF1("sigmaSkew", skewedGaussian, sigmaFitMin, sigmaFitMax, 4);
	sigmamaxSkew->SetParameters(sigmamaxFit->GetParameter(0), sigmamaxFit->GetParameter(1), sigmamaxFit->GetParameter(2), -2);
	sigmamaxSkew->SetParNames("N", "mean", "sigma", "skew");
	// sigmamax->Fit("sigmaSkew", "QIR0");
	// sigmamaxFit->Draw("same");
	// sigmamaxSkew->Draw("same");

	// Plot a and fit
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	a->SetLineWidth(2);
	TF1 *fA = new TF1("aFit", "landau", -2.5, 2.5);
	a->Fit("aFit", "QIR0");


	// Make the skewed guassian
	TF1 *fASkew = new TF1("aSkew", skewedGaussian, aFitMin, aFitMax, 4);
	fASkew->SetParameters(fA->GetParameter(0), fA->GetParameter(1), fA->GetParameter(2), -1);
	fASkew->SetParNames("N", "mean", "sigma", "skew");
	a->Fit("aSkew", "QIR0");
	a->Draw();
	fA->Draw("same");
	// fA->Draw("same");
	// fASkew->Draw("same");
	
	// Compute diffusion parameters
	double A, b;

	// set up the scaling for getting to expected sigma max
	if(expectedSigmaMax != 1){
		energyCorrection = expectedSigmaMax/(computeSkewedMean(sigmamaxSkew->GetParameter(1), sigmamaxSkew->GetParameter(2), sigmamaxSkew->GetParameter(3))/15);
		// energyCorrection = expectedSigmaMax/(sigmamaxFit->GetParameter(1)/15); 
	}
	
	cout << "a: " << computeSkewedMean(fASkew->GetParameter(1), fASkew->GetParameter(2), fASkew->GetParameter(3)) << endl;
	cout << "sigmamax: " << energyCorrection*computeSkewedMean(sigmamaxSkew->GetParameter(1), sigmamaxSkew->GetParameter(2), sigmamaxSkew->GetParameter(3))/15 << endl;
	// cout << "a: " << fA->GetParameter(1) << endl;
	// cout << "sigmamax (pixels): " << energyCorrection*sigmamaxFit->GetParameter(1)/15 << endl;

	b = TMath::Exp(-computeSkewedMean(fASkew->GetParameter(1), fASkew->GetParameter(2), fASkew->GetParameter(3))) / zd;
	A = -TMath::Power(energyCorrection*computeSkewedMean(sigmamaxSkew->GetParameter(1), sigmamaxSkew->GetParameter(2), sigmamaxSkew->GetParameter(3)), 2) / TMath::Log(1 - b * zd);

	// b = TMath::Exp(-fA->GetParameter(1)) / zd;
	// double b = TMath::Exp(-fA->GetParameter(1)) / zd;
	// double A = -TMath::Power(sigmamaxFit->GetParameter(1), 2) / TMath::Log(1 - b * zd);
	// A = -TMath::Power(energyCorrection*sigmamaxFit->GetParameter(1), 2) / TMath::Log(1 - b * zd);
	cout << "A: " << A << endl;
	cout << "b: " << b << endl;		

	// Plot comparison of diffusion parameters
	TCanvas *c3 = new TCanvas("c3", "c3", 900, 600);
	// TPad *pad1 = new TPad("p1", "", 0.0, 0.25, 1.0, 1.0);
	// TPad *pad2 = new TPad("p2", "", 0.0, 0.0, 1.0, 0.25);
	// pad1->Draw();
	// pad2->Draw();
	// pad1->cd();
	


	

	// muon data diffusion
	TF1 *fm = new TF1("fMuon", "TMath::Sqrt([0]*TMath::Log([2]-x*[1]))", 0, zd);		
	fm->SetParameter(0, -A);
	fm->SetParameter(1, b);
	fm->SetParameter(2, 1);
	fm->SetLineWidth(4);
	fm->SetLineColor(kBlue);
	fm->SetTitle("Diffusion Model");
	fm->GetXaxis()->SetTitle("Depth (#mum)");
	fm->GetYaxis()->SetTitle("#sigma_{xy} (#mum)");
	fm->Draw();

	// neutron data diffusion
	// TF1 *fn = new TF1("fNeutron", "TMath::Sqrt([0]*TMath::Log([2]-x*[1]))", 0, zd);		
	// fn->SetParameter(0, -ANeutron);
	// fn->SetParameter(1, bNeutron);
	// fn->SetParameter(2, 1);
	// fn->SetLineWidth(2);
	// fn->SetLineColor(kRed);
	// fn->Draw("same");

	// Make Legend
	TLegend *leg = new TLegend(0.35, 0.15, 0.88, 0.38);
	char fMuonLeg[200], fNeutronLeg[200];
	sprintf(fMuonLeg, "Muon Fit. A=%.2e #mum^{2}, b=%.2e #mum^{-1}", A, b);
	// sprintf(fNeutronLeg, "Iron Parameters. A=%.1e #mum^{2}, b=%.1e #mum^{-1}", ANeutron, bNeutron);
	leg->AddEntry("fMuon", fMuonLeg);
	// leg->AddEntry("fNeutron", fNeutronLeg);
	


	if(fixA){
		double Afix = ANeutron;
		double bfix = (1 - TMath::Exp(-TMath::Power(energyCorrection*computeSkewedMean(sigmamaxSkew->GetParameter(1), sigmamaxSkew->GetParameter(2), sigmamaxSkew->GetParameter(3)), 2) / Afix)) / zd;
		cout << "A: " << A << endl;
		cout << "b: " << b << endl;

		TF1 *fmFix = new TF1("fmFix", "TMath::Sqrt([0]*TMath::Log([2]-x*[1]))", 0, zd);
		fmFix->SetParameter(0, -Afix);
		fmFix->SetParameter(1, bfix);
		fmFix->SetParameter(2, 1);
		fmFix->SetLineWidth(4);
		fmFix->SetLineColor(kBlack);
		fmFix->SetLineStyle(2);
		fmFix->Draw("same");

		char fmFixLeg[200];
		sprintf(fmFixLeg, "Fixed A Value. A=%.1e #mum^{2}, b=%.1e #mum^{-1}", Afix, bfix);
		leg->AddEntry("fmFix", fmFixLeg);
	}



	leg->Draw();
	// pad2->cd();
	// TF1 *residue = new TF1("residue", "TMath::Sqrt([0]*TMath::Log([2]-x*[1])) - TMath::Sqrt([3]*TMath::Log([5]-x*[4])) ", 0, zd);
	// residue->SetTitle("");
	// residue->SetParameter(0, -ANeutron);
	// residue->SetParameter(1, bNeutron);
	// residue->SetParameter(2, 1);
	// residue->SetParameter(3, -A);
	// residue->SetParameter(4, b);
	// residue->SetParameter(5, 1);
	// residue->SetLineWidth(2);
	// residue->SetLineColor(kBlack);
	// residue->GetXaxis()->SetLabelSize(0);
	// // residue->GetYaxis()->SetTitle("#Delta#sigma_{xy} (#mum)");
	// residue->GetYaxis()->SetLabelSize(0.1);
	// residue->Draw();

	// Solve for what the depth of what the sigma max is supposed to be
	// ROOT::Math::WrappedTF1 wf1(fm);
	// ROOT::Math::BrentRootFinder brf;
	// brf.SetFunction(wf1, 0, zd)

	return;
}