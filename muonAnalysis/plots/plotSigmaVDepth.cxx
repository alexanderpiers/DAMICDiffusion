
// Plots sigma v Depth including the best fit lines
TCanvas * plotSigmaVDepth(int emin=3, int emax=4){

	const double CCDWidth = 500; // CCD width in um

	// Define constants
	const double q = 1.6e-19; // C
	const double T = 130.; // K 
	const double kb = 1.38e-23; // m^2 kg / K s^2 
	const double epsSi = 11.68 * 8.85e-12; // C / Vm 
	const double Vapp = 128.; // V
	const double rho = 3.9e17*q; // m^-3 
	const double um =1.e6;
	const double pixelsize=15; // um / pixel

	// Defining parameters
	const double neutronDataA = -215; // um^2
	const double neutronDataB = 5.6e-4;// 1 / um

	// Read files that contain the charge distribution histograms.
	TCanvas *c = new TCanvas("c", "c", 1000, 800);
	TFile *fsigma = new TFile("~/mnt/mox/DAMICDiffusion/rootfiles/sigma_depth_dedx8kev.root");

	// Read the TGraph from the 
	char *graphname = new char[100];
	char *graphtitle = new char[200];

	// sprintf(graphname, "sigmavdepth");
	sprintf(graphname, "sigma_depth_%i_%ikeV", emin, emax);
	sprintf(graphtitle, "#sigma_{x} vs. Depth for %i-%i keV events", emin, emax);

	// Get the graph and fit
	TGraph *tg = (TGraph*) fsigma->Get(graphname);
	TF1 *fMuon = tg->GetFunction("fMuon");
	

	// Scale everything back into um
	for(int i=0; i<tg->GetN(); i++){
		tg->GetX()[i] *= um;
		tg->GetY()[i] *= um;
		tg->GetEY()[i] *= um;
	}
	fMuon->SetParameter(0, fMuon->GetParameter(0)*um*um);
	fMuon->SetParameter(1, fMuon->GetParameter(1)/um);

	// Format the plots
	tg->GetXaxis()->SetLimits(0, CCDWidth);
	tg->GetYaxis()->SetRangeUser(0, tg->GetYaxis()->GetXmax()*um);
	tg->SetTitle(graphtitle);
	tg->GetXaxis()->SetTitle("Depth (#mum)");
	tg->GetYaxis()->SetTitle("#sigma_{x} (#mum)");
	tg->SetLineWidth(2);
	tg->SetLineColor(kRed + 1);
	tg->SetMarkerColor(kRed + 1);
	tg->SetMarkerStyle(21);
	tg->SetMarkerSize(0.4);
	fMuon->SetLineColor(kRed + 1);
	fMuon->SetLineColor(2);
	fMuon->SetRange(0, CCDWidth);

	// Add current neutron diffusion parameters
	TF1 *fNeutron = new TF1("fNeutron", "TMath::Sqrt( [0]*TMath::Log([2]-x*[1]))", 0, CCDWidth);
	fNeutron->SetParameter(0, neutronDataA);
	fNeutron->SetParameter(1, neutronDataB);
	fNeutron->SetParameter(2, 1);
	fNeutron->SetLineColor(kAzure - 5);
	fNeutron->SetLineWidth(2);

	tg->Draw("AP");
	fMuon->Draw("same");
	fNeutron->Draw("same");

	// Create legend
	TLegend *leg = new TLegend(0.42, 0.15, 0.88, 0.35);
	char fMuonLeg[200], fNeutronLeg[200];
	sprintf(fMuonLeg, "Muon Fit. A=%.1e #mum^{2}, b=%.1e #mum^{-1}", fMuon->GetParameter(0), fMuon->GetParameter(1));
	sprintf(fNeutronLeg, "Neutron Fit. A=%.1e #mum^{2}, b=%.1e #mum^{-1}", fNeutron->GetParameter(0), fNeutron->GetParameter(1));
	leg->AddEntry("tg", "Muon Diffusion Data", "lep");
	leg->AddEntry("fMuon", fMuonLeg);
	leg->AddEntry("fNeutron", fNeutronLeg);
	leg->Draw();

	return c;
}

