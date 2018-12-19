#include <string>

TCanvas * compareChargeDist(int dlow=100, int dhigh=110, int elow=5, int ehigh=6){
	
	// Read files that contain the charge distribution histograms.
	TFile *fDelta = new TFile("~/DAMICDiffusion/rootfile/charge_dist_deltarays_56kev.root");
	TFile *fNoDelta = new TFile("~/DAMICDiffusion/rootfile/charge_dist_nodeltarays_56kev.root");

	TFile *fDedx = new TFile("~/DAMICDiffusion/rootfile/charge_dist_dedx_6kev_56kev.root");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	// Create TCanvas and add histograms to it.
	TCanvas *c = new TCanvas("c","c", 800, 600);
	c->SetLogy();
	//c->Divide(3,1);
	double frange = 1.5;
	// Histogram name
	char histname[100];
	sprintf(histname, "depth_%i_%ium_%i_%ikev", dlow, dhigh, elow, ehigh);
	//const char * charHistname = histname.c_str();

	// Plot histograms
	TH1D * hDelta = (TH1D*)fDelta->Get(histname);
	hDelta->SetName("hDelta");
	hDelta->Scale(1./hDelta->Integral());
	hDelta->Fit("gaus", "IQR0","", -frange, frange);
	TF1 *hDeltaFit = hDelta->GetFunction("gaus");
	hDelta->GetYaxis()->SetTitle("");
	hDelta->SetLineWidth(2); hDelta->SetLineColor(38);
	hDeltaFit->SetLineWidth(3); hDeltaFit->SetLineStyle(7); hDeltaFit->SetLineColor(38);	
	hDelta->Draw("hist");
	hDeltaFit->Draw("same");
		
	// Delta ray distance filtering
	TH1D * hNoDelta = (TH1D*)fNoDelta->Get(histname);
	hNoDelta->SetName("hNoDelta");
	hNoDelta->Scale(1./hNoDelta->Integral());
	hNoDelta->Fit("gaus", "IQR0","", -frange, frange);
	TF1 *hNoDeltaFit = hNoDelta->GetFunction("gaus");
	hNoDelta->SetLineWidth(2); hNoDelta->SetLineColor(41);
	hNoDeltaFit->SetLineWidth(3); hNoDeltaFit->SetLineStyle(7); hNoDeltaFit->SetLineColor(41);	
	hNoDelta->Draw("same");
	hNoDeltaFit->Draw("same");
	
	// dEdx filtered histogram	
	TH1D * hDedx = (TH1D*)fDedx->Get(histname);
	hDedx->SetName("hDedx");
	hDedx->Scale(1./hDedx->Integral());
	hDedx->Fit("gaus", "IQR0","",  -frange, frange);
	TF1 *hDedxFit = hDedx->GetFunction("gaus");
	hDedx->SetLineWidth(2); hDedx->SetLineColor(46);
	hDedxFit->SetLineWidth(3); hDedxFit->SetLineStyle(7); hDedxFit->SetLineColor(46);	
	hDedx->Draw("same");
	hDedxFit->Draw("same");

	// Create Legend
	TLegend *leg = new TLegend(0.12, 0.7, 0.38, 0.88);
	char hDeltaLeg[200], hNoDeltaLeg[200], hDedxLeg[200];
	sprintf(hDeltaLeg, "No #delta ray Filter. #sigma=%.2f", hDeltaFit->GetParameter(2));
	sprintf(hNoDeltaLeg, "Distance #delta ray Filter. #sigma=%.2f", hNoDeltaFit->GetParameter(2));
	sprintf(hDedxLeg, "#frac{dE}{dx} Filter. #sigma=%.2f", hDedxFit->GetParameter(2));
	//leg->AddEntry("hDelta", "Test");
	//leg->AddEntry("hNoDelta", "Test");
	//leg->AddEntry("hDedx", "Test");
	leg->AddEntry("hDelta", hDeltaLeg);
	leg->AddEntry("hNoDelta", hNoDeltaLeg);
	leg->AddEntry("hDedx", hDedxLeg);
	leg->Draw();
		
}
