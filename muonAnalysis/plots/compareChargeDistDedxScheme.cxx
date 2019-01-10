#include <string>

// Plotting function for comparing the difference between the histogram charge distribution when filtering out the dedx in two methods:
// 	--Total track dedx - so summing all the energies of a segment and dividing by the segment length
//  --Individual pixel dedx - using the pixel dedx coming from the filtering by dedx length
TCanvas * compareChargeDistDedxScheme(int dlow=200, int dhigh=225, int elow=5, int ehigh=6){
	
	// Read files that contain the charge distribution histograms.
	TFile *ftrack = new TFile("~/mnt/mox/DAMICDiffusion/rootfiles/charge_dist_dedxfilt10kev.root");
	TFile *fpixel = new TFile("~/mnt/mox/DAMICDiffusion/rootfiles/charge_dist_dedxfilt10kev_newEFilter.root");
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	// Create TCanvas and add histograms to it.
	TCanvas *c = new TCanvas("c","c", 800, 600);
	c->SetLogy();
	double frange = 2.5;

	// Histogram name
	char histname[100];
	sprintf(histname, "depth_%i_%ium_%i_%ikev", dlow, dhigh, elow, ehigh);

	// Plot Histograms
	// Total Track Dedx
	TH1D * htrack = (TH1D*)ftrack->Get(histname);
	htrack->SetName("htrack");
	htrack->Scale(1./htrack->Integral());
	TF1 *htrackFit = new TF1("htrackFit", "gaus", -frange, frange);
	htrack->Fit("htrackFit", "IQR0");	
	htrack->GetYaxis()->SetTitle("");
	htrack->SetLineWidth(3); htrack->SetLineColor(38);
	htrackFit->SetLineWidth(3); htrackFit->SetLineStyle(7); htrackFit->SetLineColor(38);	
	htrack->Draw("hist");
	htrackFit->Draw("same");

	// Pixel Dedx
	TH1D * hpixel = (TH1D*)fpixel->Get(histname);
	hpixel->SetName("hpixel");
	hpixel->Scale(1./hpixel->Integral());
	TF1 *hpixelFit = new TF1("hpixelFit", "gaus", -frange, frange);
	hpixel->Fit("hpixelFit", "IQR0");	
	hpixel->GetYaxis()->SetTitle("");
	hpixel->SetLineWidth(3); hpixel->SetLineColor(41);
	hpixelFit->SetLineWidth(3); hpixelFit->SetLineStyle(7); hpixelFit->SetLineColor(41);	
	hpixel->Draw("histsame");
	hpixelFit->Draw("same");

	// Create Legend
	TLegend *leg = new TLegend(0.12, 0.7, 0.32, 0.88);
	char htrackLeg[200], hpixelLeg[200];
	sprintf(htrackLeg, "Segment #frac{dE}{dx}. #sigma=%.2f", htrackFit->GetParameter(2));
	sprintf(hpixelLeg, "Pixel #frac{dE}{dx}. #sigma=%.2f", hpixelFit->GetParameter(2));
	leg->AddEntry("htrack", htrackLeg);
	leg->AddEntry("hpixel", hpixelLeg);
	leg->Draw();

	return c; 
}