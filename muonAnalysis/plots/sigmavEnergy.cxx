void sigmavEnergy(){


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	double pixelsize = 15;

	TChain* j = new TChain("clusters_tree");
	j->Add("short2018.root");

	// Defining sigma max points from the muon data
	const int n = 2;
	double energykev[n] = {0.1, 4.15};
	double sigmamax[n] = {0.951, 0.975};

	j->SetAlias("enell","LL_eneh.fVal*2.6e-4*(  1.039*(EXTID==1) + 0.962*(EXTID==2) + 0.999*(EXTID==3) + 0.986*(EXTID==4) + 1.034*(EXTID==6) + 0.982*(EXTID==11) + 0.993*(EXTID==12) )");
    j->SetAlias("sigma","LL_sigma.fVal");

    // Plot sigma vs energy. 
    TCanvas *csigma = new TCanvas("csigma", "csigma", 800, 600);
    TH2D *sigve = new TH2D("sigve", "#sigma v Energy", 150, 0, 15, 30, 0, 1.5);
    sigve->GetXaxis()->SetTitle("Energy (keV)");
    sigve->GetYaxis()->SetTitle("#sigma (pixels)");
    j->SetMarkerStyle(20);
    j->SetMarkerSize(1); 
    sigve->SetMarkerStyle(20);
    sigve->SetMarkerSize(1);
    j->Draw("sigma:enell >> sigve", "sigma > (0.91 + 0.005*enell) && sigma < (1.0 + 0.005*enell)");
    


    // Fit the sigma maxes to a line to get the slope
    TF1 *fit = new TF1("sve", "[0]+[1]*x", 0, 15);
    fit->SetParNames("#sigma_{max}(0)", "#alpha");
    fit->SetParameters(1, 0.01);
    fit->SetLineWidth(3);
    sigve->Fit("sve", "R0", "", 0.5, 15);


    fit->SetParameter(0, sigmamax[1]-fit->GetParameter(1)*energykev[1]);
    j->Draw("sigma:enell >> sigve", "sigma < 1.2", "same");
    fit->Draw("same");

    TGraph *muonSigmaPoints = new TGraph(2, energykev, sigmamax);
    muonSigmaPoints->SetMarkerStyle(29);
    muonSigmaPoints->SetMarkerColor(kViolet);
    muonSigmaPoints->SetMarkerSize(3);
    muonSigmaPoints->Draw("sameP");
    muonSigmaPoints->SetName("msp");

    // Adding shading for between alpha [0.004, 0.005]
    // double alphaMin = 0.004;
    // double alphaMax = 0.005;
    // const int ne = 16;
    // double energy[ne];
    // double sigmaAlphaMin[ne];
    // double sigmaAlphaMax[ne];

    // // Fill array with data
    // for(int i=0; i<ne; i++){
    // 	energy[i] = i;
    // 	sigmaAlphaMin[i] = sigmamax[0] + alphaMin*i;
    // 	sigmaAlphaMax[i] = sigmamax[0] + alphaMax*i;
    // }

    // // Create the tgraphs
    // TGraph *tgAlphaMax = new TGraph(ne, energy, sigmaAlphaMax);
    // TGraph *tgAlphaMin = new TGraph(ne, energy, sigmaAlphaMin);
    // TGraph *tgAlphaFill = new TGraph(2*ne);
    // tgAlphaFill->SetName("alphafill");

    // for(int i=0; i<ne; i++){	
    // 	tgAlphaFill->SetPoint(i, energy[i], sigmaAlphaMax[i]);
    // 	tgAlphaFill->SetPoint(ne+i, energy[ne-i-1], sigmaAlphaMin[ne-i-1]);
    	
    // }
    // tgAlphaFill->SetFillStyle(1001);
    // tgAlphaFill->SetFillColorAlpha(kGreen+2, 0.6);
    // tgAlphaFill->Draw("samef");
    // tgAlphaMin->SetLineColor(kGreen+2);
    // tgAlphaMin->SetLineWidth(2);
    // tgAlphaMax->SetLineColor(kGreen+2);
    // tgAlphaMax->SetLineWidth(2);
    // tgAlphaMin->Draw("samel");
    // tgAlphaMax->Draw("samel");


    // Create legend
    TLegend *leg = new TLegend(0.12, 0.88, 0.42, 0.7);
    leg->AddEntry("sigve", "SNOLAB Background Data");
    leg->AddEntry("sve", "#sigma_{max} Energy Fit");
    leg->AddEntry("msp", "#sigma_{max} from Muons");
    // char alphaleg[100]; sprintf(alphaleg, "#alpha = [%.3f, %.3f]", alphaMin, alphaMax);
    // leg->AddEntry("alphafill", alphaleg);
    leg->Draw();


    // Plotting projection onto the y axis
    gStyle->SetOptFit(1);
    double dE = 2;
    double eRangeMin = energykev[1] - dE/2; 
    double eRangeMax = energykev[1] + dE/2;

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    TH1D *sigmaxy = sigve->ProjectionY("sigmaxy projection", sigve->GetXaxis()->FindBin(eRangeMin), sigve->GetXaxis()->FindBin(eRangeMax));
    sigmaxy->SetTitle("#sigma_{xy} of events around Muon energy scaling");
    sigmaxy->Draw("hist");
    cout << sigmaxy->GetBinCenter(sigmaxy->GetMaximumBin()) << endl;
    TF1 *sigmaxyFit = new TF1("sigmaxyFit", "gaus", 0.8, 1.1);
    sigmaxy->Fit("sigmaxyFit", "RI");
    sigmaxyFit->Draw("same");

	return;
}