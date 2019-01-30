void sigmavEnergy(){


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	TChain* j = new TChain("clusters_tree");
	j->Add("short2018.root");


	j->SetAlias("enell","LL_eneh.fVal*2.6e-4*(  1.039*(EXTID==1) + 0.962*(EXTID==2) + 0.999*(EXTID==3) + 0.986*(EXTID==4) + 1.034*(EXTID==6) + 0.982*(EXTID==11) + 0.993*(EXTID==12) )");
    j->SetAlias("sigma","LL_sigma.fVal");

    // Plot sigma vs energy. 
    TCanvas *csigma = new TCanvas("csigma", "csigma", 800, 600);
    TH2D *sigve = new TH2D("sigve", "#sigma_{max} v Energy", 150, 0, 15, 30, 0, 1.5);
    sigve->GetXaxis()->SetTitle("Energy (keV)");
    sigve->GetYaxis()->SetTitle("#sigma_{max} (pixels)");
    j->SetMarkerStyle(20);
    j->SetMarkerSize(1); 
    sigve->SetMarkerStyle(20);
    sigve->SetMarkerSize(1);
    j->Draw("sigma:enell >> sigve", "sigma > (0.91 + 0.005*enell) && sigma < (1.0 + 0.005*enell)");

    TF1 *fit = new TF1("sve", "[0]+[1]*x", 0.5, 15);
    fit->SetParNames("#sigma_{max}(0)", "#alpha");
    fit->SetParameters(1, 0.01);
    fit->SetLineWidth(3);

    sigve->Fit("sve", "QR");


	return;
}