{

gROOT->ProcessLine(".L ~/recon/Viewer.C+");
gROOT->ProcessLine(".L ~/recon/ADPmacros/muonFilter.c");

TChain* c = new TChain("clusters_tree");
c->Add("/gscratch/damic/data/uchicago/processed/D3500/SbBe_2015-04-10_FullBe/root/Image176*.root");

c->SetAlias("ene","charge_total.fVal/10300*6.4");
c->SetAlias("x","charge_mean_x.fVal");
c->SetAlias("y","charge_mean_y.fVal");

}
