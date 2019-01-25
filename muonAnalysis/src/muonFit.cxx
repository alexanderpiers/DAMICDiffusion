#include "muonFit.h"

TH2F* cluster = NULL;
Double_t sthr, sadc; //standard deviation of noise
TVector3 vfront(0,0,0); //guess position of front point
TVector3 vback(0,0,0); //guess position of end point
ofstream fout;

TH2F* ArraysToTH2F(TArrayD* x, TArrayD* y, TArrayD* val){
    
    Int_t min_x = 9999999;
    Int_t max_x = 0;
    Int_t min_y = 9999999;
    Int_t max_y = 0;
    
    for(Int_t i=0; i<x->GetSize(); i++){
        
        if(x->At(i)>max_x) max_x = (Int_t) x->At(i);
        if(x->At(i)<min_x) min_x = (Int_t) x->At(i);
        if(y->At(i)>max_y) max_y = (Int_t) y->At(i);
        if(y->At(i)<min_y) min_y = (Int_t) y->At(i);
    }
    
    TH2F* h2 = new TH2F("h2","h2",max_x-min_x+1,min_x-1,max_x,max_y-min_y+1,min_y-1,max_y);
    
    for(Int_t i=0; i<x->GetSize(); i++)
        h2->SetBinContent(((Int_t) x->At(i))-min_x+1, ((Int_t) y->At(i))-min_y+1, val->At(i));
    
    return h2;
}

Double_t maxk(Double_t obs, Double_t ne, Double_t c){
    
    //find k that offers maximum contribution
    Double_t kp = c*c*sadc*sadc*(TMath::Log(ne/2.) + TMath::EulerGamma() + obs/c/sadc/sadc);
    Int_t iters = 0;
    Int_t iiters = 0;
    while( TMath::Abs(kp/c/c/sadc/sadc + TMath::Log(kp) - TMath::Log(ne/2.) - TMath::EulerGamma() - obs/c/sadc/sadc) > 0.1 ){
        
        kp /= (kp/c/c/sadc/sadc + TMath::Log(kp))/(TMath::Log(ne/2.) + TMath::EulerGamma() + obs/c/sadc/sadc);
        iters++;
        if(iters>6){
            kp = c*c*sadc*sadc*(TMath::Log(ne/2.) + TMath::EulerGamma() + obs/c/sadc/sadc) - c*c*sadc*sadc*TMath::Log(kp);
            iters=0;
            iiters++;
        }
        //worst case
        if(iiters>2){
            kp=0;
            break;
        }
    }
    return kp;
}

Double_t sigmaxy(Double_t z, Double_t A, Double_t b, Double_t zo){
    
    Double_t s = 0;
    if(z>zo)
        s = TMath::Sqrt(-A*TMath::Log(1.-b*(z-zo)));
    
    return s;
}

Double_t fitf(Double_t z, Double_t* par){
    
    //par[0] - par[3] are startx, starty, endx, endy
    Double_t xs = par[0];
    Double_t ys = par[1];
    
    Double_t xt = par[2] - par[0];
    Double_t yt = par[3] - par[1];
    
    Double_t cfact = 0.25/(zd/TMath::Sqrt(zd*zd + yt*yt*ps*ps + xt*xt*ps*ps));
    
    //parameter is z
    Double_t x = (z/zd)*xt + xs; //in pixels
    Double_t y = (z/zd)*yt + ys; //in pixels
    
    //par[7] - par[8] are the bin numbers
    Double_t i = par[7];
    Double_t j = par[8];
    
    Double_t ne = 0;
    
    if(z>par[6]){
        
        //par[4] - par[6] are the diffusion model parameters
        //sigma comes in here
        Double_t bot = (TMath::Sqrt2()/ps)*sigmaxy(z,par[4],par[5],par[6]);
        ne = cfact*(TMath::Erf((i+0.5-x)/bot) - TMath::Erf((i-0.5-x)/bot))*(TMath::Erf((j+0.5-y)/bot) - TMath::Erf((j-0.5-y)/bot));
    }
    
    else if(x>=i-0.5 && x<i+0.5 && y>=j-0.5 && y<j+0.5){
        
        ne = cfact*4.;
    }
    
    return ne;
}

double TrackLogLikelihood(const double *xx){
    
    //we will be fitting for the start/end points and the diffusion parameters
    const Double_t fs = xx[0];
    const Double_t bs = xx[1];
    
    TVector3 vf = vfront + fs*(vback-vfront).Unit();
    TVector3 vb = vfront + bs*(vback-vfront).Unit();
    
    Double_t xs = vf.X()+0.5;
    Double_t ys = vf.Y()+0.5;
    Double_t xe = vb.X()+0.5;
    Double_t ye = vb.Y()+0.5;
    
    const Double_t smax = xx[2];
    const Double_t zo = xx[4];
    const Double_t b = (TMath::ATan(xx[3])/TMath::Pi()+0.5)/(zd-zo);
    Double_t A = -1. * smax * smax / TMath::Log(1.-b*(zd-zo));
    
    //also fit for the ADC to keV conversion factor
    const Double_t c = xx[5]*1000./3.6;
    const Double_t yield = xx[6];
    
    Int_t nx = cluster->GetNbinsX();
    Int_t ny = cluster->GetNbinsY();
    Double_t bxlow = cluster->GetXaxis()->GetBinLowEdge(1);
    Double_t bylow = cluster->GetYaxis()->GetBinLowEdge(1);
    
    //To hold information
    Double_t ll = 0;    
    Double_t par[9] = {0};
    
    //Now loop over cluster to obtain ll
    for(Int_t i=1; i<=nx; i++)
        for(Int_t j=1; j<=ny; j++){
            
            //value observed in this bin
            Double_t obs = cluster->GetBinContent(i,j);
            
            if(obs==0) continue;
            
            //update function with corresponding parameters
            par[0] = xs; par[1] = ys; par[2] = xe; par[3] = ye; par[4] = A; par[5] = b; par[6] = zo; par[7] = bxlow+i; par[8] = bylow+j;
            
            //to obtain maximum sigmaxy along z.
            Double_t sxy = sigmaxy(zd,A,b,zo) * TMath::Sqrt(zd*zd + (xe-xs)*(xe-xs)*ps*ps + (ye-ys)*(ye-ys)*ps*ps) / TMath::Sqrt((xe-xs)*(xe-xs)*ps*ps + (ye-ys)*(ye-ys)*ps*ps);
            
            //z coordinate corresponding to pixel
            Double_t zc = zd * ((xe-xs)*(i+bxlow-xs) + (ye-ys)*(j+bylow-ys)) / ((xe-xs)*(xe-xs) + (ye-ys)*(ye-ys));
            
            //set integration bounds to avoid integrating over entire range -> faster
            Double_t is = zc - 2.5*sxy;
            Double_t ie = zc + 2.5*sxy;
            
            if(is < 0.05){
                is = 0.05;
                if(is>zc) ie = is + 2.5*sxy;
            }
            
            if(ie > zd){
                
                ie = zd;
                if(zc>ie) is = zd - 2.5*sxy;
            }
            
            //round to nearest 0.1 and add 0.05
            Double_t istep = 0.1;
            
            is /= istep; is = (Int_t) is; is *= istep; is += istep/2.;
            ie /= istep; ie = (Int_t) ie; ie *= istep; ie += istep/2.;
            if(ie>zd) ie = zd;
            
            //number of electrons expected
            Double_t ne = 0;
            for(Double_t k=is; k<ie; k+=istep) ne += yield*fitf(k, par);
            ne *= istep;
            
            //incase the pixel is non-zero
            if(obs>0){
                                
                if(ne==0) {
                    
                    ll -= (obs)*(obs)/(2.*sadc*sadc) + TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
                }
                
                else{
                    
                    Double_t kp = maxk(obs,ne,c);
                                        
                    //round to integer and find maximum contribution to ll
                    Int_t ko = kp + 0.5;
                    if(ko<0) ko = 0;
                    Double_t llmax = ko*TMath::Log(ne) - ne - TMath::LnGamma(ko+1) - (obs-ko/c)*(obs-ko/c)/(2.*sadc*sadc) - TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
                    
                    //now look around
                    Int_t k = ko;
                    Double_t llk = llmax;
                    Double_t lli = llmax;
                    
                    //to the left
                    while(llmax-llk<4){
                        
                        k--;
                        if(k<0) break; //smallest value of k==0
                        llk = k*TMath::Log(ne) - ne - TMath::LnGamma(k+1) - (obs-k/c)*(obs-k/c)/(2.*sadc*sadc) - TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
                        lli += TMath::Log(1. + TMath::Exp(llk-lli));
                    }
                    
                    //to the right
                    k=ko;
                    llk = llmax;
                    while(llmax-llk<4){
                        
                        k++;
                        llk = k*TMath::Log(ne) - ne - TMath::LnGamma(k+1) - (obs-k/c)*(obs-k/c)/(2.*sadc*sadc) - TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
                        lli += TMath::Log(1. + TMath::Exp(llk-lli));
                    }
                    
                    //now increase the log likelihood
                    ll += lli;
                }
                
            }
            
            //in this case pixel is surrounding track
            else{
                
                //threshold for observation
                Double_t th = sthr*sadc;
                
                if(ne==0){
                    
                    ll += TMath::Log(1. + TMath::Erf(th/(TMath::Sqrt2()*sadc)));
                }

                else{
                    
                    //find k that offers maximum contribution
                    Double_t kp = maxk(th,ne,c);
                    
                    //round to integer and find maximum contribution to ll
                    Int_t ko = kp + 0.5;
                    if(ko<0) ko = 0;
                    Double_t y = (ko/c-th)/(TMath::Sqrt2()*sadc);
                    Double_t llmax = ko*TMath::Log(ne) - ne - TMath::LnGamma(ko+1);
                    if(y<4) llmax += TMath::Log(1. + TMath::Erf(-y));
                    else llmax -= y*y + TMath::Log(TMath::Sqrt(TMath::Pi())*y);
                    
                    //now look around
                    Int_t k = ko;
                    Double_t llk = llmax;
                    Double_t lli = llmax;
                    
                    //to the left
                    while(llmax-llk<4){
                        
                        k--;
                        if(k<0) break; //smallest value of k==0
                        y = (k/c-th)/(TMath::Sqrt2()*sadc);
                        llk = k*TMath::Log(ne) - ne - TMath::LnGamma(k+1);
                        if(y<4) llk += TMath::Log(1. + TMath::Erf(-y));
                        else llk -= y*y + TMath::Log(TMath::Sqrt(TMath::Pi())*y);
                        lli += TMath::Log(1. + TMath::Exp(llk-lli));
                    }
                    
                    //to the right
                    k=ko;
                    llk = llmax;
                    while(llmax-llk<4){
                        
                        k++;
                        y = (k/c-th)/(TMath::Sqrt2()*sadc);
                        llk = k*TMath::Log(ne) - ne - TMath::LnGamma(k+1);
                        if(y<4) llk += TMath::Log(1. + TMath::Erf(-y));
                        else llk -= y*y + TMath::Log(TMath::Sqrt(TMath::Pi())*y);
                        lli += TMath::Log(1. + TMath::Exp(llk-lli));
                    }
                    
                    //now increase the log likelihood
                    ll += lli;
                }
            }
        }
    
    return -1.*ll;
}

//pass the cluster and guesses for the start and end of the track
void muonFit(TH2F* cl, TTree *tMuon, Double_t *b, Double_t c, Double_t sadc_, Double_t sthr_){
    
    //initialize the objects
    sadc = sadc_;
    sthr = sthr_;
    
    //fit to a line
    Double_t slope = cl->GetCovariance()/cl->GetRMS(1)/cl->GetRMS(1);
    Double_t intercept = cl->GetMean(2) - slope*cl->GetMean(1);
    
    //guess front and back points
    
    //find extremes
    TVector3 exl(cl->GetXaxis()->GetBinLowEdge(1),slope*cl->GetXaxis()->GetBinLowEdge(1)+intercept,0);
    TVector3 exh(cl->GetXaxis()->GetBinUpEdge(cl->GetNbinsX()),slope*cl->GetXaxis()->GetBinUpEdge(cl->GetNbinsX())+intercept,0);
    
    TH1F* qalong = new TH1F("qalong","qalong",(Int_t) ((exh-exl).Mag()+0.5),0,(Int_t) ((exh-exl).Mag()+0.5));
    TH1F* nalong = new TH1F("nalong","nalong",(Int_t) ((exh-exl).Mag()+0.5),0,(Int_t) ((exh-exl).Mag()+0.5));
    
    for(Int_t i=1; i<=cl->GetNbinsX(); i++){
        for(Int_t j=1; j<=cl->GetNbinsY(); j++){
            
            if(cl->GetBinContent(i,j)<=0) continue;
            
            TVector3 bcv(cl->GetXaxis()->GetBinCenter(i),cl->GetYaxis()->GetBinCenter(j),0);
            Int_t bin = (bcv-exl).Dot((exh-exl).Unit())+1;
            qalong->SetBinContent(bin,qalong->GetBinContent(bin)+cl->GetBinContent(i,j));
            nalong->SetBinContent(bin,nalong->GetBinContent(bin)+1.);
        }
    }
    
    Double_t qalong_mean = qalong->Integral() / qalong->GetNbinsX();
    qalong->Smooth(); //not the best way to deal with binning effects
    nalong->Smooth();
    
    Double_t xhstart = 0;
    Double_t xhend = 0;
    Double_t xqstart = 0;
    Double_t xqend = 0;
    
    for(Int_t i=1; i<=qalong->GetNbinsX(); i++){
        
        if(qalong->GetBinContent(i)>qalong_mean/2.){
            
            xhend = qalong->GetXaxis()->GetBinUpEdge(i);
            if(xhstart==0) xhstart = qalong->GetXaxis()->GetBinLowEdge(i);
        }
        
        if(qalong->GetBinContent(i)>0){
            
            xqend = qalong->GetXaxis()->GetBinUpEdge(i);
            if(xqstart==0) xqstart = qalong->GetXaxis()->GetBinLowEdge(i);
        }
    }
    
    Double_t fmin, fmax;
    Double_t bmin, bmax;
    
    Double_t npixstart = nalong->Integral(xhstart+4,xhstart+8);
    Double_t npixend = nalong->Integral(xhend-7,xhend-3);
    
    if(npixstart < npixend){
        
        vfront = exl + xhstart*(exh-exl).Unit();
        vback = exl + xhend*(exh-exl).Unit();
        fmin = xqstart - xhstart - 1; fmax = 2.;
        bmax = xqend - xhstart + 2.; bmin = bmax - 10.;
    }
    
    else{
        
        vback = exl + xhstart*(exh-exl).Unit();
        vfront = exl + xhend*(exh-exl).Unit();
        fmin = xhend - xqend - 1; fmax = 2.;
        bmax = xhend - xqstart + 2.; bmin = bmax - 10.;
    }
    
    //want a histogram that is 1 pixel larger around
    cluster = new TH2F("cl2fit","cl2fit",cl->GetNbinsX()+2,cl->GetXaxis()->GetBinLowEdge(1)-1,cl->GetXaxis()->GetBinUpEdge(cl->GetNbinsX())+1,cl->GetNbinsY()+2,cl->GetYaxis()->GetBinLowEdge(1)-1,cl->GetYaxis()->GetBinUpEdge(cl->GetNbinsY())+1);
    
    Int_t nblanked = 0;
    
    //fill the histogram, including marking pixels "around"
    for(Int_t i=1; i<=cluster->GetNbinsX(); i++)
        for(Int_t j=1; j<=cluster->GetNbinsY(); j++){
            
            Double_t cont = 0;
            if(i>1 && j>1 && i<cluster->GetNbinsX() && j<cluster->GetNbinsY()) cont = cl->GetBinContent(i-1,j-1);
            
            if(cont!=0){
                
                //these are the pixels with collected charge
                TVector3 bcv(cluster->GetXaxis()->GetBinCenter(i),cluster->GetYaxis()->GetBinCenter(j),0);
                Int_t bin = (bcv-exl).Dot((exh-exl).Unit())+1;
                Double_t lyield = c*qalong->GetBinContent(bin)*1000./3.6*(xhend-xhstart)/TMath::Sqrt(ps*(xhend-xhstart)*ps*(xhend-xhstart) + zd*zd);
                
                if(lyield > cyield + 3.*TMath::Sqrt(cyield)) { cont = 0; nblanked++; }
                
                cluster->SetBinContent(i,j,cont);
            }
            
            else{
                
                Int_t around = 0;
                
                for(Int_t ki=-1; ki<=1; ki++)
                    for(Int_t kj=-1; kj<=1; kj++){
                        
                        Int_t bx = i-1+ki;
                        Int_t by = j-1+kj;
                        
                        if((ki==0 && kj==0) || bx<1 || by<1 || bx>cl->GetNbinsX() || by>cl->GetNbinsY()) continue;
                        if(cl->GetBinContent(bx,by)!=0) around++;
                    }
                
                //mark pixels around track to include the clustering end condition in the likelihood
                if(around>0) cluster->SetBinContent(i,j,-10000);
            }
        }
    
    cl->Draw("COLZ");
    TF1* l = new TF1("l","pol1(0)",0,5000);
    l->FixParameter(0,intercept);
    l->FixParameter(1,slope);
    // vfront.Print(); vback.Print();
    // std::cout << "Cut: " << (cyield + 2.5*TMath::Sqrt(cyield))*3.6/1000./c*TMath::Sqrt(ps*(xhend-xhstart)*ps*(xhend-xhstart) + zd*zd)/(xhend-xhstart) << std::endl;
    new TCanvas();
    cluster->Draw("COLZ");
    l->Draw("same");
    new TCanvas();
    qalong->Draw();
    new TCanvas();
    nalong->Draw();
    
    // std::cout << fmin << " " << fmax << std::endl;
    // std::cout << bmin << " " << bmax << std::endl;
    
//    double xx[7] = {6.77744e-02,4.30587e+01,7.79478e+00,5.55929e+00,0,c,7.39027e+01};
//
//    TH2F* h2 = new TH2F("h2","h2",100,10,11,100,5.65,5.75);
//    
//    for(Int_t i=1; i<=100; i++)
//        for(Int_t j=1; j<=100; j++){
//            
//            xx[4] = i*0.01 + 10 + 0.01/2.;
//            xx[5] = j*0.001 + 5.65 + 0.001/2.;
//            
//            h2->SetBinContent(i,j,TrackLogLikelihood(xx));
//        }
//    
//    h2->Draw("COLZ");
    
//    std::cout << TrackLogLikelihood(xx) << std::endl;

//    Double_t x[1000];
//    Double_t y[1000];
//    Int_t n = 0;
//    
//    for(Double_t k=40; k<45; k+=0.1){
//    
//        xx[1] = k;
//        //std::cout << k << " " << TrackLogLikelihood(xx) << std::endl;
//        x[n] = k;
//        y[n] = TrackLogLikelihood(xx);
//        n++;
//    }
//    
//    TGraph* g = new TGraph(n,x,y);
//    g->Draw("A*L");
    
    //Here is where the fit is done
    ROOT::Math::Minimizer* migrad = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
    
    // set tolerance , etc...
    migrad->SetMaxFunctionCalls(1000000);
    migrad->SetTolerance(0.0001);
    migrad->SetPrintLevel(0);
    
    // create funciton wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&TrackLogLikelihood,7);
    migrad->SetErrorDef(0.5);
    migrad->SetFunction(f);
    
    // Set the free variables to be minimized!
    migrad->SetLimitedVariable(0,"fs",0.5, 0.0001, fmin, fmax);
    migrad->SetLimitedVariable(1,"bs",(vback-vfront).Mag(), 0.001, bmin, bmax);
    migrad->SetLimitedVariable(2,"smax",7, 0.01, 1, 30);
    migrad->SetLimitedVariable(3,"a", -0.5, 0.002, -5, 5);
    //migrad->SetLimitedVariable(4,"zo", 5, 0.001, 0.001, 25);
    //migrad->SetLimitedVariable(5,"c",c,1E-6,2E-6,2E-2);
    migrad->SetLimitedVariable(6,"yield",cyield,0.1,0,500);
    
    //migrad->SetFixedVariable(0,"fs",0);
    migrad->SetFixedVariable(5,"c",c);
    migrad->SetFixedVariable(4,"zo",0);
    
    Bool_t mstat = migrad->Minimize();
    const double *mval = migrad->X();
    const double *merr = migrad->Errors();
//    
    fout << nblanked << " ";
    fout << (vfront + mval[0]*(vback-vfront).Unit()).X() << " " << (vfront + mval[0]*(vback-vfront).Unit()).Y() << " " << (vfront + mval[1]*(vback-vfront).Unit()).X() << " " << (vfront + mval[1]*(vback-vfront).Unit()).Y() << " ";
    fout << mval[2] << " " << mval[3] << " " << mval[4] << " " << mval[6] << " " << migrad->MinValue() << " " << mstat << std::endl;

    if(tMuon != nullptr){
        b[0] = mval[2];
        b[1] = mval[3];
        b[2] = merr[2];
        b[3] = merr[3];
        tMuon->Fill();
    }

    delete migrad;
    nalong->Delete();
    qalong->Delete();
    cluster->Delete();

    
    return;
}

void muonFitMany(TChain* t, TTree *tMuon, Double_t *b, Double_t c, Double_t sadc_, Double_t sthr_){
    
    //output text file
    fout.open("~/muon_fits.txt");
    
    //variables to load
    Int_t RUNID, cluster_id, EXTID;
    
    TParameter<double>* npix = NULL;
    TParameter<double>* clength = NULL;
    TParameter<double>* slength = NULL;
    TParameter<double>* corrfact = NULL;
    TParameter<double>* qtotal = NULL;
    TParameter<double>* qmax = NULL;
    TParameter<double>* mask_edge = NULL;
    TParameter<double>* track_rms = NULL;
    TParameter<double>* raw_nsat = NULL;
    
    //arrays to store object
    TArrayD* pixel_x = NULL;
    TArrayD* pixel_y = NULL;
    TArrayD* pixel_val = NULL;
    
    t->SetBranchAddress("RUNID", &RUNID);
    t->SetBranchAddress("EXTID", &EXTID);
    t->SetBranchAddress("cluster_id", &cluster_id);
    t->SetBranchAddress("size_npixels",&npix);
    t->SetBranchAddress("size_clength", &clength);
    t->SetBranchAddress("size_slength", &slength);
    t->SetBranchAddress("curve_correlation_factor", &corrfact);
    t->SetBranchAddress("charge_maximum", &qmax);
    t->SetBranchAddress("charge_total", &qtotal);
    t->SetBranchAddress("mask_edge", &mask_edge);
    t->SetBranchAddress("curve_track_rms", &track_rms);
    t->SetBranchAddress("raw_nsat", &raw_nsat);
    
    t->SetBranchAddress("pixel_x", &pixel_x);
    t->SetBranchAddress("pixel_y", &pixel_y);
    t->SetBranchAddress("pixel_val", &pixel_val);
    
    Int_t nentries = t->GetEntries();
    
    fout << "runid:extid:clid:npix:slength:nblanked:front_x:front_y:back_x:back_y:smax:a:zo:yield:ll:mstat" << std::endl;
    
    for(Int_t i=0; i<nentries; i++){
        
        Double_t prev_rid = RUNID;
        t->GetEntry(i);
        
        if(prev_rid!=RUNID)
            std::cout << "Processing run " << RUNID << std::endl;
        
        //if(npix->GetVal() > 50 && mask_edge->GetVal() == 0 && raw_nsat->GetVal() == 0 && clength->GetVal()/slength->GetVal() > 0.95 && clength->GetVal()/slength->GetVal() < 1.1 && TMath::Abs(corrfact->GetVal()) > 0.99 && qmax->GetVal()*npix->GetVal()/qtotal->GetVal() < 6. && track_rms->GetVal() < 2.){
        
    
         
        fout << RUNID << " " << EXTID << " " << cluster_id << " " << npix->GetVal() << " " << slength->GetVal() << " ";
        
        TH2F* cl = ArraysToTH2F(pixel_x,pixel_y,pixel_val);
        muonFit(cl, tMuon, b, c, sadc_, sthr_);
        if(abs(b[1]) > 9) cout << i << endl;
        cl->Delete();
        
    }
    
    fout.close();
    return;

}

void saveMuonFitToFile(const char* outfile, const char *infile, Double_t c, Double_t sadc_, Double_t sthr_){

    TFile *fMuonTracks = new TFile(infile);
    TChain *tMuonTracks = (TChain*)fMuonTracks->Get("clusters_tree");
    TFile *fMuonFit = new TFile(outfile, "RECREATE");

    double branches[4];

    // Create a tree to store data
    TTree *tMuonFit = new TTree("muonFit", "muonFit");
    tMuonFit->Branch("sigmamax", (branches));
    tMuonFit->Branch("a", (branches + 1));
    tMuonFit->Branch("sigmamaxErr", (branches + 2));
    tMuonFit->Branch("aErr", (branches + 3));

    // TH1D *histogramArray[4];
    // histogramArray[0] = new TH1D("sigmamax", "sigmamax", 50, 5, 15);
    // histogramArray[1] = new TH1D("a", "a", 100, -15, 15);
    // histogramArray[2] = new TH1D("sigmamaxErr", "sigmaMaxErr", 50, 0, 1);
    // histogramArray[3] = new TH1D("aErr", "aErr", 50, 0, 10);

    muonFitMany(tMuonTracks, tMuonFit, branches, c, sadc_, sthr_);

    fMuonFit->Write();
    delete fMuonFit;

    return;
}
