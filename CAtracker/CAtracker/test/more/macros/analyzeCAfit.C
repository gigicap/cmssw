void analyzeCAfit(TString filename, double maxval){

    double localmomentum;
    int identifier;
        
    
    gStyle->SetOptStat("");
    
    cout<<"Opening files"<<endl;
    
    TFile *f1 = new TFile(filename,"READ");
    
    cout<<"Getting tree"<<endl;
    
    TTree *t1 = f1->Get("CAtracker/triptree");
    TH1D *h123 = new TH1D("h123","Local Momentum magnitude",100,0.,maxval);
    TH1D *h234 = new TH1D("h234","Local Momentum magnitude",100,0.,maxval);
    TH1D *h345 = new TH1D("h345","Local Momentum magnitude",100,0.,maxval);
    
    h123->SetLineColor(kRed);
    h234->SetLineColor(kGreen);
    h345->SetLineColor(kBlue);
    
    TLegend *l1 = new TLegend(0.4,0.5,0.7,0.7);
    l1->AddEntry(h123,"BPix1+BPix2+BPix3","l");
    l1->AddEntry(h234,"BPix2+BPix3+TIB1","l");
    l1->AddEntry(h345,"BPix3+TIB1+TIB2","l");
    
    
    

    t1->SetBranchAddress("LocalMomentum", &localmomentum);
    t1->SetBranchAddress("Identifier", &identifier);
    
    Int_t nentries = t1->GetEntries();
    
    for(int i = 0; i<nentries; i++){
        t1->GetEntry(i);
        
        if(identifier == 123) h123->Fill(localmomentum);
        if(identifier == 234) h234->Fill(localmomentum);
        if(identifier == 345) h345->Fill(localmomentum);

    }
    
    new TCanvas;
    
    h123->Draw();
    h234->Draw("same");
    h345->Draw("same");
    l1->Draw("same");
    
    
    cout<<"setting"<<endl;
    
    return;
}
