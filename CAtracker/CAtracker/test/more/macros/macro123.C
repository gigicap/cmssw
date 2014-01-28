void macro123(){
	
gROOT->SetStyle("Plain");
gStyle->SetOptStat("");	
	
TFile *f1 = new TFile("histotrip_USE_TT.root","read");

TTree *t1 = f1->Get("CAtracker/evtree");

TH1I *h1 = new TH1I("h1","TTbar",60,0,120);
TH1I *h2 = new TH1I("h2","TTbar",60,0,120);

int nentries1 = t1->GetEntries();
    
std::cout<<"->nentries = "<<nentries1<<std::endl;
    
int prodcells123 , nseeds;

t1->SetBranchAddress("prodCells123", &prodcells123);
t1->SetBranchAddress("nSeeds", &nseeds);

for (int k = 0; k<nentries1; k++){
        t1->GetEntry(k);       
	h2->Fill(prodcells123);
	h1->Fill(nseeds);
}

	h1->SetLineColor(2);
	h2->SetLineColor(4);
	
TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetFillColor(0);

leg->AddEntry(h1,"Number of pent-uplet seeds","l");
leg->AddEntry(h2,"Number of triplet seeds","l");

new TCanvas;
	h1->Draw();
	h2->Draw("same");
	leg->Draw("same");
	
}

	
