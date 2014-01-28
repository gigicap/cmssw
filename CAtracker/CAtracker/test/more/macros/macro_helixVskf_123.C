void macro_helixVskf_123(){
	
gROOT->SetStyle("Plain");
gStyle->SetOptStat("");	
	
TFile *f_kf = new TFile("histotrip_USE_Mu.root","read");
TFile *f_helix = new TFile("histotrip.root","read");

TTree *t_kf_trip = f_kf->Get("CAtracker/triptree");
TTree *t_kf_ev = f_kf->Get("CAtracker/evtree");

TTree *t_helix_trip = f_helix->Get("CAtracker/triptree");
TTree *t_helix_ev = f_helix->Get("CAtracker/evtree");

//Histogram definition
//number of seeds
TH1I *h_ns_kf = new TH1I("h_ns_kf","Number of multiseeds",10,0,10);
TH1I *h_ns_helix = new TH1I("h_ns_helix","Number of multiseeds",10,0,10);

//Local momentum
TH1I *h_lm_kf = new TH1I("h_lm_kf","Momentum at p0",100,0,100);
TH1I *h_lm_helix = new TH1I("h_lm_helix","Momentum at p0",100,0,100);

//Local momenutm (per layer)
TH1I *h_lml_kf[3];
TH1I *h_lml_helix[3];

for(int i = 0; i<3; i++){
h_lml_kf[i] = new TH1I("h_lml_kf","Momentum at p0 (per layer configuration) [FlexibleKF]",100,0,100);
h_lml_helix[i] = new TH1I("h_lml_helix","Momentum at p0 (per layer configuration) [FastHelix]",100,0,100);
}

    std::cout<<"nentries"<<std::endl;
    
int nentries_kf_trip = t_kf_trip->GetEntries();
int nentries_kf_ev = t_kf_ev->GetEntries();
int nentries_helix_trip = t_helix_trip->GetEntries();
int nentries_helix_ev = t_helix_ev->GetEntries();    
    
int nseeds_helix , nseeds_kf;
int nprod_helix , nprod_kf;
int identifier_helix , identifier_kf;
double momentum_helix , momentum_kf;

    std::cout<<"branch address"<<std::endl;

    
t_kf_ev->SetBranchAddress("nSeeds", &nseeds_kf);
t_kf_ev->SetBranchAddress("prodCells", &nprod_kf);

t_helix_ev->SetBranchAddress("nSeeds", &nseeds_helix);
t_helix_ev->SetBranchAddress("prodCells", &nprod_helix);

t_kf_trip->SetBranchAddress("LocalMomentum", &momentum_kf);
t_kf_trip->SetBranchAddress("Identifier", &identifier_kf);

t_helix_trip->SetBranchAddress("LocalMomentum", &momentum_helix);
t_helix_trip->SetBranchAddress("Identifier", &identifier_helix);

    
    std::cout<<"loops"<<std::endl;


for (int k = 0; k<nentries_kf_ev; k++){
    t_kf_ev->GetEntry(k);       
	if(nprod_kf >2)
		h_ns_kf->Fill(nseeds_kf);
}
for (int k = 0; k<nentries_helix_ev; k++){
    t_helix_ev->GetEntry(k);       
	if(nprod_helix >2)
		h_ns_helix->Fill(nseeds_helix);
}


for (int k = 0; k<nentries_kf_trip; k++){
    t_kf_trip->GetEntry(k);       
	if(identifier_kf==123)
		h_lml_kf[0]->Fill(momentum_kf);
	if(identifier_kf==234)
		h_lml_kf[1]->Fill(momentum_kf);	
	if(identifier_kf==345)
		h_lml_kf[2]->Fill(momentum_kf);
	h_lm_kf->Fill(momentum_kf);			
}
for (int k = 0; k<nentries_helix_trip; k++){
    t_helix_trip->GetEntry(k);       
	if(identifier_helix==123)
		h_lml_helix[0]->Fill(momentum_helix);
	if(identifier_helix==234)
		h_lml_helix[1]->Fill(momentum_helix);	
	if(identifier_helix==345)
		h_lml_helix[2]->Fill(momentum_helix);
	h_lm_helix->Fill(momentum_helix);		
}

h_ns_kf->SetLineColor(2);
h_ns_helix->SetLineColor(4);

h_lm_kf->SetLineColor(2);
h_lm_helix->SetLineColor(4);

h_lml_kf[0]->SetLineColor(2);
h_lml_kf[1]->SetLineColor(4);
h_lml_kf[2]->SetLineColor(6);

h_lml_helix[0]->SetLineColor(2);
h_lml_helix[1]->SetLineColor(4);
h_lml_helix[2]->SetLineColor(6);

	
TLegend *leg1 = new TLegend(0.6,0.7,0.9,0.9);
    leg1->SetFillColor(0);
    
TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
    leg2->SetFillColor(0);    

leg1->AddEntry(h_ns_kf,"FlexibleKFFitter","l");
leg1->AddEntry(h_ns_helix,"FastHelix","l");

leg2->AddEntry(h_lml_helix[0], "BPix1 + BPix2 + BPix3","l");
leg2->AddEntry(h_lml_helix[1], "BPix2 + BPix3 + TIB1","l");
leg2->AddEntry(h_lml_helix[2], "BPix3 + TIB1 + TIB2","l");

new TCanvas;
	h_ns_kf->Draw();
	h_ns_helix->Draw("same");
	leg1->Draw("same");
	
new TCanvas;
	h_lm_kf->Draw();
	h_lm_helix->Draw("same");
	leg1->Draw("same");
	
new TCanvas;
	h_lml_kf[0]->Draw();
	h_lml_kf[1]->Draw("same");
	h_lml_kf[2]->Draw("same");	
	leg2->Draw("same");

new TCanvas;
	h_lml_helix[0]->Draw();
	h_lml_helix[1]->Draw("same");
	h_lml_helix[2]->Draw("same");	
	leg2->Draw("same");	
	
	
}

	
