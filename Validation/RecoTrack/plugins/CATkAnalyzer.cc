// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

//
// class decleration
//
class CATkAnalyzer : public edm::EDAnalyzer {
public:
  explicit CATkAnalyzer(const edm::ParameterSet&);
  ~CATkAnalyzer();  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  TFile* file;
    
  TH1D* histoPt;
  TH1D* histoChi2;
  TH1D* histoD0;
  TH1D* histoDz;
    
  TH1I* histoCharge;
  TH1I* histoNum;
  TH1I* histoNdof;
  TH1I* histoNhits;
    
  TH1I* histoLayer;  
  TH1I* histoNumLayer;
  TH2I* histoLayComb;	
  TH1I* histoLastLayer;
};

//
// constructors and destructor
//
CATkAnalyzer::CATkAnalyzer(const edm::ParameterSet& iConfig)
:trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))
{
file = new TFile("outfile.root","recreate");
  histoPt = new TH1D("Pt","Pt",100,60.,140.);
  histoNum = new TH1I("Num","Number of tracks",9,0,9);
  histoChi2 = new TH1D("Chi2","Chi2",100,0,50);
  histoCharge = new TH1I("Charge","Track charge",6,-3,3);
  histoNdof = new TH1I("Ndof","n_dof",100,0,40);
  histoD0 = new TH1D("D0","D0",100,-0.05,0.05);
  histoDz = new TH1D("Dz","Dz",100,-5,5);
  histoNhits = new TH1I("Nhits","Number of valid hits",30,0,30);
  histoLayer = new TH1I("Layer","Hit layer",20,0,20);
  histoLastLayer = new TH1I("LastLayer","Last visited layer",20,0,20);
  histoNumLayer = new TH1I("NumLayer","Number of visited layers",20,0,20);
  histoLayComb = new TH2I("LayComb","Number of visited layers vs Number of valid hits",30,0,30,30,0,30);
}

CATkAnalyzer::~CATkAnalyzer()
{
}

// ------------ method called to for each event  ------------
void CATkAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using reco::TrackCollection;
  using reco::BeamSpot;
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
//beamspot, for d0 evaluation
  Handle<BeamSpot> bs;
  iEvent.getByLabel("offlineBeamSpot",bs);

int n = 0;
  for(TrackCollection::const_iterator itTrack = tracks->begin();
      itTrack != tracks->end();                      
      ++itTrack) {
            
      
        double pt = itTrack->pt();
        double chi2 = itTrack->chi2();
        double d0 = itTrack->dxy(bs->position());
        double dz = itTrack->dz();
        int charge = itTrack->charge();
        int ndof = itTrack->ndof();
	//int nhits = itTrack->numberOfValidHits();
          int nhits = 0;

//access to hit layers
	// hit pattern of the track
	const reco::HitPattern& p = itTrack->hitPattern();
	int pass_layer = 1;
        int last_layer = 1;
	unsigned int prev_layer = 1;
	// loop over the hits of the track
	for (int i=0; i<p.numberOfHits(); i++) {
  		uint32_t hit = p.getHitPattern(i);

  		// if the hit is valid and in pixel barrel, print out the layer

  	if (p.validHitFilter(hit)){
                nhits++;
    		histoLayer->Fill(p.getLayer(hit));
                last_layer = p.getLayer(hit);
		if (p.getLayer(hit)!= prev_layer){
			prev_layer = p.getLayer(hit);
			pass_layer++;
			}
		}
	}
       
        histoNumLayer->Fill(pass_layer);
	histoLayComb->Fill(pass_layer,nhits);

 //    if(pt>0.4){
      histoPt->Fill(pt);
      histoChi2->Fill(chi2);  
      histoCharge->Fill(charge);
      histoNdof->Fill(ndof);
      histoD0->Fill(d0);
      histoDz->Fill(dz);
      histoNhits->Fill(nhits);
      histoLastLayer->Fill(last_layer);
      n++;
   //     }



	  }

    histoNum->Fill(n);
}

void CATkAnalyzer::beginJob(const edm::EventSetup&){
}

void CATkAnalyzer::endJob() {
    
    file->WriteTObject(histoNum);
    file->WriteTObject(histoPt);
    file->WriteTObject(histoChi2);
    file->WriteTObject(histoNdof);
    file->WriteTObject(histoCharge);
    file->WriteTObject(histoD0);
    file->WriteTObject(histoDz);
    file->WriteTObject(histoNhits);
    file->WriteTObject(histoLayer);
    file->WriteTObject(histoNumLayer);
    file->WriteTObject(histoLayComb);
    file->Close();
    


}

//define this as a plug-in
DEFINE_FWK_MODULE(CATkAnalyzer);

