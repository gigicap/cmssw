// -*- C++ -*-
//
// Package:    CAtracker
// Class:      CAtracker
// 
/**\class CAtracker CAtracker.cc CAtracker/CAtracker/plugins/CAtracker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gigi Cappello
//         Created:  Wed, 27 Nov 2013 10:57:35 GMT
// $Id$
//
//


// system include files
#include "CAtracker.h"  //class CAcell defined here
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <unordered_map>

//
// class declaration
//

class CAtracker : public edm::EDProducer {
   public:
      explicit CAtracker(const edm::ParameterSet&);
      ~CAtracker();
    
        //typedef std::vector<TrackingRegion*> Regions;

        void fitTripletSeeds(CAcell *trip ,
                            const edm::EventSetup& iSetup,
                            edm::ESHandle<TransientTrackingRecHitBuilder>& builder,
                            edm::ESHandle<MagneticField>& mf, int seednumber);
        void CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder>& builder, TrajectorySeed seed,  int seednumber);
        int ForwardCA();
	void BackwardCA();
	CAcell* define_used(CAcell *cell, int id, TransientTrackingRecHit::RecHitContainer multicontainer);

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

//temp storage variables for debug tree 
double tmp_FitChiSquared;
int tmp_FitNdof;
int tmp_event;
int tmp_FitFoundHits;
double tmp_LocalMomentum;
double tmp_hqbp[3];
double tmp_hdxdz[3];
double tmp_hdydz[3];
int tmp_hitid0;
int tmp_hitid1;
int tmp_hitid2;
unsigned int tmp_rawid0;
unsigned int tmp_rawid1;
unsigned int tmp_rawid2;
    double tmp_dEta;
int tmp_identifier;
int tmp_neighborly;
int tmp_CAstatus;
int tmp_isused;
    
    int tmp_procTriplets;
    int tmp_prodCells;
    int tmp_prodCells123;
    int tmp_neigCells;
    double tmp_etaCValue;
    int tmp_fitCells;
    int tmp_nSteps;
    int tmp_nSeeds; 
    

    
//maximum CA_status
int max_status; //now can be set to = 3, 

   private:
    
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
    
        //Debug switch
        bool m_debug;
 	
        //CAcell collection
        edm::InputTag m_inputTagSeeds;
        std::vector<CAcell> tripletCollection;
        std::vector<CAcell*> fittedTripletCollection;

        //TrajectorySeedCollection seedcollection;       //to be set
    
    
       // edm::ParameterSet m_regionProducerPSet;

    
        std::string m_builderName;
        std::string m_fitterName;

        edm::ESHandle<TrajectoryFitter> m_fitter;
        edm::ESHandle<TrackerGeometry> m_tracker;
    
    double dEta_cut;

	TTree *triptree;
    TTree *evtree;
	int eventcounter;
	std::vector<TransientTrackingRecHit::RecHitContainer> multiplets;
       // std::unique_ptr<TrackingRegionProducer> m_trackingRegionProducer;   //???

      
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

      	void fittreegenerator();
        void fillfittree(CAcell *t);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    
        //data structure for the neighborhood maps
        std::unordered_map<unsigned int, std::list<CAcell *>> hitUsage;


      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
/*CAtracker::CAtracker(const edm::ParameterSet& iConfig) : m_debug(iConfig.getUntrackedParameter<bool>("debug",true)), m_inputTagSeeds(
    iConfig.getUntrackedParameter <edm::InputTag> ("SeedColl")), m_regionProducerPSet(iConfig.getParameter <edm::ParameterSet> ("regionProducerConfig")),
    m_builderName(iConfig.getUntrackedParameter <std::string> ("Builder")),
    m_fitterName(iConfig.getUntrackedParameter <std::string> ("Fitter"))
{*/
    
CAtracker::CAtracker(const edm::ParameterSet& iConfig) : m_debug(iConfig.getUntrackedParameter<bool>("debug",false)), m_inputTagSeeds(iConfig.getUntrackedParameter <edm::InputTag> ("SeedColl")), m_builderName(iConfig.getUntrackedParameter <std::string> ("Builder","WithAngleAndTemplate")), m_fitterName(iConfig.getUntrackedParameter <std::string> ("Fitter","FlexibleKFFittingSmoother")), dEta_cut(iConfig.getParameter<double>("EtaCut"))
    {

    std::cout<<"Constructor starts"<<std::endl;


	//Store tree generation
	edm::Service < TFileService > fs;
	triptree = fs->make<TTree> ("triptree","");
    evtree = fs->make<TTree> ("evtree","");
    fittreegenerator();

    eventcounter = 0;

	max_status = 3;
        
    //dEta_cut = 0.1;     //To be get from the cfi
  	//dEta_cut = iConfig.getParameter<double>("EtaCut");
}


CAtracker::~CAtracker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CAtracker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

//if (m_debug)
 std::cout<<"Produce"<<std::endl;




   // Regions regions = m_trackingRegionProducer->regions(iEvent, iSetup);
	//if (m_debug)
	//	std::cout << "Created " << regions.size() << " tracking regions" << std::endl;
    
	edm::ESHandle<TrackerGeometry> tracker;
	iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
        
    
	edm::ESHandle<TransientTrackingRecHitBuilder> builder;
	iSetup.get<TransientRecHitRecord>().get(m_builderName, builder);
      
    
	iSetup.get<TrajectoryFitter::Record>().get(m_fitterName, m_fitter);

    
	iSetup.get<TrackerDigiGeometryRecord>().get(m_tracker);
             
     
	edm::ESHandle<MagneticField> theMF;
	iSetup.get<IdealMagneticFieldRecord>().get(theMF);

        edm::Handle<TrajectorySeedCollection> seedcollection;
        iEvent.getByLabel(m_inputTagSeeds, seedcollection);
	TrajectorySeedCollection scoll = *(seedcollection.product());

tripletCollection.reserve(scoll.size());
fittedTripletCollection.reserve(scoll.size());


    //DEBUG
    /*if(scoll.size() > 10000){
        std::cout << "Too many seeds! Skipping event" << std::endl;
        return;
    }*/
    
//loop on seeds
 std::cout<<"loop starts"<<std::endl;
    for (size_t is = 0; is<scoll.size(); is++) {
        if (m_debug)
            std::cout << "Starting triplet fit " << std::endl;
        //
        CAcellGenerator(builder, scoll[is], is);
        
        hitUsage[tripletCollection[is].rawId0].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId1].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId2].push_back(&tripletCollection[is]);
        
    }
    
  
  std::cout<<"collectionSize = "<<scoll.size()<<std::endl;
    
 std::cout<<"NUMBER OF TRIPLETS: "<<tripletCollection.size()<<std::endl;
    
        int  prodCells123 = 0;
    
//Loop over the triplets to build the neighborhood map
//i.e. fill the "HasNeighbor" var and IsNeighborly
int ncont = 0;
for (size_t it = 0; it<tripletCollection.size(); it++) {
    int IsNeighbor = 0;
    
    
        std::list<CAcell *> list_h0 = hitUsage[tripletCollection[it].rawId0];
        std::list<CAcell *> list_h1 = hitUsage[tripletCollection[it].rawId1];
        std::list<CAcell *> list_h2 = hitUsage[tripletCollection[it].rawId2];
    
    
        std::vector<CAcell *> j_list01 = ListIntersect(list_h0, list_h1);
        std::vector<CAcell *> j_list12 = ListIntersect(list_h1, list_h2);
    
    
    if(m_debug){
    std::cout<<"list_h0 size == "<<list_h0.size()<<std::endl;
    std::cout<<"list_h1 size == "<<list_h1.size()<<std::endl;
    std::cout<<"list_h2 size == "<<list_h2.size()<<std::endl;


   std::cout<<"j_list01 size == "<<j_list01.size()<<std::endl;
   std::cout<<"j_list12 size == "<<j_list12.size()<<std::endl;
    
    }


     for (size_t il = 0; il<j_list01.size(); il++){
         if((j_list01[il] != &tripletCollection[it]) && (j_list01[il]->tripletIdentifier != tripletCollection[it].tripletIdentifier) && fabs(j_list01[il]->dEta - tripletCollection[it].dEta)< dEta_cut)
             tripletCollection[it].left_neighbors.push_back(j_list01[il]);
     }
    
    
    for (size_t il = 0; il<j_list12.size(); il++){
        if(j_list12[il] != &tripletCollection[it] && (j_list12[il]->tripletIdentifier != tripletCollection[it].tripletIdentifier) && fabs(j_list12[il]->dEta - tripletCollection[it].dEta)< dEta_cut)
            tripletCollection[it].right_neighbors.push_back(j_list12[il]);
    }
    
    if (m_debug){ 
    std::cout<<"left_list size == "<<tripletCollection[it].left_neighbors.size()<<std::endl;
    std::cout<<"right_list size == "<<tripletCollection[it].right_neighbors.size()<<std::endl;
    }
    std::cout<<"done "<<std::endl;

    
    IsNeighbor = tripletCollection[it].left_neighbors.size()+tripletCollection[it].right_neighbors.size();
    
    tripletCollection[it].IsNeighborly = IsNeighbor;
    
    //if(m_debug)
        std::cout<<"triplet has = "<<tripletCollection[it].IsNeighborly<<" neighbors. Is in lays: "<<tripletCollection[it].tripletIdentifier<<std::endl;
    
    
    list_h0.clear();
    list_h1.clear();
    list_h2.clear();
    j_list01.clear();
    j_list12.clear();
    
    
    //for debug
    if(tripletCollection[it].tripletIdentifier == 123)
        prodCells123++;

    
    //fit only the triplets with neighbors
    if (tripletCollection[it].IsNeighborly != 0) {
        ncont++;
        fitTripletSeeds(&tripletCollection[it], iSetup, builder, theMF, it);
 	//store the neighborly and fitted triplets
    	if (tripletCollection[it].FitSuccessful) {
            tripletCollection[it].CAstatus = 1;     //initialize CAstatus
            fittedTripletCollection.push_back(&tripletCollection[it]);
        }
    }
        

}

    std::cout<<"collection clear"<<std::endl;

//Loops on fittedtriplets -> Forward CA
    int n_fiter = -1;
    n_fiter = ForwardCA();
    
    std::cout<<"Forward CA: done n_fiter =="<<n_fiter<<std::endl;
//Connect triplets into multiplets Backward CA
    BackwardCA();
    
    std::cout<<"Backward CA: done"<<std::endl;
//Make seeds from multiplets
//initialize... 
//SeedFromConsecutiveHitsCreator seedcreator;
//TrajectorySeedCollection multiSeeds;

for(size_t im = 0; im < multiplets.size(); im++){
	std::vector<TransientTrackingRecHit::ConstRecHitPointer> multiHitPointer;
	for(size_t ii = 0; ii < multiplets[im].size(); ii++)
		multiHitPointer.push_back(multiplets[im][ii]);
    	//SeedingHitSet multiset(multiHitPointer[0], multiHitPointer[1], multiHitPointer[2], multiHitPointer[3], multiHitPointer[4]);
//	seedcreator.makeSeed(multiSeeds,multiset);
}

    
std::cout<<"Multiset generated"<<std::endl;
    
//fill the debug tree (not to be kept in the final code
    for (size_t itree = 0; itree<fittedTripletCollection.size(); itree++) {
        fillfittree(fittedTripletCollection[itree]);
    }
    

    //fill event tree
    tmp_procTriplets = scoll.size();
    tmp_prodCells = tripletCollection.size();
    tmp_prodCells123 = prodCells123;
    tmp_neigCells = ncont;
    tmp_etaCValue = dEta_cut;
    tmp_fitCells = fittedTripletCollection.size();
    tmp_nSteps = n_fiter;
    tmp_nSeeds = multiplets.size();
    
    evtree->Fill();    

    
    
    if(m_debug){
    std::cout<<"================================"<<std::endl;
    std::cout<<"Processed "<<scoll.size()<<" seeds"<<std::endl;
    std::cout<<"produced  "<<tripletCollection.size()<<" CAcells, of whom:"<<std::endl;
    std::cout<<ncont<<" have neighbors"<<std::endl;
    std::cout<<fittedTripletCollection.size()<<" have been successfully fitted!"<<std::endl;
    std::cout<<"Forward iterations took n_steps = "<<n_fiter<<std::endl;

    std::cout<<"================================"<<std::endl;
    }
        
    
    scoll.clear();
    tripletCollection.clear();
    
    fittedTripletCollection.clear();
    
    hitUsage.clear();
    multiplets.clear();
    
    eventcounter++;

}


// ------------ method called once each job just before starting event loop  ------------
void 
CAtracker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CAtracker::endJob() {
}


void CAtracker::CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder> & builder, TrajectorySeed seed,  int seednumber){
    //Build the CA cells, i.e. triplets (objects containing RecHits and fit parameters
    CAcell trip;

    //HITS
    std::vector<TransientTrackingRecHit::RecHitPointer> recHitPointers;
    
    //std::cout << " B " << std::endl;
    // loop RecHits
    const TrajectorySeed::range theHitsRange = seed.recHits();
    int ip = 0;
    for( edm::OwnVector<TrackingRecHit>::const_iterator aHit = theHitsRange.first;
        aHit < theHitsRange.second; ++aHit ) {
        // this is a collection of: ReferenceCountingPointer< TransientTrackingRecHit>
        recHitPointers.push_back( builder->build( &(*aHit) ) );
		ip++;
    }
    
    
    
    //TransientTrackingRecHit hits for fits
    TransientTrackingRecHit::RecHitContainer hits;
    
    hits.push_back(recHitPointers.at(0));
    hits.push_back(recHitPointers.at(1));
    hits.push_back(recHitPointers.at(2));
    
    
    //TrackinkgRecHits for fast Id access
    SeedingHitSet sset(recHitPointers[0], recHitPointers[1], recHitPointers[2]);
    
    const TrackingRecHit* hi0 = sset[0]->hit();
    const TrackingRecHit* hi1 = sset[1]->hit();
    const TrackingRecHit* hi2 = sset[2]->hit();
    
    
    //fill triplet hits
    trip.hits = hits;
    
    
    trip.rawId0 = OmniRef(hi0);
    trip.rawId1 = OmniRef(hi1);
    trip.rawId2 = OmniRef(hi2);

    
    GlobalPoint p0 = hits[0]->globalPosition();
    GlobalPoint p2 = hits[2]->globalPosition();
    
    trip.dEta = DeltaEta(p0,p2);
    
    DetIDClassification hit0cl(hi0->geographicalId());
    DetIDClassification hit1cl(hi1->geographicalId());
    DetIDClassification hit2cl(hi2->geographicalId());
    trip.hitId0 = hit0cl.getLayer();
    trip.hitId1 = hit1cl.getLayer();
    trip.hitId2 = hit2cl.getLayer();


    trip.eventNumber = eventcounter;
    trip.tripletIdentifier = 100*trip.hitId0+10*trip.hitId1+trip.hitId2;
    
    if (m_debug) {
        std::cout << " layer 0 =  " << hit0cl.getLayer() << std::endl;
        std::cout << " layer 1 =  " << hit1cl.getLayer() << std::endl;
        std::cout << " layer 2 =  " << hit2cl.getLayer() << std::endl;
        std::cout << " identifier = "<< trip.tripletIdentifier <<std::endl;
    }
    
    
    //Fill teh hit the "where used" structure
    /*    hitUsage[trip.rawId0].push_back(&trip);
     hitUsage[trip.rawId1].push_back(&trip);
     hitUsage[trip.rawId2].push_back(&trip);*/
    
    trip.TripletNumber = seednumber;


    //Useful info for fitting
    trip.seedStartState = seed.startingState();
    
    //initialize fit variables to 0
    /*trip.FitHit0Local = traj.measurements()[2].updatedState().localParameters();
	trip.FitHit0LocalError = traj.measurements()[2].updatedState().localError();
    
	trip.FitHit1Local = traj.measurements()[1].updatedState().localParameters();
	trip.FitHit1LocalError = traj.measurements()[1].updatedState().localError();
    
	trip.FitHit2Local = traj.measurements()[0].updatedState().localParameters();
	trip.FitHit2LocalError = traj.measurements()[0].updatedState().localError();*/
    
	trip.FitSuccessful = false;
	trip.FitChiSquared = 0;
	trip.FitNdof = 0;
	trip.FitFoundHits = 0;
	trip.LocalMomentum = 0;

	//Is Used (for the final Backwards CA)
	trip.IsUsed = 0;
    
    
    tripletCollection.push_back(trip);
}



void CAtracker::fitTripletSeeds(CAcell *trip, const edm::EventSetup& es, edm::ESHandle<TransientTrackingRecHitBuilder> & builder,
                                   edm::ESHandle<MagneticField> & mf, int seednumber) {
    

TransientTrackingRecHit::RecHitContainer hits;

//std::cout << " size =  "<< recHitPointers.size() << std::endl;
hits.push_back(trip->hits[0]);
hits.push_back(trip->hits[1]);
hits.push_back(trip->hits[2]);


//reverse hits
TransientTrackingRecHit::RecHitContainer hitsContReverse;
hitsContReverse.push_back(hits[2]);
hitsContReverse.push_back(hits[1]);
hitsContReverse.push_back(hits[0]);
 
//std::cout << " E " << std::endl;
   
 
//std::cout << " F " << std::endl;
/****Define TrajectoryStateOnSurface*/
 
    
// Transform it in a TrajectoryStateOnSurface
    PTrajectoryStateOnDet seedStartState = trip->seedStartState;
    DetId seedDetId(seedStartState.detId());

//std::cout << " G " << std::endl;
    // important: the seed state will be defined at the end of the seed ( 3rd hit )
    BoundPlane const& firstHitSurface = m_tracker->idToDet(seedDetId)->surface();
//std::cout << " H " << std::endl;
    TrajectoryStateOnSurface startState =  trajectoryStateTransform::transientState(seedStartState, &firstHitSurface, &*mf );
// std::cout << " I " << std::endl;
    // generate a seed with opposite momentum
	TrajectorySeed reverseSeed(PTrajectoryStateOnDet(), TrajectorySeed::recHitContainer(), oppositeToMomentum);
    
	if (m_debug) {
		std::cout << "Starting Seed - p = " << startState.localMomentum().mag()
        << std::endl;
		std::cout << "Starting Seed - Z = " << startState.globalMomentum().z()
        << std::endl;
		std::cout << "Starting Seed - phi = "
        << startState.globalMomentum().phi() << std::endl;
		std::cout << "Starting Seed - R = "
        << startState.globalPosition().perp() << std::endl;
        
		std::cout << "TSOS from seed " << std::endl << startState << std::endl;
	}

//std::cout << " J " << std::endl;
   
	assert(startState.isValid());

	TrajectoryStateOnSurface startStateObj = TrajectoryStateOnSurface(startState.localParameters(), startState.localError(),
                                                                      startState.surface(), mf.product());
    
	startStateObj.rescaleError(100);

//std::cout << " K " << std::endl;
	
	std::vector<Trajectory> fittedTrajets = m_fitter->fit(reverseSeed, hitsContReverse, startStateObj);
    
	if (m_debug)
		std::cout << "fit done" << std::endl;
    
	if (fittedTrajets.size() != 1){
	//	std::cout << "Problem: fit returned " << fittedTrajets.size()
    //    << std::endl;
		return;
	}
    
	auto traj = fittedTrajets[0];
    
	if (!traj.isValid()) {
		std::cout << "not valid trajectory" << std::endl;
		return;
	}
    
//std::cout << " L " << std::endl;

	auto innerFit = traj.lastMeasurement();
	auto outerFit = traj.firstMeasurement();
    
	if (m_debug) {
		std::cout << "-- is valid: " << traj.isValid() << std::endl;
		std::cout << "-- found hits: " << traj.foundHits() << std::endl;
		std::cout << "-- chi squared: " << traj.chiSquared() << std::endl;
        
		std::cout << "-- local Momentum "
        << traj.firstMeasurement().forwardPredictedState().localMomentum()
        << std::endl;
        
		std::cout << "Innermost fit parameter" << std::endl;
		std::cout << innerFit.updatedState() << std::endl;
		std::cout << "Outermost fit parameter" << std::endl;
		std::cout << outerFit.updatedState() << std::endl;
	}
 
    
    
	trip->FitHit0Local = traj.measurements()[2].updatedState().localParameters();
	trip->FitHit0LocalError = traj.measurements()[2].updatedState().localError();
    
	trip->FitHit1Local = traj.measurements()[1].updatedState().localParameters();
	trip->FitHit1LocalError = traj.measurements()[1].updatedState().localError();
    
	trip->FitHit2Local = traj.measurements()[0].updatedState().localParameters();
	trip->FitHit2LocalError = traj.measurements()[0].updatedState().localError();
    
	trip->FitSuccessful = true;
	trip->FitChiSquared = traj.chiSquared();
	trip->FitNdof = traj.ndof();
	trip->FitFoundHits = traj.foundHits();
	trip->LocalMomentum = traj.firstMeasurement().forwardPredictedState().localMomentum().mag();

    
	return;
}



int CAtracker::ForwardCA(){
    int Stop = -1;
    int step_iterator = 0;
    
    int it_max = 100; //Maximum number of iteraitons
        
    while (Stop < (int)fittedTripletCollection.size() && step_iterator<it_max) {
        
        step_iterator++;
        Stop = 0;
        
        
        for (size_t t =0; t<fittedTripletCollection.size(); t++) {
            //std::cout<<" ________ "<<std::endl;
            //std::cout<<"triplet no. "<<t<<std::endl;
            if (fittedTripletCollection[t]->left_neighbors.size()!=0) {
                for (size_t il = 0; il<fittedTripletCollection[t]->left_neighbors.size(); il++) {
                    bool en_cont = false;   //useful variable
			if(m_debug){
                    	std::cout<<"neighbor list no = "<<il<<std::endl;
                    	std::cout<<"fittedTriplet status = "<<fittedTripletCollection[t]->CAstatus<<"left_neighbors status = "<<fittedTripletCollection[t]->left_neighbors[il]->CAstatus <<std::endl;
			}
                    if (fittedTripletCollection[t]->CAstatus <= fittedTripletCollection[t]->left_neighbors[il]->CAstatus && en_cont == false){
                        (fittedTripletCollection[t]->CAstatus)++;
                        //std::cout<<"NEW fittedTriplet status = "<<fittedTripletCollection[t]->CAstatus<<std::endl;
                        en_cont = true;
                    }
                if (en_cont == false)
                    Stop++;
                
                }
              
            }
            else
                Stop++;
        }
        
        //std::cout<<"Stop value = "<<Stop<<std::endl;
    }
    
    //debug
if(m_debug){
    for (size_t t =0; t<fittedTripletCollection.size(); t++) {
        std::cout<<" ._._._._._._._._ "<<std::endl;
        std::cout<<" cell number "<<t<<std::endl;
        std::cout<<" processed with status: "<<fittedTripletCollection[t]->CAstatus<<std::endl;
        
    }
    }
    
    if(step_iterator>=it_max)
        std::cout<<"WARNING: too many CA iterations!"<<std::endl;
    
    return step_iterator;
}


void CAtracker::BackwardCA(){

int multi_id = 1;

for (size_t t =0; t<fittedTripletCollection.size(); t++) {
//start only with non-used triplets with maximum CAstatus
if(fittedTripletCollection[t]->CAstatus >= max_status && fittedTripletCollection[t]->IsUsed==0){
	CAcell *current_cell = fittedTripletCollection[t];
	TransientTrackingRecHit::RecHitContainer multicontainer;
	while(current_cell){
		current_cell = define_used(current_cell, multi_id, multicontainer);
		}
		multi_id++;	
	multiplets.push_back(multicontainer);
	multicontainer.clear();
	}
}

return;
}

CAcell* CAtracker::define_used(CAcell *cell, int id, TransientTrackingRecHit::RecHitContainer multicontainer){
cell->IsUsed = id;
double pt_tmp = 99999.0;

if(cell->CAstatus==1){
	multicontainer.push_back(cell->hits[2]);
	multicontainer.push_back(cell->hits[1]);
	multicontainer.push_back(cell->hits[0]);
return NULL;	
}

else{
CAcell *next = new CAcell();
multicontainer.push_back(cell->hits[2]);
	for(size_t il = 0; il< cell->left_neighbors.size(); il++){
		if(cell->CAstatus - cell->left_neighbors[il]->CAstatus == 1){
		double pt_diff = fabs(cell->getLocalMomentum() - cell->left_neighbors[il]->getLocalMomentum());
		if(pt_diff <= pt_tmp){
			pt_tmp = pt_diff;
			next = cell->left_neighbors[il];
			} 
		}
	}
return next;
}

}


// ------------ method called when starting to processes a run  ------------

void CAtracker::beginRun(edm::Run const&, edm::EventSetup const&) {
 
  
 /*std::string regfactoryName = m_regionProducerPSet.getParameter < std::string
 > ("ComponentName");
 m_trackingRegionProducer.reset(
 TrackingRegionProducerFactory::get()->create(regfactoryName,
 m_regionProducerPSet));
 */
 //theSeedCreator->trajectorySeed(seedCollection, hits, region, es,
 //			theComparitor);
 }


// ------------ method called when ending the processing of a run  ------------

void CAtracker::endRun(edm::Run const&, edm::EventSetup const&)
{
}

 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
CAtracker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
CAtracker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 

//Useful tools to store triplets in root files
void CAtracker::fittreegenerator(){

triptree->Branch("Event",&tmp_event,"Event/I");
triptree->Branch("FitChiSquared",&tmp_FitChiSquared,"FitChiSquared/D");
triptree->Branch("FitNDof",&tmp_FitNdof,"FitNdof/I");
triptree->Branch("FitFoundHits",&tmp_FitFoundHits,"FitFoundHits/I");
triptree->Branch("LocalMomentum",&tmp_LocalMomentum,"LocalMomentum/D");

triptree->Branch("hqbp",tmp_hqbp,"hqbp[3]/D");
triptree->Branch("hdxdz",tmp_hdxdz,"hdxdz[3]/D");
triptree->Branch("hdydz",tmp_hdydz,"hdydz[3]/D");

triptree->Branch("hitId0",&tmp_hitid0,"hitId0/I");
triptree->Branch("hitId1",&tmp_hitid1,"hitId1/I");
triptree->Branch("hitId2",&tmp_hitid2,"hitId2/I");
    
triptree->Branch("rawId0",&tmp_rawid0,"rawId0/I");
triptree->Branch("rawId1",&tmp_rawid1,"rawId1/I");
triptree->Branch("rawId2",&tmp_rawid2,"rawId2/I");
    
triptree->Branch("Eta",&tmp_dEta,"Eta/D");
    
    
triptree->Branch("Neighborly",&tmp_neighborly,"Neighborly/I");
triptree->Branch("CAstatus",&tmp_CAstatus,"CAstatus/I");
    
    
triptree->Branch("Identifier",&tmp_identifier,"Identifier/I");

triptree->Branch("IsUsed",&tmp_isused,"IsUsed/I");
    
    
    //event tree
    evtree->Branch("procTriplets",&tmp_procTriplets,"procTriplets/I");
    evtree->Branch("prodCells",&tmp_prodCells,"prodCells/I");
    evtree->Branch("prodCells123",&tmp_prodCells123,"prodCells123/I");
    evtree->Branch("neigCells",&tmp_neigCells,"neigCells/I");
    evtree->Branch("etaCValue",&tmp_etaCValue,"etaCValue/D");
    evtree->Branch("fitCells",&tmp_fitCells,"fitCells/I");
    evtree->Branch("nSteps",&tmp_nSteps,"nSteps/I");
    evtree->Branch("nSeeds",&tmp_nSeeds,"nSeeds/I");

    

return;
}

void CAtracker::fillfittree(CAcell *t){

tmp_FitChiSquared = t->getFitChiSquared();
tmp_FitNdof = t->getFitNdof();
tmp_FitFoundHits = t->getFitFoundHits();
tmp_LocalMomentum = t->getLocalMomentum();


tmp_hqbp[0] = t->geth0qbp();
tmp_hqbp[1] = t->geth1qbp();
tmp_hqbp[2] = t->geth2qbp();

tmp_hdxdz[0] = t->geth0dxdz();
tmp_hdxdz[1] = t->geth1dxdz();
tmp_hdxdz[2] = t->geth2dxdz();

tmp_hdydz[0] = t->geth0dydz();
tmp_hdydz[1] = t->geth1dydz();
tmp_hdydz[2] = t->geth2dydz();

tmp_event = t->geteventNumber();

tmp_hitid0 = t->gethitId0();
tmp_hitid1 = t->gethitId1();
tmp_hitid2 = t->gethitId2();
    
tmp_rawid0 = t->rawId0;
tmp_rawid1 = t->rawId1;
tmp_rawid2 = t->rawId2;
    
tmp_dEta = t->dEta;
    
tmp_neighborly = t->IsNeighborly;
tmp_identifier = t->tripletIdentifier;

tmp_CAstatus = t->CAstatus;

tmp_isused = t->IsUsed;
    
triptree->Fill();
    

return;
}




// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CAtracker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CAtracker);
