// -*- C++ -*-
//
// Package:    CAHitsGenerator
// Class:      CAHitsGenerator
// 
/**\class CAHitsGenerator CAHitsGenerator.cc CAtracker/CAtracker/plugins/CAHitsGenerator.cc*/
//
// Original Author:  Gigi Cappello
//         Created:  Wed, 27 Nov 2013 10:57:35 GMT
//

// system include files
#include "CAHitsGenerator.h"  //class CAcell defined here
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace ctfseeding;
    

CAHitsGenerator::CAHitsGenerator(const edm::ParameterSet& cfg) :
m_debug(cfg.getUntrackedParameter<bool>("debug")),
m_builderName(cfg.getUntrackedParameter <std::string> ("Builder")),
dEta_cut(cfg.getParameter<double>("EtaCut"))
    {

   // aPset = cfg.getParameter<edm::ParameterSet>("CAHitsGenerator");
     
   // m_debug = aPset.getUntrackedParameter<bool>("debug");
   // m_builderName = aPset.getUntrackedParameter <std::string> ("Builder");
   // dEta_cut = aPset.getParameter<double>("EtaCut");
    std::cout<<"Constructor starts"<<std::endl;
        
    //dEta_cut = 0.0256;

    //triplet producer definition (???)
    oPset = cfg.getParameter<edm::ParameterSet>("GeneratorPSet");
    oName = oPset.getParameter<std::string>("ComponentName");

	std::cout<<"generating a : "<<oName<<std::endl;

	//Store tree generation
	edm::Service < TFileService > fs;
	triptree = fs->make<TTree> ("triptree","");
	evtree = fs->make<TTree> ("evtree","");
    fittreegenerator();

    std::cout<<"Constructor done"<<std::endl;
        
        
        std::cout<<"setting seedinglayer names"<<std::endl;

	//seedinglayersets
 	layerBuilderName = cfg.getParameter<std::string>("SeedingLayers");




        
        //generator = HitTripletGeneratorFromPairAndLayersFactory::get()->create(oName,oPset);
        //generator = OrderedHitsGeneratorFactory::get()->create(oName,oPset);
        //generator = PixelTripletHLTGenerator::get()->create(oName,oPset);
        //generator = new PixelTripletHLTGenerator(oPset);
        
        //	scoll = *(tripletcollection.product());
        
        std::cout<<"seedinglayerset Has been set"<<std::endl;
        
        
	//region producer

/*	edm::ParameterSet regfactoryPSet = cfg.getParameter<edm::ParameterSet>("RegionFactoryPSet");
	std::string regfactoryName = regfactoryPSet.getParameter<std::string>("ComponentName");
	//Regions = TrackingRegionProducerFactory::get()->create(regfactoryName,regfactoryPSet,consumesCollector());
	Regions = TrackingRegionProducerFactory::get()->create(regfactoryName,regfactoryPSet);
        std::cout<<"Regions set"<<std::endl;*/


    eventcounter = 0;
    max_status = 3;
	
    initialised = false;
}


CAHitsGenerator::~CAHitsGenerator()
{
  GeneratorContainer::const_iterator it;
   for (it = theGenerators.begin(); it!= theGenerators.end(); it++) {
     delete (*it);
   }
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
//just a dummy slice of code...
void CAHitsGenerator::init(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

    //if (m_debug)
    std::cout<<"Begin init---"<<std::endl;
    
    std::cout<<"getting all parameters"<<std::endl;  

	   //!!!!!!
    	//generating regions


    
  /*  
    // Regions regions = m_trackingRegionProducer->regions(iEvent, iSetup);
	//if (m_debug)
	//	std::cout << "Created " << regions.size() << " tracking regions" << std::endl;
    
    
	iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
 	iSetup.get<TransientRecHitRecord>().get(m_builderName, builder);
	iSetup.get<TrackerDigiGeometryRecord>().get(m_tracker);
	iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    
    
    //edm::Handle<OrderedHitTriplets> tripletcollection;
    //iEvent.getByLabel(m_inputTagTriplets, tripletcollection);
    
    generator = HitTripletGeneratorFromPairAndLayersFactory::get()->create(oName,oPset);
    
//	scoll = *(tripletcollection.product());
    
    std::cout<<"Generator Has been set"<<std::endl;
    
//	tripletCollection.reserve(scoll.size());
//	fittedTripletCollection.reserve(scoll.size());
//	result.clear();

    */
    
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  iSetup.get<TransientRecHitRecord>().get(m_builderName, builder);
 // iSetup.get<TrackerDigiGeometryRecord>().get(m_tracker);
  iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    
   iSetup.get<TrackerDigiGeometryRecord>().get(layerBuilderName, layerBuilder);
   SeedingLayerSets layerSets  =  layerBuilder->layers(iSetup); 
 
   std::vector<LayerTriplets::LayerPairAndLayers>::const_iterator it;
   std::vector<LayerTriplets::LayerPairAndLayers> trilayers=LayerTriplets(layerSets).layers();
 
   for (it = trilayers.begin(); it != trilayers.end(); it++) {
     SeedingLayer first = (*it).first.first;
     SeedingLayer second = (*it).first.second;
     std::vector<SeedingLayer> thirds = (*it).second;
 
 	//building generators
	 std::cout<<"Generator Has to be defined"<<std::endl;
     HitTripletGeneratorFromPairAndLayers * aGen = HitTripletGeneratorFromPairAndLayersFactory::get()->create(oName,oPset);
     aGen->init( HitPairGeneratorFromLayerPair( first, second, &theLayerCache),thirds, &theLayerCache);
 	 std::cout<<"Generator Has been set"<<std::endl;
     theGenerators.push_back( aGen);
     }
 
initialised = true;
    
return;
}


void CAHitsGenerator::hitSets(const TrackingRegion& region, OrderedMultiHits& result, const edm::Event& iEvent, const edm::EventSetup& iSetup)
//void CAHitsGenerator::hitSets(OrderedMultiHits& result, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
    
if(!initialised) init(iEvent,iSetup);

    //From Init
    //if (m_debug)
 //   std::cout<<"Begin init---"<<std::endl;
    
 //   std::cout<<"getting all parameters"<<std::endl;
    
    
    // Regions regions = m_trackingRegionProducer->regions(iEvent, iSetup);
	//if (m_debug)
	//	std::cout << "Created " << regions.size() << " tracking regions" << std::endl;
    

    
    //edm::Handle<OrderedHitTriplets> tripletcollection;
    //iEvent.getByLabel(m_inputTagTriplets, tripletcollection);
    
 
/*	std::cout<<"Generating regions"<<std::endl;
    	typedef std::vector<TrackingRegion* > theRegions;
	theRegions regions = Regions->regions(iEvent,iSetup);
	TrackingRegion& region1 = *regions[0];  
        
    
    */
    
    std::cout<<"Generating triplets"<<std::endl;

    //produce triplets
    //const OrderedSeedingHits & scoll1 = generator->run(region1,iEvent,iSetup);
    GeneratorContainer::const_iterator i;
   int gensize = 0;
   for (i=theGenerators.begin(); i!=theGenerators.end(); i++) {
     (**i).hitTriplets(region, scoll, iEvent, iSetup);
	gensize++;
   }
   theLayerCache.clear();
	//scoll = generator->run(region,iEvent,iSetup);

    std::cout<<"generator size = "<<gensize<<std::endl;
    std::cout<<scoll.size()<<std::endl;

    std::cout<<"Generated!"<<std::endl;

    
    tripletCollection.reserve(scoll.size());
    fittedTripletCollection.reserve(scoll.size());
    

    for (size_t is = 0; is<scoll.size(); is++) {
        if (m_debug)
            std::cout << "Starting hit usage " << std::endl;
     
        	
   	//SeedingHitSet sset(scoll[is].inner(), scoll[is].middle(), scoll[is].outer());
	SeedingHitSet sset = scoll[is];

        CAcellGenerator(builder, sset, is);
        
        hitUsage[tripletCollection[is].rawId0].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId1].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId2].push_back(&tripletCollection[is]);
        
    }
    
  
//  std::cout<<"collectionSize = "<<scoll.size()<<std::endl;
    
// std::cout<<"NUMBER OF TRIPLETS: "<<tripletCollection.size()<<std::endl;
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
 //   std::cout<<"done "<<std::endl;

    
    IsNeighbor = tripletCollection[it].left_neighbors.size()+tripletCollection[it].right_neighbors.size();
    
    tripletCollection[it].IsNeighborly = IsNeighbor;
    
    if(m_debug)
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
    // fast helix fitter
        fitTripletSeeds(&tripletCollection[it], theMF);
 	//store the neighborly and fitted triplets
    	if (tripletCollection[it].FitSuccessful) {
            tripletCollection[it].CAstatus = 1;     //initialize CAstatus
            fittedTripletCollection.push_back(&tripletCollection[it]);
        }
    }
        

}


//Loops on fittedtriplets -> Forward CA
    int n_fiter = -1;
    n_fiter = ForwardCA();
//Connect triplets into multiplets Backward CA
    BackwardCA();
    
//Make seeds from multiplets
//initialize... 
//SeedFromConsecutiveHitsCreator seedcreator;
//TrajectorySeedCollection multiSeeds;

for(size_t im = 0; im < multiplets.size(); im++){
	std::vector<TransientTrackingRecHit::ConstRecHitPointer> multiHitPointer;
	for(size_t ii = 0; ii < multiplets[im].size(); ii++)
		multiHitPointer.push_back(multiplets[im][ii]);
	if(multiHitPointer.size()!=5)
		std::cout<<"WARNING: multiset size = "<<multiHitPointer.size()<<std::endl;
	else{
    		SeedingHitSet multiset(multiHitPointer[0], multiHitPointer[1], multiHitPointer[2], multiHitPointer[3], multiHitPointer[4]);
		result.push_back(multiset);
	}
	//seedcreator.makeSeed(multiSeeds,multiset);
}

    
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
tmp_nSeeds = result.size(); 

evtree->Fill();    
    
    if(m_debug){
    std::cout<<"================================"<<std::endl;
    std::cout<<"Processed "<<scoll.size()<<" seeds"<<std::endl;
    std::cout<<"produced  "<<tripletCollection.size()<<" CAcells, of whom:"<<std::endl;
    std::cout<<ncont<<" have neighbors"<<std::endl;
    std::cout<<fittedTripletCollection.size()<<" have been successfully fitted!"<<std::endl;
    std::cout<<"Forward iterations took n_steps = "<<n_fiter<<std::endl;
    std::cout<<"Final number of multiSeeds = "<<result.size()<<std::endl;
    std::cout<<"================================"<<std::endl;
    }
/*********************************************/        
    
    scoll.clear();
    tripletCollection.clear();
    
    fittedTripletCollection.clear();
    
    hitUsage.clear();
    multiplets.clear();
    
    result.clear();

    eventcounter++;

}


void CAHitsGenerator::CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder> & builder, SeedingHitSet sset,  int seednumber){
    //Build the CA cells, i.e. triplets (objects containing RecHits and fit parameters
    CAcell trip;
    
    //translate hits into a SeedingHitSet
	trip.hits = sset;
    
    //TrackinkgRecHits for fast Id access
    
    const TrackingRecHit* hi0 = sset[0]->hit();
    const TrackingRecHit* hi1 = sset[1]->hit();
    const TrackingRecHit* hi2 = sset[2]->hit();
    
 
    
    trip.rawId0 = OmniRef(hi0);
    trip.rawId1 = OmniRef(hi1);
    trip.rawId2 = OmniRef(hi2);


    GlobalPoint p0 = sset[0]->globalPosition();
    GlobalPoint p1 = sset[1]->globalPosition();
    GlobalPoint p2 = sset[2]->globalPosition();


    trip.p0 = p0;	
    trip.p1 = p1;	   
    trip.p2 = p2;

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
    //trip.seedStartState = seed.startingState();
    
    //initialize fit variables to 0
    /*trip.FitHit0Local = traj.measurements()[2].updatedState().localParameters();
	trip.FitHit0LocalError = traj.measurements()[2].updatedState().localError();
    
	trip.FitHit1Local = traj.measurements()[1].updatedState().localParameters();
	trip.FitHit1LocalError = traj.measurements()[1].updatedState().localError();
    
	trip.FitHit2Local = traj.measurements()[0].updatedState().localParameters();
	trip.FitHit2LocalError = traj.measurements()[0].updatedState().localError();*/
    
	trip.FitSuccessful = false;
	//trip.FitChiSquared = 0;
	//trip.FitNdof = 0;
	//trip.FitFoundHits = 0;
	trip.LocalMomentum = 0;

	//Is Used (for the final Backwards CA)
	trip.IsUsed = 0;
    
    
    tripletCollection.push_back(trip);
}


void CAHitsGenerator::fitTripletSeeds(CAcell *trip, edm::ESHandle<MagneticField> & mf) {

 
    
    double nomField = mf->nominalValue();

    FastHelix fit(trip->p2, trip->p1, trip->p0, nomField , mf.product());
    
    GlobalTrajectoryParameters params = fit.stateAtVertex();
    
    if (fit.isValid()){
    trip->FitSuccessful = true;
	trip->LocalMomentum = params.momentum().mag();
    }

return;
}



int CAHitsGenerator::ForwardCA(){
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


void CAHitsGenerator::BackwardCA(){

int multi_id = 1;

for (size_t t =0; t<fittedTripletCollection.size(); t++) {
//start only with non-used triplets with maximum CAstatus
if(fittedTripletCollection[t]->CAstatus >= max_status && fittedTripletCollection[t]->IsUsed==0){
	CAcell *current_cell = fittedTripletCollection[t];
	TransientTrackingRecHit::RecHitContainer multicontainer;
	while(current_cell){
		current_cell = define_used(current_cell, multi_id, &multicontainer);
		}
		multi_id++;	
	multiplets.push_back(multicontainer);
	multicontainer.clear();
	}
}

return;
}

CAcell* CAHitsGenerator::define_used(CAcell *cell, int id, TransientTrackingRecHit::RecHitContainer *multicontainer){
cell->IsUsed = id;
double pt_tmp = 99999.0;

if(cell->CAstatus==1){
	multicontainer->push_back(cell->hits[2]);
	multicontainer->push_back(cell->hits[1]);
	multicontainer->push_back(cell->hits[0]);
return NULL;	
}

else{
CAcell *next = new CAcell();
multicontainer->push_back(cell->hits[2]);
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



//Useful tools to store triplets in root files
void CAHitsGenerator::fittreegenerator(){

triptree->Branch("Event",&tmp_event,"Event/I");
//triptree->Branch("FitChiSquared",&tmp_FitChiSquared,"FitChiSquared/D");
//triptree->Branch("FitNDof",&tmp_FitNdof,"FitNdof/I");
//triptree->Branch("FitFoundHits",&tmp_FitFoundHits,"FitFoundHits/I");
triptree->Branch("LocalMomentum",&tmp_LocalMomentum,"LocalMomentum/D");

//triptree->Branch("hqbp",tmp_hqbp,"hqbp[3]/D");
//triptree->Branch("hdxdz",tmp_hdxdz,"hdxdz[3]/D");
//triptree->Branch("hdydz",tmp_hdydz,"hdydz[3]/D");

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


void CAHitsGenerator::fillfittree(CAcell *t){

//tmp_FitChiSquared = t->getFitChiSquared();
//tmp_FitNdof = t->getFitNdof();
//tmp_FitFoundHits = t->getFitFoundHits();
tmp_LocalMomentum = t->getLocalMomentum();


/*tmp_hqbp[0] = t->geth0qbp();
tmp_hqbp[1] = t->geth1qbp();
tmp_hqbp[2] = t->geth2qbp();

tmp_hdxdz[0] = t->geth0dxdz();
tmp_hdxdz[1] = t->geth1dxdz();
tmp_hdxdz[2] = t->geth2dxdz();

tmp_hdydz[0] = t->geth0dydz();
tmp_hdydz[1] = t->geth1dydz();
tmp_hdydz[2] = t->geth2dydz();*/

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


//Other functions

std::vector<CAcell *> CAHitsGenerator::ListIntersect(std::list<CAcell *> list_h1, std::list<CAcell *> list_h2){
    
    std::vector<CAcell *> temp;
    
    std::list<CAcell *>::const_iterator L1itr = list_h1.begin();
    std::list<CAcell *>::const_iterator L2itr = list_h2.begin();
    
    
    while(L1itr != list_h1.end() && L2itr != list_h2.end()){
        
        //check to see if they're equal, add to our temp list
        if(*L1itr == *L2itr){
            temp.push_back(*L1itr);
            L1itr++;
            L2itr++;
        }
        else if(*L1itr < *L2itr){
            L1itr++;
        }
        else{
            L2itr++;
        }
    }
    return temp;
}



unsigned int CAHitsGenerator::OmniRef(const TrackingRecHit* rhit){

    unsigned int ocref = 0;
    
    int subdetid = rhit->geographicalId().subdetId();
    if (subdetid==PixelSubdetector::PixelBarrel||subdetid==PixelSubdetector::PixelEndcap) {
        const SiPixelRecHit* pRHit = dynamic_cast<const SiPixelRecHit*>(rhit);
        if (!pRHit->cluster().isNonnull())
            edm::LogError("TrackAssociator") << ">>> RecHit does not have an associated cluster!" << " file: " << __FILE__ << " line: " << __LINE__;
        ocref = pRHit->omniClusterRef().rawIndex();
    }
    
    else if (subdetid==SiStripDetId::TIB||subdetid==SiStripDetId::TOB||subdetid==SiStripDetId::TID||subdetid==SiStripDetId::TEC) {
        const std::type_info &tid = typeid(*rhit);
        
        if (tid == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D* sMatchedRHit = dynamic_cast<const SiStripMatchedRecHit2D*>(rhit);
            if (!sMatchedRHit->monoHit().cluster().isNonnull() || !sMatchedRHit->stereoHit().cluster().isNonnull())
                edm::LogError("TrackAssociator") << ">>> RecHit does not have an associated cluster!" << " file: " << __FILE__ << " line: " << __LINE__;
            ocref = sMatchedRHit->monoClusterRef().rawIndex() + sMatchedRHit->stereoClusterRef().rawIndex();
        }
        else if (tid == typeid(SiStripRecHit2D)) {
            const SiStripRecHit2D* sRHit = dynamic_cast<const SiStripRecHit2D*>(rhit);
            if (!sRHit->cluster().isNonnull())
                edm::LogError("TrackAssociator") << ">>> RecHit does not have an associated cluster!" << " file: " << __FILE__ << " line: " << __LINE__;
            ocref = sRHit->omniClusterRef().rawIndex();
        }
        else if (tid == typeid(SiStripRecHit1D)) {
            const SiStripRecHit1D* sRHit = dynamic_cast<const SiStripRecHit1D*>(rhit);
            if (!sRHit->cluster().isNonnull())
                edm::LogError("TrackAssociator") << ">>> RecHit does not have an associated cluster!" << " file: " << __FILE__ << " line: " << __LINE__;
            ocref = sRHit->omniClusterRef().rawIndex();
        }
        
        else {
            edm::LogError("TrackAssociator") << ">>> getMatchedClusters: TrackingRecHit not associated to any SiStripCluster! subdetid = " << subdetid;
        }
        
    }
    else {
        edm::LogError("TrackAssociator") << ">>> getMatchedClusters: TrackingRecHit not associated to any cluster! subdetid = " << subdetid;
    }
    
    return ocref;
}


double CAHitsGenerator::DeltaEta(GlobalPoint p1, GlobalPoint p2){
    
    GlobalVector gV(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
    double eta = gV.eta();
    
    return eta;
}
