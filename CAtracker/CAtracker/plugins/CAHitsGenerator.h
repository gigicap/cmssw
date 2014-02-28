#ifndef CAHitsGenerator_h
#define CAHitsGenerator_h

// cmssw includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/ParameterSet/interface/ParameterDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"


#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGenerator.h"
#include "RecoTracker/TkSeedingLayers/interface/OrderedSeedingHits.h"
#include "RecoTracker/TkSeedGenerator/interface/MultiHitGenerator.h"
#include "RecoTracker/TkSeedingLayers/interface/OrderedMultiHits.h"

#include "RecoTracker/TkSeedGenerator/interface/SeedCreatorFactory.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducerFactory.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGeneratorFromPairAndLayersFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGeneratorFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"        
#include "RecoTracker/TkSeedGenerator/plugins/SeedFromConsecutiveHitsCreator.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSets.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPair.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGeneratorFromPairAndLayers.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGeneratorFromPairAndLayersFactory.h"
#include "RecoPixelVertexing/PixelTriplets/interface/LayerTriplets.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitTriplets.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/PixelTripletHLTGenerator.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducerFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

//root
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TVector3.h"

#include <unordered_map>


#include "DetIdClassification.h"
#include "CACell.h"

//
// CAHitsGenerator class declaration
//


class CAHitsGenerator : public MultiHitGenerator {
   public:
      CAHitsGenerator(const edm::ParameterSet&);
      ~CAHitsGenerator();
    
      void init(const edm::Event& ev, const edm::EventSetup& es);
      void hitSets(const TrackingRegion& , OrderedMultiHits& , const edm::Event& , const edm::EventSetup& );
	void hitSets(OrderedMultiHits& , const edm::Event& , const edm::EventSetup& );
      //virtual const OrderedSeedingHits & run(const TrackingRegion& , edm::Event& , const edm::EventSetup&);

    
        //typedef std::vector<TrackingRegion*> Regions;

    	void fitTripletSeeds(CAcell *trip, edm::ESHandle<MagneticField> & mf);
        void CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder>& builder, SeedingHitSet seed,  int seednumber);
        int ForwardCA();
		void BackwardCA();
		CAcell* define_used(CAcell *cell, int id, TransientTrackingRecHit::RecHitContainer* multicontainer);
    
        std::vector<CAcell *> ListIntersect(std::list<CAcell *>, std::list<CAcell *>);
        unsigned int OmniRef(const TrackingRecHit*);
        double DeltaEta(GlobalPoint, GlobalPoint);

    
        //once external functions

      //  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

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
    
    
        
        //void hitSets(OrderedSeedingHits , edm::Event& , const edm::EventSetup& );
    
        //Debug switch
        bool m_debug;
 	
        //CAcell collection
        //edm::InputTag m_inputTagTriplets;
        std::vector<CAcell> tripletCollection;
        std::vector<CAcell*> fittedTripletCollection;

        //TrajectorySeedCollection seedcollection;       //to be set
    
    
       // edm::ParameterSet m_regionProducerPSet;

    
        std::string m_builderName;
        std::string m_fitterName;

        edm::ESHandle<TrajectoryFitter> m_fitter;
        edm::ESHandle<TrackerGeometry> m_tracker;
    
    //triplet generator stuff
    //edm::ParameterSet aPset;
    edm::ParameterSet oPset;
    std::string oName;
    
    double dEta_cut;

	TTree *triptree;
	TTree *evtree;
	int eventcounter;
	std::vector<TransientTrackingRecHit::RecHitContainer> multiplets;
       // std::unique_ptr<TrackingRegionProducer> m_trackingRegionProducer;   //???

 
      	void fittreegenerator();
        void fillfittree(CAcell *t);
 
    
        //data structure for the neighborhood maps
        std::unordered_map<unsigned int, std::list<CAcell *>> hitUsage;

	//input collection
	OrderedHitTriplets  scoll;
	//const OrderedSeedingHits & scoll;
	//final seedinghitset collection
	OrderedMultiHits result;
    //OrderedSeedingHits result;

    //PixelTripletHLTGenerator *generator;
    typedef std::vector<HitTripletGeneratorFromPairAndLayers* > GeneratorContainer;
    GeneratorContainer        theGenerators;

    	//OrderedHitsGenerator *generator;
	//TrackingRegionProducer* Regions;
	edm::ESHandle<TrackerGeometry> tracker;
	edm::ESHandle<TransientTrackingRecHitBuilder> builder;
	edm::ESHandle<MagneticField> theMF;

	std::string layerBuilderName;
	edm::ESHandle<SeedingLayerSetsBuilder> layerBuilder;
	
	bool initialised;
	typedef LayerHitMapCache  LayerCacheType;
 	LayerCacheType            theLayerCache;
      // ----------member data ---------------------------
};




#endif
