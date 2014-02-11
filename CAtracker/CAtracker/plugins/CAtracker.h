#ifndef CAtracker_h
#define CAtracker_h

// system include files
#include <memory>
#include <vector>
#include <map>
#include <deque>
#include <queue>
#include <iostream>
#include <set>
#include <limits>
#include <array>
#include <algorithm>

#include <math.h>

// cmssw
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterDescription.h"
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


#include "RecoTracker/TkSeedGenerator/interface/SeedCreatorFactory.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducerFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
#include "RecoTracker/TkSeedGenerator/plugins/SeedFromConsecutiveHitsCreator.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

#include "DetIdClassification.h"


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

#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

//root
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TVector3.h"


//Class CAcell definition
class CAcell{
    
public:
	//CAcell hits (from seeds)
    TransientTrackingRecHit::RecHitContainer hits;
	
	//fit parameters
	 LocalTrajectoryParameters FitHit0Local;
	 LocalTrajectoryError FitHit0LocalError;
    
	 LocalTrajectoryParameters FitHit1Local;
	 LocalTrajectoryError FitHit1LocalError;
    
	 LocalTrajectoryParameters FitHit2Local;
	 LocalTrajectoryError FitHit2LocalError;
    
	bool FitSuccessful;
	double FitChiSquared;
	int FitNdof;
	int FitFoundHits;
	double LocalMomentum;

	int eventNumber;
    
    //triplet identifier e.g. (1,2,3) =  123
    int tripletIdentifier;

	int hitId0;
	int hitId1;
	int hitId2;
    
    double dEta;
    
    unsigned int rawId0;
    unsigned int rawId1;
    unsigned int rawId2;
    
    //neighborMAP
    std::vector<CAcell *> left_neighbors;       //list of the cells sharing hits 0-1 with the cell
    std::vector<CAcell *> right_neighbors;      //list of the cells sharing hits 1-2 with the cell

    //How many neighbors does the cell have?
    int IsNeighborly;
    
    
    //just for convenience and debug purposes
    int TripletNumber;

    //for the backwardCA
    int IsUsed;
    
    //for fitting
    PTrajectoryStateOnDet seedStartState;
    
     CAcell();
    
double getFitChiSquared(){ return FitChiSquared;}
int getFitNdof(){return FitNdof;}
int getFitFoundHits(){return FitFoundHits;}
double getLocalMomentum(){ return LocalMomentum;}

double geth0qbp(){ return FitHit0Local.qbp();}
double geth0dxdz(){ return FitHit0Local.dxdz();}
double geth0dydz(){ return FitHit0Local.dydz();}

double geth1qbp(){ return FitHit0Local.qbp();}
double geth1dxdz(){ return FitHit0Local.dxdz();}
double geth1dydz(){ return FitHit0Local.dydz();}

double geth2qbp(){ return FitHit2Local.qbp();}
double geth2dxdz(){ return FitHit2Local.dxdz();}
double geth2dydz(){ return FitHit2Local.dydz();}

int geteventNumber(){return eventNumber;}
    
int gethitId0(){return hitId0;}
int gethitId1(){return hitId1;}
int gethitId2(){return hitId2;}

    
//CAstuff
    int CAstatus;
    
    
};


CAcell::CAcell(){
    
    TripletNumber = 0;
    CAstatus = 0;
    tripletIdentifier = 0;
    
}


//Other functions

std::vector<CAcell *> ListIntersect(std::list<CAcell *> list_h1, std::list<CAcell *> list_h2){
    
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



unsigned int OmniRef(const TrackingRecHit* rhit){
	
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


double DeltaEta(GlobalPoint p1, GlobalPoint p2){

    GlobalVector gV(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
    double eta = gV.eta();
    
    return eta;
}


#endif
