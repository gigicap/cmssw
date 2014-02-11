#ifndef TripGenerator_H
#define TripGenerator_H

/** A HitTripletGenerator from HitPairGenerator and vector of
    Layers. The HitPairGenerator provides a set of hit pairs.
    For each pair the search for compatible hit(s) is done among
    provided Layers
 */

#include "RecoTracker/TkHitPairs/interface/HitPairGenerator.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGenerator.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/CombinedHitTripletGenerator.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGeneratorFromPairAndLayers.h"

//Root tvector include for TIPfilter
//#include <ROOT.h>
#include "TVector3.h"
#include "TVector2.h"

#include <utility>
#include <vector>


class SeedComparitor;

class TripGenerator : public HitTripletGeneratorFromPairAndLayers {

typedef CombinedHitTripletGenerator::LayerCacheType       LayerCacheType;

public:
  TripGenerator( const edm::ParameterSet& cfg); 

  virtual ~TripGenerator();

  virtual void init( const HitPairGenerator & pairs,
      const std::vector<ctfseeding::SeedingLayer> & layers, LayerCacheType* layerCache);

  virtual void hitTriplets( const TrackingRegion& region, OrderedHitTriplets & trs, 
      const edm::Event & ev, const edm::EventSetup& es);

  const HitPairGenerator & pairGenerator() const { return *thePairGenerator; }
  const std::vector<ctfseeding::SeedingLayer> & thirdLayers() const { return theLayers; }
  
      
    float minzl[50];
    float maxzl[50];
    
    const float Bz = 3.8112; //Tesla
    const float minPt = 1.00; //GeV


private:
    bool TIPFilter(OrderedHitTriplet *hittriplet , float deltaTIP) const;
    bool PhiFilter(OrderedHitTriplet *hittriplet , float deltaPhi) const;
    bool ThetaFilter(OrderedHitTriplet *hittriplet , float deltaTheta) const;
    float minRadiusCurvature(float minPt, float Bz) const;
    float wraparound(float angle) const;


private:
  HitPairGenerator * thePairGenerator;
  std::vector<ctfseeding::SeedingLayer> theLayers;
  LayerCacheType * theLayerCache;

  bool useFixedPreFiltering;
  float extraHitRZtolerance;
  float extraHitRPhitolerance;
    std::vector<double> dphis;
    std::vector<double> dthetas;
    std::vector<double> dTIPs;
  SeedComparitor * theComparitor;

};
#endif

