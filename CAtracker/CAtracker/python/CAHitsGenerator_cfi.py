import FWCore.ParameterSet.Config as cms

from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import *

CAHitsGenerator = cms.PSet (
   debug = cms.bool(True),
   Builder = cms.untracked.string('WithAngleAndTemplate'),
   EtaCut = cms.double(0.0256),									#Can be changed in the _cfg
   ComponentName = cms.string('CAHitsGenerator'),
   OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitTripletGenerator'),
        SeedingLayers = cms.string('PixelLayerTriplets'),
        GeneratorPSet = cms.PSet(PixelTripletHLTGenerator)
   )

)

#RegionFactoryPSet = cms.PSet(
#RegionPsetFomBeamSpotBlock,
#omponentName = cms.string('GlobalRegionProducerFromBeamSpot')
#)
