import FWCore.ParameterSet.Config as cms

from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import *


# seeding layers (123 , 234 , 345)
BarrelPentuplets = cms.ESProducer("SeedingLayersESProducer",
                                  ComponentName = cms.string('BarrelPentuplets'),
                                  layerList = cms.vstring('BPix1+BPix2+BPix3','BPix2+BPix3+TIB1','BPix3+TIB1+TIB2'),
                                  BPix = cms.PSet(
                                                  useErrorsFromParam = cms.bool(True),
                                                  hitErrorRPhi = cms.double(0.0027),
                                                  TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
                                                  HitProducer = cms.string('siPixelRecHits'),
                                                  hitErrorRZ = cms.double(0.006)
                                                  ),
                                  TIB = cms.PSet(
                                                 TTRHBuilder = cms.string('WithTrackAngle'),
                                                 matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
                                                 #skipClusters = cms.InputTag('pixelLessStepClusters')
                                                 )
                                  )




CAHitsGenerator = cms.PSet (
   debug = cms.untracked.bool(True),
   Builder = cms.untracked.string('WithAngleAndTemplate'),
   EtaCut = cms.double(0.0256),									#Can be changed in the _cfg
                            #   ComponentName = cms.string('CAHitsGenerator'),
    #OrderedHitsFactoryPSet = cms.PSet(
    TripletFactoryForCAPSet = cms.PSet(
        #ComponentName = cms.string('StandardHitTripletGenerator'),
                                       #ComponentName = cms.string('PixelTripletHLTGenerator'),
                                       #SeedingLayers = cms.string('BarrelPentuplets'),
                                       #GeneratorPSet = cms.PSet(PixelTripletHLTGenerator)
        PixelTripletHLTGenerator
   )
)

#RegionFactoryPSet = cms.PSet(
#RegionPsetFomBeamSpotBlock,
#omponentName = cms.string('GlobalRegionProducerFromBeamSpot')
#)
 