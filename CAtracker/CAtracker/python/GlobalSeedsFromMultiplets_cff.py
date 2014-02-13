import FWCore.ParameterSet.Config as cms

#standard imports (from GlobalSeedsFromTriplets_cff.py)
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
from RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi import *
from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *

#import triplet producer an cahitsgenerator
from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
from CAtracker.CAtracker.CAHitsGenerator_cfi import CAHitsGenerator as CellularAutomaton

#import Multiplet seedinglayer (now just defined here)
BarrelPentuplets = cms.ESProducer("SeedingLayersESProducer",
    ComponentName = cms.string('BarrelPentuplets'),
    layerList = cms.vstring('BPix1+BPix2+BPix3+TIB1+TIB2'),
    BPix = cms.PSet(
           useErrorsFromParam = cms.bool(True),
           hitErrorRPhi = cms.double(0.0027),
           TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
           HitProducer = cms.string('siPixelRecHits'),
          hitErrorRZ = cms.double(0.006)
        ),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
            #    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
      )
)


#globalSeedsFromMultiplets
import RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi
GlobalSeedsFromMultiplets = RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi.seedGeneratorFromRegionHitsEDProducer.clone(
OrderedHitsFactoryPSet = cms.PSet(
          CellularAutomaton,
          ComponentName = cms.string('CAHitsGenerator'), #right?
          SeedingLayers = cms.string('BarrelPentuplets')
                                  #GeneratorPSet = cms.PSet(CellularAutomaton)
                            
     )
)
