import FWCore.ParameterSet.Config as cms
### Build seeds from CA pentuplets (Gigi's Code)### 
### STEP 0 ###
 
# hit building
from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *
from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *

#import triplet producer an cahitsgenerator
from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
from CAtracker.CAtracker.CAHitsGenerator_cfi import CAHitsGenerator as CellularAutomaton


#import Multiplet seedinglayer (now just defined here)
#BarrelPentuplets = cms.ESProducer("SeedingLayersESProducer",
#    ComponentName = cms.string('BarrelPentuplets'),
#    layerList = cms.vstring('BPix1+BPix2+BPix3+TIB1+TIB2'),
#    BPix = cms.PSet(
#           useErrorsFromParam = cms.bool(True),
#           hitErrorRPhi = cms.double(0.0027),
#           TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
#           HitProducer = cms.string('siPixelRecHits'),
#          hitErrorRZ = cms.double(0.006)
#        ),
#    TIB = cms.PSet(
#        TTRHBuilder = cms.string('WithTrackAngle'),
          #    matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
#      )
#)

BarrelPentuplets = cms.ESProducer("SeedingLayersESProducer",
    ComponentName = cms.string('BarrelPentuplets'),
    layerList = cms.vstring('BPix1+BPix2+BPix3','BPix2+BPix3+TIB1','BPix3+TIB1+TIB2'),
    #layerList = cms.vstring('BPix2+BPix3+TIB1'),
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
      )
)


from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock

#globalSeedsFromMultiplets
import RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi
initialStepSeedsCA = RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi.seedGeneratorFromRegionHitsEDProducer.clone(
OrderedHitsFactoryPSet = cms.PSet(
          CellularAutomaton,
          ComponentName = cms.string('CAHitsGenerator'), #right?
          SeedingLayers = cms.string('BarrelPentuplets')
                                  #GeneratorPSet = cms.PSet(CellularAutomaton)                            
     ),
RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
       ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
       RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
       		ptMin = 0.6,
       		originRadius = 0.02,
       		nSigmaZ = 4.0
       		)
       )#,
#SeedCreatorPSet = cms.PSet( 
#       ComponentName = cms.string( "SeedFromConsecutiveHitsCreator" ),
#       propagator = cms.string( "PropagatorWithMaterial" )
#     )
)


from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
#initialStepSeedsCA.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'

#initialStepSeedsCA.OrderedHitsFactoryPSet.SeedingLayers = 'BarrelPentuplets'


# building
import TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi
initialStepTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.clone(
     ComponentName = 'initialStepTrajectoryFilter',
     filterPset = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.filterPset.clone(
     minimumNumberOfHits = 5,  #???
     minPt = 0.2
     )
     )
 
import TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi
initialStepChi2Est = TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi.Chi2MeasurementEstimator.clone(
     ComponentName = cms.string('initialStepChi2Est'),
     nSigma = cms.double(3.0),
     MaxChi2 = cms.double(30.0)
 )
 
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi
initialStepTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi.GroupedCkfTrajectoryBuilder.clone(
     ComponentName = 'initialStepTrajectoryBuilder',
     trajectoryFilterName = 'initialStepTrajectoryFilter',
     alwaysUseInvalidHits = True,
     maxCand = 6,
     estimator = cms.string('initialStepChi2Est'),
     maxDPhiForLooperReconstruction = cms.double(2.0),
     maxPtForLooperReconstruction = cms.double(0.7) 
     )

import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
initialStepTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
     src = cms.InputTag('initialStepSeedsCA'),
     ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
     numHitsForSeedCleaner = cms.int32(50),
     #onlyPixelHitsForSeedCleaner = cms.bool(True),
 	onlyPixelHitsForSeedCleaner = cms.bool(False),
     TrajectoryBuilder = 'initialStepTrajectoryBuilder',
     doSeedingRegionRebuilding = True,
     useHitsSplitting = True
  )
 
 # fitting

import RecoTracker.TrackProducer.TrackProducer_cfi
initialStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
     src = 'initialStepTrackCandidates',
     AlgorithmName = cms.string('iter0'),
     Fitter = cms.string('FlexibleKFFittingSmoother')
     )
 
 # Final selection
import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
initialStepSelector = RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.multiTrackSelector.clone(
     src='initialStepTracks',
     useAnyMVA = cms.bool(True),
     GBRForestLabel = cms.string('MVASelectorIter0'),
     trackSelectors= cms.VPSet(
         RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
             name = 'initialStepLoose',
             ), #end of pset
         RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.tightMTS.clone(
             name = 'initialStepTight',
             preFilterName = 'initialStepLoose',
             ),
         RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.highpurityMTS.clone(
             name = 'initialStep',
             preFilterName = 'initialStepTight',
             ),
         ) #end of vpset
     ) #end of clone
 
 # Final sequence
InitialStepCAPentuplets = cms.Sequence(initialStepSeedsCA*
                            initialStepTrackCandidates*
                            initialStepTracks*
                            initialStepSelector)
