import FWCore.ParameterSet.Config as cms
#to change parameters
import sys

process = cms.Process("CATRACKER")
#process.load("Configuration/StandardSequences/GeometryPilot2_cff")
process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


# source
readFiles = cms.untracked.vstring()

source = cms.Source ("PoolSource",fileNames = readFiles)


#singleMuPt1 (25k)
readFiles.extend( [
                   '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/543A06E4-02E1-E211-8DFA-02163E008E6A.root'
               ] );

#singleMuPt10 (25k)
#readFiles.extend( [
# '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/4E3340B2-02E1-E211-8BCA-003048F0E39C.root'
#] );

#singleMuPt100 (9k)
#readFiles.extend( [
#                '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt100/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/A2495409-00E1-E211-B3E9-00304896B904.root'
#           ] );

#singleMuPt1000 (9k)
    #readFiles.extend( [
    #               '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1000/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/FE430C9B-05E1-E211-A0A6-0025B3203BA4.root'
#               ] );

#TTbar
    #readFiles.extend( [
    #'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/2A456B65-E6E0-E211-94FE-003048D374EC.root',
#'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/2CED6BEF-E8E0-E211-A0D8-003048F1DB62.root',
    #'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/E6E9A0E9-E4E0-E211-AD5C-003048FEAE6C.root'
#] );

#TTbar-PU
#readFiles.extend( [
    #               '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/0C80F82F-DCE4-E211-976F-002481E94B9C.root',
 #                  '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/22444AF4-D7E4-E211-9C7C-003048CFC64E.root',
  #                 '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/6085F0D4-D6E4-E211-B3F1-001E67398A34.root',
  #                 '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/6C160D24-D3E4-E211-B8BF-00215AEDFCCC.root',
  #                 '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/782B3885-E4E4-E211-BF39-003048F23D78.root',
  #                 '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/A469A789-D4E4-E211-9DBC-001E6739823C.root',
  #             '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/A8A55915-E3E4-E211-9A4F-003048FF3AE2.root'
#] );

#readFiles.extend(['root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v1/00000/0C80F82F-DCE4-E211-976F-002481E94B9C.root'])

process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'


# hit building
from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *
from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *



#from Triplets.TripProducer.TripGenerator_cfi import *

#seeding layers (Pixel only)
#initialStepSeedLayers = cms.ESProducer("SeedingLayersESProducer",
#    ComponentName = cms.string('initialStepSeedLayers'),
#    layerList = cms.vstring('BPix1+BPix2+BPix3'),
#    BPix = cms.PSet(
#           useErrorsFromParam = cms.bool(True),
#           hitErrorRPhi = cms.double(0.0027),
#           TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
#           HitProducer = cms.string('siPixelRecHits'),
#           hitErrorRZ = cms.double(0.006)
#           )
#)

# seeding layers (123 , 234 , 345) 
process.initialStepSeedLayers2 = cms.ESProducer("SeedingLayersESProducer",
    ComponentName = cms.string('initialStepSeedLayers2'),
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


# seeding
from RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff import *
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock

process.initialStepSeedsGigi = RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff.globalSeedsFromTriplets.clone(
       RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
       ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
       RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
       ptMin = 0.6,
       originRadius = 0.02,
       nSigmaZ = 4.0
       )
       )
       )

from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
process.initialStepSeedsGigi.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'

process.initialStepSeedsGigi.OrderedHitsFactoryPSet.SeedingLayers = 'initialStepSeedLayers2'
#process.initialStepSeedsGigi.OrderedHitsFactoryPSet.SeedingLayers = 'MixedLayerTriplets'
#process.initialStepSeedsGigi.OrderedHitsFactoryPSet.SeedingLayers = 'PixelLayerTriplets'

#SWITCHER
#process.initialStepSeedsGigi.OrderedHitsFactoryPSet.GeneratorPSet = cms.PSet(TripGenerator)
process.initialStepSeedsGigi.OrderedHitsFactoryPSet.GeneratorPSet = cms.PSet(PixelTripletHLTGenerator)

from CAtracker.CAtracker.CAtracker_cfi import *

val_EtaCut=float(sys.argv[2])
process.CAtracker = cms.EDProducer('CAtracker',SeedColl=cms.untracked.InputTag('initialStepSeedsGigi'),EtaCut=cms.double(val_EtaCut))


process.evtInfo = cms.OutputModule("AsciiOutputModule")

process.TFileService = cms.Service("TFileService", fileName = cms.string("histotrip.root") )

process.s = cms.Sequence(process.siPixelRecHits + process.siStripMatchedRecHits + process.initialStepSeedsGigi + process.CAtracker)
#process.s = cms.Sequence(process.siPixelRecHits + process.initialStepSeedsGigi)

process.p = cms.Path(process.s)

process.ep = cms.EndPath(process.evtInfo)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
