import FWCore.ParameterSet.Config as cms
#to change parameters
import sys

process = cms.Process("CASEEDER")

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


# source
readFiles = cms.untracked.vstring()

source = cms.Source ("PoolSource",fileNames = readFiles)

###TTbar
#readFiles.extend( [
#     '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/START70_V2-v4/00000/0E58E7AC-9A5D-E311-9C8B-0025905A48B2.root',
#     '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/START70_V2-v4/00000/5265FB10-9A5D-E311-938F-0025905A60F2.root',
#     '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/START70_V2-v4/00000/6CA502DF-975D-E311-9B99-0025905A60CA.root' ] );


    #readFiles.extend( [
    #   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/124EA7D9-6D5D-E311-A660-0025905A48F0.root',
    #  '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/5EFC5EDC-6D5D-E311-BBFB-0025905A60DE.root',
    #   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/681894EA-6D5D-E311-B58F-0025905A6080.root',
    #   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/7C894AE1-6D5D-E311-8F2B-0025905A6138.root',
    #   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/A0C9DDE4-6D5D-E311-BFEC-0025905A613C.root',
    #   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/C8C3EC2C-775D-E311-9C9D-0025905A60AA.root',
    #   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/CCD3F5DC-6D5D-E311-90AB-0025905A60A8.root',
    #   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/F42724E2-6D5D-E311-B8E9-0025905A60FE.root',
#   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V2-v4/00000/F8F07AE6-6D5D-E311-ACBE-0025905A60B4.root' ] );

###SingleMuPt100
#readFiles.extend( [
#'/store/relval/CMSSW_7_0_0_pre9/RelValSingleMuPt100/GEN-SIM-RECO/START70_V2-v4/00000/24AE5597-B25D-E311-8609-0025905A60D2.root' ] );

###TTbar+PU
readFiles.extend( [
                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/3A1952FE-5B5C-E311-95ED-BCAEC518FF7E.root',
                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/6CC96A39-5D5C-E311-BEC6-02163E009E6C.root',
                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/7AD3DB9A-5E5C-E311-A04E-02163E008D92.root',
                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/82B5BEFD-605C-E311-987E-02163E008EDD.root',
                   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/9C6545EE-655C-E311-AF91-003048FEB8F2.root',
                '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/E862A659-5F5C-E311-B330-02163E00A0F5.root' ] );


process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START70_V2::All'


# hit building

process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
process.load("RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff")

#from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *
#from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *

#triplets otf
process.load("RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi")

# seeding
process.load("CAtracker.CAtracker.GlobalSeedsFromMultiplets_cff")


#from CAtracker.CAtracker.GlobalSeedsFromMultiplets_cff import *
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock

#from RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi import *
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock

#process.initialStepSeedsGigi = CAtracker.CAtracker.GlobalSeedsFromMultiplets_cff.globalSeedsFromMultiplets.clone()


#from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
#process.initialStepSeedsGigi.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'

process.BarrelPentuplets = cms.ESProducer("SeedingLayersESProducer",
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



#process.CASeedingStep = CAtracker.CAtracker.GlobalSeedsFromMultiplets_cff.GlobalSeedsFromMultiplets.clone(
process.CASeedingStep = process.GlobalSeedsFromMultiplets.clone(
  RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
       ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
       RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
       ptMin = 0.6,
       originRadius = 0.02,
       nSigmaZ = 4.0
       )
       )
)

val_EtaCut=float(sys.argv[2])
process.CASeedingStep.OrderedHitsFactoryPSet.EtaCut=cms.double(val_EtaCut)

process.CASeedingStep.OrderedHitsFactoryPSet.SeedingLayers = 'BarrelPentuplets'

process.evtInfo = cms.OutputModule("AsciiOutputModule")

process.TFileService = cms.Service("TFileService", fileName = cms.string("histotrip.root") )

process.s = cms.Sequence(process.siPixelRecHits + process.siStripMatchedRecHits + process.CASeedingStep)
#process.s = cms.Sequence(process.siPixelRecHits + process.initialStepSeedsGigi)

process.p = cms.Path(process.s)

process.ep = cms.EndPath(process.evtInfo)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
