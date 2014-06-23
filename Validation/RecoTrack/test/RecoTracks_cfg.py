import FWCore.ParameterSet.Config as cms
process = cms.Process("RP")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("Configuration.EventContent.EventContent_cff")





# source
readFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles)
###TTbar
#readFiles.extend( [
#     '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/START70_V2-v4/00000/0E58E7AC-9A5D-E311-9C8B-0025905A48B2.root',
#                       '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/START70_V2-v4/00000/5265FB10-9A5D-E311-938F-0025905A60F2.root',
#               '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/START70_V2-v4/00000/6CA502DF-975D-E311-9B99-0025905A60CA.root' ] );


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
readFiles.extend( [
 '/store/relval/CMSSW_7_0_0_pre9/RelValSingleMuPt100/GEN-SIM-RECO/START70_V2-v4/00000/24AE5597-B25D-E311-8609-0025905A60D2.root' ] );

###TTbar+PU
#readFiles.extend( [
#                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/3A1952FE-5B5C-E311-95ED-BCAEC518FF7E.root',
#                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/6CC96A39-5D5C-E311-BEC6-02163E009E6C.root',
#                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/7AD3DB9A-5E5C-E311-A04E-02163E008D92.root',
#                    '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/82B5BEFD-605C-E311-987E-02163E008EDD.root',
#                   '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/9C6545EE-655C-E311-AF91-003048FEB8F2.root',
#                '/store/relval/CMSSW_7_0_0_pre9/RelValTTbar/GEN-SIM-RECO/PU_START70_V2-v4/00000/E862A659-5F5C-E311-B330-02163E00A0F5.root' ] );


process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")
#process.load("RecoTracker.IterativeTracking.InitialStep_cff")


process.source = source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


process.GlobalTag.globaltag = 'START70_V2::All'
#process.GlobalTag.globaltag = 'GR_H_V33::All'
#process.RECO = cms.OutputModule("PoolOutputModule",
#process.RECOSIMEventContent,
#fileName = cms.untracked.string('recoTracks.root')
#)

process.evtInfo = cms.OutputModule("AsciiOutputModule")


process.demo = process.demo = cms.EDAnalyzer('CATkAnalyzer', tracks = cms.untracked.InputTag("initialStepTracks"))

process.p1 = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits*process.InitialStepCAPentuplets*process.demo
                  )

    #process.p1 = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits*process.InitialStep*process.demo
#              )


process.ep = cms.EndPath(process.evtInfo)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
