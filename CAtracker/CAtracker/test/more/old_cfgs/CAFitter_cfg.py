import FWCore.ParameterSet.Config as cms

process = cms.Process("CATRACKER")
#process.load("Configuration/StandardSequences/GeometryPilot2_cff")
process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


# source
readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring()


#singleMuPt1 (25k)
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
#readFiles.extend( [
#                   '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/543A06E4-02E1-E211-8DFA-02163E008E6A.root'
#               ] );

#secFiles.extend( [
# '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/1E32E771-D2E0-E211-8E6C-003048CFAF86.root',
# '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/2EC1BA85-ECE0-E211-9A25-003048D37670.root',
# '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/7A356E90-D2E0-E211-A6FC-0025B32035DA.root'
#              ] );



#singleMuPt10 (25k)
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
#readFiles.extend( [
# '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/4E3340B2-02E1-E211-8BCA-003048F0E39C.root'
#] );

#secFiles.extend( [    
#'/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/06B5F324-E9E0-E211-A7E5-02163E008D7B.root',    
#'/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/4E98DBCA-E9E0-E211-8FDA-001D09F24353.root',    
#'/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/5A72C6C0-E1E0-E211-9EBE-003048D3750C.root'
           #] );


#singleMuPt100 (9k)
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( [
                '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt100/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/A2495409-00E1-E211-B3E9-00304896B904.root'
           ] );

    #secFiles.extend( [
    #          '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt100/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/1232ED18-E6E0-E211-8BB9-003048FE9B78.root'
#       ] );

#singleMuPt1000 (9k)
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
    #readFiles.extend( [
    # '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1000/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/FE430C9B-05E1-E211-A0A6-0025B3203BA4.root'
#          ] );

    #secFiles.extend( [
    #  '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/54B5D0ED-E5E0-E211-A571-001D09F2906A.root',
    # '/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/CC86B8A7-D2E0-E211-A4BB-00304896B910.root'
#    ] );


#TTbar
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
    #readFiles.extend( [
    # '/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/2A456B65-E6E0-E211-94FE-003048D374EC.root',
#'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/2CED6BEF-E8E0-E211-A0D8-003048F1DB62.root',
    #'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/E6E9A0E9-E4E0-E211-AD5C-003048FEAE6C.root'
# ] );

#secFiles.extend( [
#'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/08103FCF-D2E0-E211-AA58-003048F16F46.root',
#'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/1CFC80F6-D2E0-E211-901F-003048976D36.root',
    #'/store/relval/CMSSW_6_2_0_pre8/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/48D4C698-D2E0-E211-AB9F-003048F1C766.root'
# ] );



process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

#process.MessageLogger = cms.Service("MessageLogger",
# default = cms.untracked.PSet( limit = cms.untracked.int32(50) )
#)


process.load("Triplets.TripProducer.InitialStepGigiTrip_cff")

#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
#process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
#process.load("Validation.RecoTrack.cuts_cff")
#process.initialStepSeedClusterMask.trajectories = cms.InputTag("initialStepSeedsGigi")
#process.newCombinedSeeds.seedCollections[0] = cms.InputTag("initialStepSeedsGigi")

process.CAtracker = cms.EDProducer('CAtracker',SeedColl=cms.untracked.InputTag('initialStepSeedsGigi'))


process.evtInfo = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.initialStepSeedsGigi*process.CAtracker)

process.ep = cms.EndPath(process.evtInfo)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
