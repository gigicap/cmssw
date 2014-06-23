import FWCore.ParameterSet.Config as cms

process = cms.Process("RP")

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.EventContent.EventContent_cff")

process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")


### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'

# source
#readFiles = cms.untracked.vstring()


process.source = cms.Source("PoolSource",
    catalog = cms.untracked.string('PoolFileCatalog.xml'),
    fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_0_0_pre7/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v2/00000/BE70D2E4-2A46-E311-8E7C-0026189438D7.root')
)


#singleMuPt10 (25k)
#source = cms.Source ("PoolSource",fileNames = readFiles)
#readFiles.extend( [
#              '/store/relval/CMSSW_7_0_0_pre7/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v2/00000/BE70D2E4-2A46-E311-8E7C-0026189438D7.root'     
#              ] );


#process.source = source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))



process.RECO = cms.OutputModule("PoolOutputModule",
    process.RECOSIMEventContent,
    fileName = cms.untracked.string('recoTracks.root')
)

process.p1 = cms.Path(process.trackerlocalreco*process.offlineBeamSpot*process.ckftracks)
process.outpath = cms.EndPath(process.RECO)
