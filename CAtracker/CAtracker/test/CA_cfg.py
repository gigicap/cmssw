import FWCore.ParameterSet.Config as cms
process = cms.Process("CATest")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

# source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
'/store/relval/CMSSW_6_2_0_pre8/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/4E3340B2-02E1-E211-8BCA-003048F0E39C.root'
))


process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet( threshold = cms.untracked.string('INFO'))
)

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.TrackerSimData.trackerSimGeometryXML_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
#process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")

#from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *
#process.siPixelClusters.src = cms.InputTag('simSiPixelDigis')


from RecoPixelVertexing.PixelTrackFitting.PixelTracks_cfi import *
from RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi import *

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'
# hit building
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
process.load("RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff")

#triplets otf
from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *



# seeding
#process.load("CAtracker.CAtracker.CAHitsGenerator_cfi")
from CAtracker.CAtracker.CAHitsGenerator_cfi import CAHitsGenerator as CellularAutomaton
#CellularAutomaton.debug = cms.untracked.bool(True)
#CellularAutomaton.Builder = cms.untracked.string('WithAngleAndTemplate')
#CellularAutomaton.EtaCut = cms.double(0.0256)
#CellularAutomaton.OrderedHitsFactoryPSet = cms.PSet(PixelTripletHLTGenerator)

#process.load("CAtracker.CAtracker.CAHitsGenerator_cfi")
#process.CellularAutomaton = process.CAHitsGenerator.clone(
#	debug = cms.untracked.bool(True),
#	Builder = cms.untracked.string('WithAngleAndTemplate'),
#	EtaCut = cms.double(0.0256),
#	OrderedHitsFactoryPSet = cms.PSet(PixelTripletHLTGenerator)	
#)

#print CellularAutomaton.EtaCut

# seeding layers (123 , 234 , 345)
process.BarrelPentuplets = cms.ESProducer("SeedingLayersESProducer",
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



process.CApentuplets = cms.EDAnalyzer("HitCAProducer",
  #OrderedHitsFactoryPSet = cms.PSet(
  #  ComponentName = cms.string("CAHitsGenerator"),
  #  SeedingLayers = cms.string("BarrelPentuplets"),
  #  GeneratorPSet = cms.PSet(CellularAutomaton)
  #),
  OrderedHitsFactoryPSet = cms.PSet(
    ComponentName = cms.string("CAHitsGenerator"),
    debug = cms.untracked.bool(True),
    Builder = cms.untracked.string('WithAngleAndTemplate'),
    EtaCut = cms.double(0.0256),
    #ComponentName = cms.string('CAHitsGenerator'),
    GeneratorPSet = cms.PSet(PixelTripletHLTGenerator),
    SeedingLayers = cms.string("BarrelPentuplets")
  ),
    RegionFactoryPSet = cms.PSet(
        RegionPSetBlock,
        ComponentName = cms.string('GlobalRegionProducer')
  )
)

#process.CApentuplets.OrderedHitsFactoryPSet.GeneratorPSet.debug = cms.untracked.bool(True)
#process.CApentuplets.OrderedHitsFactoryPSet.GeneratorPSet.Builder = cms.untracked.string('WithAngleAndTemplate')
#process.CApentuplets.OrderedHitsFactoryPSet.GeneratorPSet.EtaCut = cms.double(0.0256)
#process.CApentuplets.OrderedHitsFactoryPSet.GeneratorPSet.OrderedHitsFactoryPSet = cms.PSet(PixelTripletHLTGenerator)

process.p = cms.Path(process.siPixelRecHits + process.siStripMatchedRecHits + process.CApentuplets)
#process.p = cms.Path(process.CApentuplets)
