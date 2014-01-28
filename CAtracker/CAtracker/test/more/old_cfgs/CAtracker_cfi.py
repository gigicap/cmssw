import FWCore.ParameterSet.Config as cms
 
CAtracker = cms.PSet(
    maxElement = cms.uint32(100000),
     ComponentName = cms.string('CAtracker'),
                         
	#cut on hit etas
	#Layer dependent?		
    EtaCut = cms.double(0),

    #weight to be assigned for backward CA
    #Layer dependent?
    alpha_qovp = cms.double(0),
    alpha_kink = cms.double(0),
    alpha_pt = cms.double(0)


 )
