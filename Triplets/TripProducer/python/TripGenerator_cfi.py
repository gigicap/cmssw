import FWCore.ParameterSet.Config as cms
 
TripGenerator = cms.PSet(
    maxElement = cms.uint32(100000),
     ComponentName = cms.string('TripGenerator'),
                         
	#cut vectors: array index == seeding layer index 
	# 0 = BPIX1 1 = BPIX2 (not used) 2 = BPIX2 (1_2_3 triplets) 3 = TIB1 (2_3_4 triplets)		
                         phicuts = cms.vdouble(0, 0, 0.075, 0.2),
                         thetacuts = cms.vdouble(0, 0 ,0.01, 0.2),
                         TIPcuts = cms.vdouble(0, 0, 0.5, 0.5),
	

     SeedComparitorPSet = cms.PSet(
      ComponentName = cms.string('none')
      )
 )
