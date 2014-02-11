import FWCore.ParameterSet.Config as cms

# Iterative steps
from RecoTracker.IterativeTracking.iterativeTk_cff import *
from RecoTracker.IterativeTracking.ElectronSeeds_cff import *
# My iterative step
from Triplets.TripProducer.InitialStepGigiTrip_cff import *

#import copy

#dEdX reconstruction
#from RecoTracker.DeDx.dedxEstimators_cff import *

#BeamHalo tracking
#from RecoTracker.Configuration.RecoTrackerBHM_cff import *


#special sequences, such as pixel-less
#from RecoTracker.Configuration.RecoTrackerNotStandard_cff import *

#InitialStepGigiTrip
MyiterTracking = cms.Sequence(InitialStep*
                            LowPtTripletStep*
                            PixelPairStep*
                            DetachedTripletStep*
                            MixedTripletStep*
                            PixelLessStep*
                            TobTecStep*
                            earlyGeneralTracks*
                            muonSeededStep*
                            preDuplicateMergingGeneralTracks*
                            generalTracksSequence*
                            ConvStep*
                            conversionStepTracks
                            )

Myckftracks = cms.Sequence(MyiterTracking*electronSeedsSeq)


#ckftracks_woBH = cms.Sequence(MyiterTracking*electronSeedsSeq*doAlldEdXEstimators)
#Myckftracks = ckftracks_woBH.copy() #+ beamhaloTracksSeq) # temporarily out, takes too much resources

#ckftracks_wodEdX = Myckftracks.copy()
#ckftracks_wodEdX.remove(doAlldEdXEstimators)


#Myckftracks_plus_pixelless = cms.Sequence(Myckftracks*ctfTracksPixelLess)


#from RecoJets.JetAssociationProducers.trackExtrapolator_cfi import *
#trackingGlobalReco = cms.Sequence(Myckftracks*trackExtrapolator)