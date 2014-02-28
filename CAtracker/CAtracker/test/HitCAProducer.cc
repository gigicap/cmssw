#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGeneratorFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGenerator.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducerFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

//#include "FWCore/Framework/interface/ESWatcher.h"
//#include "UserCode/konec/test/R2DTimerObserver.h"
#include "TH1D.h"
#include "TFile.h"


class HitCAProducer : public edm::EDAnalyzer {
public:
  explicit HitCAProducer(const edm::ParameterSet& conf);
  ~HitCAProducer();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
private:
  void init(const edm::EventSetup& es);
  edm::ParameterSet theConfig;
  OrderedHitsGenerator * theGenerator;
  TrackingRegionProducer* theRegionProducer;
  TH1D *hCPU, *hNum;
};

HitCAProducer::HitCAProducer(const edm::ParameterSet& conf)
  : theConfig(conf), theGenerator(0)
{
  edm::LogInfo("HitCAProducer")<<" CTOR";
  hCPU = new TH1D ("hCPU","hCPU",140,0.,0.070);
  hNum = new TH1D ("hNum","hNum",250,0.,500.);
}

HitCAProducer::~HitCAProducer()
{
  edm::LogInfo("HitCAProducer")<<" DTOR";
  delete theGenerator;

  TFile rootFile("analysis.root", "RECREATE", "my histograms");
  hCPU->Write();
  hNum->Write();
  rootFile.Close();
}

void HitCAProducer::init(const edm::EventSetup& es)
{

  std::cout << "INIT called" << std::endl;
  edm::ParameterSet orderedPSet =
  theConfig.getParameter<edm::ParameterSet>("OrderedHitsFactoryPSet");
  //theConfig.getParameter<edm::ParameterSet>("CAHitsGeneratorPSet");
  std::string orderedName = orderedPSet.getParameter<std::string>("ComponentName");
  theGenerator = OrderedHitsGeneratorFactory::get()->create( orderedName, orderedPSet);


  edm::ParameterSet regfactoryPSet =
      theConfig.getParameter<edm::ParameterSet>("RegionFactoryPSet");
  std::string regfactoryName = regfactoryPSet.getParameter<std::string>("ComponentName");
  //theRegionProducer = TrackingRegionProducerFactory::get()->create(regfactoryName,regfactoryPSet,consumesCollector());
  theRegionProducer = TrackingRegionProducerFactory::get()->create(regfactoryName,regfactoryPSet);
}

void HitCAProducer::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  static unsigned int lastRun=0;
  if (ev.id().run() != lastRun) { lastRun=ev.id().run(); init(es); }

// static edm::ESWatcher<TrackerDigiGeometryRecord> recordWatcher;
// if (recordWatcher.check(es)) init(es);

// GlobalTrackingRegion region;

  typedef std::vector<TrackingRegion* > Regions;
  Regions regions = theRegionProducer->regions(ev,es);
  const TrackingRegion & region = *regions[0];

// static R2DTimerObserver timer("**** MY TIMING REPORT ***");
// timer.start();
  edm::LogInfo("HitCAProducer") << "call CA producer! ";

//sostituire con il run...
  const OrderedSeedingHits & pentuplets = theGenerator->run(region,ev,es);

// timer.stop();
// hCPU->Fill( timer.lastMeasurement().real() );
  hNum->Fill(pentuplets.size());
  edm::LogInfo("HitCAProducer") << "size of multiplets: "<<pentuplets.size();

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HitCAProducer);
