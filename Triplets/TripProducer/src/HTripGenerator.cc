#include "../interface/HTripGenerator.h"
 
HTripGenerator::HTripGenerator(unsigned int nSize)
 {
   theTriplets.reserve(nSize);
 }
 
const OrderedHitTriplets & HTripGenerator::run(
     const TrackingRegion& region, const edm::Event & ev, const edm::EventSetup& es)
 {
   theTriplets.clear();
   hitTriplets(region, theTriplets, ev, es);
   return theTriplets;
 }
 
void HTripGenerator::clear() 
 {
   theTriplets.clear();
 } 