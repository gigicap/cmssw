#include "TripGenerator.h"

#include <cmath>
#include "RecoPixelVertexing/PixelTriplets/plugins/ThirdHitPredictionFromInvParabola.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/ThirdHitRZPrediction.h"
#include "RecoTracker/TkMSParametrization/interface/PixelRecoUtilities.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/ThirdHitCorrection.h"
#include "RecoTracker/TkHitPairs/interface/RecHitsSortedInPhi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>

#include "RecoTracker/TkSeedingLayers/interface/SeedComparitorFactory.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedComparitor.h"

#include "DataFormats/GeometryVector/interface/Pi.h"
//#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerAlgo.h"
//#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerTools.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerAlgo.h" //amend to point at your copy...
#include "RecoPixelVertexing/PixelTriplets/plugins/KDTreeLinkerTools.h"

#include<cstdio>

using pixelrecoutilities::LongitudinalBendingCorrection;
typedef PixelRecoRange<float> Range;

using namespace std;
using namespace ctfseeding;

TripGenerator:: TripGenerator(const edm::ParameterSet& cfg)
  : thePairGenerator(0),
    theLayerCache(0)
    //extraHitRZtolerance(cfg.getParameter<double>("extraHitRZtolerance")),
    //extraHitRPhitolerance(cfg.getParameter<double>("extraHitRPhitolerance"))
{
	theMaxElement=cfg.getParameter<unsigned int>("maxElement");
    
    //theta, phi, tip cuts from cfg
    dphis = cfg.getParameter<std::vector<double>>("phicuts");
    dthetas = cfg.getParameter<std::vector<double>>("thetacuts");
    dTIPs = cfg.getParameter<std::vector<double>>("TIPcuts");
        
    edm::ParameterSet comparitorPSet =
    cfg.getParameter<edm::ParameterSet>("SeedComparitorPSet");
    std::string comparitorName = comparitorPSet.getParameter<std::string>("ComponentName");
    theComparitor = (comparitorName == "none") ?
    0 :  SeedComparitorFactory::get()->create( comparitorName, comparitorPSet);
}

TripGenerator::~TripGenerator() { 
  delete thePairGenerator;
  delete theComparitor;
}

void TripGenerator::init( const HitPairGenerator & pairs,
				     const std::vector<SeedingLayer> & layers,
				     LayerCacheType* layerCache)
{
  thePairGenerator = pairs.clone();
  theLayers = layers;
  theLayerCache = layerCache;
}

void TripGenerator::hitTriplets(const TrackingRegion& region, 
					   OrderedHitTriplets & result,
					   const edm::Event & ev,
					   const edm::EventSetup& es)
{
  //R_min (for phi calculation. pt_min = 1GeV/c)
    float R_MIN=minRadiusCurvature(minPt, Bz);

  if (theComparitor) theComparitor->init(es);
  
  auto const & doublets = thePairGenerator->doublets(region,ev,es);
  
  if (doublets.empty()) return;

  //auto outSeq =  doublets.detLayer(HitDoublets::outer)->seqNum();


  // std::cout << "pairs " << doublets.size() << std::endl;
  
  float regOffset = region.origin().perp(); //try to take account of non-centrality (?)
  int size = theLayers.size();
  
   
  const RecHitsSortedInPhi * thirdHitMap[size];
  typedef RecHitsSortedInPhi::Hit Hit;

  using NodeInfo = KDTreeNodeInfo<unsigned int>;
  std::vector<NodeInfo > layerTree; // re-used throughout

  KDTreeLinkerAlgo<unsigned int> hitTree[size];
 // float rzError[size]; //save maximum errors
  float maxphi = Geom::ftwoPi(), minphi = -maxphi; // increase to cater for any range
  
  // fill the prediction vector
  for (int il=0; il!=size; ++il) {
    thirdHitMap[il] = &(*theLayerCache)(&theLayers[il], region, ev, es);
    auto const & hits = *thirdHitMap[il];
    //pred.initLayer(theLayers[il].detLayer());
    //pred.initTolerance(extraHitRZtolerance);
    
    layerTree.clear();
    float minv=999999.0, maxv= -999999.0; // Initialise to extreme values in case no hits
    //float maxErr=0.0f;
    for (unsigned int i=0; i!=hits.size(); ++i) {
      auto angle = hits.phi(i);
      auto v =  hits.gv(i);
      //use (phi,r) for endcaps rather than (phi,z)
      minv = std::min(minv,v);  maxv = std::max(maxv,v);
      //float myerr = hits.dv[i];
      //maxErr = std::max(maxErr,myerr);
      layerTree.emplace_back(i, angle, v); // save it
      if (angle < 0)  // wrap all points in phi
	{ layerTree.emplace_back(i, angle+Geom::ftwoPi(), v);}
      else
	{ layerTree.emplace_back(i, angle-Geom::ftwoPi(), v);}
    }
    KDTreeBox phiZ(minphi, maxphi, minv-0.01f, maxv+0.01f);  // declare our bounds
    //add fudge factors in case only one hit and also for floating-point inaccuracy
    hitTree[il].build(layerTree, phiZ); // make KDtree
	//rzError[il] = maxErr; //save error
    // std::cout << "layer " << theLayers[il].detLayer()->seqNum() << " " << layerTree.size() << std::endl; 
  }
  
   
  for (std::size_t ip =0;  ip!=doublets.size(); ip++) {
    auto xi = doublets.x(ip,HitDoublets::inner);
    auto yi = doublets.y(ip,HitDoublets::inner);
    auto zi = doublets.z(ip,HitDoublets::inner);
    //auto rvi = doublets.rv(ip,HitDoublets::inner);
    auto xo = doublets.x(ip,HitDoublets::outer);
    auto yo = doublets.y(ip,HitDoublets::outer);
    auto zo = doublets.z(ip,HitDoublets::outer);
    auto rvo = doublets.rv(ip,HitDoublets::outer);
    
        //extra based
        //float d02 = sqrt(((yi*(xo-xi)-xi*(yo-yi))*(yi*(xo-xi)-xi*(yo-yi)))/((xo-xi)*(xo-xi)+(yo-yi)*(yo-yi)));
        //float r0 = sqrt(xi*xi+yi*yi-d02*d02);
        //float cotTheta = (zo-zi)/sqrt(xo*xo+yo*yo-r0*r0);
        
                
        float phi2 = atan2(yo,xo);
        float phi1 = atan2(yi,xi);
        //Phi of the 0-1 line
        
        float dphi = fabs(phi2-phi1);
        
        float dy = yo-yi;
        float dx = xo-xi;
        float dz = zo-zi;
      
        float d = sqrt(dx*dx+dy*dy);
      
        //For angle-based theta prediction.
        int s = (yo>0) - (yo<0);
        float r = s*d;
        float Theta = atan2(r,dz);
        float Thetacorr = fabs(atan(r/dz)); 

    // std::cout << ip << ": " << point1.r() << ","<< point1.z() << " " 
    //                        << point2.r() << ","<< point2.z() <<std::endl;

    for (int il=0; il!=size; ++il) {
      if (hitTree[il].empty()) continue; // Don't bother if no hits
      
      auto const & hits = *thirdHitMap[il];
      
      const DetLayer * layer = theLayers[il].detLayer();
      //auto barrelLayer = layer->isBarrel();
      
            //layer index (lix) Is it necessary or is lix == il?
            int lix = theLayers[il].detLayer()->seqNum();
            
            //Get the correct cut values for the selected layer
            float dtheta = dthetas[lix];
            float dphiss = dphis[lix];
            float dTIP = dTIPs[lix];
            
                            double dthetac = dtheta;

              //  if(fabs(-log(tan(Theta/2)))>1.0)
	     //	dthetac = dtheta*6*exp(fabs(-log(tan(Theta/2))))/2.71;
	
           
           // cout<<"Layer n. "<<il<<" is ID = "<<lix<<". MaxError = "<<rzError[il]<<endl;

            float zpos[4] = {0,0,0,0};
            float zmin, zmax;
            //float minz, maxz;
            
            //static float nSigmaRZ = std::sqrt(12.f); // ...and continue as before
            
            float dphi1 = 0;
            float dphi2 = 0;
            //float sigmaZ = 0;
        
            float thetamin = 0;
            float thetamax = 0;
            float deltarmin = 0;
            float deltarmax = 0;
        
            if (layer->location() == GeomDetEnumerators::barrel) {
                const BarrelDetLayer& bl = dynamic_cast<const BarrelDetLayer&>(*layer);
                float halfThickness  = bl.surface().bounds().thickness()/2;
                float radius = bl.specificSurface().radius();
             
                float rlmin = radius-halfThickness;
                float rlmax = radius+halfThickness;
                
                dphi1 = fabs(acos(d/(2*R_MIN))-acos(rlmin/(2*R_MIN)));
                
                if(dphi1 > Geom::pi()) dphi1 = dphi1 - Geom::twoPi();
                if(dphi1 < -Geom::pi()) dphi1 = dphi1 + Geom::twoPi();
                
                dphi2 = fabs(acos(d/(2*R_MIN))-acos(rlmax/(2*R_MIN)));
                
                if(dphi2 > Geom::pi()) dphi2 = dphi2 - Geom::twoPi();
                if(dphi2 < -Geom::pi()) dphi2 = dphi2 + Geom::twoPi();
            
                //extra based r prediction
                //sigmaZ = nSigmaRZ*rzError[il];
            
          
                //zpos[0] = zi + sqrt(rlmin*rlmin-d02*d02)*cotTheta;
                //zpos[1] = zi + sqrt(rlmax*rlmax-d02*d02)*cotTheta;
                

                
                //angle based r prediction
				thetamin = Theta - Thetacorr*dthetac;
                //thetamin = (1-dthetac)*Theta;
				//thetamin = Theta - dthetac;
                thetamin = wraparound(thetamin);
				thetamax = Theta + Thetacorr*dthetac;
                
                //thetamax = (1+dthetac)*Theta;
				//thetamax = Theta + dthetac;
                thetamax = wraparound(thetamax);
                
                float ro = rvo*s;
                
                deltarmin = s*rlmin-ro;
                deltarmax = s*rlmax-ro;
                zpos[0] = (zo +deltarmin*1/tan(thetamin)) < (zo +deltarmin*1/tan(thetamax)) ? (zo +deltarmin*1/tan(thetamin)) : (zo +deltarmin*1/tan(thetamax));
                zpos[1] = (zo +deltarmin*1/tan(thetamin)) > (zo +deltarmin*1/tan(thetamax)) ? (zo +deltarmin*1/tan(thetamin)) : (zo +deltarmin*1/tan(thetamax));
                
                zpos[2] = (zo +deltarmax*1/tan(thetamin)) < (zo +deltarmax*1/tan(thetamax)) ? (zo +deltarmax*1/tan(thetamin)) : (zo +deltarmax*1/tan(thetamax));
                zpos[3] = (zo +deltarmax*1/tan(thetamin)) > (zo +deltarmax*1/tan(thetamax)) ? (zo +deltarmax*1/tan(thetamin)) : (zo +deltarmax*1/tan(thetamax));

  
            }
            // nothing for the endcaps now
        
            //angle based
            zmin = zpos[0]<=zpos[2] ? zpos[0] : zpos[2];
		    //zmin = zmin - sigmaZ - regOffset;
		    zmin = zmin - regOffset;
            
            zmax = zpos[1]>=zpos[3] ? zpos[1] : zpos[3];
		    //zmax = zmax + sigmaZ + regOffset;
		    zmax = zmax + regOffset;
             
            

            
            dphi1 = dphi1>=dphi2 ? dphi1 : dphi2;
            dphi = dphi >= dphi1 ? dphi : dphi1;
            
            
            float prmin = phi2 - dphi;
            float prmax = phi2 + dphi;
            
            
            //Still to keep?
            //static float nSigmaPhi = 3.f;
            
            
            layerTree.clear(); // Now recover hits in bounding box...
            
            if ((prmax-prmin) > Geom::twoPi())
            {
                prmax=Geom::pi(); prmin = -Geom::pi();
            }
            else
            {
                while (prmax>maxphi) { prmin -= Geom::twoPi(); prmax -= Geom::twoPi();}
                while (prmin<minphi) { prmin += Geom::twoPi(); prmax += Geom::twoPi();}
                // This needs range -twoPi to +twoPi to work
            }
            //End of phirange definition

            KDTreeBox phiZ(prmin, prmax, zmin, zmax);
            hitTree[il].search(phiZ, layerTree);

        
            //float eta = -log(tan(Theta/2));
//cout<<phi2<<"\t"<<prmin<<"\t"<<prmax<<"\t"<<prmax-prmin<<"\t"<<Theta<<"\t"<<thetamin<<"\t"<<thetamax<<"\t"<<thetamax-thetamin<<"\t"<<zmin<<"\t"<<zmax<<"\t"<<zmax-zmin<<endl;

      // std::cout << ip << ": " << theLayers[il].detLayer()->seqNum() << " " << layerTree.size() << " " << prmin << " " << prmax << std::endl;


  for (auto const & ih : layerTree) {
	
	if (theMaxElement!=0 && result.size() >= theMaxElement){
	  result.clear();
	  edm::LogError("TooManyTriplets")<<" number of triples exceeds maximum. no triplets produced.";
	  return;
	}

	auto KDdata = ih.data;


	    // insert here check with comparitor
	    OrderedHitTriplet hittriplet( doublets.hit(ip,HitDoublets::inner), doublets.hit(ip,HitDoublets::outer), hits.theHits[KDdata].hit());
       if (TIPFilter(&hittriplet, dTIP) && PhiFilter(&hittriplet, dphiss) && ThetaFilter(&hittriplet, dthetac)) {
      //if (TIPFilter(&hittriplet, dTIP) && PhiFilter(&hittriplet, dphiss)) {
	  	result.push_back( hittriplet );
	    } else {
	      LogDebug("RejectedTriplet") << "rejected triplet from comparitor ";
	    }
	 	
      }
    }
  }
  // std::cout << "triplets " << result.size() << std::endl;
}

bool TripGenerator::ThetaFilter(OrderedHitTriplet *hittriplet, float deltaTheta) const
{
    
    float h0x, h0y, h0r, h0z;
    float h1x, h1y, h1r, h1z;
    float h2x, h2y, h2r, h2z;
    
    h0x = hittriplet->inner()->globalPosition().x();
    h0y = hittriplet->inner()->globalPosition().y();
    h0z = hittriplet->inner()->globalPosition().z();
    h0r = 	sqrt(h0x*h0x + h0y*h0y);
    
    h1x = hittriplet->middle()->globalPosition().x();
    h1y = hittriplet->middle()->globalPosition().y();
    h1z = hittriplet->middle()->globalPosition().z();
    h1r = 	sqrt(h1x*h1x + h1y*h1y);
    
    h2x = hittriplet->outer()->globalPosition().x();
    h2y = hittriplet->outer()->globalPosition().y();
    h2z = hittriplet->outer()->globalPosition().z();
    h2r = 	sqrt(h2x*h2x + h2y*h2y);

    float dr1 = h1r - h0r;
    float dz1 = h1z - h0z;
    float dr2 = h2r - h1r;
    float dz2 = h2z - h1z;
    
    //float theta1 = atan2(dr1,dz1);
    //float theta2 = atan2(dr2,dz2);
	
	float theta1 = atan(dr1/dz1);
	float theta2 = atan(dr2/dz2);

    
    //if(theta1 < 0 ) theta1 = -theta1;
    //if(theta2 < 0 ) theta2 = -theta2;
    
    float dtheta = theta2/theta1;
    //float dtheta = fabs(theta2-theta1);

    if(fabs(dtheta -1)<= deltaTheta) return true;
    //if(dtheta<= deltaTheta) return true;
    else
     return false;
    //return true;
}

bool TripGenerator::PhiFilter(OrderedHitTriplet *hittriplet, float deltaPhi) const
{
    
    float h0x, h0y;
    float h1x, h1y;
    float h2x, h2y;
    
    h0x = hittriplet->inner()->globalPosition().x();
    h0y = hittriplet->inner()->globalPosition().y();
    
    
    h1x = hittriplet->middle()->globalPosition().x();
    h1y = hittriplet->middle()->globalPosition().y();
    
    h2x = hittriplet->outer()->globalPosition().x();
    h2y = hittriplet->outer()->globalPosition().y();
    
    float dx1 = h1x - h0x;
    float dy1 = h1y - h0y;
    float dx2 = h2x - h1x;
    float dy2 = h2y - h1y;
    
    float phi1 = atan2(dy1,dx1);
    float phi2 = atan2(dy2,dx2);
    
    float dphi = fabs(phi2-phi1);
    if (dphi > Geom::pi()) dphi = dphi - Geom::twoPi();
    if (dphi < -Geom::pi()) dphi  = dphi +Geom::twoPi();
 
    if(dphi<= deltaPhi)	return true;
    else
    return false;
    //return true;
}



bool TripGenerator::TIPFilter(OrderedHitTriplet *hittriplet, float deltaTIP) const
{
    //evaluation of the Transverse impact parameter using the formula from Dtrandlie et al.
	
	
    float h0x, h0y, h0r;
    float h1x, h1y, h1r;
    float h2x, h2y, h2r;
    
    h0x = hittriplet->inner()->globalPosition().x();
    h0y = hittriplet->inner()->globalPosition().y();
    h0r = 	h0x*h0x + h0y*h0y;
    
    h1x = hittriplet->middle()->globalPosition().x();
    h1y = hittriplet->middle()->globalPosition().y();
    h1r = 	h1x*h1x + h1y*h1y;
    
    h2x = hittriplet->outer()->globalPosition().x();
    h2y = hittriplet->outer()->globalPosition().y();
    h2r = 	h2x*h2x + h2y*h2y;
    
    
    TVector3 gp0(h0x, h0y, h0r);
    TVector3 gp1(h1x, h1y, h1r);
    TVector3 gp2(h2x, h2y, h2r);
    
    TVector3  a = gp1 - gp0;
    TVector3  b = gp2 - gp1;
    TVector3  n = a.Cross(b);
    n = n.Unit();
    
    TVector2  cOrigin((-n.X())/(2*n.Z()) , (-n.Y())/(2*n.Z()));
    float c = -(n.X()*gp0.X() + n.Y()*gp0.Y() + n.Z()*gp0.Z()); 
    
    //float sqrArg = (1.0 - n.Z()*n.Z() - 4.0*c*n.Z())/(4.0*n.Z()*n.Z());  //Unused?
    float cR = sqrt((1.0 - n.Z()*n.Z() - 4.0*c*n.Z())/(4.0*n.Z()*n.Z()));
    
    TVector2  v(-cOrigin.X() , -cOrigin.Y());
    v = v.Unit();
    
    TVector2  pCA(cOrigin.X() + cR*v.X(), cOrigin.Y() + cR*v.Y());
    
    float TIP = sqrt(pCA.X()*pCA.X() + pCA.Y()*pCA.Y());
    
    if(TIP<= deltaTIP)	return true;
    else
        return false;
    
	
}


float TripGenerator::wraparound(float angle) const{
    int w = 1;
    if(fabs(angle)>Geom::pi()) w = -1;
    
    if(angle < - Geom::pi()) angle = angle+Geom::twoPi();
    if(angle > Geom::pi()) angle = angle-Geom::twoPi();
    
    return angle*w;
}


float TripGenerator::minRadiusCurvature(float minPt, float Bz)
const {
	// e = 1.602177×10^-19 C  (coulombs)
	const double Q = 1.602177E-19;
	// 1 GeV/c = 5.344286×10^-19 J s/m  (joule seconds per meter)
	const double GEV_C = 5.344286E-19;
    
	return minPt * GEV_C / (Bz * Q) * 1E2;
}


