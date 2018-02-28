#ifndef HcalTupleMaker_HGCRecHits_h
#define HcalTupleMaker_HGCRecHits_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

// HGCAL objects
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

//
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

// Used:
// http://cmslxr.fnal.gov/source/AnalysisAlgos/SiStripClusterInfoProducer/plugins/SiStripProcessedRawDigiProducer.cc
// as an example

class HcalTupleMaker_HGCRecHits : public edm::EDProducer {
 protected:

  std::vector<edm::InputTag> m_HGCRecHitsTags;
  std::vector<edm::EDGetTokenT<HGCRecHitCollection>> m_HGCRecHitsTokens;

  std::vector<std::string> m_geometrySource;

  const std::string     m_prefix;
  const std::string     m_suffix;

  bool debug=false;
  
  void produce( edm::Event & iEvent, const edm::EventSetup & iSetup) { 
    
    //-----------------------------------------------------
    // Prepare to put things into event
    //-----------------------------------------------------
    
    loadAlgo();

    //-----------------------------------------------------
    // edm::Handles
    //-----------------------------------------------------
    
    edm::Handle<HGCRecHitCollection> HGCRecHits;

    //
    // Loop over HGCRecHitCollections
    //
    for( typename std::vector<edm::EDGetTokenT<HGCRecHitCollection> >::const_iterator
	   token = m_HGCRecHitsTokens.begin(); token != m_HGCRecHitsTokens.end(); ++token ) {
      unsigned index(token - m_HGCRecHitsTokens.begin());
      HGCRecHits.clear();

      //
      // getByToken
      // 
      iEvent.getByToken(*token, HGCRecHits);
      if( HGCRecHits.isValid() && !HGCRecHits->empty() ) {
	edm::LogInfo("Input found") << m_HGCRecHitsTags.at(index);
      } else {
	edm::LogInfo("Input not found") << m_HGCRecHitsTags.at(index);
	return;
      }

      //
      // Geometry & looping over rechits
      //
      std::string nameDetector_ = m_geometrySource[index];
      if (debug) std::cout << nameDetector_ << std::endl;

      if (nameDetector_ == "HCal") {

	edm::ESHandle<CaloGeometry> geom;
	iSetup.get<CaloGeometryRecord>().get(geom);
	if (!geom.isValid()) {
	  edm::LogWarning("HcalTupleMaker_HGCRecHits") 
	    << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
	  return;
	}
	const CaloGeometry* geom0 = geom.product();

	for (const auto & it : *(HGCRecHits.product())) {
	  
	  DetId detId = it.id();
	  int ilayer = HcalDetId(detId).depth();
	  
	  run(detId, ilayer, index, geom0, &it);

	}

      } else {

	edm::ESHandle<HGCalGeometry> geom;
	iSetup.get<IdealGeometryRecord>().get(nameDetector_, geom);
	if (!geom.isValid()) {
	  edm::LogWarning("HcalTupleMaker_HGCRecHits") 
	    << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
	  return;	
	}
	const HGCalGeometry* geom0 = geom.product();

	for (const auto & it : *(HGCRecHits.product())) {
	  
	  DetId detId = it.id();
	  int ilayer   = HGCalDetId(detId).layer();

	  run(detId, ilayer, index, geom0, &it);

	}

      } // Hcal or HGcal for geometry

    } // Looping over different rechit collections

    //-----------------------------------------------------
    // Put things into the event
    //-----------------------------------------------------
    
    dumpAlgo(iEvent);

    }
  
 public:
  
 HcalTupleMaker_HGCRecHits(const edm::ParameterSet& iConfig) :
  m_HGCRecHitsTags (iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("source")),
    m_HGCRecHitsTokens(edm::vector_transform(m_HGCRecHitsTags, [this](edm::InputTag const & tag){
	  return consumes<HGCRecHitCollection>(tag);})),
    m_geometrySource (iConfig.getUntrackedParameter<std::vector<std::string> >("geometrySource")),
    m_prefix         (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
    m_suffix         (iConfig.getUntrackedParameter<std::string>  ("Suffix")) {

    produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Layer"  + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Energy" + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posx"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posy"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posz"   + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Subdet" + m_suffix );
    //produces<std::vector<int>   > ( m_prefix + "Flags"  + m_suffix );
    //produces<std::vector<int>   > ( m_prefix + "Aux"    + m_suffix );
    //produces<std::vector<float> > ( m_prefix + "Time"   + m_suffix );
  }

  std::unique_ptr<std::vector<float> > eta;
  std::unique_ptr<std::vector<float> > phi;
  std::unique_ptr<std::vector<int  > > layer;
  std::unique_ptr<std::vector<float> > energy;
  std::unique_ptr<std::vector<float> > posx;
  std::unique_ptr<std::vector<float> > posy;
  std::unique_ptr<std::vector<float> > posz;
  std::unique_ptr<std::vector<int>   > subdet;

 protected:

  void loadAlgo(){
    eta    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    phi    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    layer  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    energy = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    posx   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    posy   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    posz   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    subdet = std::unique_ptr<std::vector<int>   > ( new std::vector<int  > ());
  }
  
  void dumpAlgo( edm::Event & iEvent ){
    iEvent.put( move(eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(phi    ), m_prefix + "Phi"    + m_suffix );
    iEvent.put( move(layer  ), m_prefix + "Layer"  + m_suffix );
    iEvent.put( move(energy ), m_prefix + "Energy" + m_suffix );
    iEvent.put( move(posx   ), m_prefix + "Posx"   + m_suffix );
    iEvent.put( move(posy   ), m_prefix + "Posy"   + m_suffix );
    iEvent.put( move(posz   ), m_prefix + "Posz"   + m_suffix );
    iEvent.put( move(subdet ), m_prefix + "Subdet" + m_suffix );
  }  

  template<class T1, class T2>
  void run(DetId & detId, int ilayer, int index,
				      const T1* geom, T2 it) {

    GlobalPoint global = geom->getPosition(detId);
    //double      energy = it->energy();

    /*
    float globalx = global.x();
    float globaly = global.y();
    float globalz = global.z();
    */

    if (debug) std::cout << ilayer << " " << global.z() << std::endl;

    eta    -> push_back ( global.eta() );
    phi    -> push_back ( global.phi() );
    layer  -> push_back ( ilayer       );
    energy -> push_back ( it->energy() );
    posx   -> push_back ( global.x()   );
    posy   -> push_back ( global.y()   );
    posz   -> push_back ( global.z()   );
    subdet -> push_back ( index        );

  }
    
};

#endif
