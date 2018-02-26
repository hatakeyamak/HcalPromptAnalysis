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

//#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_HGCRecHitAlgorithm.h"

// Used:
// http://cmslxr.fnal.gov/source/AnalysisAlgos/SiStripClusterInfoProducer/plugins/SiStripProcessedRawDigiProducer.cc
// as an example

class HcalTupleMaker_HGCRecHits : public edm::EDProducer {
 protected:

  /*
  const edm::InputTag   m_HGCRecHitsTag;
  //edm::EDGetTokenT<RecHitCollection> m_HGCRecHitsToken;
  edm::EDGetTokenT<HGCRecHitCollection> m_HGCRecHitsToken;
  */

  std::vector<edm::InputTag> m_HGCRecHitsTags;
  std::vector<edm::EDGetTokenT<HGCRecHitCollection>> m_HGCRecHitsTokens;

  std::vector<std::string> m_geometrySource;

  const std::string     m_prefix;
  const std::string     m_suffix;

  //HcalTupleMaker_HGCRecHitAlgorithm algo;
  
  void produce( edm::Event & iEvent, const edm::EventSetup & iSetup) { 
    
    //-----------------------------------------------------
    // Prepare to put things into event
    //-----------------------------------------------------
    
    loadAlgo();
    /*
    std::unique_ptr<std::vector<double> >            eta     ( new std::vector<double>           ());
    std::unique_ptr<std::vector<double> >            phi     ( new std::vector<double>           ());
    std::unique_ptr<std::vector<int   > >            layer   ( new std::vector<int>              ());
    std::unique_ptr<std::vector<double> >            energy  ( new std::vector<double>           ());
    */

    //-----------------------------------------------------
    // edm::Handles
    //-----------------------------------------------------
    
    //bool run_algo = true;

    //std::vector<edm::Handle<HGCRecHitCollection> > HGCRecHitCollections;
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
	//run_algo = true;
      } else {
	edm::LogInfo("Input not found") << m_HGCRecHitsTags.at(index);
	return;
	//run_algo = false;
      }

      //
      // Geometry & looping over rechits
      //
      std::string nameDetector_ = m_geometrySource[index];
      std::cout << nameDetector_ << std::endl;

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
	  
	  run(detId, ilayer, geom0, &it);

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

	  run(detId, ilayer, geom0, &it);
	  std::cout << ilayer << std::endl;

	}

      } // Hcal or HGcal for geometry

    } // Looping over different rechit collections

      /*
      edm::ESHandle<HGCalGeometry> geom;
      iSetup.get<IdealGeometryRecord>().get(nameDetector, geom);
      if (!geom.isValid()) 
	edm::LogWarning("HGCalValidation") << "Cannot get valid HGCalGeometry Object for " << nameDetector;
      const HGCalGeometry* geometry = geom.product();
      */

      //
      // Loop over rechits
      //
      //if ( run_algo ) algo.run ( *HGCRecHits, *geom0 );

	
      //layer = HGCalDetId(detId).layer();
	
      //}


      //}

    // Loop starts ----------

    //std::string label = findInput(m_HGCRecHitsTags, m_HGCRecHitsTokens, iEvent);

    /*
    bool gotHGCRecHits = iEvent.getByToken(m_HGCRecHitsToken, HGCRecHits);
    if (!gotHGCRecHits ) {
      std::cout << "Could not find HGCRecHits with tag " << m_HGCRecHitsTag << std::endl;
      run_algo = false;
    }
    */

    //-----------------------------------------------------
    // edm::ESHandles
    //-----------------------------------------------------

    /*

    std::string nameDetector = "HGCalHESiliconSensitive";

    edm::ESHandle<HGCalGeometry> geom;
    iSetup.get<IdealGeometryRecord>().get(nameDetector, geom);
    if (!geom.isValid()) 
      edm::LogWarning("HGCalValidation") << "Cannot get valid HGCalGeometry Object for " << nameDetector;
    const HGCalGeometry* geometry = geom.product();

    */
    
    //-----------------------------------------------------
    // Run the algorithm
    //-----------------------------------------------------

    /*
    if ( run_algo ) algo.run ( *HGCRecHits, *geometry );
    */

    // Loop may end ----------

    //-----------------------------------------------------
    // Put things into the event
    //-----------------------------------------------------
    
    dumpAlgo(iEvent);
    /*
    iEvent.put(move(eta)    , m_prefix + "Eta"           + m_suffix );
    iEvent.put(move(phi)    , m_prefix + "Phi"           + m_suffix );
    iEvent.put(move(layer)  , m_prefix + "Layer"         + m_suffix );
    iEvent.put(move(energy) , m_prefix + "Energy"        + m_suffix );
    */

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
    //produces<std::vector<int>   > ( m_prefix + "Subdet" + m_suffix );
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

 protected:

  void loadAlgo(){
    eta    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    phi    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    layer  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    energy = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    posx   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    posy   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    posz   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
  }
  
  void dumpAlgo( edm::Event & iEvent ){
    iEvent.put( move(eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(phi    ), m_prefix + "Phi"    + m_suffix );
    iEvent.put( move(layer  ), m_prefix + "Layer"  + m_suffix );
    iEvent.put( move(energy ), m_prefix + "Energy" + m_suffix );
    iEvent.put( move(posx   ), m_prefix + "Posx"   + m_suffix );
    iEvent.put( move(posy   ), m_prefix + "Posy"   + m_suffix );
    iEvent.put( move(posz   ), m_prefix + "Posz"   + m_suffix );
  }  

  template<class T1, class T2>
  void run(DetId & detId, int ilayer, 
				      const T1* geom, T2 it) {

    GlobalPoint global = geom->getPosition(detId);
    //double      energy = it->energy();

    /*
    float globalx = global.x();
    float globaly = global.y();
    float globalz = global.z();
    */

    eta    -> push_back ( global.eta() );
    phi    -> push_back ( global.phi() );
    layer  -> push_back ( ilayer       );
    energy -> push_back ( it->energy() );
    posx   -> push_back ( global.x()   );
    posy   -> push_back ( global.y()   );
    posz   -> push_back ( global.z()   );

  }
    
};

#endif
