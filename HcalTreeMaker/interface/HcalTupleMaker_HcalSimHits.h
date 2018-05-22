#ifndef HcalTupleMaker_HcalSimHits_h
#define HcalTupleMaker_HcalSimHits_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/transform.h"

// PCaloHits objects
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

// HGCAL & HCAL Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"

#include "DataFormats/HcalDetId/interface/HcalTestNumbering.h"

#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"

//
// Class definition
// 
class HcalTupleMaker_HcalSimHits : public edm::EDProducer {
 protected:

  std::vector<edm::InputTag> m_PCaloHitsTags;
  std::vector<edm::EDGetTokenT<edm::PCaloHitContainer>> m_PCaloHitsTokens;

  std::vector<std::string> m_geometrySource;

  const std::string     m_prefix;
  const std::string     m_suffix;

  bool debug=false;
  bool testNumber_=true;
  
  //HCal Geometry                                                                                                     
  const HcalDDDSimConstants*            hcCons_;
  const HcalDDDRecConstants*            hcConr_;
  const CaloSubdetectorGeometry*        hcGeometry_;

  edm::ESHandle<CaloGeometry> geometry ;
  
  std::map<uint32_t, HepGeom::Transform3D> transMap_;
    
  void produce( edm::Event & iEvent, const edm::EventSetup & iSetup ) { 

    if (debug) std::cout << "HcalTupleMaker_HcalSimHits: produce starts" << std::endl;
    
    //-----------------------------------------------------
    // Prepare to put things into event
    //-----------------------------------------------------
    
    loadAlgo();

    //-----------------------------------------------------
    // edm::Handles
    //-----------------------------------------------------

    edm::Handle<edm::PCaloHitContainer> PCaloHits;

    //
    // Loop over PCaloHit containers
    //
    for( typename std::vector<edm::EDGetTokenT<edm::PCaloHitContainer> >::const_iterator
	   token = m_PCaloHitsTokens.begin(); token != m_PCaloHitsTokens.end(); ++token ) {
      unsigned index(token - m_PCaloHitsTokens.begin());
      PCaloHits.clear();

      //
      // Geometry & looping over rechits
      //                                                                                                              
      std::string nameDetector_ = m_geometrySource[index];
      if (debug) std::cout << nameDetector_ << std::endl;

      //
      // getByToken
      // 
      iEvent.getByToken(*token, PCaloHits);
      if( PCaloHits.isValid() && !PCaloHits->empty() ) {
	if (debug) std::cout << "Input found" << std::endl;
	edm::LogInfo("Input found") << m_PCaloHitsTags.at(index);
      } else {
	if (debug) std::cout << "Input not found" << std::endl;
	edm::LogInfo("Input not found") << m_PCaloHitsTags.at(index);
	continue;
      }
      if (debug) std::cout << nameDetector_ << " " << PCaloHits->size() << std::endl;

      //
      // Loop over PCaloHits
      //
      for (const auto & it : *(PCaloHits.product())) {	  

	int    subdet(0);
	HepGeom::Point3D<float> gcoord;

	unsigned int id_ = it.id();
	HcalDetId detId;

	// 
	if (nameDetector_ == "HCal") {

	  //
	  // Direct HcalTestNumbering method:
	  // - Validation/HGCalValidation/plugins/HGCalHitValidation.cc
          // - Validation/HGCalValidation/plugins/HGCGeometryValidation.cc
	  //

	  /*
	  int z, depth, ieta, iphi, lay;
	  
	  HcalTestNumbering::unpackHcalIndex(it.id(), subdet, z, depth, ieta, iphi, lay);
	  //KH if (subdet != static_cast<int>(HcalEndcap)) continue;

	  
	  HcalCellType::HcalCell hccell = hcCons_->cell(subdet, z, lay, ieta, iphi);
	  //double zp  = hccell.rz/10*tanh(hccell.eta);  // mm -> cm
	  double zp  = hccell.rz/10;  // mm -> cm, rz is actually Z?
	  int sign = (z==0)?(-1):(1);
	  zp      *= sign;
	  double rho = zp*tan(2.0*atan(exp(-hccell.eta)));
	  double xp  = rho * cos(hccell.phi); //cm
	  double yp  = rho * sin(hccell.phi); //cm
	  */

	  //
	  // HcalHitRelabeller method:
	  // - Validation/HGCalValidation/plugins/HGCalSimHitValidation.cc
	  //
	  if (testNumber_) detId = HcalHitRelabeller::relabel(id_,hcConr_);
	  else detId = HcalDetId(id_);

	  subdet           = detId.subdet();
	  //KH if (subdet != static_cast<int>(HcalEndcap)) continue;
	  //cell             = detId.ietaAbs();
	  //sector           = detId.iphi();
	  //subsector        = 1;
	  //depth              = detId.depth();
	  //zside            = detId.zside();

	  if (debug) std::cout << it.energy() << " " << subdet << std::endl;

	  if (debug)
	  if (it.energy()>0.5) std::cout << "HcalTupleMaker_HcalSimHits: " 
					 << it.energy() << " " 
					 << nameDetector_ << " " 
					 << subdet << std::endl;

	  /*
	  std::pair<double,double> etaphi = hcConr_->getEtaPhi(subdet,zside*cell,sector);
	  double rz = hcConr_->getRZ(subdet,zside*cell,layer);	  // This is actually Z?
	  
	  gcoord = HepGeom::Point3D<float>(rz*cos(etaphi.second)/cosh(etaphi.first)/tanh(etaphi.first),
					   rz*sin(etaphi.second)/cosh(etaphi.first)/tanh(etaphi.first),
					   rz);
	  */

	  //
	  // Use CaloCellGeometry getPosition
	  // 
	  iSetup.get<CaloGeometryRecord>().get (geometry);
	  auto cellGeometry = geometry->getSubdetectorGeometry(detId)->getGeometry(detId);
	  
	  if (debug) 
	  std::cout << "HCAL geom comparison: "
	    //<< "(" << xp         << ", " << yp         << ", " << zp         << ") "  
	    //<< rho << " "
	    //<< "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		    << "(" << cellGeometry->getPosition().x() << ", " << cellGeometry->getPosition().y() << ", " << cellGeometry->getPosition().z() << ") "  
		    << std::endl;


	  //
	  // Use CaloCellGeometry getPosition() method at the end
	  // 
	  gcoord = HepGeom::Point3D<float>(cellGeometry->getPosition().x(),
					   cellGeometry->getPosition().y(),
					   cellGeometry->getPosition().z());


	} else {

	}  // if nameDetector_ 

	double tof = (gcoord.mag()*CLHEP::mm)/CLHEP::c_light; 

	if (debug) std::cout << "HcalTupleMaker_HcalSimHits: " << nameDetector_ << " "
		  << it.time() << " "
		  << it.time()-tof << " "
		  << subdet << " "
   	          << detId.depth() << " "
   	          << detId.ieta()  << " "
   	          << detId.iphi()  << " "
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << gcoord.getEta() << " "
		  << gcoord.getPhi() << " "
		  << std::endl;
	
	v_energy -> push_back ( it.energy() );
	v_time   -> push_back ( it.time() );	
	v_time_tof -> push_back ( it.time()-tof );	
	v_subdet -> push_back ( subdet );
	v_ieta   -> push_back ( detId.ieta() );
	v_iphi   -> push_back ( detId.iphi() );
	v_depth  -> push_back ( detId.depth() );
	v_index  -> push_back ( index );
	v_eta    -> push_back ( gcoord.getEta() );
	v_phi    -> push_back ( gcoord.getPhi() );
	v_eta    -> push_back ( gcoord.getEta() );
	v_phi    -> push_back ( gcoord.getPhi() );
	v_posx   -> push_back ( gcoord.x() );
	v_posy   -> push_back ( gcoord.y() );
	v_posz   -> push_back ( gcoord.z() );

      } // for-loop of PCaloHits

    } // Looping over different PCaloHit collections

    //-----------------------------------------------------
    // Put things into the event
    //-----------------------------------------------------
    
    dumpAlgo(iEvent);

  }
  
 public:
  
 HcalTupleMaker_HcalSimHits(const edm::ParameterSet& iConfig) :
  m_PCaloHitsTags (iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("source")),
    m_PCaloHitsTokens(edm::vector_transform(m_PCaloHitsTags, [this](edm::InputTag const & tag){
	  return consumes<edm::PCaloHitContainer>(tag);})),
    m_geometrySource (iConfig.getUntrackedParameter<std::vector<std::string> >("geometrySource")),
    m_prefix         (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
    m_suffix         (iConfig.getUntrackedParameter<std::string>  ("Suffix")) {

    produces<std::vector<float> > ( m_prefix + "Energy" + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Time"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "TimeTOF"   + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Subdet" + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Ieta"   + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Iphi"   + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Depth"  + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Index"  + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posx"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posy"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posz"   + m_suffix );

  }

  std::unique_ptr<std::vector<float> > v_energy;
  std::unique_ptr<std::vector<float> > v_time;
  std::unique_ptr<std::vector<float> > v_time_tof;
  std::unique_ptr<std::vector<int>   > v_id;
  std::unique_ptr<std::vector<int>   > v_index; // index for different input collections
  std::unique_ptr<std::vector<int>   > v_subdet;

  std::unique_ptr<std::vector<int  > > v_ieta;
  std::unique_ptr<std::vector<int  > > v_iphi;
  std::unique_ptr<std::vector<int  > > v_depth;
  std::unique_ptr<std::vector<int  > > v_zside;

  std::unique_ptr<std::vector<float> > v_posx;
  std::unique_ptr<std::vector<float> > v_posy;
  std::unique_ptr<std::vector<float> > v_posz;
  std::unique_ptr<std::vector<float> > v_eta;
  std::unique_ptr<std::vector<float> > v_phi;

 private:

  void beginRun(const edm::Run&, const edm::EventSetup& iSetup){

    //initiating HGC Geometry
    for (size_t i=0; i<m_geometrySource.size(); i++) {
      
      // HCAL 
      if (m_geometrySource[i].find("HCal") != std::string::npos) {

	//---HcalDDDSimConstants
	/*
	edm::ESHandle<HcalDDDSimConstants> pHSNDC;
	iSetup.get<HcalSimNumberingRecord>().get(pHSNDC);
	if (pHSNDC.isValid()) {
	  hcCons_ = pHSNDC.product();
	  //KH hgcCons_.push_back(0);
	} else {
	  edm::LogWarning("HcalTupleMaker_HcalSimHits") << "Cannot initiate HcalDDDSimConstants: "
					<< m_geometrySource[i] << std::endl;
	}
	*/

	//---HcalDDDRecConstants
	edm::ESHandle<HcalDDDRecConstants> pHRNDC;
	iSetup.get<HcalRecNumberingRecord>().get(pHRNDC);
	if (pHRNDC.isValid()) {
	  hcConr_ = pHRNDC.product();
	} else {
	  edm::LogWarning("HcalTupleMaker_HcalSimHits") << "Cannot initiate HcalDDDRecConstants: "
					<< m_geometrySource[i] << std::endl;
	}

	//---CaloGeometryRecord
	edm::ESHandle<CaloGeometry> caloG;
	iSetup.get<CaloGeometryRecord>().get(caloG);
	if (caloG.isValid()) {
	  const CaloGeometry* geo = caloG.product();
	  hcGeometry_ = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
	  //KH hgcGeometry_.push_back(0);
	} else {
	  edm::LogWarning("HcalTupleMaker_HcalSimHits") << "Cannot initiate HcalGeometry for "
					<< m_geometrySource[i] << std::endl;
	}

      }
      // Not Hcal ?
      else {
      }
      
    }

  }

 protected:

  void loadAlgo(){
    v_energy = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_time   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_time_tof   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_subdet = std::unique_ptr<std::vector<int>   > ( new std::vector<int  > ());
    v_ieta   = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    v_iphi   = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    v_depth  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    v_index  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    v_eta    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_phi    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posx   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posy   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posz   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());

  }
  
  void dumpAlgo( edm::Event & iEvent ){
    iEvent.put( move(v_energy ), m_prefix + "Energy" + m_suffix );
    iEvent.put( move(v_time   ), m_prefix + "Time"   + m_suffix );
    iEvent.put( move(v_time_tof   ), m_prefix + "TimeTOF"   + m_suffix );
    iEvent.put( move(v_subdet ), m_prefix + "Subdet" + m_suffix );
    iEvent.put( move(v_ieta   ), m_prefix + "Ieta"   + m_suffix );
    iEvent.put( move(v_iphi   ), m_prefix + "Iphi"   + m_suffix );
    iEvent.put( move(v_depth  ), m_prefix + "Depth"  + m_suffix );
    iEvent.put( move(v_index  ), m_prefix + "Index"  + m_suffix );
    iEvent.put( move(v_eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(v_phi    ), m_prefix + "Phi"    + m_suffix );
    iEvent.put( move(v_posx   ), m_prefix + "Posx"   + m_suffix );
    iEvent.put( move(v_posy   ), m_prefix + "Posy"   + m_suffix );
    iEvent.put( move(v_posz   ), m_prefix + "Posz"   + m_suffix );

  }  
    
};

#endif
