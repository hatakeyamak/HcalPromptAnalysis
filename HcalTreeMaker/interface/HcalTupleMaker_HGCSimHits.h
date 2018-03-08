#ifndef HcalTupleMaker_HGCSimHits_h
#define HcalTupleMaker_HGCSimHits_h

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
class HcalTupleMaker_HGCSimHits : public edm::EDProducer {
 protected:

  std::vector<edm::InputTag> m_PCaloHitsTags;
  std::vector<edm::EDGetTokenT<edm::PCaloHitContainer>> m_PCaloHitsTokens;

  std::vector<std::string> m_geometrySource;

  const std::string     m_prefix;
  const std::string     m_suffix;

  bool debug=false;

  //const HGCalDDDConstants   *hgccons_;
  //const HcalDDDRecConstants *hcalcons_;

  //HGC Geometry                                                                                                     
  std::vector<const HGCalDDDConstants*> hgcCons_;
  std::vector<const HGCalGeometry*>     hgcGeometry_;
  const HcalDDDSimConstants*            hcCons_;
  const HcalDDDRecConstants*            hcConr_;
  const CaloSubdetectorGeometry*        hcGeometry_;

  std::map<uint32_t, HepGeom::Transform3D> transMap_;
    
  void produce( edm::Event & iEvent, const edm::EventSetup & iSetup ) { 

    if (debug) std::cout << "HcalTupleMaker_HGCSimHits: produce starts" << std::endl;
    
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

	int    cell, sector, subsector, layer, zside;
	int    subdet(0);
	HepGeom::Point3D<float> gcoord;

	unsigned int id_ = it.id();

	// 
	if (nameDetector_ == "HCal") {

	  //
	  // Newer
	  // - Validation/HGCalValidation/plugins/HGCalHitValidation.cc
          // - Validation/HGCalValidation/plugins/HGCGeometryValidation.cc
	  //
	  int z, depth, eta, phi, lay;
	  HcalTestNumbering::unpackHcalIndex(it.id(), subdet, z, depth, eta, phi, lay);
	  if (subdet != static_cast<int>(HcalEndcap)) continue;

	  HcalCellType::HcalCell hccell = hcCons_->cell(subdet, z, lay, eta, phi);
	  double zp  = hccell.rz/10;  // mm -> cm
	  int sign = (z==0)?(-1):(1);
	  zp      *= sign;
	  double rho = zp*tan(2.0*atan(exp(-hccell.eta)));
	  double xp  = rho * cos(hccell.phi); //cm
	  double yp  = rho * sin(hccell.phi); //cm

	  //
	  // Older 
	  // - Validation/HGCalValidation/plugins/HGCalSimHitValidation.cc
	  //
	  HcalDetId detId = HcalHitRelabeller::relabel(id_,hcConr_);
	  subdet           = detId.subdet();
	  if (subdet != static_cast<int>(HcalEndcap)) continue;
	  cell             = detId.ietaAbs();
	  sector           = detId.iphi();
	  subsector        = 1;
	  layer            = detId.depth();
	  zside            = detId.zside();

	  if (debug) std::cout << it.energy() << " " << subdet << std::endl;

	  if (it.energy()>0.5) std::cout << "HcalTupleMaker_HGCSimHits: " 
					 << it.energy() << " " 
					 << nameDetector_ << " " 
					 << subdet << " " << cell << " " << sector << std::endl;

	  std::pair<double,double> etaphi = hcConr_->getEtaPhi(subdet,zside*cell,sector);
	  double rz = hcConr_->getRZ(subdet,zside*cell,layer);	  
	  
	  gcoord = HepGeom::Point3D<float>(rz*cos(etaphi.second)/cosh(etaphi.first),
					   rz*sin(etaphi.second)/cosh(etaphi.first),
					   rz*tanh(etaphi.first));

	  if (debug) std::cout << "HCAL geom comparison: "
		  << "(" << xp         << ", " << yp         << ", " << zp         << ") "  
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << std::endl;
	  //
	  // Use HcalTestNumbering::unpackHcalIndex based method at the end
	  // 
	  gcoord = HepGeom::Point3D<float>(xp,yp,zp);

	  if (debug) std::cout << "HCAL geom comparison: "
		  << "(" << xp         << ", " << yp         << ", " << zp         << ") "  
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << std::endl;

	} else {

	  if (it.energy()>0.5) std::cout << "HcalTupleMaker_HGCSimHits: " 
					 << it.energy() << " " 
					 << nameDetector_ << " " 
					 << subdet << " " << cell << " " << sector << std::endl;
	  
	  if (debug) std::cout << "HcalTupleMaker_HGCSimHits: " << hgcCons_[index]->geomMode() << std::endl;

	  if (hgcCons_[index]->geomMode() == HGCalGeometryMode::Square) {

	    if (debug) std::cout << "HcalTupleMaker_HGCSimHits: in the square mode." << std::endl;

	    HGCalTestNumbering::unpackSquareIndex(id_, zside, layer, sector, subsector, cell);
	    std::pair<float,float> xy = hgcCons_[index]->locateCell(cell,layer,subsector,false);
	    const HepGeom::Point3D<float> lcoord(xy.first,xy.second,0);
	    bool symmDet_=true;
	    int subs = (symmDet_ ? 0 : subsector);
	    id_      = HGCalTestNumbering::packSquareIndex(zside,layer,sector,subs,0);
	    gcoord   = (transMap_[id_]*lcoord); // 

	  } else {

	    if (debug) std::cout << "HcalTupleMaker_HGCSimHits: in the non-square mode." << std::endl;

	    HGCalTestNumbering::unpackHexagonIndex(id_, subdet, zside, layer, sector, subsector, cell);
	    std::pair<float,float> xy = hgcCons_[index]->locateCell(cell,layer,sector,false);
	    double zp = hgcCons_[index]->waferZ(layer,false);
	    if (zside < 0) zp = -zp;
	    float  xp = (zp < 0) ? -xy.first/10 : xy.first/10; // mm->cm
	    float  yp = xy.second/10; //mm->cm
	    gcoord = HepGeom::Point3D<float>(xp,yp,zp); // 

	  }

	  //
	  // 
	  //  
	  HGCalTestNumbering::unpackHexagonIndex(id_, subdet, zside, layer, sector, subsector, cell);      
	  // sector: wafer
	  // subsector: celltype
	  std::pair<float, float> xy;
	  std::pair<int,float> layerIdx;
	  double zp, xp, yp;
          xy = hgcCons_[index]->locateCell(cell,layer,sector,false); //mm
          zp = hgcCons_[index]->waferZ(layer,false); //cm 
          if (zside < 0) zp = -zp;
          xp = (zp<0) ? -xy.first/10 : xy.first/10; //mm
          yp = xy.second/10; //mm	  

	  if (debug) std::cout << "HGC geom comparison: "
		  << "(" << xp         << ", " << yp         << ", " << zp         << ") "  
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << std::endl;

	}  // if nameDetector_ 

	double tof = (gcoord.mag()*CLHEP::mm)/CLHEP::c_light; 
	
	v_energy -> push_back ( it.energy() );

	if (debug) std::cout << "HcalTupleMaker_HGCSimHits: " << nameDetector_ << " "
		  << it.time() << " "
		  << it.time()-tof << " "
		  << subdet << " "
		  << layer << " "
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << gcoord.getEta() << " "
		  << gcoord.getPhi() << " "
		  << std::endl;

	v_time   -> push_back ( it.time() );	
	//v_time   -> push_back ( it.time()-tof );	
	v_subdet -> push_back ( subdet );
	v_layer  -> push_back ( layer );
	v_index  -> push_back ( index );
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
  
 HcalTupleMaker_HGCSimHits(const edm::ParameterSet& iConfig) :
  m_PCaloHitsTags (iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("source")),
    m_PCaloHitsTokens(edm::vector_transform(m_PCaloHitsTags, [this](edm::InputTag const & tag){
	  return consumes<edm::PCaloHitContainer>(tag);})),
    m_geometrySource (iConfig.getUntrackedParameter<std::vector<std::string> >("geometrySource")),
    m_prefix         (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
    m_suffix         (iConfig.getUntrackedParameter<std::string>  ("Suffix")) {

    produces<std::vector<float> > ( m_prefix + "Energy" + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Time"   + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Subdet" + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Layer"  + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Index"  + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posx"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posy"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posz"   + m_suffix );

  }

  std::unique_ptr<std::vector<float> > v_energy;
  std::unique_ptr<std::vector<float> > v_energyem;
  std::unique_ptr<std::vector<float> > v_energyhad;
  std::unique_ptr<std::vector<float> > v_time;
  std::unique_ptr<std::vector<int>   > v_id;
  std::unique_ptr<std::vector<int>   > v_index; // index for different input collections
  std::unique_ptr<std::vector<int>   > v_subdet;

  std::unique_ptr<std::vector<int  > > v_cell; // 
  std::unique_ptr<std::vector<int  > > v_sector; // wafer
  std::unique_ptr<std::vector<int  > > v_subsector; // type
  std::unique_ptr<std::vector<int  > > v_layer;
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
      
      // HCAL for BH/HEB
      if (m_geometrySource[i].find("HCal") != std::string::npos) {
	edm::ESHandle<HcalDDDSimConstants> pHSNDC;
	iSetup.get<HcalSimNumberingRecord>().get(pHSNDC);
	if (pHSNDC.isValid()) {
	  hcCons_ = pHSNDC.product();
	  hgcCons_.push_back(0);
	} else {
	  edm::LogWarning("HGCalValid") << "Cannot initiate HcalDDDSimConstants: "
					<< m_geometrySource[i] << std::endl;
	}
	edm::ESHandle<HcalDDDRecConstants> pHRNDC;
	iSetup.get<HcalRecNumberingRecord>().get(pHRNDC);
	if (pHRNDC.isValid()) {
	  hcConr_ = pHRNDC.product();
	} else {
	  edm::LogWarning("HGCalValid") << "Cannot initiate HcalDDDRecConstants: "
					<< m_geometrySource[i] << std::endl;
	}
	edm::ESHandle<CaloGeometry> caloG;
	iSetup.get<CaloGeometryRecord>().get(caloG);
	if (caloG.isValid()) {
	  const CaloGeometry* geo = caloG.product();
	  hcGeometry_ = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
	  hgcGeometry_.push_back(0);
	} else {
	  edm::LogWarning("HGCalValid") << "Cannot initiate HcalGeometry for "
					<< m_geometrySource[i] << std::endl;
	}
      }
      // HGC for EE & HEF
      else {
	edm::ESHandle<HGCalDDDConstants> hgcCons;
	iSetup.get<IdealGeometryRecord>().get(m_geometrySource[i],hgcCons);
	if (hgcCons.isValid()) {
	  hgcCons_.push_back(hgcCons.product());
	} else {
	  edm::LogWarning("HGCalValid") << "Cannot initiate HGCalDDDConstants for "
					<< m_geometrySource[i] << std::endl;
	}
	edm::ESHandle<HGCalGeometry> hgcGeom;
	iSetup.get<IdealGeometryRecord>().get(m_geometrySource[i],hgcGeom);	
	if(hgcGeom.isValid()) {
	  hgcGeometry_.push_back(hgcGeom.product());	
	} else {
	  edm::LogWarning("HGCalValid") << "Cannot initiate HGCalGeometry for "
					<< m_geometrySource[i] << std::endl;
	}
      }
    }

  }

 protected:

  void loadAlgo(){
    v_energy = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_time   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_subdet = std::unique_ptr<std::vector<int>   > ( new std::vector<int  > ());
    v_layer  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
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
    iEvent.put( move(v_subdet ), m_prefix + "Subdet" + m_suffix );
    iEvent.put( move(v_layer  ), m_prefix + "Layer"  + m_suffix );
    iEvent.put( move(v_index  ), m_prefix + "Index"  + m_suffix );
    iEvent.put( move(v_eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(v_phi    ), m_prefix + "Phi"    + m_suffix );
    iEvent.put( move(v_posx   ), m_prefix + "Posx"   + m_suffix );
    iEvent.put( move(v_posy   ), m_prefix + "Posy"   + m_suffix );
    iEvent.put( move(v_posz   ), m_prefix + "Posz"   + m_suffix );

  }  
    
};

#endif
