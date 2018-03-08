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

// test
#include <cmath>
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DetectorDescription/Core/interface/DDExpandedView.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include "Geometry/HGCalCommonData/interface/HGCalGeometryMode.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include <CLHEP/Geometry/Transform3D.h>

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

  const HGCalDDDConstants   *hgccons_;
  const HcalDDDRecConstants *hcalcons_;

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
	bool symmDet_=true;
	if (nameDetector_ == "HCal") {

	  //
	  // Newer
	  // - Validation/HGCalValidation/plugins/HGCalHitValidation.cc
          // - Validation/HGCalValidation/plugins/HGCGeometryValidation.cc
	  //
	  std::cout << "Hcal geom check starts" << std::endl;
	  int z, depth, eta, phi, lay;
	  HcalTestNumbering::unpackHcalIndex(it.id(), subdet, z, depth, eta, phi, lay);
	  if (subdet != static_cast<int>(HcalEndcap)) continue;

	  std::cout << "Hcal geom check starts in endcap" << std::endl;
	  HcalCellType::HcalCell hccell = hcCons_->cell(subdet, z, lay, eta, phi);
	  double zp  = hccell.rz/10;  // mm -> cm
	  int sign = (z==0)?(-1):(1);
	  zp      *= sign;
	  double rho = zp*tan(2.0*atan(exp(-hccell.eta)));
	  double xp  = rho * cos(hccell.phi); //cm
	  double yp  = rho * sin(hccell.phi); //cm

	  //KH HcalDDDRecConstants::HcalID idx = hcConr_->getHCID(subdet,eta,phi,lay,depth);
	  //KH HcalDetId id = HcalDetId(HcalEndcap,sign*idx.eta,idx.phi,idx.depth);

	  //
	  // Older 
	  // - Validation/HGCalValidation/plugins/HGCalSimHitValidation.cc
	  //
	  std::cout << "Hcal geom check starts - older method" << std::endl;
	  HcalDetId detId = HcalHitRelabeller::relabel(id_,hcalcons_);
	  subdet           = detId.subdet();
	  if (subdet != static_cast<int>(HcalEndcap)) continue;
	  cell             = detId.ietaAbs();
	  sector           = detId.iphi();
	  subsector        = 1;
	  layer            = detId.depth();
	  zside            = detId.zside();

	  if (debug) std::cout << it.energy() << " " << subdet << std::endl;

	  if (it.energy()>0.5) std::cout << it.energy() << " " 
					 << nameDetector_ << " " 
					 << subdet << " " << cell << " " << sector << std::endl;

	  std::pair<double,double> etaphi = hcalcons_->getEtaPhi(subdet,zside*cell,sector);
	  double rz = hcalcons_->getRZ(subdet,zside*cell,layer);	  
	  
	  gcoord = HepGeom::Point3D<float>(rz*cos(etaphi.second)/cosh(etaphi.first),
					   rz*sin(etaphi.second)/cosh(etaphi.first),
					   rz*tanh(etaphi.first));

	  std::cout << "HCAL geom comparison: "
		  << "(" << xp         << ", " << yp         << ", " << zp         << ") "  
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << std::endl;

	  gcoord = HepGeom::Point3D<float>(xp,yp,zp);

	  std::cout << "HCAL geom comparison: "
		  << "(" << xp         << ", " << yp         << ", " << zp         << ") "  
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << std::endl;

	} else {

	  if (it.energy()>0.5) std::cout << it.energy() << " " 
					 << nameDetector_ << " " 
					 << subdet << " " << cell << " " << sector << std::endl;

	  edm::ESHandle<HGCalDDDConstants>  pHGDC;
	  iSetup.get<IdealGeometryRecord>().get(nameDetector_, pHGDC);
	  hgccons_ = &(*pHGDC);

	  //get geometry
	  /*
	  edm::ESHandle<CaloGeometry> geom;
	  iSetup.get<CaloGeometryRecord>().get(geom);

	  const HGCalGeometry* gHGCal_;

	  if(nameDetector_=="HGCalEESensitive") 
	    gHGCal_ = dynamic_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::Forward, HGCEE));
	  else 
	    gHGCal_ = dynamic_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::Forward, HGCHEF));

	  const auto& topo     = gHGCal_->topology();
	  const auto* dddConst = topo.dddConstants();
	  
	  std::cout << hgccons_->geomMode() << " " << dddConst.geomMode() << std::endl;
	  */
	  std::cout << "HcalTupleMaker_HGCSimHits: " << hgccons_->geomMode() << std::endl;

	  if (hgccons_->geomMode() == HGCalGeometryMode::Square) {

	    std::cout << "HcalTupleMaker_HGCSimHits: in the square mode." << std::endl;

	    HGCalTestNumbering::unpackSquareIndex(id_, zside, layer, sector, subsector, cell);
	    std::pair<float,float> xy = hgccons_->locateCell(cell,layer,subsector,false);
	    const HepGeom::Point3D<float> lcoord(xy.first,xy.second,0);
	    int subs = (symmDet_ ? 0 : subsector);
	    id_      = HGCalTestNumbering::packSquareIndex(zside,layer,sector,subs,0);
	    gcoord   = (transMap_[id_]*lcoord);

	  } else {

	    std::cout << "HcalTupleMaker_HGCSimHits: in the non-square mode." << std::endl;

	    HGCalTestNumbering::unpackHexagonIndex(id_, subdet, zside, layer, sector, subsector, cell);
	    std::pair<float,float> xy = hgccons_->locateCell(cell,layer,sector,false);
	    double zp = hgccons_->waferZ(layer,false);
	    if (zside < 0) zp = -zp;
	    float  xp = (zp < 0) ? -xy.first : xy.first;
	    gcoord = HepGeom::Point3D<float>(xp,xy.second,zp);

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

	  GlobalPoint xyz = hgcGeometry_[index]->getPosition(id_);

	  std::cout << "HGC geom comparison: "
		  << "(" << xp         << ", " << yp         << ", " << zp         << ") "  
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << "(" << xyz.x()    << ", " << xyz.y()    << ", " << xyz.z()    << ") "  
		  << std::endl;

	}  // if nameDetector_ 
	double tof = (gcoord.mag()*CLHEP::mm)/CLHEP::c_light; 
	
	v_energy -> push_back ( it.energy() );

	std::cout << "HcalTupleMaker_HGCSimHits: " << nameDetector_ << " "
		  << it.time() << " "
		  << it.time()-tof << " "
		  << subdet << " "
		  << layer << " "
		  << "(" << gcoord.x() << ", " << gcoord.y() << ", " << gcoord.z() << ") "  
		  << gcoord.getEta() << " "
		  << gcoord.getPhi() << " "
		  << std::endl;
	/*
	v_time   -> push_back ( it.time()-tof );	
	v_subdet -> push_back ( subdet );
	v_layer  -> push_back ( layer );
	v_eta    -> push_back ( gcoord.getEta() );
	v_phi    -> push_back ( gcoord.getPhi() );
	v_posx   -> push_back ( gcoord.x() );
	v_posy   -> push_back ( gcoord.y() );
	v_posz   -> push_back ( gcoord.z() );
	*/

	/* uint32_t id       = hits[i].id(); */
	/* 0195     int    cell, sector, subsector, layer, zside; */
	/* 0196     int    subdet(0); */
	/* 0197     if (heRebuild_) { */
	/*   0198       HcalDetId detId  = HcalDetId(id_); */
	/*   0199       subdet           = detId.subdet(); */
	/*   0200       if (subdet != static_cast<int>(HcalEndcap)) continue; */
	/*   0201       cell             = detId.ietaAbs(); */
	/*   0202       sector           = detId.iphi(); */
	/*   0203       subsector        = 1; */
	/*   0204       layer            = detId.depth(); */
	/*   0205       zside            = detId.zside(); */
	/*   0206     } else { */
	/*   0207       if (hgcons_->geomMode() == HGCalGeometryMode::Square) { */
	/*     0208     HGCalTestNumbering::unpackSquareIndex(id_, zside, layer, sector, subsector, cell); */
	/*     0209       } else { */
	/*     0210     HGCalTestNumbering::unpackHexagonIndex(id_, subdet, zside, layer, sector, subsector, cell); */
	/*     0211       } */
	/*   0212     } */



	/* DetId detId = it.id(); */
	/* int ilayer = HcalDetId(detId).depth();	  */
	//run(detId, ilayer, index, geom0, &it); 
      } // for-loop

      /* // */
      /* // Geometry & looping over rechits */
      /* // */
      /* std::string nameDetector_ = m_geometrySource[index]; */
      /* if (debug) std::cout << nameDetector_ << std::endl; */

      /* if (nameDetector_ == "HCal") { */

      /* 	edm::ESHandle<CaloGeometry> geom; */
      /* 	iSetup.get<CaloGeometryRecord>().get(geom); */
      /* 	if (!geom.isValid()) { */
      /* 	  edm::LogWarning("HcalTupleMaker_HGCSimHits")  */
      /* 	    << "Cannot get valid HGCalGeometry Object for " << nameDetector_; */
      /* 	  return; */
      /* 	} */
      /* 	const CaloGeometry* geom0 = geom.product(); */

      /* 	for (const auto & it : *(HGCSimHits.product())) { */
	  
      /* 	  DetId detId = it.id(); */
      /* 	  int ilayer = HcalDetId(detId).depth(); */
	  
      /* 	  run(detId, ilayer, index, geom0, &it); */

      /* 	} */

      /* } else { */

      /* 	edm::ESHandle<HGCalGeometry> geom; */
      /* 	iSetup.get<IdealGeometryRecord>().get(nameDetector_, geom); */
      /* 	if (!geom.isValid()) { */
      /* 	  edm::LogWarning("HcalTupleMaker_HGCSimHits")  */
      /* 	    << "Cannot get valid HGCalGeometry Object for " << nameDetector_; */
      /* 	  return;	 */
      /* 	} */
      /* 	const HGCalGeometry* geom0 = geom.product(); */

      /* 	for (const auto & it : *(HGCSimHits.product())) { */
	  
      /* 	  DetId detId = it.id(); */
      /* 	  int ilayer   = HGCalDetId(detId).layer(); */

      /* 	  run(detId, ilayer, index, geom0, &it); */

      /* 	} */

      /* } // Hcal or HGcal for geometry */

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
    /*
    produces<std::vector<int>   > ( m_prefix + "Subdet" + m_suffix );
    produces<std::vector<int>   > ( m_prefix + "Layer"  + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Eta"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Phi"    + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posx"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posy"   + m_suffix );
    produces<std::vector<float> > ( m_prefix + "Posz"   + m_suffix );
    */

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

    // -----

    std::cout << "a1" << std::endl; 
    edm::ESHandle<HcalDDDRecConstants> pHRNDC;
    iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
    hcalcons_  = &(*pHRNDC);
    
    std::cout << "a2" << std::endl; 
    std::string nameDetector_ = "HGCalEESensitive";
    //edm::ESHandle<HGCalDDDConstants>  pHGDC;
    //iSetup.get<IdealGeometryRecord>().get(nameDetector_, pHGDC);
    //hgccons_ = &(*pHGDC);
    //
    //edm::ESTransientHandle<DDCompactView> pDD;
    //iSetup.get<IdealGeometryRecord>().get( pDD );
    //defineGeometry(pDD);
    
    std::cout << "a3" << std::endl; 

  }

 protected:

  void loadAlgo(){
    v_energy = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    /*
    v_subdet = std::unique_ptr<std::vector<int>   > ( new std::vector<int  > ());
    v_layer  = std::unique_ptr<std::vector<int  > > ( new std::vector<int  > ());
    v_eta    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_phi    = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posx   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posy   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    v_posz   = std::unique_ptr<std::vector<float> > ( new std::vector<float> ());
    */
  }
  
  void dumpAlgo( edm::Event & iEvent ){
    iEvent.put( move(v_energy ), m_prefix + "Energy" + m_suffix );
    /*
    iEvent.put( move(v_subdet ), m_prefix + "Subdet" + m_suffix );
    iEvent.put( move(v_layer  ), m_prefix + "Layer"  + m_suffix );
    iEvent.put( move(v_eta    ), m_prefix + "Eta"    + m_suffix );
    iEvent.put( move(v_phi    ), m_prefix + "Phi"    + m_suffix );
    iEvent.put( move(v_posx   ), m_prefix + "Posx"   + m_suffix );
    iEvent.put( move(v_posy   ), m_prefix + "Posy"   + m_suffix );
    iEvent.put( move(v_posz   ), m_prefix + "Posz"   + m_suffix );
    */
  }  

  /* template<class T1, class T2> */
  /* void run(DetId & detId, int ilayer, int index, */
  /* 				      const T1* geom, T2 it) { */

  /*   GlobalPoint global = geom->getPosition(detId); */

  /*   if (debug) std::cout << ilayer << " " << global.z() << std::endl; */

  /*   energy -> push_back ( it->energy() ); */
  /*   subdet -> push_back ( index        ); */
  /*   layer  -> push_back ( ilayer       ); */
  /*   eta    -> push_back ( global.eta() ); */
  /*   phi    -> push_back ( global.phi() ); */
  /*   posx   -> push_back ( global.x()   ); */
  /*   posy   -> push_back ( global.y()   ); */
  /*   posz   -> push_back ( global.z()   ); */

  /* } */
    
};

#endif
