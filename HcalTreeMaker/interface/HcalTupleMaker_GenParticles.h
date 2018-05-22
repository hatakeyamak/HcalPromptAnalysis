#ifndef HcalTupleMaker_GenParticles_h
#define HcalTupleMaker_GenParticles_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

class HcalTupleMaker_GenParticles : public edm::EDProducer {
 public:
  explicit HcalTupleMaker_GenParticles(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genCollectionToken_;
  edm::EDGetTokenT<edm::HepMCProduct> HepMCToken_;
  bool debug;
  bool bool_HepMCProduct;

};

#endif
