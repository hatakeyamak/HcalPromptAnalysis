#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_GenParticles.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/ESHandle.h"

HcalTupleMaker_GenParticles::HcalTupleMaker_GenParticles(const edm::ParameterSet& iConfig):
  inputTag    (iConfig.getUntrackedParameter<edm::InputTag>("source")),
  prefix      (iConfig.getUntrackedParameter<std::string>  ("Prefix")),
  suffix      (iConfig.getUntrackedParameter<std::string>  ("Suffix")),
  bool_HepMCProduct (iConfig.getUntrackedParameter<bool>("HepMCProduct"))
{

  produces< std::vector< double > >(prefix + "Pt"  + suffix );
  produces< std::vector< double > >(prefix + "Eta" + suffix );
  produces< std::vector< double > >(prefix + "Phi" + suffix );
  produces< std::vector< double > >(prefix + "M"   + suffix );
  produces< std::vector< int > >(prefix + "PdgId" + suffix );
  produces< std::vector< int > >(prefix + "Status" + suffix );
  //produces< std::vector< int > >(prefix + "Parent" + suffix );
  //produces< std::vector< int > >(prefix + "ParentId" + suffix );

  genCollectionToken_ = consumes<edm::View<reco::GenParticle>>(inputTag);
  HepMCToken_ = consumes<edm::HepMCProduct>(inputTag);

}

void HcalTupleMaker_GenParticles::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<std::vector<double> >            pt                ( new std::vector<double>           ());
  std::unique_ptr<std::vector<double> >            eta               ( new std::vector<double>           ());
  std::unique_ptr<std::vector<double> >            phi               ( new std::vector<double>           ());
  std::unique_ptr<std::vector<double> >            mass              ( new std::vector<double>           ());

  std::unique_ptr<std::vector<int   > >            pdgid             ( new std::vector<int>              ());
  std::unique_ptr<std::vector<int   > >            status            ( new std::vector<int>              ());
  //std::unique_ptr<std::vector<int   > >            parent            ( new std::vector<int>              ());
  //std::unique_ptr<std::vector<int   > >            parentid          ( new std::vector<int>              ());

  if (!bool_HepMCProduct) {
    
    edm::Handle< edm::View<reco::GenParticle> > genPartCands;
    iEvent.getByToken(genCollectionToken_, genPartCands);

    for(edm::View<reco::GenParticle>::const_iterator iPart = genPartCands->begin(); iPart != genPartCands->end(); ++iPart)
      {
	
	pt->push_back(iPart->pt());
	eta->push_back(iPart->eta());
	phi->push_back(iPart->phi());
	mass->push_back(iPart->mass());

	pdgid->push_back(iPart->pdgId());
	status->push_back(iPart->status());
	
      }

  } else {

    edm::Handle<edm::HepMCProduct> evtMC;
    iEvent.getByToken(HepMCToken_,evtMC);

    const HepMC::GenEvent * myGenEvent = evtMC->GetEvent();
    for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();
	  p != myGenEvent->particles_end(); ++p ) {

      pt->push_back((*p)->momentum().perp());
      eta->push_back((*p)->momentum().eta());
      phi->push_back((*p)->momentum().phi());
      mass->push_back((*p)->momentum().m());
      pdgid->push_back((*p)->pdg_id());
      status->push_back((*p)->status());
      
    }

  }
    
  iEvent.put(move( pt              ) , prefix + "Pt"            + suffix );
  iEvent.put(move( eta             ) , prefix + "Eta"           + suffix );
  iEvent.put(move( phi             ) , prefix + "Phi"           + suffix );
  iEvent.put(move( mass            ) , prefix + "M"             + suffix );

  iEvent.put(move( pdgid           ) , prefix + "PdgId"         + suffix );
  iEvent.put(move( status          ) , prefix + "Status"        + suffix );
  //iEvent.put(move( parent          ) , prefix + "Parent"        + suffix );
  //iEvent.put(move( parentid        ) , prefix + "ParentId"      + suffix );

}
