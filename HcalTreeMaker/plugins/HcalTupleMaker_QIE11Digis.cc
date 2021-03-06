#include <iostream>
#include <ostream>
#include <string>
#include "HcalPromptAnalysis/HcalTreeMaker/interface/HcalTupleMaker_QIE11Digis.h"
// #include "DataFormats/HcalDigi/interface/HcalLaserDigi.h"
// #include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"


// NEEDS UPDATING
double adc2fC_QIE11[256]={
  // - - - - - - - range 0 - - - - - - - -
  //subrange0
  1.58, 4.73, 7.88, 11.0, 14.2, 17.3, 20.5, 23.6,
  26.8, 29.9, 33.1, 36.2, 39.4, 42.5, 45.7, 48.8,
  //subrange1
  53.6, 60.1, 66.6, 73.0, 79.5, 86.0, 92.5, 98.9,
  105, 112, 118, 125, 131, 138, 144, 151,
  //subrange2
  157, 164, 170, 177, 186, 199, 212, 225,
  238, 251, 264, 277, 289, 302, 315, 328,
  //subrange3
  341, 354, 367, 380, 393, 406, 418, 431,
  444, 464, 490, 516, 542, 568, 594, 620,

  // - - - - - - - range 1 - - - - - - - -
  //subrange0
  569, 594, 619, 645, 670, 695, 720, 745,
  771, 796, 821, 846, 871, 897, 922, 947,
  //subrange1
  960, 1010, 1060, 1120, 1170, 1220, 1270, 1320,
  1370, 1430, 1480, 1530, 1580, 1630, 1690, 1740,
  //subrange2
  1790, 1840, 1890, 1940,  2020, 2120, 2230, 2330,
  2430, 2540, 2640, 2740, 2850, 2950, 3050, 3150,
  //subrange3
  3260, 3360, 3460, 3570, 3670, 3770, 3880, 3980,
  4080, 4240, 4450, 4650, 4860, 5070, 5280, 5490,

  // - - - - - - - range 2 - - - - - - - -
  //subrange0
  5080, 5280, 5480, 5680, 5880, 6080, 6280, 6480,
  6680, 6890, 7090, 7290, 7490, 7690, 7890, 8090,
  //subrange1
  8400, 8810, 9220, 9630, 10000, 10400, 10900, 11300,
  11700, 12100, 12500, 12900, 13300, 13700, 14100, 14500,
  //subrange2
  15000, 15400, 15800, 16200, 16800, 17600, 18400, 19300,
  20100, 20900, 21700, 22500, 23400, 24200, 25000, 25800,
  //subrange3
  26600, 27500, 28300, 29100, 29900, 30700, 31600, 32400,
  33200, 34400, 36100, 37700, 39400, 41000, 42700, 44300,

  // - - - - - - - range 3 - - - - - - - - -
  //subrange0
  41100, 42700, 44300, 45900, 47600, 49200, 50800, 52500,
  54100, 55700, 57400, 59000, 60600, 62200, 63900, 65500,
  //subrange1
  68000, 71300, 74700, 78000, 81400, 84700, 88000, 91400,
  94700, 98100, 101000, 105000, 108000, 111000, 115000, 118000,
  //subrange2
  121000, 125000, 128000, 131000, 137000, 145000, 152000, 160000,
  168000, 176000, 183000, 191000, 199000, 206000, 214000, 222000,
  //subrange3
  230000, 237000, 245000, 253000, 261000, 268000, 276000, 284000,
  291000, 302000, 316000, 329000, 343000, 356000, 370000, 384000

};

HcalTupleMaker_QIE11Digis::HcalTupleMaker_QIE11Digis(const edm::ParameterSet& iConfig):
  prefix          (iConfig.getUntrackedParameter<std::string>("Prefix")),
  suffix          (iConfig.getUntrackedParameter<std::string>("Suffix")),
  storelaser      (iConfig.getUntrackedParameter<bool>("StoreLaser")),
  _taguMNio       (iConfig.getUntrackedParameter<edm::InputTag>("taguMNio",edm::InputTag("hcalDigis"))),
  m_qie11DigisTag (iConfig.getUntrackedParameter<edm::InputTag>("tagQIE11", edm::InputTag("hcalDigis")))
{ 
  qie11digisToken_ = consumes<HcalDataFrameContainer<QIE11DataFrame> >(m_qie11DigisTag);
 
  if (storelaser) {
    std::cout << "Storing uMNio laser informaiton" << std::endl;
    _tokuMNio = consumes<HcalUMNioDigi>(_taguMNio);
  }
  
  produces<std::vector<int>   >                  ( "QIE11DigiIEta"       );
  produces<std::vector<int>   >                  ( "QIE11DigiIPhi"       );
  produces<std::vector<int>   >                  ( "QIE11DigiSubdet"     );
  produces<std::vector<int>   >                  ( "QIE11DigiDepth"      );
  produces<std::vector<int>   >                  ( "QIE11DigiRawID"      );
  produces<std::vector<int>   >                  ( "QIE11DigiLinkError"  );
  produces<std::vector<int>   >                  ( "QIE11DigiCapIDError" );
  produces<std::vector<int>   >                  ( "QIE11DigiFlags"      );
  produces<std::vector<int>   >                  ( "QIE11DigiSOI"        );
  produces<std::vector<std::vector<int>   > >    ( "QIE11DigiADC"        );
  produces<std::vector<std::vector<double>   > > ( "QIE11DigiFC"         );
  produces<std::vector<std::vector<double>   > > ( "QIE11DigiRawFC"      );
  produces<std::vector<std::vector<double>   > > ( "QIE11DigiRC"         );
  produces<std::vector<std::vector<double>   > > ( "QIE11DigiGain"       );
  produces<std::vector<std::vector<int>   > >    ( "QIE11DigiTDC"        );
  produces<std::vector<std::vector<int>   > >    ( "QIE11DigiCapID"      );
  produces <int>                                 ( "laserType"           );
}

void HcalTupleMaker_QIE11Digis::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::unique_ptr<std::vector<int> >                    ieta    ( new std::vector<int>   ());
  std::unique_ptr<std::vector<int> >                    iphi    ( new std::vector<int>   ());
  std::unique_ptr<std::vector<int> >                    subdet  ( new std::vector<int>   ());
  std::unique_ptr<std::vector<int> >                    depth   ( new std::vector<int>   ());
  std::unique_ptr<std::vector<int> >                    rawId   ( new std::vector<int>   ());
  std::unique_ptr<std::vector<int> >                    linkEr  ( new std::vector<int>   ());
  std::unique_ptr<std::vector<int> >                    capidEr ( new std::vector<int>   ());
  std::unique_ptr<std::vector<int> >                    flags   ( new std::vector<int>   ());
  // std::unique_ptr<int>                                  lasertype (new int() );
  std::unique_ptr<std::vector<int> >                    soi     ( new std::vector<int>   ());
  std::unique_ptr<std::vector<std::vector<int  > > >    adc     ( new std::vector<std::vector<int  > >    ());
  std::unique_ptr<std::vector<std::vector<double  > > > rawfc   ( new std::vector<std::vector<double  > > ());
  std::unique_ptr<std::vector<std::vector<double  > > > fc      ( new std::vector<std::vector<double  > > ());
  std::unique_ptr<std::vector<std::vector<double  > > > gain    ( new std::vector<std::vector<double  > > ());
  std::unique_ptr<std::vector<std::vector<double  > > > respcorr( new std::vector<std::vector<double  > > ());
  std::unique_ptr<std::vector<std::vector<int  > > >    tdc     ( new std::vector<std::vector<int  > >    ());
  std::unique_ptr<std::vector<std::vector<int  > > >    capid   ( new std::vector<std::vector<int  > >    ());
  
  bool use_event = true;
  
  edm::Handle<HcalDataFrameContainer<QIE11DataFrame> >  qie11Digis;
  bool gotqie11digis = iEvent.getByToken(qie11digisToken_, qie11Digis);
  
  HcalCalibrations calibrations;
  CaloSamples tool;

  iSetup.get<HcalDbRecord > ().get(conditions);
  
  if (!gotqie11digis ) {
	std::cout << "Could not find QIE11 digis " <<  m_qie11DigisTag << std::endl;
	use_event = false;
  }
   
  if (use_event) {

    //
    // Laser related
    // 
    if (storelaser) {
      edm::Handle<HcalUMNioDigi> cumnio;
      std::cout << "Only using laser events" << std::endl;
      bool gotuMNio = iEvent.getByToken(_tokuMNio,cumnio);                           
      if (!gotuMNio ) {
	std::cout << "Could not find uMNio " << _taguMNio << std::endl;
	use_event = false;
      }
      std::unique_ptr<int> lasertype (new int(cumnio -> valueUserWord(0)));
      iEvent.put(move( lasertype )         , "laserType"      );
      //std::cout << "Laser type is " << lasertype << std::endl;
    }
    else {
      std::unique_ptr<int> lasertype (new int());
      iEvent.put(move( lasertype )         , "laserType"      );
    }

    //
    // QIE11
    //    
    for (uint32_t i=0; i<qie11Digis->size(); i++) {
      // From: https://github.com/awhitbeck/HFcommissioningAnalysis/blob/b3456c9fe66ef9bcc6c54773d60f768c269a5c74/src/HFanalyzer.cc#L429
      QIE11DataFrame qie11df = (*qie11Digis)[i];
      // QIE10 structure: static_cast<QIE11DataFrame>((*qie11Digis)[i]);

      // Extract info on detector location
      DetId detid = qie11df.detid();
      HcalDetId hcaldetid = HcalDetId(detid);

      int sub = hcaldetid.subdet();
      if (sub!=1 && sub!=2) continue;
      //KH std::cout << "check1: " << sub << std::endl;
      HcalCalibrations calibrations = conditions->getHcalCalibrations(detid);
      const HcalQIECoder* channelCoder = conditions->getHcalCoder(detid);
      const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
      HcalCoderDb coder(*channelCoder, *shape);
      coder.adc2fC(qie11df, tool);
      int soi_digi = tool.presamples();
      //int lastbin = tool.size() - 1;
      //KH std::cout << "check2" << std::endl;

      ieta    -> push_back ( hcaldetid.ieta()        );
      iphi    -> push_back ( hcaldetid.iphi()        );
      subdet  -> push_back ( 8/*hcaldetid.subdet()*/ );
      depth   -> push_back ( hcaldetid.depth()       );
      rawId   -> push_back ( hcaldetid.rawId()       );
      linkEr  -> push_back ( qie11df.linkError()     );
      capidEr -> push_back ( qie11df.capidError()    );
      flags   -> push_back ( qie11df.flags()         );
      soi     -> push_back ( soi_digi );
      
      if (0) {
	std::cout << "Printing raw dataframe" << std::endl;
	std::cout << qie11df << std::endl;
	std::cout << "Printing content of samples() method" << std::endl;
	std::cout << qie11df.samples() << std::endl;
      }
      
      //KH soi   -> push_back ( std::vector<int  >   () ) ;
      adc   -> push_back ( std::vector<int  >   () ) ;
      rawfc -> push_back ( std::vector<double  >() ) ;
      fc    -> push_back ( std::vector<double  >() ) ;
      gain  -> push_back ( std::vector<double  >() ) ;
      respcorr -> push_back ( std::vector<double  >() ) ;
      tdc   -> push_back ( std::vector<int  >   () ) ;
      capid -> push_back ( std::vector<int  >   () ) ;
      size_t last_entry = adc -> size() - 1;

      // TS
      int nTS = qie11df.samples();
      
      //KH const bool saveEffectivePeds = channelInfo->hasEffectivePedestals();
      
      for (int its=0; its<nTS; ++its) {
	//KH (*soi  )[last_entry].push_back ( qie11df[its].soi()               ); // soi is a bool, but stored as an int
	(*adc  )[last_entry].push_back ( qie11df[its].adc()               );
	//KH (*fc   )[last_entry].push_back ( adc2fC_QIE11[qie11df[its].adc()] );

	(*rawfc)[last_entry].push_back ( tool[its] );
	(*fc   )[last_entry].push_back ( tool[its]-calibrations.effpedestal(qie11df[its].capid()) );

	double respcorrgain = calibrations.respcorrgain(qie11df[its].capid());
	double rawgain = calibrations.rawgain(qie11df[its].capid());
	double rc = respcorrgain;
	if (rawgain!=0.) rc /= rawgain;
	//double rc = 0.;
	//double rawgain = 0.;
	(*respcorr)[last_entry].push_back ( rc  );
	(*gain)[last_entry].push_back ( rawgain );

	(*tdc  )[last_entry].push_back ( qie11df[its].tdc()                );
	(*capid)[last_entry].push_back ( qie11df[its].capid()              );	

	/*
	std::cout << tool[its] << " "
		  << calibrations.effpedestal(qie11df[its].capid()) << " "
		  << calibrations.pedestal(qie11df[its].capid()) << " "
		  << rawgain << " " << rc << std::endl;
	*/

      }
      
    }
    
    iEvent.put(move( ieta   )       , "QIE11DigiIEta"      );
    iEvent.put(move( iphi   )       , "QIE11DigiIPhi"      );
    iEvent.put(move( subdet )       , "QIE11DigiSubdet"    );
    iEvent.put(move( depth  )       , "QIE11DigiDepth"     );
    iEvent.put(move( rawId  )       , "QIE11DigiRawID"     );
    iEvent.put(move( linkEr )       , "QIE11DigiLinkError" );
    iEvent.put(move( capidEr)       , "QIE11DigiCapIDError");
    iEvent.put(move( flags  )       , "QIE11DigiFlags"     );
    iEvent.put(move( soi    )       , "QIE11DigiSOI"       );
    iEvent.put(move( adc    )       , "QIE11DigiADC"       );
    iEvent.put(move( rawfc  )       , "QIE11DigiRawFC"     );
    iEvent.put(move( fc     )       , "QIE11DigiFC"        );
    iEvent.put(move( gain   )       , "QIE11DigiGain"      );
    iEvent.put(move( respcorr )     , "QIE11DigiRC"        );
    iEvent.put(move( tdc    )       , "QIE11DigiTDC"       );
    iEvent.put(move( capid  )       , "QIE11DigiCapID"     );
  }
}
