#include "TruthTracker.h"
#include <iostream>

#include <vector>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>

#include <IMPL/TrackImpl.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"

#include <ILDDetectorIDs.h>

using namespace lcio ;
using namespace marlin ;


bool TrackerHitSortPredicate(const TrackerHit* hit1, const TrackerHit* hit2)
{
  // TrackerHits must be checked beforehand that they have one and only one SimTrackerHit
  SimTrackerHit * simHit1 = dynamic_cast<SimTrackerHit*>(hit1->getRawHits().at(0)) ;
  SimTrackerHit * simHit2 = dynamic_cast<SimTrackerHit*>(hit2->getRawHits().at(0)) ;

  if( simHit1->getMCParticle() == simHit2->getMCParticle() ) {
    return simHit1->getTime() < simHit2->getTime() ;
  }
  else { 
    return simHit1->getMCParticle() < simHit2->getMCParticle() ;
  }
}


TruthTracker aTruthTracker ;


TruthTracker::TruthTracker() : Processor("TruthTracker") {
  
  // modify processor description
  _description = "Creates Track Collection from MC Truth" ;
  
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::TRACKERHIT,
			   "VTXCollectionName" , 
			   "Name of the collection of VTX tracker hits"  ,
			   _colNameTrkHitsVTX ,
			   std::string("VTXTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
			   "FTDCollectionName" , 
			   "Name of the collection of FTD tracker hits"  ,
			   _colNameTrkHitsFTD ,
			   std::string("FTDTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
			   "SITCollectionName" , 
			   "Name of the collection of SIT tracker hits"  ,
			   _colNameTrkHitsSIT ,
			   std::string("SITTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
			   "TPCCollectionName" , 
			   "Name of the collection of TPC tracker hits"  ,
			   _colNameTrkHitsTPC ,
			   std::string("TPCTrackerHits") ) ;

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "VTXSimTrkhitCollectionName" ,
			   "Name of the VTX sim trk hit collection" ,
			   _colNameSimTrkHitsVTX ,
			   std::string("VXDCollection") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "SITSimTrkhitCollectionName" ,
			   "Name of the SIT sim trk hit collection" ,
			   _colNameSimTrkHitsSIT ,
			   std::string("SITCollection") ) ;

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "FTDSimTrkhitCollectionName" ,
			   "Name of the FTD sim trk hit collection" ,
			   _colNameSimTrkHitsFTD ,
			   std::string("FTDCollection") ) ;

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "TPCSimTrkhitCollectionName" , 
			   "Name of the TPC sim trk hit collection"  ,
			   _colNameSimTrkHitsTPC,
			   std::string("TPCCollection") ) ;

  registerOutputCollection( LCIO::TRACK,
			   "OutputTrackCollectionName" , 
			   "Name of the output track collection"  ,
			   _output_track_col_name ,
			   std::string("TruthTracks") ) ;

  registerOutputCollection( LCIO::LCRELATION,
			   "OutputTrackRelCollection" , 
			   "Name of the MCParticle-Track Relations collection for output tracks"  ,
			   _output_track_rel_name ,
			   std::string("TruthTracksMCP") ) ;

 registerProcessorParameter( "nEventPrintout",
			      "Print out progress every N events ",
			      _nEventPrintout,
			      int(1000));

 registerProcessorParameter( "MCpThreshold",
			      "Momentume Threshold MC particles which will produce tracks GeV",
			      _MCpThreshold,
			      float(0.5));


 _n_run = 0 ;
 _n_evt = 0 ;
 
  _colMCP = NULL ;
  _nMCP = 0 ;
  
 _VTXSimHits = NULL ;
 _SITSimHits = NULL ;
 _FTDSimHits = NULL ;
 _TPCSimHits = NULL ;
 
 _nVTXSimHits = 0 ;
 _nSITSimHits = 0 ;
 _nFTDSimHits = 0 ;
 _nTPCSimHits = 0 ;
 
 _VTXTrkHits = NULL ;
 _SITTrkHits = NULL ;
 _FTDTrkHits = NULL ;
 _TPCTrkHits = NULL ;
 
 _nVTXTrkHits = 0 ;
 _nSITTrkHits = 0 ;
 _nFTDTrkHits = 0 ;
 _nTPCTrkHits = 0 ;

}


void TruthTracker::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
    
  // usually a good idea to
  printParameters() ;

 

}

void TruthTracker::processRunHeader( LCRunHeader* run) { 

  ++_n_run ;
} 

void TruthTracker::processEvent( LCEvent * evt ) { 

    
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  if( _n_evt % _nEventPrintout == 0 ) streamlog_out(DEBUG) << "   processing event: " << _n_evt 
		       << std::endl ;

  // get the input collections and set the memeber variable pointers to them
  this->SetupInputCollections(evt) ;


  // establish the track collection that will be created 
  LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    
  
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackVec->setFlag( trkFlag.getFlag()  ) ;
  
  // establish the track relations collection that will be created 
  LCCollectionVec* trackRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;


  // loop over all tracker hits and check which MCPartilce created them 
  // check each TrackerHit to see that it has only one SimTrackerHit
  std::vector<TrackerHit*> trkhits; 

  for(int i=0; i< _nVTXTrkHits ; ++i)
    {
      TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_VTXTrkHits->getElementAt( i ));      
      if( trkhit->getRawHits().size()==1 ) {
	trkhits.push_back(trkhit);
      }
    }

  for(int i=0; i< _nSITTrkHits ; ++i)
    {
      TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_SITTrkHits->getElementAt( i ));      
      if( trkhit->getRawHits().size()==1 ) {
	trkhits.push_back(trkhit);
      }
    }

  for(int i=0; i< _nFTDTrkHits ; ++i)
    {
      TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_FTDTrkHits->getElementAt( i ));      
      if( trkhit->getRawHits().size()==1 ) {
	trkhits.push_back(trkhit);
      }
    }
  
  for(int i=0; i< _nTPCTrkHits ; ++i)
    {
      TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_TPCTrkHits->getElementAt( i ));      
      if( trkhit->getRawHits().size()==1 ) {
	trkhits.push_back(trkhit);
      }
    }
  

  streamlog_out( DEBUG4 ) << "Number of Tracker hits = " << trkhits.size() << std::endl;     

  for(unsigned int i=0; i< trkhits.size() ; ++i)
    {
      SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*>(trkhits[i]->getRawHits().at(0)) ;
      streamlog_out( DEBUG1 ) << "Tracker hit: = " << trkhits[1] << "  mcp = " << simHit->getMCParticle() << " time = " << simHit->getTime() << std::endl;     
    }

  // now order the hits by MCParticle and then in time  
  std::sort(trkhits.begin(), trkhits.end(), TrackerHitSortPredicate);

  int iTrack(0) ;
  
  streamlog_out( DEBUG4 ) << "Add Tracker hits to Tracks" << std::endl;
  
  if( trkhits.size() > 0) {
    MCParticle* mcplast = NULL;
    
    std::vector<TrackerHit*> hit_list;
    
    for(unsigned int i=0; i< trkhits.size() ; ++i)
      {

	MCParticle* mcp = dynamic_cast<SimTrackerHit*>(trkhits.at(i)->getRawHits().at(0))->getMCParticle();
	double const* p    = mcp->getMomentum() ;
	float  const pmag2 = p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ; 

	if( (mcp != mcplast || i == (trkhits.size()-1)) && hit_list.size() > 2 )
	  { // create track from vector of hits
	    TrackImpl* Track = new TrackImpl ; 
	    
	    std::map<int, int> hitNumbers; 
	    
	    const int theVXD = ILDDetectorIDs::DetID::VXD;
	    const int theSIT = ILDDetectorIDs::DetID::SIT;
	    const int theFTD = ILDDetectorIDs::DetID::FTD;
	    const int theTPC = ILDDetectorIDs::DetID::TPC;
	    
	    hitNumbers[theVXD] = 0;
	    hitNumbers[theSIT] = 0;
	    hitNumbers[theFTD] = 0;
	    hitNumbers[theTPC] = 0;
	    
	    streamlog_out( DEBUG4 ) << "Set the hitcounter for VXD code = " << theVXD << std::endl;     
	    streamlog_out( DEBUG4 ) << "Set the hitcounter for SIT code = " << theSIT << std::endl;     
	    streamlog_out( DEBUG4 ) << "Set the hitcounter for FTD code = " << theFTD<< std::endl;     
	    streamlog_out( DEBUG4 ) << "Set the hitcounter for TPC code = " << theTPC << std::endl;     
	    
	    for(unsigned int j=0; j<hit_list.size(); ++j)
	      {
		Track->addHit(hit_list.at(j)) ;
		unsigned int layerID =  hit_list.at(j)->ext<ILDDetectorIDs::HitInfo>()->layerID;
		++hitNumbers[(layerID / ILDDetectorIDs::DetID::Factor)];
		streamlog_out( DEBUG1 ) << "Hit from Detector " << (layerID / ILDDetectorIDs::DetID::Factor) << std::endl;     
	      }
	    
	    Track->subdetectorHitNumbers().resize(12);
	    Track->subdetectorHitNumbers()[0] = hitNumbers[theVXD];
	    Track->subdetectorHitNumbers()[1] = hitNumbers[theFTD];
	    Track->subdetectorHitNumbers()[2] = hitNumbers[theSIT];
	    Track->subdetectorHitNumbers()[3] = hitNumbers[theTPC];
	    Track->subdetectorHitNumbers()[4] = int(0);
	    Track->subdetectorHitNumbers()[5] = int(0);
	    Track->subdetectorHitNumbers()[6] = hitNumbers[theVXD];
	    Track->subdetectorHitNumbers()[7] = hitNumbers[theFTD];
	    Track->subdetectorHitNumbers()[8] = hitNumbers[theSIT];
	    Track->subdetectorHitNumbers()[9] = hitNumbers[theTPC];
	    Track->subdetectorHitNumbers()[10] = int(0);
	    Track->subdetectorHitNumbers()[11] = int(0);
	    
	    hit_list.clear();	
	    streamlog_out( DEBUG4 ) << "Add Track " << iTrack++ << "  processed " << i << " hits" << std::endl;     
	    trackVec->addElement( Track );
	  }

	else
	  {
	    if( pmag2  > (_MCpThreshold*_MCpThreshold) )
	      {
		hit_list.push_back(trkhits.at(i)) ;
		mcplast = mcp;
	      }
	  }
      }
  }    
    
  evt->addCollection( trackVec , _output_track_col_name) ;
  evt->addCollection( trackRelVec , _output_track_rel_name) ;

  
  ++_n_evt ;
}



void TruthTracker::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TruthTracker::end(){ 
  
  streamlog_out(DEBUG) << "TruthTracker::end()  " << name() 
 	    << " processed " << _n_evt << " events in " << _n_run << " runs "
 	    << std::endl ;

}

LCCollection* TruthTracker::GetCollection(  LCEvent * evt, std::string colName ){

  LCCollection* col = NULL;
  
  int nElements = 0;

  try{
    col = evt->getCollection( colName.c_str() ) ;
    nElements = col->getNumberOfElements()  ;
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() << " track collection found in event = " << col << " number of elements " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() <<  " collection absent in event" << std::endl;     
  }

  return col; 

}

LCRelationNavigator* TruthTracker::GetRelations( EVENT::LCEvent * evt , std::string RelName ) {

  LCRelationNavigator* nav = NULL ;

  try{
    nav = new LCRelationNavigator(evt->getCollection( RelName.c_str() ));
    streamlog_out( DEBUG2 ) << "TruthTracker --> " << RelName << " track relation collection in event = " << nav << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG2 ) << "TruthTracker --> " << RelName.c_str() << " track relation collection absent in event" << std::endl;     
  }
  
  return nav;
  
}


void TruthTracker::SetupInputCollections( LCEvent * evt ){

  _colMCP = GetCollection( evt, _colNameMC) ;
  streamlog_out( DEBUG2 ) << "number of MCParticles: " << _nMCP << std::endl ;

  _VTXSimHits = GetCollection( evt, _colNameSimTrkHitsVTX ) ;
  _SITSimHits = GetCollection( evt, _colNameSimTrkHitsSIT ) ;
  _FTDSimHits = GetCollection( evt, _colNameSimTrkHitsFTD ) ;
  _TPCSimHits = GetCollection( evt, _colNameSimTrkHitsTPC ) ;
  
  if(_VTXSimHits) _nVTXSimHits = _VTXSimHits->getNumberOfElements() ; 
  if(_SITSimHits) _nSITSimHits = _SITSimHits->getNumberOfElements() ;
  if(_FTDSimHits) _nFTDSimHits = _FTDSimHits->getNumberOfElements() ;
  if(_TPCSimHits) _nTPCSimHits = _TPCSimHits->getNumberOfElements() ;
  
  streamlog_out( DEBUG2 ) << "number of VTXSimHits: " << _nVTXSimHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of SITSimHits: " << _nSITSimHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of FTDSimHits: " << _nFTDSimHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of TPCSimHits: " << _nTPCSimHits << std::endl ;

  _VTXTrkHits = GetCollection( evt, _colNameTrkHitsVTX ) ;
  _SITTrkHits = GetCollection( evt, _colNameTrkHitsSIT ) ;
  _FTDTrkHits = GetCollection( evt, _colNameTrkHitsFTD ) ;
  _TPCTrkHits = GetCollection( evt, _colNameTrkHitsTPC ) ;
  
  if(_VTXTrkHits) _nVTXTrkHits = _VTXTrkHits->getNumberOfElements() ; 
  if(_SITTrkHits) _nSITTrkHits = _SITTrkHits->getNumberOfElements() ;
  if(_FTDTrkHits) _nFTDTrkHits = _FTDTrkHits->getNumberOfElements() ;
  if(_TPCTrkHits) _nTPCTrkHits = _TPCTrkHits->getNumberOfElements() ;
  
  streamlog_out( DEBUG2 ) << "number of VTXTrkHits: " << _nVTXTrkHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of SITTrkHits: " << _nSITTrkHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of FTDTrkHits: " << _nFTDTrkHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of TPCTrkHits: " << _nTPCTrkHits << std::endl ;


  if(_nVTXTrkHits) { 
    // check to whether the ID have been set 
    if( ! dynamic_cast<TrackerHit*>( _VTXTrkHits->getElementAt( 0 ))->ext<ILDDetectorIDs::HitInfo>() ) {      
    }
  }

  if(_nSITTrkHits) { 
    // check to whether the ID have been set 
    if( ! dynamic_cast<TrackerHit*>( _SITTrkHits->getElementAt( 0 ))->ext<ILDDetectorIDs::HitInfo>() ) {      
    }
  }

  if(_nFTDTrkHits) { 
    // check to whether the ID have been set 
    if( ! dynamic_cast<TrackerHit*>( _FTDTrkHits->getElementAt( 0 ))->ext<ILDDetectorIDs::HitInfo>() ) {      
    }
  }

  if(_nTPCTrkHits) { 
    // check to whether the ID have been set 
    if( ! dynamic_cast<TrackerHit*>( _TPCTrkHits->getElementAt( 0 ))->ext<ILDDetectorIDs::HitInfo>() ) {      
    }
  }
}
