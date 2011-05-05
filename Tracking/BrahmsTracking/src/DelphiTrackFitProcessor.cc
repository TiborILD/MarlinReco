#include "DelphiTrackFitProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <EVENT/Track.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"

#include <ILDDetectorIDs.h>

#include "MarlinDelphiTrack.h"
#include "MaterialDB_F77.hh"



using namespace lcio ;
using namespace marlin ;
using namespace marlin_delphiF77 ;

// needed for setting verbosity of F77 Delphi code
extern "C" {
  extern struct {
    int ideb ;
    int ihis ;
  } fkdebug_;
}

DelphiTrackFitProcessor aDelphiTrackFitProcessor ;


DelphiTrackFitProcessor::DelphiTrackFitProcessor() : Processor("DelphiTrackFitProcessor") {
  
  // modify processor description
  _description = "DelphiTrackFitProcessor refits an input track collection, producing a new collection of tracks." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::TRACK,
			   "InputTrackCollectionName" , 
			   "Name of the input track collection"  ,
			   _input_track_col_name ,
			   std::string("LDCTracks") ) ;

  registerOutputCollection( LCIO::TRACK,
			   "OutputTrackCollectionName" , 
			   "Name of the output track collection"  ,
			   _output_track_col_name ,
			   std::string("RefittedLDCTracks") ) ;

 registerProcessorParameter("MultipleScatteringOn",
			     "Use MultipleScattering in Fit",
			     _MSOn,
			     bool(true));

registerProcessorParameter("F77DebugLevel",
			     "Set the level of verbosity for the Delphi F77 code",
			     _F77DebugLevel,
			     int(0));

}


void DelphiTrackFitProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
    
  // usually a good idea to
  printParameters() ;

  // set up the geometery needed by F77
  MaterialDB_F77::Instance()->switchONMaterial();

  _n_run = 0 ;
  _n_evt = 0 ;

  fkdebug_.ideb = _F77DebugLevel ;
  fkdebug_.ihis = 0 ;
  
}

void DelphiTrackFitProcessor::processRunHeader( LCRunHeader* run) { 

  ++_n_run ;
} 

void DelphiTrackFitProcessor::processEvent( LCEvent * evt ) { 

    
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << _n_evt 
		       << std::endl ;
   
  // get input collection and relations 
  LCCollection* input_track_col = this->GetCollection( evt, _input_track_col_name ) ;

  if( input_track_col != 0 ){

    // establish the track collection that will be created 
    LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    

    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    trackVec->setFlag( trkFlag.getFlag()  ) ;

    int nTracks = input_track_col->getNumberOfElements()  ;

    // loop over the input tacks and refit using KalTest    
    for(int i=0; i< nTracks ; ++i)
      {
      
	Track* track = dynamic_cast<Track*>( input_track_col->getElementAt( i ) ) ;
		
	IMarlinTrack* marlin_trk = new MarlinDelphiTrack(track) ;
	
	bool fit_success = marlin_trk->fit( false ) ; // SJA:FIXME: false means from out to in here i.e. backwards. This would be better if had a more meaningful name perhaps fit_fwd and fit_rev

	if( fit_success ){ 

	  TrackImpl* refittedTrack = marlin_trk->getIPFit() ;

	  float point[3] ;
	  //	  point[0] = point[1] = point[2] = 0.0 ; // IP
	  point[0] = -1806. ;
	  point[1] = -13.0 ;
	  point[2] = 63.0 ; 

	  float r_cylinder = 500.0 ;

	  //	  TrackImpl* nearestFit = marlin_trk->getNearestFitToPoint(point) ;	
	  TrackImpl* nearestFit = marlin_trk->getNearestFitToCylinder(r_cylinder) ;	

	  if(nearestFit){
	    
	    //  get ref and cov matrix for printing 
	    const float* ref = nearestFit->getReferencePoint() ;
	    EVENT::FloatVec cov = nearestFit->getCovMatrix() ;
	    //	    streamlog_out( DEBUG ) << " MarlinDelphi track parameters: nearest to " << point[0] << " " << point[1] <<  " " << point[2] << " : "
	    streamlog_out( DEBUG ) << " MarlinDelphi track parameters: nearest to cylinder " << r_cylinder << " : "
				   << " chi2/ndf " <<  nearestFit->getChi2() /  nearestFit->getNdf()  
				   << " chi2 " <<  nearestFit->getChi2() << std::endl 
	      
				   << "\t D0 "          <<  nearestFit->getD0()         <<  "[+/-" << sqrt( cov[0] ) << "] " 
				   << "\t Phi :"        <<  nearestFit->getPhi()        <<  "[+/-" << sqrt( cov[2] ) << "] " 
				   << "\t Omega "       <<  nearestFit->getOmega()      <<  "[+/-" << sqrt( cov[5] ) << "] " 
				   << "\t Z0 "          <<  nearestFit->getZ0()         <<  "[+/-" << sqrt( cov[9] ) << "] " 
				   << "\t tan(Lambda) " <<  nearestFit->getTanLambda()  <<  "[+/-" << sqrt( cov[14]) << "] " 
				   << "\t ref : [" << ref[0] << ", " << ref[1] << ", "  << ref[2] << "]"     
				   << std::endl ;

	    // create NEW track to be returned. The responsiblitiy for deletion lies with the caller.
	    IMPL::TrackImpl* copyTrk = new IMPL::TrackImpl ;
	    
	    copyTrk->setD0( nearestFit->getD0() ) ;  
	    copyTrk->setPhi( nearestFit->getPhi() ) ;   
	    copyTrk->setOmega( nearestFit->getOmega() ) ;
	    copyTrk->setZ0( nearestFit->getZ0() ) ;  
	    copyTrk->setTanLambda( nearestFit->getTanLambda() ) ;  
	    
	    copyTrk->setChi2( nearestFit->getChi2() ) ;
	    copyTrk->setNdf( nearestFit->getNdf() ) ;
	    
	    copyTrk->setReferencePoint( const_cast<float*>(nearestFit->getReferencePoint()) ) ;
	    copyTrk->setIsReferencePointPCA(false) ;
	    
	    copyTrk->setCovMatrix( nearestFit->getCovMatrix() ) ;
	    
	    copyTrk->subdetectorHitNumbers().resize(12);
	    for ( unsigned int detIndex = 0 ;  detIndex < copyTrk->subdetectorHitNumbers().size() ; detIndex++ ) 
	      {
		copyTrk->subdetectorHitNumbers()[detIndex] = track->getSubdetectorHitNumbers()[detIndex] ;
	      }

	    // add the hits. Currently all hits are added and no bookeeping is made of which hits failed to be included in the fit.
	    EVENT::TrackerHitVec lcioHits = track->getTrackerHits();
	    EVENT::TrackerHitVec::iterator it = lcioHits.begin();
	    
	    for( it = lcioHits.begin() ; it != lcioHits.end() ; ++it )
	      { 
		copyTrk->addHit( *it ) ;
	      }
	    
	    //	    bool addtrack = true;
	    bool addtrack = true;
	    if( addtrack ) {
	      trackVec->addElement( copyTrk ) ; // add it to the event ... 
	    } else {
	      delete copyTrk ; // or delete it 
	    }

	  }
	  else {
	    streamlog_out( DEBUG ) << " MarlinDelphi track parameters not present for : " << point[0] << " " << point[1] <<  " " << point[2] << std::endl ;
	  }


	  //	//SJA:FIXME: This has to go away. The use of hardcoded number here is completely error prone ...
	  refittedTrack->subdetectorHitNumbers().resize(12);
	  for ( unsigned int detIndex = 0 ;  detIndex < refittedTrack->subdetectorHitNumbers().size() ; detIndex++ ) 
	    {
	      refittedTrack->subdetectorHitNumbers()[detIndex] = track->getSubdetectorHitNumbers()[detIndex] ;
	    }
	  
	  //	  trackVec->addElement( refittedTrack );
	  delete refittedTrack;
	}
	
	delete marlin_trk;
	
      }
    evt->addCollection( trackVec , _output_track_col_name) ;

  }  
  ++_n_evt ;
}



void DelphiTrackFitProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void DelphiTrackFitProcessor::end(){ 
  
  streamlog_out(DEBUG) << "DelphiTrackFitProcessor::end()  " << name() 
 	    << " processed " << _n_evt << " events in " << _n_run << " runs "
 	    << std::endl ;

}

LCCollection* DelphiTrackFitProcessor::GetCollection( LCEvent * evt, std::string colName ){

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



