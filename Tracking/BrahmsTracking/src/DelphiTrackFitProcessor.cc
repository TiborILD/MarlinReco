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
  MaterialDB_F77::Instance().initialise(_MSOn);
  //  _kaltest = new MarlinKalTest( *marlin::Global::GEAR ) ;

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

	  MarlinDelphiTrack* delphi_trk = static_cast<MarlinDelphiTrack*>(marlin_trk);

	  std::vector<TanagraFit const*> fits ;

	  delphi_trk->getTanagraFits(fits) ;

	  for( unsigned int iext=0; iext < fits.size(); ++iext) {

	    float cov[15] ;
	    fits[iext]->get_cov(cov);
	    
	    streamlog_out( MESSAGE ) << "Track Parameters for surface " << iext << " at R = " << fits[iext]->get_r() 
				   << "    rphi "        <<  fits[iext]->get_rphi()  <<  "[+/-" << sqrt( cov[0] ) << "] " 
				   << "    z "           <<  fits[iext]->get_z()  <<  "[+/-" << sqrt( cov[2] ) << "] " 
				   << "    theta "       <<  fits[iext]->get_theta()  <<  "[+/-" << sqrt( cov[5] ) << "] " 
				   << "    phi "         <<  fits[iext]->get_phi()  <<  "[+/-" << sqrt( cov[9]) << "] " 
				   << "    invP "        <<  fits[iext]->get_invp()  <<  "[+/-" << sqrt( cov[14]) << "] " 
				   <<  std::endl ;
	  }
	


	  //	//SJA:FIXME: This has to go away. The use of hardcoded number here is completely error prone ...
	  refittedTrack->subdetectorHitNumbers().resize(12);
	  for ( unsigned int detIndex = 0 ;  detIndex < refittedTrack->subdetectorHitNumbers().size() ; detIndex++ ) 
	    {
	      refittedTrack->subdetectorHitNumbers()[detIndex] = track->getSubdetectorHitNumbers()[detIndex] ;
	    }
	  
	  trackVec->addElement( refittedTrack );

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



