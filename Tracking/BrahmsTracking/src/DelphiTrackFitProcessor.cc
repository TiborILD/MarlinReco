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
#include <gear/GearParameters.h>

#include <ILDDetectorIDs.h>

#include "MarlinDelphiTrack.h"
#include "MaterialDB_F77.hh"

#include <TrackPropagators.h>


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

  streamlog_out(DEBUG4) << "   processing event: " << _n_evt 
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
		
	IMarlinTrack* marlin_trk = new MarlinDelphiTrack(track) ; // here you can supply any track which has hits attached
	
	bool fit_success = marlin_trk->fit( false ) ; // SJA:FIXME: false means from out to in here i.e. backwards. This would be better if had a more meaningful name perhaps fit_fwd and fit_rev

	if( fit_success ){ // if the fit converged 

	  TrackImpl* refittedTrack = marlin_trk->getIPFit() ;

	  float point[3] ;
	  point[0] = point[1] = point[2] = 0.0 ; // IP

	  float r_cylinder = 300.0 ;

	  // SJA:FIXME: these need to be const so that no body can change them ... 
	  //	  TrackImpl* nearestFit = marlin_trk->getNearestFitToPoint(point) ;	
	  TrackImpl* nearestFit = marlin_trk->getNearestFitToCylinder(r_cylinder) ;	

	  if(nearestFit){
	   
	    //  get ref and cov matrix for printing 
	    const float* ref = nearestFit->getReferencePoint() ;
	    EVENT::FloatVec cov = nearestFit->getCovMatrix() ;
	    //	    streamlog_out( DEBUG ) << " MarlinDelphi track parameters: nearest to " << point[0] << " " << point[1] <<  " " << point[2] << " : "
	    streamlog_out( DEBUG4 ) << " MarlinDelphi track parameters: nearest to cylinder r = " << r_cylinder << " : "
				   << " chi2/ndf " <<  nearestFit->getChi2() /  nearestFit->getNdf()  
				   << " chi2 " <<  nearestFit->getChi2() << std::endl 
	      
				   << "\t D0 "          <<  nearestFit->getD0()         <<  "[+/-" << sqrt( cov[0] ) << "] " 
				   << "\t Phi :"        <<  nearestFit->getPhi()        <<  "[+/-" << sqrt( cov[2] ) << "] " 
				   << "\t Omega "       <<  nearestFit->getOmega()      <<  "[+/-" << sqrt( cov[5] ) << "] " 
				   << "\t Z0 "          <<  nearestFit->getZ0()         <<  "[+/-" << sqrt( cov[9] ) << "] " 
				   << "\t tan(Lambda) " <<  nearestFit->getTanLambda()  <<  "[+/-" << sqrt( cov[14]) << "] " 
				    << "\t ref : [" << ref[0] << ", " << ref[1] << ", "  << ref[2] << "] radius ref = " << sqrt(ref[0]*ref[0]+ref[1]*ref[1])
				   << std::endl ;
	    
	    streamlog_out( DEBUG4 ) << "Move Track to SIT layers 2 and 1" << std::endl ;

	    const gear::GearParameters& pSITDet = Global::GEAR->getGearParameters("SIT") ; 

	    int nSITR = int(pSITDet.getDoubleVals("SITLayerRadius").size()) ;
	    int nLayersSIT = 0 ;

	    nLayersSIT = nSITR ;

	    for (int iL=nLayersSIT-1;iL>-1;--iL) {

	      double layer_rad = float(pSITDet.getDoubleVals("SITLayerRadius")[iL]) ;

	      // Now propagate the fit to the SIT layers using the TrackPropagators
	      IMPL::TrackImpl* trk_at_SIT = TrackPropagators::PropagateLCIOToCylinder(nearestFit, layer_rad, 0.0, 0.0, 0) ;

	      EVENT::FloatVec cov_at_SIT = trk_at_SIT->getCovMatrix() ;

	      if(trk_at_SIT){
		streamlog_out( DEBUG4 ) << "Track At SIT layer "<< iL << std::endl
					<< "\t D0 "          <<  trk_at_SIT->getD0()           <<  "[+/-" << sqrt( cov_at_SIT[0] ) << "] "             
					<< "\t Phi :"        <<  trk_at_SIT->getPhi()          <<  "[+/-" << sqrt( cov_at_SIT[2] ) << "] "  
					<< "\t Omega "       <<  trk_at_SIT->getOmega()        <<  "[+/-" << sqrt( cov_at_SIT[5] ) << "] "  
					<< "\t Z0 "          <<  trk_at_SIT->getZ0()           <<  "[+/-" << sqrt( cov_at_SIT[9] ) << "] "  
					<< "\t tan(Lambda) " <<  trk_at_SIT->getTanLambda()    <<  "[+/-" << sqrt( cov_at_SIT[14]) << "] "  
					<< "\t ref : [" << trk_at_SIT->getReferencePoint()[0] << ", " << trk_at_SIT->getReferencePoint()[1] << ", "  << trk_at_SIT->getReferencePoint()[2] 
					<< "] radius ref = " << sqrt(trk_at_SIT->getReferencePoint()[0]*trk_at_SIT->getReferencePoint()[0]+trk_at_SIT->getReferencePoint()[1]*trk_at_SIT->getReferencePoint()[1])
					<< std::endl ;
		delete trk_at_SIT ;
	      }
	      else{ 
		streamlog_out( DEBUG4 ) << "Track At SIT layer "<< iL << ": No intersection found" << std::endl; 
	      } 
	    }

	    // You can also propagate to a vtx ladder i.e. a plane parallel to the z-axis using
	    // Propagate track to a new reference point taken as its crossing point with a plane parallel to the z axis, containing points x1,x2 and y1,y2. Tolerance for intersection determined by epsilon.
	    // For direction== 0  the closest crossing point will be taken
	    // For direction== 1  the first crossing traversing in positive s will be taken
	    // For direction==-1  the first crossing traversing in negative s will be taken
	    // IMPL::TrackImpl* trk_at_VXD = TrackPropagators::PropagateLCIOToPlaneParralelToZ(nearestFit, x1, y1, x2, y2, 0) ;
	    
	    // if we want to check the distance to a hit and use the errors from the covariance matrix we will need to convert the d0 parameter into x and y with errors ... 
	    // then we would get the hit and make the comparison ... 

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

	    // here we could add more hits to the track and create another track using 
	    // IMarlinTrack* marlin_trk = new MarlinDelphiTrack(copyTrk) ;
	    // it would be nice to extend IMarlinTrack so that it had an add hit and fit method, but that will take a bit of time ... 
	    // so for now we should just create a new track each time we add a hit ... 
	    
	    //	    bool addtrack = false;
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
	  
	  //	    bool addtrack = true;
	  bool addtrack = false;
	  if( addtrack ) {
	    trackVec->addElement( refittedTrack ) ; // add it to the event ... 
	  } else {
	    delete refittedTrack ; // or delete it 
	  } 
	  
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



