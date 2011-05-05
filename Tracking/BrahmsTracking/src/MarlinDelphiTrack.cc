
#include "MarlinDelphiTrack.h"

#include <cmath>
#include <float.h>

#include "MarlinTrackFit.h"

#include "TanagraTrack.h"

#include "TrackPropagators.h"

#include <lcio.h>
#include <IMPL/TrackImpl.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/LCCollection.h>

#include"tktkbank.h"

#include <ILDDetectorIDs.h>
#include "MaterialDB_F77.hh" 

#include "marlin/Global.h"
#include "gear/GEAR.h"
#include <gear/BField.h>

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "streamlog/streamlog.h"

namespace marlin_delphiF77{

  struct compare_r {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ; 
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ; 
      return ( r_a_sqd < r_b_sqd ) ; 
    }
  } ; 


  // constructor 
  MarlinDelphiTrack::MarlinDelphiTrack( EVENT::Track * lcTrk) 
    : _initialLCTrack(lcTrk), _ipPropagatedLCTrack(NULL), _fit_done(false)
  {

    Tk_Tk_Bank::Instance().clear();

    _lcioHits = lcTrk->getTrackerHits() ;
  
    // sort hits in R
    sort(_lcioHits.begin(), _lcioHits.end(), compare_r() );
  
    // add hits from the track supplied lcio

    EVENT::TrackerHitVec::iterator it = _lcioHits.begin();  
    for( it = _lcioHits.begin() ; it != _lcioHits.end() ; ++it )
      {

	EVENT::TrackerHit* trkhit = (*it);
	//      int layerID = trkhit->ext<ILDDetectorIDs::HitInfo>()->layerID ;
	int layerID = 0;
	if( trkhit->ext<ILDDetectorIDs::HitInfo>() ) {
	  layerID = trkhit->ext<ILDDetectorIDs::HitInfo>()->layerID ;
	}

	streamlog_out(DEBUG1) << "hit " << it - _lcioHits.begin() 
			      << " has type " << trkhit->getType() 
			      << " and layer "
			      << layerID
			      << std::endl ;
      

	if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::VXD ){	

	  streamlog_out(DEBUG1) << "MarlinDelphiTrack::MarlinDelphiTrack add VXD hit with layerID = " << layerID << std::endl ; 

	}
	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::SIT ){	

	  streamlog_out(DEBUG1) << "MarlinDelphiTrack::MarlinDelphiTrack add SIT hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::TPC ){	

	  streamlog_out(DEBUG1) << "MarlinDelphiTrack::MarlinDelphiTrack add TPC hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::FTD ){	

	  streamlog_out(DEBUG1) << "MarlinDelphiTrack::MarlinDelphiTrack add FTD hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::SET ){	

	  streamlog_out(DEBUG1) << "MarlinDelphiTrack::MarlinDelphiTrack add SET hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::ETD ){	

	  streamlog_out(DEBUG1) << "MarlinDelphiTrack::MarlinDelphiTrack add ETD hit with layerID = " << layerID << std::endl ; 

	}
	
      }


    streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack track created with" 
			  << "number of hits = " << _lcioHits.size() 
			  << std::endl ;
  } 


  MarlinDelphiTrack::~MarlinDelphiTrack(){

    for( unsigned int i=0; i < _tanagra_fits.size(); ++i) {
      delete _tanagra_fits[i] ;
    }

    for( unsigned int i=0; i < _lcio_tracks.size(); ++i) {
      delete _lcio_tracks[i] ;
    }

  }


  bool MarlinDelphiTrack::fit( bool fitDirection ) {

    streamlog_out(DEBUG2) << "MarlinDelphiTrack::fit() called " << std::endl ;

    MaterialDB_F77::Instance()->switchONMaterial();

    if (_lcioHits.size() < 3) {
    
      streamlog_out( ERROR) << "<<<<<< MarlinDelphiTrack::fitTrack(): Shortage of Hits! nhits = "  
			    << _lcioHits.size() << " >>>>>>>" << std::endl;
      return false;

    }

    if ( ! _fit_done ) {

      // inputs for MarlinTrackFit
      int useExtraPoint = 0 ;
      int fitOpt = 3 ;
      float chi2PrefitCut = 1000.0 ;
      
      float bField = float(marlin::Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z()) ;
      
      int nhits = _lcioHits.size() ;
      
      float * xhit = new float[nhits];
      float * yhit = new float[nhits];
      float * zhit = new float[nhits];
      float * rphireso = new float[nhits];
      float * zreso = new float[nhits];
      int   * idet = new int[nhits];
      int   * itype = new int[nhits];

      for( int i=0 ; i < nhits ; ++i )
	{ 

	  EVENT::TrackerHit * trkhit = _lcioHits[i];

	  xhit[i] = trkhit->getPosition()[0];
	  yhit[i] = trkhit->getPosition()[1];
	  zhit[i] = trkhit->getPosition()[2];
	
	  rphireso[i] = sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]) ;
	  zreso[i]    = sqrt(trkhit->getCovMatrix()[5]) ; 

	  // set detector type and coordinate type
	  int layerID = trkhit->ext<ILDDetectorIDs::HitInfo>()->layerID;

	  if(      (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::VXD )
	    {
	      idet[i]  = int(3);
	      itype[i] = int(3);
	    }
	  else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::SIT )
	    {
	      idet[i]  = int(7);
	      itype[i] = int(3);
	    }
	  else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::FTD )
	    {
	      idet[i]  = int(2);
	      itype[i] = int(2);
	    }
	  else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::TPC )
	    {
	      idet[i]  = int(1);
	      itype[i] = int(3);
	    }
	  else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::SET )
	    {
	      idet[i]  = int(8);
	      itype[i] = int(3);
	    }
	  else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::ETD )
	    {
	      idet[i]  = int(9);
	      itype[i] = int(2);
	    }

	}

      // outputs
      float param[5]; 
      float eparam[15];
      float PCA[3];
      float chi2;
      int   ndf;
      float chi2rphi;
      float chi2z;
    
      int * lhits = new int[nhits];

      MarlinTrackFit mtrk;

      const int storeExtraFits = 1 ;

      //SJA:FIXME: dofitting is only used at the moment to passed the inputs to the F77 code. The output fits are retreived via the tktkbanks.

      int error = mtrk.DoFitting( useExtraPoint, fitOpt, 
				  nhits, bField, idet, itype, 
				  chi2PrefitCut, 
				  xhit, yhit, zhit, rphireso, zreso, 
				  param, eparam, PCA, chi2, ndf,
				  chi2rphi, chi2z, lhits, storeExtraFits) ;

      if( error != 0 ){
	streamlog_out( DEBUG4 ) << "MarlinDelphiTrack::fit(): DoFitting() returned error code " << error << std::endl ;
	return false ;
      }

      if( Tk_Tk_Bank::Instance().size() == 0 ){
	streamlog_out( DEBUG4 ) << "MarlinDelphiTrack::fit(): Tk_Tk_Bank has size 0 " << std::endl ;
	return false ;
      }

      streamlog_out( DEBUG3 ) << "MarlinDelphiTrack::fit(): track parameters returned by DoFitting(): "
			     << " chi2/ndf " <<  chi2 /  ndf  
			     << " chi2 " <<  chi2 << std::endl 
    
			     << "    D0 "          <<  param[3]  <<  "[+/-" << sqrt( eparam[0] ) << "] " 
			     << "    Phi :"        <<  param[2]  <<  "[+/-" << sqrt( eparam[2] ) << "] " 
			     << "    Omega "       <<  param[0]  <<  "[+/-" << sqrt( eparam[5] ) << "] " 
			     << "    Z0 "          <<  param[4]  <<  "[+/-" << sqrt( eparam[9] ) << "] " 
			     << "    tan(Lambda) " <<  param[1]  <<  "[+/-" << sqrt( eparam[14]) << "] " 
			     << " - error code : " << error 
			     << std::endl 
			     << "\t\t PCA: " << PCA[0] << "\t" << PCA[1] << "\t" <<  PCA[2] 
			     << std::endl ;

      streamlog_out( DEBUG4 ) << "MarlinDelphiTrack::fit(): "
			     << "number of extrapolation surfaces = " << Tk_Tk_Bank::Instance().size()
			     << std::endl ;    

      IMPL::TrackImpl* trkClosestToIP = NULL ;
      double rMinSurf = DBL_MAX ;
      
      for(int iext=0; iext < Tk_Tk_Bank::Instance().size(); ++iext) 
	{

	  float ex_tanagra_param[6]; 
	  float ex_tanagra_eparam[15];
	
	  ex_tanagra_param[0] = Tk_Tk_Bank::Instance().getCoord1_of_ref_point(iext);
	  ex_tanagra_param[1] = Tk_Tk_Bank::Instance().getCoord2_of_ref_point(iext);
	  ex_tanagra_param[2] = Tk_Tk_Bank::Instance().getCoord3_of_ref_point(iext);
	  ex_tanagra_param[3] = Tk_Tk_Bank::Instance().getTheta(iext);
	  ex_tanagra_param[4] = Tk_Tk_Bank::Instance().getPhi(iext);
	  ex_tanagra_param[5] = Tk_Tk_Bank::Instance().getInvp(iext);

	  int fit_code      = Tk_Tk_Bank::Instance().getMeasurement_code(iext);
	  int ndf_surf      = Tk_Tk_Bank::Instance().getNdf(iext);
	  float chi2_surf   = Tk_Tk_Bank::Instance().getChi2(iext);

	  ex_tanagra_eparam[0] = Tk_Tk_Bank::Instance().getCovmatrix1(iext) ;
	  ex_tanagra_eparam[1] = Tk_Tk_Bank::Instance().getCovmatrix2(iext) ;
	  ex_tanagra_eparam[2] = Tk_Tk_Bank::Instance().getCovmatrix3(iext) ;
	  ex_tanagra_eparam[3] = Tk_Tk_Bank::Instance().getCovmatrix4(iext) ;
	  ex_tanagra_eparam[4] = Tk_Tk_Bank::Instance().getCovmatrix5(iext) ;
	  ex_tanagra_eparam[5] = Tk_Tk_Bank::Instance().getCovmatrix6(iext) ;
	  ex_tanagra_eparam[6] = Tk_Tk_Bank::Instance().getCovmatrix7(iext) ;
	  ex_tanagra_eparam[7] = Tk_Tk_Bank::Instance().getCovmatrix8(iext) ;
	  ex_tanagra_eparam[8] = Tk_Tk_Bank::Instance().getCovmatrix9(iext) ;
	  ex_tanagra_eparam[9] = Tk_Tk_Bank::Instance().getCovmatrix10(iext) ;
	  ex_tanagra_eparam[10] = Tk_Tk_Bank::Instance().getCovmatrix11(iext) ;
	  ex_tanagra_eparam[11] = Tk_Tk_Bank::Instance().getCovmatrix12(iext) ;
	  ex_tanagra_eparam[12] = Tk_Tk_Bank::Instance().getCovmatrix13(iext) ;
	  ex_tanagra_eparam[13] = Tk_Tk_Bank::Instance().getCovmatrix14(iext) ;
	  ex_tanagra_eparam[14] = Tk_Tk_Bank::Instance().getCovmatrix15(iext) ;
	
	  float tcov[15];
	  for(int icov = 0; icov < 15; ++icov){
	    tcov[icov] =  ex_tanagra_eparam[icov];
	  }

	  TanagraTrack* tf = new TanagraTrack(  ex_tanagra_param[0], ex_tanagra_param[1], ex_tanagra_param[2], ex_tanagra_param[3], ex_tanagra_param[4], ex_tanagra_param[5], ndf_surf, chi2_surf, fit_code, tcov);
	  
	  _tanagra_fits.push_back(tf);

	  IMPL::TrackImpl* trkAtSurf = TrackPropagators::convertTanagraToLCIO(tf) ;

	  _lcio_tracks.push_back(trkAtSurf) ;
	  
	  //  get cov matrix for printing 
	  EVENT::FloatVec cov = trkAtSurf->getCovMatrix() ;
	  
	  const float* ref = trkAtSurf->getReferencePoint() ;
	  
	  // add the hits. Currently all hits are added and no bookeeping is made of which hits failed to be included in the fit.
	  EVENT::TrackerHitVec::iterator it = _lcioHits.begin();
	  
	  for( it = _lcioHits.begin() ; it != _lcioHits.end() ; ++it )
	    { 
	      trkAtSurf->addHit( *it ) ;
	    }
	  
	  if( ! trkClosestToIP ) {
	    trkClosestToIP = trkAtSurf ;
	    rMinSurf = ex_tanagra_param[0] ;
	  } 
	  else if( ex_tanagra_param[0] < rMinSurf ) {
	    trkClosestToIP = trkAtSurf ;
	    rMinSurf = ex_tanagra_param[0] ;
	  }

	  streamlog_out( DEBUG2 ) << " MarlinDelphi track parameters: for surface " << iext << " at R = " << ex_tanagra_param[0] * 10.0 << " mm "  
				 << " chi2/ndf " <<  trkAtSurf->getChi2() /  trkAtSurf->getNdf()  
				 << " chi2 " <<  trkAtSurf->getChi2() << std::endl 
	    
				 << "\t D0 "          <<  trkAtSurf->getD0()         <<  "[+/-" << sqrt( cov[0] ) << "] " 
				 << "\t Phi :"        <<  trkAtSurf->getPhi()        <<  "[+/-" << sqrt( cov[2] ) << "] " 
				 << "\t Omega "       <<  trkAtSurf->getOmega()      <<  "[+/-" << sqrt( cov[5] ) << "] " 
				 << "\t Z0 "          <<  trkAtSurf->getZ0()         <<  "[+/-" << sqrt( cov[9] ) << "] " 
				 << "\t tan(Lambda) " <<  trkAtSurf->getTanLambda()  <<  "[+/-" << sqrt( cov[14]) << "] " 
	    
				 << "\t ref : [" << ref[0] << ", " << ref[1] << ", "  << ref[2] 
				 << " - r: " << std::sqrt( ref[0]*ref[0]+ref[1]*ref[1] ) << "]" 
	
				 << "    rphi "        <<  ex_tanagra_param[1]  <<  "[+/-" << sqrt( ex_tanagra_eparam[0] ) << "] " 
				 << "    z "           <<  ex_tanagra_param[2]  <<  "[+/-" << sqrt( ex_tanagra_eparam[2] ) << "] " 
				 << "    theta "       <<  ex_tanagra_param[3]  <<  "[+/-" << sqrt( ex_tanagra_eparam[5] ) << "] " 
				 << "    phi "         <<  ex_tanagra_param[4]  <<  "[+/-" << sqrt( ex_tanagra_eparam[9]) << "] " 
				 << "    invP "        <<  ex_tanagra_param[5]  <<  "[+/-" << sqrt( ex_tanagra_eparam[14]) << "] " 
	    
				 << std::endl ;
	  
	}

      
      _ipPropagatedLCTrack = TrackPropagators::PropagateLCIOToNewRef(trkClosestToIP, 0.0, 0.0, 0.0) ;

      EVENT::FloatVec covAtIP = _ipPropagatedLCTrack->getCovMatrix() ;
      
      const float* refClosestToIP = trkClosestToIP->getReferencePoint() ;

      streamlog_out( DEBUG3 ) << " MarlinDelphi track parameters propagated to IP "  
			      << " chi2/ndf " <<  _ipPropagatedLCTrack->getChi2() / _ipPropagatedLCTrack->getNdf()  
			      << " chi2 " <<  _ipPropagatedLCTrack->getChi2() << std::endl 
	
			      << "\t D0 "          <<  _ipPropagatedLCTrack->getD0()         <<  "[+/-" << sqrt( covAtIP[0] ) << "] " 
			      << "\t Phi :"        <<  _ipPropagatedLCTrack->getPhi()        <<  "[+/-" << sqrt( covAtIP[2] ) << "] " 
			      << "\t Omega "       <<  _ipPropagatedLCTrack->getOmega()      <<  "[+/-" << sqrt( covAtIP[5] ) << "] " 
			      << "\t Z0 "          <<  _ipPropagatedLCTrack->getZ0()         <<  "[+/-" << sqrt( covAtIP[9] ) << "] " 
			      << "\t tan(Lambda) " <<  _ipPropagatedLCTrack->getTanLambda()  <<  "[+/-" << sqrt( covAtIP[14]) << "] " 
	
			      << "\t from ref : [" << refClosestToIP[0] << ", " << refClosestToIP[1] << ", "  << refClosestToIP[2] 
			      << " - r: " << std::sqrt( refClosestToIP[0]*refClosestToIP[0]+refClosestToIP[1]*refClosestToIP[1] ) << "]" 
	
			      << std::endl ;

      
      _lcio_tracks.push_back(_ipPropagatedLCTrack) ;
      
      streamlog_out(DEBUG3) << "MarlinDelphiTrack::fit()  _ipPropagatedLCTrack = " << _ipPropagatedLCTrack << std::endl ;

      _fit_done = true ;

    } 

    return true;

  }


  IMPL::TrackImpl* MarlinDelphiTrack::getIPFit(){

    streamlog_out(DEBUG3) << "MarlinDelphiTrack::getIPFit() called " << std::endl ;

    if ( _fit_done == false ) this->fit(false) ;

    streamlog_out(DEBUG3) << "MarlinDelphiTrack::getIPFit()  _ipPropagatedLCTrack = " << _ipPropagatedLCTrack << std::endl ;

    if ( _ipPropagatedLCTrack == NULL ) return NULL ;

    // create NEW track to be returned. The responsiblitiy for deletion lies with the caller.
    IMPL::TrackImpl* trkAtIP = new IMPL::TrackImpl ;

    trkAtIP->setD0( _ipPropagatedLCTrack->getD0() ) ;  
    trkAtIP->setPhi( _ipPropagatedLCTrack->getPhi() ) ;   
    trkAtIP->setOmega( _ipPropagatedLCTrack->getOmega() ) ;
    trkAtIP->setZ0( _ipPropagatedLCTrack->getZ0() ) ;  
    trkAtIP->setTanLambda( _ipPropagatedLCTrack->getTanLambda() ) ;  
    
    trkAtIP->setChi2( _ipPropagatedLCTrack->getChi2() ) ;
    trkAtIP->setNdf( _ipPropagatedLCTrack->getNdf() ) ;

    trkAtIP->setReferencePoint( const_cast<float*>(_ipPropagatedLCTrack->getReferencePoint()) ) ;
    trkAtIP->setIsReferencePointPCA(false) ;

    trkAtIP->setCovMatrix( _ipPropagatedLCTrack->getCovMatrix() ) ;

    //  get cov matrix for printing 
    EVENT::FloatVec cov = trkAtIP->getCovMatrix() ;

    // add the hits. Currently all hits are added and no bookeeping is made of which hits failed to be included in the fit.
    EVENT::TrackerHitVec::iterator it = _lcioHits.begin();
  
    for( it = _lcioHits.begin() ; it != _lcioHits.end() ; ++it )
      { 
	trkAtIP->addHit( *it ) ;
      }
  
    streamlog_out( DEBUG3 ) << " MarlinDelphi track parameters at IP: "
			   << " chi2/ndf " <<  trkAtIP->getChi2() /  trkAtIP->getNdf()  
			   << " chi2 " <<  trkAtIP->getChi2() << std::endl 
    
			   << "\t D0 "          <<  trkAtIP->getD0()         <<  "[+/-" << sqrt( cov[0] ) << "] " 
			   << "\t Phi :"        <<  trkAtIP->getPhi()        <<  "[+/-" << sqrt( cov[2] ) << "] " 
			   << "\t Omega "       <<  trkAtIP->getOmega()      <<  "[+/-" << sqrt( cov[5] ) << "] " 
			   << "\t Z0 "          <<  trkAtIP->getZ0()         <<  "[+/-" << sqrt( cov[9] ) << "] " 
			   << "\t tan(Lambda) " <<  trkAtIP->getTanLambda()  <<  "[+/-" << sqrt( cov[14]) << "] " 
    
			   << std::endl ;
    
    return trkAtIP ;

  }


  IMPL::TrackImpl* MarlinDelphiTrack::getNearestFitToPoint(float* point){

    if ( _fit_done == false ) this->fit(false) ;

    IMPL::TrackImpl* trk = NULL ;
    double mindist2 = DBL_MAX ;
    
    for( unsigned int i=0; i < _lcio_tracks.size(); ++i) {
      
      const float* ref = _lcio_tracks[i]->getReferencePoint() ;
      
      const double deltaX = ref[0] - point[0] ;
      const double deltaY = ref[1] - point[1] ;
      const double deltaZ = ref[2] - point[2] ;
      
      const double dist2 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ ; 

      if( mindist2 > dist2 ) {
	trk = _lcio_tracks[i] ;
	mindist2 = dist2 ;
      }
      
    }
    
    return trk ; 

  } 
  
  IMPL::TrackImpl* MarlinDelphiTrack::getNearestFitToCylinder(float r) {
    
    if ( _fit_done == false ) this->fit(false) ;

    IMPL::TrackImpl* trk = NULL ;
    double mindist2 = DBL_MAX ;
    
    for( unsigned int i=0; i < _lcio_tracks.size(); ++i) {
      
      const float* ref = _lcio_tracks[i]->getReferencePoint() ;

      const double r2_point = ref[0] * ref[0] + ref[1] * ref[1] ;
      
      const double delta_r2 =  r2_point - r*r ;

      const double dist2 = delta_r2 * delta_r2;
      
      if( mindist2 >  dist2 ) {
	trk = _lcio_tracks[i] ;
	mindist2 = dist2 ;
      }
      
    }
    
    return trk ; 

  }


  IMPL::TrackImpl* MarlinDelphiTrack::getNearestFitToZPlane(float z) {
    
    if ( _fit_done == false ) this->fit(false) ;
    
    IMPL::TrackImpl* trk = NULL ;
    double mindist2 = DBL_MAX ;
    
    for( unsigned int i=0; i < _lcio_tracks.size(); ++i) {
      
      const float* ref = _lcio_tracks[i]->getReferencePoint() ;
      
      const double deltaZ = ref[2] - z ;
      
      const double dist2 = deltaZ*deltaZ ; 
      
      if( mindist2 > dist2 ) {
	trk = _lcio_tracks[i] ;
	mindist2 = dist2 ;
      }
      
    }
    
    return trk ; 

  }




} // end of marlin_delphiF77 namespace 
