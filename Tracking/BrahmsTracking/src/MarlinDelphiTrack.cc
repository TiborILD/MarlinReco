
#include "MarlinDelphiTrack.h"

#include <cmath>

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
    : _initialLCTrack(lcTrk), _currentLCTrack(NULL), _fit_done(false)
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

	streamlog_out(DEBUG3) << "hit " << it - _lcioHits.begin() 
			      << " has type " << trkhit->getType() 
			      << " and layer "
			      << trkhit->ext<ILDDetectorIDs::HitInfo>()->layerID
			      << std::endl ;
      

	if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::VXD ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add VXD hit with layerID = " << layerID << std::endl ; 

	}
	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::SIT ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add SIT hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::TPC ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add TPC hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::FTD ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add FTD hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::SET ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add SET hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::ETD ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add ETD hit with layerID = " << layerID << std::endl ; 

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

  }


  bool MarlinDelphiTrack::fit( bool fitDirection ) {

    streamlog_out(DEBUG) << "MarlinDelphiTrack::fit() called " << std::endl ;

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

      int error = mtrk.DoFitting( useExtraPoint, fitOpt, 
				  nhits, bField, idet, itype, 
				  chi2PrefitCut, 
				  xhit, yhit, zhit, rphireso, zreso, 
				  param, eparam, PCA, chi2, ndf,
				  chi2rphi, chi2z, lhits) ;

      streamlog_out( DEBUG ) << "MarlinDelphiTrack::fit(): track parameters returned by DoFitting(): "
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

      streamlog_out( DEBUG ) << "MarlinDelphiTrack::fit(): "
			     << "number of extrapolation surfaces = " << Tk_Tk_Bank::Instance().size()
			     << std::endl ;    

      IMPL::TrackImpl* trkClosestToIP(NULL) ;
      double rMinSurf = 0.0 ;
      
      for(int iext=0; iext < Tk_Tk_Bank::Instance().size(); ++iext) 
	{

	  float ex_tanagra_param[6]; 
	  float ex_tanagra_eparam[15];

	  float ex_lcio_param[5]; 
	  float ex_lcio_eparam[15];
	
	  ex_tanagra_param[0] = Tk_Tk_Bank::Instance().getCoord1_of_ref_point(iext);
	  ex_tanagra_param[1] = Tk_Tk_Bank::Instance().getCoord2_of_ref_point(iext);
	  ex_tanagra_param[2] = Tk_Tk_Bank::Instance().getCoord3_of_ref_point(iext);
	  ex_tanagra_param[3] = Tk_Tk_Bank::Instance().getTheta(iext);
	  ex_tanagra_param[4] = Tk_Tk_Bank::Instance().getPhi(iext);
	  ex_tanagra_param[5] = Tk_Tk_Bank::Instance().getInvp(iext);

	  int fit_code = Tk_Tk_Bank::Instance().getMeasurement_code(iext);
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

	  //	  IMPL::TrackImpl* trkAtSurf = tf->convertTANAGRAtoLC() ;
	  IMPL::TrackImpl* trkAtSurf = TrackPropagators::convertTanagraToLCIO(tf) ;
	  
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

	  streamlog_out( DEBUG ) << " MarlinDelphi track parameters: for surface " << iext << " at R = " << ex_tanagra_param[0] * 10.0 << " mm "  
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
      

      //      IMPL::TrackImpl* trkAtIP = this->PropegateLCToNewRef(trkClosestToIP, 0.0, 0.0, 0.0) ;
      IMPL::TrackImpl* trkAtIP = TrackPropagators::PropagateLCIOToNewRef(trkClosestToIP, 0.0, 0.0, 0.0) ;

      EVENT::FloatVec covAtIP = trkAtIP->getCovMatrix() ;
      
      const float* refClosestToIP = trkClosestToIP->getReferencePoint() ;

      streamlog_out( DEBUG ) << " MarlinDelphi track parameters propagated to IP "  
			     << " chi2/ndf " <<  trkAtIP->getChi2() / trkAtIP->getNdf()  
			     << " chi2 " <<  trkAtIP->getChi2() << std::endl 
	
			     << "\t D0 "          <<  trkAtIP->getD0()         <<  "[+/-" << sqrt( covAtIP[0] ) << "] " 
			     << "\t Phi :"        <<  trkAtIP->getPhi()        <<  "[+/-" << sqrt( covAtIP[2] ) << "] " 
			     << "\t Omega "       <<  trkAtIP->getOmega()      <<  "[+/-" << sqrt( covAtIP[5] ) << "] " 
			     << "\t Z0 "          <<  trkAtIP->getZ0()         <<  "[+/-" << sqrt( covAtIP[9] ) << "] " 
			     << "\t tan(Lambda) " <<  trkAtIP->getTanLambda()  <<  "[+/-" << sqrt( covAtIP[14]) << "] " 
	
			     << "\t ref : [" << refClosestToIP[0] << ", " << refClosestToIP[1] << ", "  << refClosestToIP[2] 
			     << " - r: " << std::sqrt( refClosestToIP[0]*refClosestToIP[0]+refClosestToIP[1]*refClosestToIP[1] ) << "]" 
	
			     << std::endl ;


      // create track to be returned
      _currentLCTrack = new IMPL::TrackImpl();
      
      //  this is for incomming tracks ...
      double omega     =  param[0] ;            
      double tanLambda =  param[1] ;
      double phi0      =  param[2] ;
      double d0        =  param[3] ;
      double z0        =  param[4] ;
      
      _currentLCTrack->setD0( d0 ) ;  
      _currentLCTrack->setPhi( phi0  ) ;   
      _currentLCTrack->setOmega( omega  ) ;
      _currentLCTrack->setZ0( z0  ) ;  
      _currentLCTrack->setTanLambda( tanLambda ) ;  
      
      _currentLCTrack->setChi2( chi2 ) ;
      
      float pivot[3] ;
      pivot[0] =  PCA[0] + d0 * sin(phi0) ;
      pivot[1] =  PCA[1] - d0 * cos(phi0) ;
      pivot[2] =  PCA[2]  ;
      
      _currentLCTrack->setReferencePoint( pivot ) ;
      _currentLCTrack->setIsReferencePointPCA(false);

      EVENT::FloatVec cov( 15 )  ; 
      cov[ 0] =  eparam[ 0] ; //   d0,   d0
		        
      cov[ 1] =  eparam[ 1] ; //   phi0, d0
      cov[ 2] =  eparam[ 2] ; //   phi0, phi
		        
      cov[ 3] =  eparam[ 3] ; //   omega, d0
      cov[ 4] =  eparam[ 4] ; //   omega, phi
      cov[ 5] =  eparam[ 5] ; //   omega, omega
		        
      cov[ 6] =  eparam[ 6] ; //   z0  , d0
      cov[ 7] =  eparam[ 7] ; //   z0  , phi
      cov[ 8] =  eparam[ 8] ; //   z0  , omega
      cov[ 9] =  eparam[ 9] ; //   z0  , z0
		        
      cov[10] =  eparam[10] ; //   tanl, d0
      cov[11] =  eparam[11] ; //   tanl, phi
      cov[12] =  eparam[12] ; //   tanl, omega
      cov[13] =  eparam[13] ; //   tanl, z0
      cov[14] =  eparam[14] ; //   tanl, tanl

      _currentLCTrack->setCovMatrix( cov ) ;

      streamlog_out(DEBUG) << "MarlinDelphiTrack::fit()  _currentLCTrack = " << _currentLCTrack << std::endl ;

      _fit_done = true ;

    } 

    return true;

  }



  IMPL::TrackImpl* MarlinDelphiTrack::getIPFit(){

    streamlog_out(DEBUG) << "MarlinDelphiTrack::getIPFit() called " << std::endl ;

    streamlog_out(DEBUG) << "MarlinDelphiTrack::getIPFit()  _currentLCTrack = " << _currentLCTrack << std::endl ;

    IMPL::TrackImpl* trk = _currentLCTrack ;

//  get cov matrix for printing 
    EVENT::FloatVec cov = trk->getCovMatrix() ;

    float ref[3] ;

    if(  trk->isReferencePointPCA() == true ) {
      ref[0] = ref[1] = ref[2] = 0.0 ;
    }
    else {
      const float* pivot = trk->getReferencePoint() ;
      ref[0] = pivot[0] ;
      ref[1] = pivot[1] ;
      ref[2] = pivot[2] ;
    }

    // add the hits. Currently all hits are added and no bookeeping is made of which hits failed to be included in the fit.
    EVENT::TrackerHitVec::iterator it = _lcioHits.begin();
  
    for( it = _lcioHits.begin() ; it != _lcioHits.end() ; ++it )
      { 
	trk->addHit( *it ) ;
      }
  
    streamlog_out( DEBUG ) << " MarlinDelphi track parameters: "
			   << " chi2/ndf " <<  trk->getChi2() /  trk->getNdf()  
			   << " chi2 " <<  trk->getChi2() << std::endl 
    
			   << "\t D0 "          <<  trk->getD0()         <<  "[+/-" << sqrt( cov[0] ) << "] " 
			   << "\t Phi :"        <<  trk->getPhi()        <<  "[+/-" << sqrt( cov[2] ) << "] " 
			   << "\t Omega "       <<  trk->getOmega()      <<  "[+/-" << sqrt( cov[5] ) << "] " 
			   << "\t Z0 "          <<  trk->getZ0()         <<  "[+/-" << sqrt( cov[9] ) << "] " 
			   << "\t tan(Lambda) " <<  trk->getTanLambda()  <<  "[+/-" << sqrt( cov[14]) << "] " 
    
			   << "\t ref : [" << ref[0] << ", " << ref[1] << ", "  << ref[2] 
			   << " - r: " << std::sqrt( ref[0]*ref[0]+ref[1]*ref[1] ) << "]" 
			   << std::endl ;
  

  
    return trk;

  }



IMPL::TrackImpl* MarlinDelphiTrack::getNearestFit(float* point){
  //


  return NULL ; // not yet implemented
} 



} // end of marlin_delphiF77 namespace 
