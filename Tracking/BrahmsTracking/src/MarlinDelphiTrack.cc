
#include "MarlinDelphiTrack.h"

//SJA:try using MarlinTrackFit for now
#include "MarlinTrackFit.h"

#include <cmath>

#include <lcio.h>
#include <IMPL/TrackImpl.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/LCCollection.h>

#include"tktkbank.h"

#include <ILDDetectorIDs.h>
#include "MaterialDB_F77.hh" 

//---- GEAR ---- // this is only need for the bfield to pass to MarlinTrackFit it should not be needed and should be removed
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include <gear/BField.h>

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

    // use global TkTkBank pointer for now. The TkTkBank should be made into a singleton
    //    if( TkTkBank != NULL) throw ;
    //    TkTkBank = new Tk_Tk_Bank ;  
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

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add VXD hit " << std::endl ; 

	}
	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::SIT ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add SIT hit with layerID = " << layerID << std::endl ; 

	}

	else if( (layerID / ILDDetectorIDs::DetID::Factor) == ILDDetectorIDs::DetID::TPC ){	

	  streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack add TPC hit " << std::endl ; 

	}
	
      }


    streamlog_out(DEBUG3) << "MarlinDelphiTrack::MarlinDelphiTrack track created with" 
			  << "number of hits = " << _lcioHits.size() 
			  << std::endl ;


  } 

  MarlinDelphiTrack::~MarlinDelphiTrack(){
    //    delete TkTkBank;
    //    TkTkBank = NULL;
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
      // inputs
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
      float RefPoint[3];
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
				  param, eparam, RefPoint, chi2, ndf,
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
			     << "\t\t RefPoint: " << RefPoint[0] << "\t" << RefPoint[1] << "\t" <<  RefPoint[2] 
			     << std::endl ;

      streamlog_out( DEBUG ) << "MarlinDelphiTrack::fit(): "
			     << "number of extrapolation surfaces = " << Tk_Tk_Bank::Instance().size()
			     << std::endl ;    

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
	
	  ConvertTANAGRAtoLC( ex_tanagra_param, bField, fit_code, ex_lcio_param );
	  ConvertCovMatrix( ex_tanagra_param, ex_tanagra_eparam, ex_lcio_param, bField, fit_code, ex_lcio_eparam) ;	

	  float Rref;
	  if( fit_code == 1){
	    Rref = ex_tanagra_param[0]; 
	  } else {
	    Rref = sqrt( ex_tanagra_param[0] * ex_tanagra_param[0] + ex_tanagra_param[1] * ex_tanagra_param[1]) ;
	  }

	  // convert to mm 
	  Rref *= 10.0 ;

	  ex_lcio_param[0] = 0.1*ex_lcio_param[0];
	  ex_lcio_param[3] = 10.*ex_lcio_param[3];
	  ex_lcio_param[4] = 10.*ex_lcio_param[4];

	  // convert from cm to mm
	  float r     = ex_tanagra_param[0] * 10.0 ;
	  float rphi  = ex_tanagra_param[1] * 10.0 ;
	  float z     = ex_tanagra_param[2] * 10.0 ;
	  float theta = ex_tanagra_param[3] ;
	  float phi   = ex_tanagra_param[4] ;
	  float invp  = ex_tanagra_param[5] ;

	  float tcov[15];
	  for(int icov = 0; icov < 15; ++icov){
	    tcov[icov] =  ex_tanagra_eparam[icov];
	  }

	  // convert from cm to mm
	  tcov[0] =  tcov[0] * 10.0 * 10.0 ;  // (R-Phi,R-Phi)
	  tcov[1] =  tcov[1] * 10.0 * 10.0 ;  // (R-Phi,ZR) 
	  tcov[2] =  tcov[2] * 10.0 * 10.0 ;  // (ZR,ZR)
	  tcov[3] =  tcov[3] * 10.0 ;         // (R-Phi,Theta)
	  tcov[4] =  tcov[4] * 10.0 ;         // (ZR,Theta) 
//	  tcov[5] =  tcov[5] ;                // (Theta,Theta)
	  tcov[6] =  tcov[6] * 10.0 ;         // (R-Phi,PhiR)
	  tcov[7] =  tcov[7] * 10.0 ;         // (ZR,PhiR) 
//	  tcov[8] =  tcov[8] ;                // (Theta,PhiR)
//	  tcov[9] =  tcov[9] ;                // (PhiR,PhiR)
	  tcov[10] =  tcov[10] * 10.0 ;       // (R-Phi,1/p)
	  tcov[11] =  tcov[11] * 10.0 ;       // (ZR,1/p)
//	  tcov[12] =  tcov[12] ;              // (Theta,1/p) 
//	  tcov[13] =  tcov[13] ;              // (PhiR,1/p) 
//	  tcov[14] =  tcov[14] ;              // (1/p,1/p)

	  TanagraFit* tf = new TanagraFit(  r,  rphi,  z,  theta,  phi,  invp,  tcov);

	  _tanagra_fits.push_back(tf);

	  streamlog_out( DEBUG ) << "Track Parameters for surface " << iext << " at R = " << Rref 
				 << "    D0 "          <<  ex_lcio_param[3]  <<  "[+/-" << sqrt( ex_lcio_eparam[0] ) << "] " 
				 << "    Phi "         <<  ex_lcio_param[2]  <<  "[+/-" << sqrt( ex_lcio_eparam[2] ) << "] " 
				 << "    Omega "       <<  ex_lcio_param[0]  <<  "[+/-" << sqrt( ex_lcio_eparam[5] ) << "] " 
				 << "    Z0 "          <<  ex_lcio_param[4]  <<  "[+/-" << sqrt( ex_lcio_eparam[9] ) << "] " 
				 << "    tan(Lambda) " <<  ex_lcio_param[1]  <<  "[+/-" << sqrt( ex_lcio_eparam[14]) << "] " 
				 << std::endl
				 << "    rphi "        <<  ex_tanagra_param[1]  <<  "[+/-" << sqrt( ex_tanagra_eparam[0] ) << "] " 
				 << "    z "           <<  ex_tanagra_param[2]  <<  "[+/-" << sqrt( ex_tanagra_eparam[2] ) << "] " 
				 << "    theta "       <<  ex_tanagra_param[3]  <<  "[+/-" << sqrt( ex_tanagra_eparam[5] ) << "] " 
				 << "    phi "         <<  ex_tanagra_param[4]  <<  "[+/-" << sqrt( ex_tanagra_eparam[9]) << "] " 
				 << "    invP "        <<  ex_tanagra_param[5]  <<  "[+/-" << sqrt( ex_tanagra_eparam[14]) << "] " 
				 <<  std::endl ;
	}

      // create track to be returned
      _currentLCTrack = new IMPL::TrackImpl();

      //  this is for incomming tracks ...
      double omega     =  param[0] ;            
      double tanLambda =  param[1] ;
      double phi       =  param[2] ;
      double d0        =  param[3] ;
      double z0        =  param[4] ;

      _currentLCTrack->setD0( d0 ) ;  
      _currentLCTrack->setPhi( phi  ) ; // fi0  - M_PI/2.  ) ;  
      _currentLCTrack->setOmega( omega  ) ;
      _currentLCTrack->setZ0( z0  ) ;  
      _currentLCTrack->setTanLambda( tanLambda ) ;  
  
      _currentLCTrack->setChi2( chi2 ) ;
  
      float pivot[3] ;
      pivot[0] =  RefPoint[0]  ;
      pivot[1] =  RefPoint[1]  ;
      pivot[2] =  RefPoint[2]  ;
  
      _currentLCTrack->setReferencePoint( pivot ) ;
      _currentLCTrack->setIsReferencePointPCA(true);

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
    const float * pivot = trk->getReferencePoint() ;


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
    
			   << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", "  << pivot[2] 
			   << " - r: " << std::sqrt( pivot[0]*pivot[0]+pivot[1]*pivot[1] ) << "]" 
			   << std::endl ;
  

  
    return trk;

  }

  void MarlinDelphiTrack::getTanagraFits( std::vector<TanagraFit const*>& fits ){

    if( fits.size()!=0 ) { // vector must be empty
      streamlog_out( DEBUG ) << "MarlinDelphiTrack::getTanagraFits vector passed not empty" << std::endl ;
      throw ; 
    }
    else {
      for( unsigned int i=0; i < _tanagra_fits.size(); ++i) {
	fits.push_back(_tanagra_fits[i]) ;
      }
    }
  }


void MarlinDelphiTrack::ConvertTANAGRAtoLC(float * rfit, float & bField, int & fitCode, float * param) {

//    // Calculation of helix parameters 
//    // in canonical form
//    // measurement code (R,R-PHI,Z)
//    float xx = rfit[0]*cos(rfit[1]/rfit[0]);
//    float yy = rfit[0]*sin(rfit[1]/rfit[0]);
//    float pos[3];
//    // conversion to mm
//    // measurement code (R,R-PHI,Z)
//    pos[0] = 10*xx;
//    pos[1] = 10*yy;
//    pos[2] = 10*rfit[2];
//    if (fitCode == 0) { // measurement code (X,Y,Z)
//      pos[0] = 10*rfit[0];
//      pos[1] = 10*rfit[1];
//    }
//    float mom[3];
//    mom[0] = sin(rfit[3])*cos(rfit[4])/fabs(rfit[5]);
//    mom[1] = sin(rfit[3])*sin(rfit[4])/fabs(rfit[5]);
//    mom[2] = cos(rfit[3])/fabs(rfit[5]);
//    float charge = -rfit[5]/fabs(rfit[5]);
//    
//    HelixClass helix;
//    helix.Initialize_VP(pos,mom,charge,bField);
//    
//    // conversion to cm
//    param[0] = 10.*helix.getOmega();
//    param[1] = helix.getTanLambda();
//    param[2] = helix.getPhi0();
//    param[3] = 0.1*helix.getD0();
//    param[4] = 0.1*helix.getZ0();    
//
//
//
    // Calculation of helix parameters 
    // in canonical form
    // measurement code (R,R-PHI,Z)

  // note everything still done in cm at this point

  const double FCT = 2.997924e-3 ;

  double r_ref     = rfit[0] ;
  double rphi_ref  = rfit[1] ;
  double phi_ref   = rphi_ref / r_ref ;
  double z_ref     = rfit[2] ;
  double theta     = rfit[3];
  double phi       = rfit[4];
  double inv_p     = fabs(rfit[5]);
  double charge    = -rfit[5]/fabs(rfit[5]);

  double omega = charge * (FCT*double(bField)*inv_p) / sin(theta) ;
  
  double tanLambda = 1.0 / tan(theta) ;

  double x_ref = r_ref*cos(phi_ref) ;
  double y_ref = r_ref*sin(phi_ref) ;
  
  if (fitCode == 0) { // measurement code (X,Y,Z)
    x_ref =  rfit[0] ;
    y_ref =  rfit[1] ;
  }
  
  double radius = fabs(1.0/omega) ;

  double x_centre = x_ref + radius*cos( phi - 0.5*M_PI * charge) ;
  double y_centre = y_ref + radius*sin( phi - 0.5*M_PI * charge) ;
  double r_centre = sqrt(x_centre*x_centre+y_centre*y_centre) ;
  
  double phi_of_pca = atan2(y_centre,x_centre) ;

  double phi0 = phi_of_pca + (charge * 0.5*M_PI) ;

  while ( phi0 < 0 )  phi0 += 2.0*M_PI ;
  while ( phi0 >= 2.0*M_PI ) phi0 -= 2.0*M_PI ;

  double x_pca = x_centre - radius*cos( phi_of_pca );
  double y_pca = y_centre - radius*sin( phi_of_pca );

  double delta_x = x_pca - x_ref ;
  double delta_y = y_pca - y_ref ;
  
  double alpha = -omega * delta_x * cos(phi) - omega * delta_y * sin(phi)  ;
  double beta  = 1.0 - omega * delta_x * sin(phi) + omega * delta_y * cos(phi)  ;
  
  double s = atan2(-alpha,beta) / omega ;
  
  double z0 = z_ref + s * tanLambda ;
    
  double d0;

  if (charge > 0) {
    d0 = charge*radius - r_centre ;
  }
  else {
    d0 = charge*radius + r_centre ;
  }

  param[0] = omega ;
  param[1] = tanLambda ;
  param[2] = phi0;
  param[3] = d0 ;
  param[4] = z0 ;    

//  std::cout << " atan2(-alpha,beta) = " << atan2(-alpha,beta) << std::endl;
//  std::cout << " x_pca = " << x_pca << std::endl;
//  std::cout << " y_pca = " << y_pca << std::endl;
//  std::cout << " delta_x = " << delta_x << std::endl;
//  std::cout << " delta_y = " << delta_y << std::endl;
//  std::cout << " alpha = "   << alpha << std::endl;
//  std::cout << " beta  = "   << beta << std::endl;
//  std::cout << " x ref = "  << x_ref << std::endl;
//  std::cout << " y ref = "  << y_ref << std::endl;
//  std::cout << " z ref = "  << z_ref << std::endl;
//  std::cout << " tanLambda = " << tanLambda << std::endl;
//  std::cout << " s = " <<  s << std::endl;

    
}


void MarlinDelphiTrack::ConvertCovMatrix(float * rfit, float * rfite, float * param, float & bField, int fitCode, float * eparam) {
 // Calculation of cov matrix in canonical form
  // Subroutine converts covariance matrix from the 
  // TANAGRA format 
  // (R-Phi,R-Phi)
  // (R-Phi,ZR)     (ZR,ZR)
  // (R-Phi,Theta)  (ZR,Theta)  (Theta,Theta)
  // (R-Phi,phi)    (ZR,phi)    (Theta,phi)    (phi,phi)
  // (R-Phi,1/p)    (ZR,1/p)    (Theta,1/p)    (phi,1/p)  (1/p,1/p)
  // into canonical format
  // (Omega,Omega)
  // (Omega,TanLambda) (TanLambda,TanLambda)
  // (Omega,Phi0)      (TanLambda,Phi0)       (Phi0,Phi0)
  // (Omega,D0)        (TanLambda,D0)         (Phi0,D0)    (D0,D0)
  // (Omega,Z0)        (TanLambda,Z0)         (Phi0,Z0)    (D0,Z0)  (Z0,Z0)
  // TANAGRA   Vector Eta = (R-Phi,ZR,Theta,phi,1/p)
  // Canonical Vector X = (Omega,TanLambda,Phi,D0,Z0)
  // UX[i,j] = VEta[l,m]*dXdEta[i,l]*dXdEta[j,m]
  //
  // Output of the DELPHI fit routine
  // rfit[0] = fixed radius of the reference point (rho, eta0)
  // rfit[1] = R*Phi at the reference point 
  // rfit[2] = Theta angle of the momentum vector @ the reference point
  // rfit[3] = phi angle @ the reference point
  // declaration of cov matricies and Jacobian
  double VEta[5][5];
  double UX[5][5];
  double dXdEta[5][5];

   // Fill array of initial covariance matrix
  int counter = 0;
  for (int i=0;i<5;++i) {
    for (int j=0;j<=i;++j) {
      VEta[i][j] = rfite[counter];
      VEta[j][i] = VEta[i][j];
      counter++;
    }
  }

//   std::cout << "Counter = " << counter << std::endl;

//   for (int i=0;i<5;++i) {
//     printf("%19.15f %19.15f %19.15f %19.15f %19.15f\n",
//  	   VEta[i][0],VEta[i][1],VEta[i][2],VEta[i][3],VEta[i][4]);
//   }

  // Calculation of Jacobian 

  // some variables needed  for 
  // calculations of Jacobian matrix elements
  // Variables are set in cm !
  double cot_eta3 = 1.0/tan(rfit[3]);
  double sin_eta3 = sin(rfit[3]);
  double sin_eta4 = sin(rfit[4]);
  double cos_eta4 = cos(rfit[4]);
  double R = -100.0*sin_eta3/(0.299792458*bField*rfit[5]); 
  double rho = rfit[0];
  // Coordinates of the Circle central point (XC,YC)
  // XC = xRef + R*sin(phi)
  // YC = yRef - R*cos(phi)
  double XC = double(rho*cos(rfit[1]/rho) + R*sin_eta4);
  double YC = double(rho*sin(rfit[1]/rho) - R*cos_eta4);
  if (fitCode == 0) {
      XC = double(rfit[0] + R*sin_eta4);
      YC = double(rfit[1] - R*cos_eta4);
  }
  double RC2 = double(XC*XC+YC*YC);

  // X0 = Omega
  // Omega = -FCT*BField/(p*sin(Theta))
  dXdEta[0][0] = 0;
  dXdEta[0][1] = 0;
  dXdEta[0][2] = -cot_eta3/R;
  dXdEta[0][3] = 0;
  dXdEta[0][4] = 1./(R*rfit[5]);
    
  // X1 = TanLambda
  // TanLambda = cot(Theta)
  dXdEta[1][0] = 0;
  dXdEta[1][1] = 0;
  dXdEta[1][2] = -1/(sin_eta3*sin_eta3);
  dXdEta[1][3] = 0;
  dXdEta[1][4] = 0;
  
  // X2 = Phi0
  // Phi = atan2(-YC,-XC) - pi*q/2,
  // where q is the charge of particles
  double dXCdEta[5];
  double dYCdEta[5];
  
  dXCdEta[0] = -sin(rfit[1]/rho);
  dXCdEta[1] = double(0);
  dXCdEta[2] = R*cot_eta3*sin_eta4;
  dXCdEta[3] = R*cos_eta4;
  dXCdEta[4] = -R*sin_eta4/rfit[5];
  
  dYCdEta[0] = cos(rfit[1]/rho);
  dYCdEta[1] = double(0);
  dYCdEta[2] = -R*cot_eta3*cos_eta4;
  dYCdEta[3] = R*sin_eta4;
  dYCdEta[4] = R*cos_eta4/rfit[5];

  if (fitCode == 0) {
      dXCdEta[0] = double(1);
      dXCdEta[1] = double(0);
      dYCdEta[0] = double(0);
      dYCdEta[1] = double(1);      
  }

  for (int i=0;i<5;++i) 
    dXdEta[2][i] = (dYCdEta[i]*XC-dXCdEta[i]*YC)/RC2;
  
  
  // X3 = D0
  // D0 = R - XC/sin(Phi0) or D0 = R + YC/cos(Phi0)
  double dRdEta[5];
  dRdEta[0] = 0.0;
  dRdEta[1] = 0.0;
  dRdEta[2] = R*cot_eta3;
  dRdEta[3] = 0.0;
  dRdEta[4] = -R/rfit[5];
  
  // method using formulae XC = (R-d0)*sin(Phi0)
  //                       YC = (d0-R)*cos(Phi0) 


  double sinPhi0=sin(param[2]);
  double cosPhi0=cos(param[2]);
  double sin2Phi0=sinPhi0*sinPhi0;
  double cos2Phi0=cosPhi0*cosPhi0;
  
  if (fabs(cosPhi0)>fabs(sinPhi0)) {
    for (int i=0;i<5;++i) {
      dXdEta[3][i] = dRdEta[i] + (dYCdEta[i]*cosPhi0+YC*sinPhi0*dXdEta[2][i])/cos2Phi0;
    }
  }
  else {
    for (int i=0;i<5;++i) {
      dXdEta[3][i] = dRdEta[i] - (dXCdEta[i]*sinPhi0-XC*cosPhi0*dXdEta[2][i])/sin2Phi0;
    }    
  }
  

  // method using formulae : d0 = R-sqrt(XC**2+YC**2) , R > 0
  //                         d0 = R+sqrt(XC**2+YC**2) , R < 0
  // 
  //  double rXCRC = double(XC/RC);
  //  double rYCRC = double(YC/RC);
  //  double RC  = double(sqrt(RC2));

  //   for (int i=0;i<5;++i) {
  //     if (R > 0)
  //       dXdEta[3][i] =  dRdEta[i] - rXCRC*dXCdEta[i] - rYCRC*dYCdEta[i];
  //     else 
  //       dXdEta[3][i] =  dRdEta[i] + rXCRC*dXCdEta[i] + rYCRC*dYCdEta[i];
  //   }

  // X4 = Z0
  // Z0 = ZR + (Phi-PhiR)*sign(Omega)*TanLambda
  dXdEta[4][0] = -param[1]*R*dXdEta[2][0];
  dXdEta[4][1] = 1;
  dXdEta[4][2] = -dXdEta[1][2]*R*param[2] -param[1]*R*dXdEta[2][2] -param[1]*param[2]*dRdEta[2]
    +dXdEta[1][2]*R*rfit[4] +param[1]*dRdEta[2]*rfit[4];
  dXdEta[4][3] = -param[1]*R*dXdEta[2][3]+param[1]*R;
  dXdEta[4][4] = -param[1]*dRdEta[4]*param[2]-param[1]*R*dXdEta[2][4]+param[1]*dRdEta[4]*rfit[4];
  

  // Transformation of Cov matrix -->
  for (int i=0;i<5;++i) {
    for (int j=0;j<5;++j) {
      UX[i][j] = 0.0;
      for (int k=0;k<5;++k) {
	for (int l=0;l<5;++l)
	  UX[i][j] += VEta[k][l]*dXdEta[i][k]*dXdEta[j][l];
      }
    }
  }

//   for (int iC=0;iC<15;++iC) {
//     std::cout << eparam[iC] << " " ;
//   }
//   std::cout << std::endl;
  
//   for (int i=0;i<5;++i) {
//     printf("%19.15f %19.15f %19.15f %19.15f %19.15f\n",
//  	   UX[i][0],UX[i][1],UX[i][2],UX[i][3],UX[i][4]);
//   }

//   std::cout << std::endl;


  //.. (omega,tanl,phi,d0,z0) ->
  //.. (d0,phi0,omega,z0,tanl)

  // reordering errors according to the ILC Convention
  // T.Kraemer,  LC-DET-2006-004
  // Convertion from cm to mm           COV MATRIX
  //                                             ELEMENT          DIMENSION
//   eparam[0] = float(0.01*UX[0][0]);        // (Omega,Omega)           L-2
  
//   eparam[1] = float(0.1*UX[1][0]);         // (Omega,TanLambda)       L-1
//   eparam[2] = float(UX[1][1]);             // (TanLambda,TanLambda)   L0
    
//   eparam[3] = float(0.1*UX[2][0]);         // (Omega,Phi)             L-1
//   eparam[4] = float(UX[2][1]);             // (TanLambda,Phi)         L0
//   eparam[5] = float(UX[2][2]);             // (Phi,Phi)               L0
    
//   eparam[6] = float(UX[3][0]);             // (Omega,D0)              L0
//   eparam[7] = float(10*UX[3][1]);          // (TanLambda,D0)          L1
//   eparam[8] = float(10*UX[3][2]);          // (Phi,D0)                L1
//   eparam[9] = float(100*UX[3][3]);         // (D0,D0)                 L2
  
//   eparam[10] = float(UX[4][0]);            // (Omega,Z0)              L0
//   eparam[11] = float(10*UX[4][1]);         // (TanLambda,Z0)          L1
//   eparam[12] = float(10*UX[4][2]);         // (Phi,Z0)                L1
//   eparam[13] = float(100*UX[4][3]);        // (D0,Z0)                 L2
//   eparam[14] = float(100*UX[4][4]);        // (Z0,Z0)                 L2

  eparam[0] = float(100*UX[3][3]);         // (D0,D0)                 L2
  
  eparam[1] = float(10*UX[2][3]);          // (Phi,D0)                L1
  eparam[2] = float(UX[2][2]);             // (Phi,Phi)               L0
    
  eparam[3] = float(UX[0][3]);             // (Omega,D0)              L0
  eparam[4] = float(0.1*UX[0][2]);         // (Omega,Phi)             L-1
  eparam[5] = float(0.01*UX[0][0]);        // (Omega,Omega)           L-2
    
  eparam[6] = float(100*UX[4][3]);         // (Z0,D0)                 L2
  eparam[7] = float(10*UX[4][2]);          // (Z0,Phi)                L1
  eparam[8] = float(UX[4][0]);             // (Z0,Omega)              L0
  eparam[9] = float(100*UX[4][4]);         // (Z0,Z0)                 L2
  
  eparam[10] = float(10*UX[1][3]);         // (TanLambda,D0)          L1
  eparam[11] = float(UX[1][2]);            // (TanLambda,Phi)         L0
  eparam[12] = float(0.1*UX[1][0]);        // (TanLambda,Omega)       L-1
  eparam[13] = float(10*UX[1][4]);         // (TanLambda,Z0)          L1
  eparam[14] = float(UX[1][1]);            // (TanLambda,TanLambda)   L0

  //  std::cout << "sigma(D0)=" << sqrt(eparam[9]) << std::endl;
  //  std::cout << std::endl;

}



} // end of marlin_delphiF77 namespace 
