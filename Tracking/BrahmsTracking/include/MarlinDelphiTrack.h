#ifndef MarlinDelphiTrack_h
#define MarlinDelphiTrack_h

#include "IMarlinTrack.h"

#include <vector>

namespace EVENT{
  class TrackerHit ;
  class Track ;
}

namespace IMPL{
  class TrackImpl ;
}

class TanagraFit {
 public:
  
  TanagraFit( float r, float rphi, float z, float theta, float phi, float invp, float* cov)
    : _r(r), _rphi(rphi), _z(z), _theta(theta), _phi(phi), _invp(invp)
  {    for(int i=0; i<15; ++i) { _cov[i] = cov[i] ; } } ;

  float  get_r() const { return _r ; } ;	
  float  get_rphi() const { return _rphi ; } ;	
  float  get_z() const { return _z ; } ;	
  float  get_theta() const { return _theta ; } ;
  float  get_phi() const { return _phi ; } ;	
  float  get_invp() const { return _invp ; } ; 
  void   get_cov( float* cov ) const { for(int i=0; i<15; ++i) { cov[i] = _cov[i] ; } } ; 

 private:
  float _r ;
  float _rphi ;
  float _z ;
  float _theta ;
  float _phi ;
  float _invp ;  
  float _cov[15] ;

} ;


namespace marlin_delphiF77{
  class MarlinDelphiTrack : public IMarlinTrack {

  public:

    MarlinDelphiTrack( EVENT::Track* trk ) ;

    ~MarlinDelphiTrack() ;

    void getTanagraFits( std::vector<TanagraFit const*>& fits) ;

  protected:
  
  private:

    // make member functions private to force use through interface
    bool fit( bool fitDirection );
  
    // returns an LCIO Track with fit parameters determinded at the IP. SJA:FIXME: who will delete this?
    IMPL::TrackImpl* getIPFit();
  
    void ConvertCovMatrix(float * rfit, float * rfite, float * param, float & bField, int fitCode, float * eparam) ;
    void ConvertTANAGRAtoLC(float * rfit, float & bField, int & fitCode, float * param) ;

    // memeber variables 
    EVENT::Track*    _initialLCTrack ;
    IMPL::TrackImpl* _currentLCTrack ;


    EVENT::TrackerHitVec _lcioHits ; 

    std::vector<TanagraFit *> _tanagra_fits ;

    bool _fit_done ;

  } ;

} // end of marlin_kaltest namespace 

#endif
