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

namespace marlin_delphiF77{
  class TanagraTrack ;
}

namespace marlin_delphiF77{

  class MarlinDelphiTrack : public IMarlinTrack {

  public:

    MarlinDelphiTrack( EVENT::Track* trk ) ;

    ~MarlinDelphiTrack() ;

  protected:
  
  private:

    // make member functions private to force use through interface
    bool fit( bool fitDirection );
  
    // returns an LCIO Track with fit parameters determinded at the IP. SJA:FIXME: who will delete this?
    IMPL::TrackImpl* getIPFit();

    // returns an LCIO Track whose referece point is the closest to the specified point. SJA:FIXME: who will delete this?
    IMPL::TrackImpl* getNearestFit(float* point) ;

    //    void ConvertCovMatrix(float* rfit, float* rfite, float* param, float& bField, int fitCode, float* eparam);

    //    void ConvertTANAGRAtoLC(float* rfit, float & bField, int & fitCode, float* param) ;

    //    IMPL::TrackImpl* PropegateLCToNewRef( EVENT::Track* trk, float xref, float yref, float zref ) ;

    // memeber variables 
    EVENT::Track*    _initialLCTrack ;
    IMPL::TrackImpl* _currentLCTrack ;


    EVENT::TrackerHitVec _lcioHits ; 

    std::vector<TanagraTrack *> _tanagra_fits ;

    bool _fit_done ;

  } ;

} // end of marlin_delphiF77 namespace 

#endif
