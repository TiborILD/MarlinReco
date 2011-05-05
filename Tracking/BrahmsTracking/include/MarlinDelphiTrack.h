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
  
  // returns a pointer to a New LCIO Track with fit parameters determinded at the IP. The responsiblitiy for deletion lies with the caller.
    IMPL::TrackImpl* getIPFit() ;

    // returns a pointer to an LCIO Track whose referece point is the closest to the specified point.
    IMPL::TrackImpl* getNearestFitToPoint(float* point) ;

    // returns a pointer to an LCIO Track whose referece point is the closest to a cylinder of radius r which is centered at the origin parallel to the z axis.
    IMPL::TrackImpl* getNearestFitToCylinder(float r) ;
    
    // returns a pointer to an LCIO Track whose referece point is the closest to a plane normal to the z axis.
    IMPL::TrackImpl* getNearestFitToZPlane(float z) ;

    // memeber variables 
    EVENT::Track*    _initialLCTrack ;
    IMPL::TrackImpl* _ipPropagatedLCTrack ;


    EVENT::TrackerHitVec _lcioHits ; 

    std::vector< TanagraTrack* > _tanagra_fits ;

    std::vector< IMPL::TrackImpl* > _lcio_tracks ;

    bool _fit_done ;

  } ;

} // end of marlin_delphiF77 namespace 

#endif
