#ifndef TruthTracker_h
#define TruthTracker_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

using namespace lcio ;
using namespace marlin ;

namespace EVENT{
  class MCParticle ;
  class TrackerHit ;
  class Track ;
}

namespace UTIL{
  class LCRelationNavigator ;
}

/**  Track creation based on MC truth. 
 * 

 *  <h4>Input - Prerequisites</h4>
 *  Needs a collections of LCIO TrackerHits. 
 *
 *  <h4>Output</h4> 
 *  LCIO Track Collection
 * 
 * @param InputTrackerHitCollectionName Name of the TrackerHit collections 
 * @param OutputTrackCollectionName Name of the output Track collection
 * 
 * @author S. J. Aplin, DESY
 */

class TruthTracker : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TruthTracker ; }
  
  
  TruthTracker() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /* helper function to get collection using try catch block */
  LCCollection* GetCollection(  LCEvent * evt, std::string colName ) ;
  
  /* helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations( LCEvent * evt, std::string RelName ) ;

  /* helper function to get collection using try catch block */
  void SetupInputCollections( LCEvent * evt ) ;


  /** input MCParticle collection name.
   */
  std::string _colNameMC ;

  LCCollection* _colMCP;

  /** input SimTrackerHit collections
   */
  std::string _colNameSimTrkHitsVTX;
  std::string _colNameSimTrkHitsSIT;
  std::string _colNameSimTrkHitsFTD;
  std::string _colNameSimTrkHitsTPC;

  LCCollection* _VTXSimHits;
  LCCollection* _SITSimHits;
  LCCollection* _FTDSimHits;
  LCCollection* _TPCSimHits;

  /** input TrackerHit collections
   */
  std::string _colNameTrkHitsVTX ;
  std::string _colNameTrkHitsFTD ;
  std::string _colNameTrkHitsSIT ;
  std::string _colNameTrkHitsTPC ;

  LCCollection* _VTXTrkHits;
  LCCollection* _SITTrkHits;
  LCCollection* _FTDTrkHits;
  LCCollection* _TPCTrkHits;

  /** output track collection 
   */
  std::string _output_track_col_name ;
  
 /** Output track relations
  */
  std::string _output_track_rel_name ;


  std::map< MCParticle*, std::vector<TrackerHit*> > _MCParticleTrkHitMap;

  int _nMCP;

  int _nVTXSimHits ;
  int _nSITSimHits ;
  int _nFTDSimHits ;
  int _nTPCSimHits ;
  
  int _nVTXTrkHits ;
  int _nSITTrkHits ;
  int _nFTDTrkHits ;
  int _nTPCTrkHits ;

  int _nEventPrintout ;
  int _n_run ;
  int _n_evt ;

  float _MCpThreshold ;

} ;

#endif



