#ifndef DelphiTrackFitProcessor_h
#define DelphiTrackFitProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <UTIL/LCRelationNavigator.h>

using namespace lcio ;
using namespace marlin ;


/**  Track Refitter processor for marlin. Refits an input track collection using the Delphi Kalman Filter implemented in F77, producing a new collection of tracks
 * 

 *  <h4>Input - Prerequisites</h4>
 *  Needs a collection of LCIO Tracks. 
 *
 *  <h4>Output</h4> 
 *  Refitted LCIO Track Collection
 * 
 * @param InputTrackCollectionName Name of the Track collection to be refitted
 * @param OutputTrackCollectionName Name of the refitted Track collection to be refitted
 * 
 * @author S. J. Aplin, DESY
 */

class DelphiTrackFitProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DelphiTrackFitProcessor ; }
  
  
  DelphiTrackFitProcessor() ;
  
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
  LCCollection* GetCollection( LCEvent * evt, std::string colName ) ;
  
  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name ;

  /** refitted track collection name.
   */
  std::string _output_track_col_name ;

  int  _F77DebugLevel ;
  bool _MSOn ;

  int _n_run ;
  int _n_evt ;


} ;

#endif



