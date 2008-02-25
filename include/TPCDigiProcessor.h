/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef TPCDigiProcessor_h
#define TPCDigiProcessor_h 1

#include <marlin/Processor.h>
#include <lcio.h>
#include <string>
#include <gsl/gsl_rng.h>

#ifdef MARLIN_USE_AIDA

#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>


//#define EXPERTCHECKPLOTS 


#ifdef EXPERTCHECKPLOTS
// includes all AIDA header files
#include <AIDA/AIDA.h>
#endif

#endif


using namespace lcio ;
using namespace marlin ;
#ifdef MARLIN_USE_AIDA
using namespace AIDA ;
#endif


/** ====== TPCDigiProcessor ====== <br>
 *
 * This Processor depends on Circle.h from MarlinUtil
 * 
 * Caution: This digitiser presently does not process space-point like SimTrackerHits which have been flagged with CellIDs set to the negetive row number. This must be implemented in future. 
 *Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in r-phi and z. 
 * Double hits are identified but are currently not added to the collection. This may be change 
 * at a later date when criteria for their seperation is defined. The resolutions are defined in 
 * the GEAR stearing file. 
 *
 * Resolution in r-phi is calculated according to the formular <br>
 * sigma_{point}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}
 * Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
 * Cd = 25 ( microns / cm^(1/2) )
 * (this is for B=4T, h is the pad height = pad-row pitch in mm,
 * theta is the polar angle)       
 *
 * At the moment resolution in z assumed to be independent of drift length. <br>
 *
 * The type of TPC TrackerHit is set to 500 via method TrackerHitImpl::setType(int type) <br>
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits in TPC <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in TPC <br>
 * @param CollectionName The name of input SimTrackerHit collection <br>
 * (default name STpc01_TPC)
 * @param RejectCellID0 Whether or not to reject SimTrackerHits with Cell ID 0. Mokka drivers
 * TPC00-TPC03 encoded the pad row number in the cell ID, which should always be non-zero anyway.
 * Drivers TPC04 and TPC05 do not simulate pad rows and thus have the cell ID set to zero for all hits.
 * You will need to set RejectCellID0 to 0 in order to use this processor with these drivers, but note
 * that the implications for track reconstruction are not strictly defined. Mokka driver TPC06 uses
 * a mixed approach with one hit per pad row having non-zero cell ID, extra hits having 0 cell ID.
 * Typically, unless you use TPC04 or TPC05, you should not touch this parameter. <br>
 * (default value 1)
 * @param TPCTrackerHitsCol The name of output collection of TrackerHits <br>
 * (default name TPCTrackerHits) <br>
 * <br>
 * @authors S. Aplin, DESY and A.Raspereza, MPI
 *
 * Changed 7/9/07 so that the const and diffusion resolution terms are taken as processor parameters rather than the gear file.
 * The parameters _pixZ and pixRP were also changed from gear parameters to processor parameters
 * clare.lynch@bristol.ac.uk
 *
 */
class TPCDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new TPCDigiProcessor ; }
  
  
  TPCDigiProcessor() ;
  
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

  /** Input collection name.
   */
  std::string _colName ;
  std::string _TPCTrackerHitsCol ;

  int _rejectCellID0;

  int _nRun ;
  int _nEvt ;

  // gsl random number generator
  gsl_rng * _random ;



  float _pointResoRPhi0; // Coefficient for RPhi point res independant of drift length 
  float _pointResoPadPhi; // Coefficient for the point res dependance on relative phi angle to the pad verticle 
  float _diffRPhi; // Coefficient for the rphi point res dependance on diffusion 
  int   _nEff; // number of effective electrons 


  float _pointResoZ0; // Coefficient Z point res independant of drift length 
  float _diffZ; // Coefficient for the Z point res dependance on diffusion 


  float _pixZ;
  float _pixRP;

#ifdef EXPERTCHECKPLOTS
  IAnalysisFactory * AF;
  ITreeFactory * TRF;
  ITree * TREE;
  IHistogramFactory * HF;
  IHistogram1D * phiDiffHisto;
  IHistogram1D * thetaDiffHisto;
  IHistogram1D * phiRelHisto;
  IHistogram1D * thetaRelHisto;

  IHistogram1D * phiDistHisto;
  IHistogram1D * phiPullHisto;
  IHistogram1D * rDiffHisto;
#endif

  //FIXME: Cathode is hard coded
  const static double _cathode;

} ;

  //FIXME: Cathode is hard coded
const double TPCDigiProcessor::_cathode=5.0/2.0; // cathode is 5mm thick 

#endif



