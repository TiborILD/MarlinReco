#ifndef PIDvarPDF_hh
#define PIDvarPDF_hh 1

#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "PIDVariables.hh"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

/*************************************************************
 * PIDvarPDF processor
 *
 * Produces histograms of variables such as used by the
 * LikelihoodPIDProcessor.
 * TODO: Adapt histogram names for actual use with
 * LikelihoodPIDProcessor.
 *
 * Fills a ROOT tree with relevant information on particles
 * found in the MCParticle and ReconstructedParticle
 * collections -- for a study of the variable distributions.
 *
 * Produces dummy histograms with associated TF1 functions
 * containing the Bethe-Bloch curve for each particle
 * (as currently implemented in the PIDParticles class) --
 * for visualisation and comparison to dEdx vs. p scatter
 * plots.
 *
 * Adapted from PIDTree processor by J. List.
 *
 * Strahinja Lukić, Jan. 2016
 *
**************************************************************/


class PIDvarPDF : public Processor {
  
 public:
 

  virtual Processor*  newProcessor() { return new PIDvarPDF ; }
  
  
  PIDvarPDF() ;
  ~PIDvarPDF() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

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

  typedef PIDParticles::particleType ParticleType;
  typedef PIDParticles::ParticleMap ParticleMap;
  typedef PIDParticles::ParticleMap::iterator particle_iterator;
  typedef PIDParticles::ParticleMap::const_iterator particle_c_iterator;
  typedef PIDVariables::VarMap variableMap;
  typedef PIDVariables::VarMap::iterator variable_iterator;
  typedef PIDVariables::VarMap::const_iterator variable_c_iterator;


 protected:

  /** Input collection name.
   */
  std::string _trueToReco ;
  std::string _recoToTrue ;
  std::string _mcParticleCollectionName ;
  std::string _trackColName ;
  std::string _pandoraPFOs ;

  int nEvt ;

 private:

  // declaration of trees
  TTree *varTree ;
  int    nMCParticles;
  
  int nMCPtot, nRec, nTrkCaloMismatch;

  PIDVariables pidVars;
  vector<double> sensitiveVars[PIDVariables::N_VarTypes];
  TH1 *sensVarHistos[PIDParticles::nParticleTypes][PIDVariables::N_VarTypes];
  TH1F *bbHistos[PIDParticles::nParticleTypes];
  TF1 *bbFunction[PIDParticles::nParticleTypes];

  vector<double> trueP;
  vector<double> truePt ;
  vector<double> trueTheta ;
  vector<double> truePhi ;
  vector<double> trueCharge;
  vector<double> trued0;
  vector<double> truez0;
  vector<int>    truePDG;
  vector<int>    trueMother;
  
  vector<bool> isReconstructed;
  vector<double> isSeen;  // store max weight of relation here!
  vector<double> seenP;
  vector<double> seenPt ;
  vector<double> seenTheta ;
  vector<double> seenPhi ;
  vector<double> seenDEdx ;
  vector<double> seenCharge ;


} ;


#endif // PIDvarPDF_hh



