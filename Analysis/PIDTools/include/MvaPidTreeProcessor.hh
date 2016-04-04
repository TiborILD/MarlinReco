#ifndef MvaPidTreeProcessor_hh

#define MvaPidTreeProcessor_hh 1

#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>

#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include <PIDVariables.hh>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

/*************************************************************
 * MvaPidTreeProcessor
 *
 * Produces a tree with MvaPid variables for particles
 * Useful to study MvaPid variables and their ability
 * to discriminate particle types.
 *
 * Adapted from PIDTree processor by J. List.
 *
 * Strahinja Lukić, April 2016
 *
**************************************************************/


class MvaPidTreeProcessor : public Processor {
  
 public:
 

  virtual Processor*  newProcessor() { return new MvaPidTreeProcessor ; }
  
  
  MvaPidTreeProcessor() ;
//  ~MvaPidTreeProcessor() ;
  
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
  typedef PIDVariables_base::VarVec VarVec;
  typedef PIDVariables_base::VarVec::iterator variable_iterator;
  typedef PIDVariables_base::VarVec::const_iterator variable_c_iterator;


 protected:

  /** LCIO collection names.
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

  PIDVariables_MvaPid _variables;

  vector<vector<double>*> _treeVars;

  vector<double> trueP;
  vector<double> truePt ;
  vector<double> trueTheta ;
  vector<double> truePhi ;
  vector<double> trueCharge;
  vector<int>    truePDG;
  
  vector<bool> isReconstructed;
  vector<double> isSeen;  // store max weight of relation here!
  vector<double> seenP;
  vector<double> seenPt ;
  vector<double> seenTheta ;
  vector<double> seenPhi ;
  vector<double> seenCharge ;


} ;


#endif // MvaPidTreeProcessor_hh



