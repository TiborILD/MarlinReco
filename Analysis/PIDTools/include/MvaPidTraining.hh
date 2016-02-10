/*
 * MvaPidTraining.hh
 *
 * Training processor for MvaPid
 *
 *  Created on: Feb 8, 2016
 *      Author: S. Lukic
 */

#ifndef MVAPIDTRAINING_HH_
#define MVAPIDTRAINING_HH_


/*
 * MvaPidProcessor.hh
 *
 * Distinguish particle type using MVA with sensitive variables formed
 * from calorimetric deposits, shower shapes and specific energy
 * loss and momentum measured in the tracker
 *
 * FIXME: Does this actually need to be a Marlin processor?
 * It does not read any lcio files, only root files. It could be a
 * standalone utility.
 *
 *  Created on: Feb 4, 2016
 *      Author: Strahinja Lukic
 */

#ifndef MVAPIDPROCESSOR_HH_
#define MVAPIDPROCESSOR_HH_

#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>

#include "TMVA/Factory.h"

#include "PIDParticles.hh"
#include "PIDVariables.hh"

using namespace lcio ;
using namespace marlin ;

using PIDParticles::MVAHypothesesMap;
using PIDParticles::particleType;

class LowMomentumMuPiSeparationPID_BDTG;

class MvaPidTraining : public Processor{
public:
  virtual Processor*  newProcessor() { return new MvaPidTraining ; }
  MvaPidTraining();
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();

  typedef MVAHypothesesMap::iterator hypotheses_iterator;
  typedef MVAHypothesesMap::const_iterator hypotheses_c_iterator;

  typedef PIDVariables::VarMap VariableMap;
  typedef PIDVariables::VarMap::const_iterator variable_c_iterator;
  typedef PIDVariables::VarMap::iterator variable_iterator;

protected:

 /** LCIO collection names - steerable
  */
 std::string _trueToReco ;
 std::string _recoToTrue ;
 std::string _mcParticleCollectionName ;
 std::string _trackColName ;
 std::string _pandoraPFOs ;



private:

  std::string _description;

  TMVA::Factory *_factory;
  TFile* _rootfile;
  TTree* _tree;

  PIDVariables _variables;
  float _seenP;
  int _truePDG;
  bool _isReconstructed;

  // Steerables:
  // MVA method used
  int _signalPDG;
  std::string _rootFileName;
  std::string _mvaMethod;
  std::string _weightFileName;

  // Counters
  unsigned int _nEvt;
  unsigned int _nMCPtot, _nRec, _nTrkCaloMismatch;
};




#endif /* MVAPIDTRAINING_HH_ */
