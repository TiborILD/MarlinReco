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


#include <marlin/Processor.h>



#include "PIDParticles.hh"
#include "PIDVariables.hh"

using namespace lcio ;
using namespace marlin ;

using PIDParticles::MVAHypothesesMap;
using PIDParticles::particleType;

class LowMomentumMuPiSeparationPID_BDTG;
class TTree;
class TH1F;


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

  typedef PIDVariables::varType variableType;
  typedef PIDVariables::VarMap VariableMap;
  typedef PIDVariables::VarMap::const_iterator variable_c_iterator;
  typedef PIDVariables::VarMap::iterator variable_iterator;

protected:

 /** LCIO collection names - steerable
  */
 std::string _trueToReco ;
 std::string _recoToTrue ;
 std::string _mcParticleCollectionName ;



private:

  std::string _description;

  TTree* _tree;

  std::map<variableType, float> _trainingVars;
  PIDVariables _variables;
  float _seenP;
  int _truePDG;
  bool _isReconstructed, _hasClusters, _hasdEdx;

  // Steerables:
  int _signalPDG;
  std::string _mvaResponseFileName;
  std::string _mvaMethod;
  std::string _mvaMethodOptions;
  std::string _weightFileName;

  // Counters
  unsigned int _nEvt;
  unsigned int _nMCPtot, _nRec, _nTrkCaloMismatch;
  unsigned int _nEmptyClusters, _nEmptyTracks, _nEmptyShapes, _nZerodEdx;

  // Number of channels for the histogram of Q-statistic
  static const int nChanQ = 1000;
};




#endif /* MVAPIDTRAINING_HH_ */
