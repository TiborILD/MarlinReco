/*
 * MvaPidProcessor.hh
 *
 * Distinguish particle type using MVA with sensitive variables formed
 * from calorimetric deposits, shower shapes and specific energy
 * loss and momentum measured in the tracker
 *
 *  Created on: Feb 4, 2016
 *      Author: Strahinja Lukic
 */

#ifndef MVAPIDPROCESSOR_HH_
#define MVAPIDPROCESSOR_HH_

#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>
#include <UTIL/PIDHandler.h>
// Moving the reader to PIDParticles
// #include "TMVA/Reader.h"

#include "PIDParticles.hh"
#include "PIDVariables.hh"

using UTIL::PIDHandler;

using namespace lcio ;
using namespace marlin ;

using PIDParticles::MVAHypothesesMap;
using PIDParticles::particleType;

class LowMomentumMuPiSeparationPID_BDTG;

class MvaPidProcessor : public Processor{
public:
  virtual Processor*  newProcessor() { return new MvaPidProcessor ; }
  MvaPidProcessor();
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

  static const char *algoName;

private:

  // Updates _variables, sets _mvaVars, evaluates MVA, selects best hypothesis
  // Fills _pidPars
  void Identify(ReconstructedParticle*);
  // Test statistic for making decision if multiple hypotheses make the MVA cut
  double GetQ(hypotheses_c_iterator ith) const { return ith->second.GetQ(); } ;

  hypotheses_c_iterator _bestHypothesis; // Found by Identify()

// Moving this to the hypotheses map
//  typedef std::map<particleType, TMVA::Reader*> ReaderMap;
//  ReaderMap _readerMap;
  PIDVariables *_variables;
  // Temporary copy of the variable values for the MVA reader
  // Workaround because the TMVA::Reader::AddVariable(...)
  // does not take const pointers
  std::map<const char*, float> _mvaVars;
  MVAHypothesesMap *_hypotheses;

  std::string _description;

  LCCollection* _pfoCol;
  // Parameters to write to LCIO (should these be const?)
  // Would it be better if there was an overloaded PIDHandler::setParticleID()
  // that took a std::map<std::string, float>?
  FloatVec _pidPars;
  StringVec _pidParNames;

  PIDHandler *_pidh;

  // Steerables:
  // MVA method used
  std::string _mvaMethod;
  std::string _inputPFOsCollection;
  std::vector<std::string> _weightFileNames;

  // mu-pi separation
  std::vector<std::string> _muPiWeightFileNames; // steerable
  LowMomentumMuPiSeparationPID_BDTG *_mupiPID;

  unsigned int _nEvt, _nPFO, _nUnidentified, _nDecisionQ;
  std::map<particleType, unsigned int> _mapNDecisionTot;
  std::map<particleType, unsigned int> _mapNDecisionQ;
};




#endif /* MVAPIDPROCESSOR_HH_ */
