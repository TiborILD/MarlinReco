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
#include "TMVA/Reader.h"

#include "PIDParticles.hh"
#include "PIDVariables.hh"

using namespace lcio ;
using namespace marlin ;

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

  typedef PIDParticles::MVAHypothesesMap ParticleMap;
  typedef PIDParticles::MVAHypothesesMap::iterator particle_iterator;
  typedef PIDParticles::MVAHypothesesMap::const_iterator particle_c_iterator;

  typedef PIDVariables::VarMap::const_iterator variable_c_iterator;
  typedef PIDVariables::VarMap::iterator variable_iterator;


private:

  // Test statistic for making decision if multiple hypotheses make the MVA cut
  double GetQ(particle_iterator);

  TMVA::Reader *_reader;
  PIDVariables *_variables;
  ParticleMap *_hypotheses;

  std::string _description;

  LCCollection* _pfoCol;
  // Parameters to write to LCIO (should this be static const?)
  std::map<std::string, float> _pidPars;

  // Steerables:
  // MVA method used
  std::string _mvaMethod;
  std::string _inputPFOsCollection;
  std::string _weightFileName;

  // mu-pi separation
  std::vector<std::string> _muPiWeightFileNames; // steerable
  EVENT::FloatVec _energyBoundary; // built in
  LowMomentumMuPiSeparationPID_BDTG *_mupiPID;

  unsigned int _nEvt;
};




#endif /* MVAPIDPROCESSOR_HH_ */
