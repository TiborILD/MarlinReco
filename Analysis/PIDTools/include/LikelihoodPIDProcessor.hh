#ifndef LikelihoodPIDProcessor_hh
#define LikelihoodPIDProcessor_hh 1

#include <string>
#include <vector>
#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>

#include "PIDParticles.hh"
#include "PIDVariables.hh"

using namespace lcio ;
using namespace marlin ;

class LikelihoodPID;
class LowMomentumMuPiSeparationPID_BDTG;

class LikelihoodPIDProcessor : public Processor{
public:
  virtual Processor*  newProcessor() { return new LikelihoodPIDProcessor ; }
  LikelihoodPIDProcessor();
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();

  typedef PIDParticles::LLHypothesesMap ParticleMap;
  typedef PIDParticles::LLHypothesesMap::iterator particle_iterator;
  typedef PIDParticles::LLHypothesesMap::const_iterator particle_c_iterator;
  typedef PIDVariables::VarMap::const_iterator variable_c_iterator;
  typedef PIDVariables::VarMap::iterator variable_iterator;

 
private:
  void createParticleIDClass(ReconstructedParticle *part, PIDHandler &pidh, int algoID, float MVAoutput);
  
  LikelihoodPID *_myPID;
  std::string _description;
  std::string _inputPFOsCollection;
  std::string _PDFName;
  std::vector<std::string> _weightFileName;
  EVENT::FloatVec _energyBoundary;
  LCCollection* _pfoCol;
  std::vector<std::string> _parNames;
  std::string _algoName;

  LowMomentumMuPiSeparationPID_BDTG *_mupiPID;

  std::vector<float> _particlePriors;

  bool _basicFlg;
  bool _dEdxFlg;
  bool _showerShapesFlg;

  unsigned int _nEvt;
};

#endif 
