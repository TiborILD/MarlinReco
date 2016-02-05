/*
 * PIDParticles.hh
 *
 * PIDParticle
 * Helper class that contains and manages parameters for a
 * particle type, such as constant inherent particle properties
 * and variable parameters used by PIDTools processors
 *
 * PIDParticles
 * Helper struct (effectively no more than a namespace)
 * to create a map of PIDParticle objects for use by PIDTools
 * processors etc.
 *
 *  Created on: Dec 23, 2015
 *      Author: S. Lukic
 */

#ifndef PIDPARTICLES_HH_
#define PIDPARTICLES_HH_

#include "TMath.h"
#include <map>


namespace PIDParticles
{

class PIDParticle_base {
public:
  PIDParticle_base (const char *name, int _pdg, double _mass, const double* BBpars) :
    pdg(_pdg), mass(_mass), _name(name)
  {
    for (int i=0; i<5; i++) _BBpars[i] = BBpars[i];
  }
  PIDParticle_base (const PIDParticle_base &base) :
    pdg(base.pdg), mass(base.mass), _name(base.Name())
  {
    for (int i=0; i<5; i++) _BBpars[i] = base.GetBBpars()[i];
  }
  ~PIDParticle_base() {}

  const int pdg;
  const double mass;

  const double * GetBBpars() const { return _BBpars; }
  const char* Name() const {return _name;}

private:
  // There is no setter for BBpars - they are set in the constructor and do not change!
  double _BBpars[5];

  const char *_name;
};


/*** Predefined intrinsic particle properties -- all constants in one place ***/

enum particleType {electron, muon, pion, kaon, proton, lowEmuon, nParticleTypes};

static const double BBparsElectron[5] = {-2.40638e-03, 7.10337e-01, 2.87718e-01, -1.71591e+00, 0.0};
static const double BBparsMuon[5] = {8.11408e-02, 9.92207e-01, 7.58509e+05, -1.70167e-01, 4.63670e-04};
static const double BBparsPion[5] = {8.10756e-02, -1.45051e+06, -3.09843e+04, 2.84056e-01, 3.38131e-04};
static const double BBparsKaon[5] = {7.96117e-02, 4.13335e+03, 1.13577e+06, 1.80555e-01, -3.15083e-04};
static const double BBparsProton[5] = {7.78772e-02, 4.49300e+04, 9.13778e+04, 1.50088e-01, -6.64184e-04};

static const PIDParticle_base electronProperties("electron", 11, .000510998, BBparsElectron);
static const PIDParticle_base muonProperties("muon", 13, .105658, BBparsMuon);
static const PIDParticle_base pionProperties("pion", 211, .139570, BBparsPion);
static const PIDParticle_base kaonProperties("kaon", 321, .493677, BBparsKaon);
static const PIDParticle_base protonProperties("proton", 2212, .938272, BBparsProton);


/*** Derived classes ***/

// PIDParticle with parameters for LikelihoodPID

class LLPIDHypothesis : public PIDParticle_base {
public:

  LLPIDHypothesis (const char *_name, int _pdg, double _mass,
               float _prior, const double* _BBpars) :
    PIDParticle_base(_name, _pdg, _mass, _BBpars),
    prior(_prior), _posterior(0), _logL(0), _threshold(0)
  {  }

  LLPIDHypothesis (const PIDParticle_base &base, float _prior) :
    PIDParticle_base(base),
    prior(_prior), _posterior(0), _logL(0), _threshold(0)
  {  }

  ~LLPIDHypothesis() {}

  const float prior;

  double Posterior() const { return _posterior; }
  double LogL() const { return _logL; }
  double Threshold() const { return _threshold; }

  void SetPosterior(double posterior) { _posterior = posterior; }
  void SetThreshold(double thr) { _threshold = thr; }

  void AddLogL(double logLpartial) { _logL += logLpartial; }

  void ResetLogL() { _logL=0.; }

private:
  double _posterior, _logL;
  double _threshold;
};



// PIDParticle with parameters for MVA PID

class MVAPIDHypothesis : public PIDParticle_base {
public:

  MVAPIDHypothesis (const char *_name, int _pdg, double _mass, const double* _BBpars, const float mvaCut=0.) :
    PIDParticle_base(_name, _pdg, _mass, _BBpars),
    _mva(0), _q(0), _mvaCut(mvaCut)
  {  }

  MVAPIDHypothesis (const PIDParticle_base &base, const float mvaCut=0.) :
    PIDParticle_base(base),
    _mva(0), _q(0), _mvaCut(mvaCut)
  {  }

  ~MVAPIDHypothesis() {};

  double GetMVAout() const { return _mva; }
  double GetQ() const { return _q; }
  const float GetMVAcut() const { return _mvaCut; }

  void SetMVAout(double mva) { _mva = mva; }
  void SetQ(double q) { _q = q; }

private:
  double _mva, _q;
  const float _mvaCut;
};




typedef std::map<particleType, PIDParticle_base> ParticleMap;
typedef std::map<particleType, LLPIDHypothesis> LLHypothesesMap;
typedef std::map<particleType, MVAPIDHypothesis> MVAHypothesesMap;

ParticleMap* CreateParticleMap();
LLHypothesesMap* CreateLLPIDMap(std::vector<float> priors);
MVAHypothesesMap* CreateMVAPIDMap();

}


#endif /* PIDPARTICLES_HH_ */