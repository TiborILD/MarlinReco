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

/**
*/

class PIDParticle {
public:
  PIDParticle (const char *_name, int _pdg, double _mass,
               double _prior, const double* _BBpars) :
    pdg(_pdg), mass(_mass), prior(_prior),
    posterior(0), logL(0), threshold(0), dEdxDist(0), name(_name)
  {
    for (int i=0; i<5; i++) BBpars[i] = _BBpars[i];
  };
  ~PIDParticle() {};

  const int pdg;
  const double mass;
  const double prior;

  const double * GetBBpars() const { return BBpars; };
  double Posterior() const { return posterior; }
  double LogL() const { return logL; }
  double Threshold() const { return threshold; }
  double DEdxDist() const { return dEdxDist; }

  void SetPosterior(double _posterior) { posterior = _posterior; };
  void SetThreshold(double _threshold) { threshold = _threshold; };
  void SetDEdxDist(double _dEdxDist) { dEdxDist = _dEdxDist; } ;

  void AddLogL(double logLpartial) { logL += logLpartial; };

  void ResetLogL() { logL=0.; };

  const char* Name() const {return name;};

private:
  // There is no setter for BBpars - they are set in the constructor and do not change!
  double BBpars[5];
  double posterior, logL;
  double threshold;
  double dEdxDist;

  const char *name;
};



struct PIDParticles {

  enum particleType {electron, muon, pion, kaon, proton, lowEmuon, nParticleTypes};

  typedef std::map<particleType, PIDParticle> ParameterMap;

  static ParameterMap* CreateMap(std::vector<float> priors);

};


#endif /* PIDPARTICLES_HH_ */
