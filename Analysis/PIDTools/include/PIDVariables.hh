/*
 * PIDVariables.hh
 *
 * PIDVariable
 * Helper class that contains and manages a variable
 * used by PIDTools processors.
 *
 * PIDVariables
 * Helper class that creates and manages a map of PIDVariable
 * objects, Update()s it from a reconstructed particle etc.
 * This class takes care and standardises the calculation of
 * variables used for PID.
 *
 *  Created on: Dec 23, 2015
 *      Author: S. Lukic
 */

#ifndef PIDVARIABLES_H_
#define PIDVARIABLES_H_ 1

#include "TVector3.h"
#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "EVENT/ReconstructedParticle.h"

#include "PIDParticles.hh"

// test comment

class PIDVariable {
public:
  PIDVariable(const char*_name, const char *_axisTitle, const double _loLim, const double _hiLim) :
    value(0.), name(_name), axisTitle(_axisTitle), loLim(_loLim), hiLim(_hiLim) {};
  ~PIDVariable() {};

  double Value() const { return value; };
  const char *Name() const { return name; };
  const char *AxisTitle() const { return axisTitle; };
  double LoLim() const { return loLim; };
  double HiLim() const { return hiLim; };

  void SetValue(double val) { value = val; };

private:
  double value;
  const char *name;
  const char *axisTitle;
  const double loLim, hiLim;
};


class PIDVariables {

public:
  PIDVariables();
  PIDVariables(EVENT::ReconstructedParticle*);
  ~PIDVariables() { delete particlePars; };

  // WARNING: If editing enum varType, don't forget
  // to update the First-last vars for algorithms (below)!
  // And update the PIDVariables default consturctor
  enum varType {
    CALO_Total,
    CALO_EFrac,
    CALO_MuSys,
    CLUSHAPE_Chi2,
    CLUSHAPE_DiscrL,
    CLUSHAPE_DiscrT,
    CLUSHAPE_xl20,
    DEDX_Chi2electron,
    DEDX_Chi2muon,
    DEDX_Chi2pion,
    DEDX_Chi2kaon,
    DEDX_Chi2proton,
    N_VarTypes
  };

  typedef std::map<varType, PIDVariable> VarMap;
  typedef PIDParticles::ParameterMap ParticleMap;

  // First and last variables used in different PID algorithms
  static const varType basic_first;
  static const varType calo_beyond;
  static const varType basic_beyond;
  static const varType clushape_first;
  static const varType clushape_beyond;
  static const varType dEdx_first;
  static const varType dEdx_beyond;

  // Various cuts used when calculating the variables
  static const double caloCut;
  static const double muSysCut;
  static const double muSysPCut;

  static const double dEdx_MIP;


  double GetVariable(varType i) const { return varMap.at(i).Value(); };
  const VarMap* GetMap() const { return &varMap; };
  const ParticleMap* GetParticleMap() const { return particlePars; };
  double GetDEdx() const { return dEdx; };
  double GetP() const { return p; };

  void Update(EVENT::ReconstructedParticle*);
  void Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p);

protected:
  VarMap varMap;
//  double values[N_VarTypes]; This doesn't work because the compiler does not know whether varType are contiguous int

private:
  void PopulateMap();
  double dEdx, p;

  double get_dEdxChi2(PIDParticle* hypothesis) const;
  double BetheBloch(PIDParticle* hypothesis) const;

  // This class maintains its own separate particle parameter map with
  // just the basic properties
  ParticleMap *particlePars;

};



#endif // PIDVARIABLES_H_
