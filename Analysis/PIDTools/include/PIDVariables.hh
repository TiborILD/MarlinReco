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
  PIDVariable(const char*_name, const char *_axisTitle,
      const int _nBinsP, const double *_pBins,
      const double _loLim, const double _hiLim, const int _nBinsV=50) :
    value(0.), name(_name), axisTitle(_axisTitle),
    nBinsV(_nBinsV), loLim(_loLim), hiLim(_hiLim),
    nBinsP(_nBinsP), pBins(_pBins) {};
  ~PIDVariable() {};

  float Value() const { return value; };
  const char *Name() const { return name; };
  const float &Reference() const { return value; };

  const char *AxisTitle() const { return axisTitle; };
  int NBinsV() const { return nBinsV; };
  double LoLim() const { return loLim; };
  double HiLim() const { return hiLim; };
  int NBinsP() const { return nBinsP; };
  const double *PBins() const { return pBins; };

  void SetValue(double val) { value = val; };

private:
  float value;
  const char *name;
  const char *axisTitle;
  // No. of bins and limits of the variable distribution histogram
  const int nBinsV;
  const double loLim, hiLim;
  // Slicing in reconstructed momentum
  const int nBinsP;
  const double *pBins;
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
  typedef PIDParticles::ParticleMap ParticleMap;

  // First and last variables used in different PID algorithms
  static const varType basic_first;
  static const varType calo_beyond;
  static const varType basic_beyond;
  static const varType clushape_first;
  static const varType clushape_beyond;
  static const varType dEdx_first;
  static const varType dEdx_beyond;

  // Various cuts used when calculating the variables
  static const double ptCut; // Avoid division by zero for "CaloTotal"
  static const double caloCut;
  static const double muSysCut;
  static const double muSysPCut;

  static const double dEdx_MIP;


  double GetValue(varType i) const { return varMap.at(i).Value(); };
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

  double get_dEdxChi2(PIDParticles::PIDParticle_base* hypothesis) const;
  double BetheBloch(PIDParticles::PIDParticle_base* hypothesis) const;

  // This class maintains its own separate particle parameter map with
  // just the basic properties
  ParticleMap *particlePars;

};



#endif // PIDVARIABLES_H_
