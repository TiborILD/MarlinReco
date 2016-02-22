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
#include "TString.h"
#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include "EVENT/ReconstructedParticle.h"

#include "PIDParticles.hh"

class TRandom3;

class PIDVariable {
public:
  PIDVariable(const char*name, const char *description, const char *unit,
      const int nBinsP, const double *pBins,
      const double loLim, const double hiLim, const int nBinsV=50) :
    value(0.), _name(name), _description(description), _unit(unit),
    _nBinsV(nBinsV), _loLim(loLim), _hiLim(hiLim),
    _nBinsP(nBinsP), _pBins(pBins) {};
  ~PIDVariable() {};

  float Value() const { return value; };
  const char *Name() const { return _name; };
  const float *Address() const { return &value; };

  const char *AxisTitle() const { if(_unit[0] == '\0') return _description;
                                  else return Form("%s (%s)", _description, _unit); };
  const char *Description() const { return _description; };
  const char *Unit() const { return _unit; } ;
  int NBinsV() const { return _nBinsV; };
  double LoLim() const { return _loLim; };
  double HiLim() const { return _hiLim; };
  int NBinsP() const { return _nBinsP; };
  const double *PBins() const { return _pBins; };

  void SetValue(double val) { value = val; };

protected:
  float value;
  const char *_name;
  const char *_description;
  const char *_unit;
  // No. of bins and limits of the variable distribution histogram
  const int _nBinsV;
  const double _loLim, _hiLim;
  // Slicing in reconstructed momentum
  const int _nBinsP;
  const double *_pBins;
};


class PIDVariables {

public:
  PIDVariables();
  PIDVariables(EVENT::ReconstructedParticle*);
  ~PIDVariables();


  static const short MASK_EmptyClusters  = 1 ;
  static const short MASK_EmptyTracks   = 1 << 1;
  static const short MASK_EmptyShapes = 1 << 2;
  static const short MASK_ZerodEdx   = 1 << 3;

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


  float GetValue(varType i) const { return varMap.at(i).Value(); };
  const VarMap* GetMap() const { return &varMap; };
  const ParticleMap* GetParticleMap() const { return particlePars; };
  float GetDEdx() const { return dEdx; };
  float GetP() const { return p; };


  int Update(EVENT::ReconstructedParticle*);
  int Update(const EVENT::ClusterVec, const EVENT::TrackVec, const TVector3 p);
  void SetOutOfRange();

protected:
  VarMap varMap;

private:
  void PopulateMap();
  float dEdx, p;

  float get_dEdxChi2(PIDParticles::PIDParticle_base* hypothesis) const;
  float get_dEdxSignedLogChi2(PIDParticles::PIDParticle_base* hypothesis) const;
  double BetheBloch(PIDParticles::PIDParticle_base* hypothesis) const;

  // This class maintains its own separate particle parameter map with
  // just the basic properties
  ParticleMap *particlePars;

  TRandom3 *_rand;
};



#endif // PIDVARIABLES_H_
