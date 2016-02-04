#ifndef LIKELIHOODPID_hh
#define LIKELIHOODPID_hh 1

#include <string>
#include <sstream>
#include "TLorentzVector.h"
#include "TVectorT.h"
#include "EVENT/Cluster.h"
#include "EVENT/Track.h"
#include <UTIL/PIDHandler.h>
#include "TFile.h"
#include "TH1.h"

#include "PIDParticles.hh"
#include "PIDVariables.hh"



class LikelihoodPID{
public:
  LikelihoodPID(std::string fname, std::vector<float> priors);
  ~LikelihoodPID() ;
   
  // Flag masks for the user to tell which algorithms to use
  static const short MASK_Basic;
  static const short MASK_dEdx;
  static const short MASK_Shapes;

  // logL value that you get if L=0 or sensitive variable is out of bounds
  static const double lowestLogL;

  typedef PIDParticles::particleType parType;
  typedef PIDParticles::LLHypothesesMap ParticleMap;
  typedef PIDParticles::LLHypothesesMap::iterator particle_iterator;
  typedef PIDParticles::LLHypothesesMap::const_iterator particle_c_iterator;
  typedef PIDVariables::varType varType;
  typedef PIDVariables::VarMap::const_iterator variable_c_iterator;
  typedef PIDVariables::VarMap::iterator variable_iterator;

  Int_t Classification(IMPL::ReconstructedParticleImpl*);
  Int_t Classification();

  const ParticleMap * GetParticlePars() const { return particlePars; };
  int GetBestType() const;
  int GetBestPDG() const;
  double GetBestLikelihood() const;
  double GetBestProbability() const;
  // This function returns the prior probability, which is a property of the hypothesis
  float GetPrior(parType hypothesis) const { return particlePars->at(hypothesis).prior; };
  // This function returns the stored posterior probability of the hypothesis, calculated
  // in the CalcPosteriors() function
  double GetPosterior(parType hypothesis) const { return particlePars->at(hypothesis).Posterior(); };
  // This function returns the stored cumulative Log(L) of the hypothesis, calculated
  // in the CalcPosteriors() function
  double GetLogL(parType hypothesis) const { return particlePars->at(hypothesis).LogL(); };
  Double_t getCorrEnergy(TLorentzVector pp, EVENT::Track* trk, EVENT::ClusterVec& cluvec);
  Double_t getCorrEnergy(TLorentzVector pp, parType);

  //for leptonIDlikelihood
  Double_t get_dEdxDist(parType);

  void setBestParticle(parType best) {
    if (best==PIDParticles::lowEmuon) best=PIDParticles::muon;
    bestParticle = particlePars->find(best);
  };

  void   setBasicFlg(Bool_t flg)        { if(flg) {_algoFlags |=  MASK_Basic; }
                                          else { _algoFlags &= ~MASK_Basic; } }
  void   setdEdxFlg(Bool_t flg)         { if(flg) {_algoFlags |=  MASK_dEdx;  }
                                          else { _algoFlags &= ~MASK_dEdx;  } }
  void setShowerShapesFlg(Bool_t flg)   { if(flg) {_algoFlags |=  MASK_Shapes;}
                                          else {_algoFlags &= ~MASK_Shapes;} }

protected:

  void CalcPosteriors();
  // This function returns Log(L) for the hypothesis, based on the value "value" of the variable "valtype"
  const Double_t LogL(parType, varType valtype);
  Double_t getPenalty(Int_t ptype, Int_t hypothesis, Double_t p);

  // Particle to process
  IMPL::ReconstructedParticleImpl* _particle;

//  Double_t par[5][5];
  TFile* fpdf;
  TH1* pdf[PIDParticles::nParticleTypes][PIDVariables::N_VarTypes];

  // particle properties
  ParticleMap *particlePars;
  particle_c_iterator bestParticle;

  PIDVariables variables;

  //for shower profile
  EVENT::FloatVec shapes;

  short _algoFlags;
};

#endif 
