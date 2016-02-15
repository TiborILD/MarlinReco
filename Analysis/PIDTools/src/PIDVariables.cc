#include "PIDVariables.hh"
//#include "PIDParticles.hh"
#include "LikelihoodPID.hh"

#include "TRandom3.h"

/*******************************************************
 *
 *   Implementation of the PIDVariables class
 *
 ******************************************************/

const PIDVariables::varType PIDVariables::basic_first = PIDVariables::CALO_Total;
const PIDVariables::varType PIDVariables::calo_beyond = PIDVariables::CALO_MuSys;
const PIDVariables::varType PIDVariables::basic_beyond = PIDVariables::CLUSHAPE_Chi2;
const PIDVariables::varType PIDVariables::clushape_first = PIDVariables::CLUSHAPE_Chi2;
const PIDVariables::varType PIDVariables::clushape_beyond = PIDVariables::DEDX_Chi2electron;
const PIDVariables::varType PIDVariables::dEdx_first = PIDVariables::DEDX_Chi2electron;
const PIDVariables::varType PIDVariables::dEdx_beyond = PIDVariables::N_VarTypes;

// TODO: Optimize these cuts
const double PIDVariables::ptCut = 0.01; // 1/10 of the minimum reconstructible pT at ILD
const double PIDVariables::caloCut = 0.05; // Minimum ECAL+HCAL to calculate ECAL/(ECAL+HCAL)
const double PIDVariables::muSysCut = 0.01;
const double PIDVariables::muSysPCut = 5.;

const double PIDVariables::dEdx_MIP = 1.35e-7;

PIDVariables::PIDVariables() : dEdx(0), p(0) {
  // Create map of particle properties.
  particlePars = PIDParticles::CreateParticleMap();
  _rand = new TRandom3;

  PopulateMap();
}

PIDVariables::PIDVariables(EVENT::ReconstructedParticle* _particle) {
  // Create map of particle properties. This must be done first
  particlePars = PIDParticles::CreateParticleMap();
  _rand = new TRandom3;

  PopulateMap();
  Update(_particle);
}

void PIDVariables::PopulateMap() {
  const int nBinsCommon = 8;
  double * binsCommon = new double[nBinsCommon+1];
  binsCommon[0] =  0.;
  binsCommon[1] =  1.;
  binsCommon[2] =  2.;
  binsCommon[3] =  5.;
  binsCommon[4] = 10.;
  binsCommon[5] = 20.;
  binsCommon[6] = 50.;
  binsCommon[7] = 100.;
  binsCommon[8] = 5000.;
  varMap.insert(std::pair<varType, PIDVariable>(CALO_Total, PIDVariable("CaloTotal", "(ECAL+HCAL)/p", "", nBinsCommon, binsCommon, 0., 2., 40)));
  varMap.insert(std::pair<varType, PIDVariable>(CALO_EFrac, PIDVariable("CaloEFrac", "ECAL/(ECAL+HCAL)", "", nBinsCommon, binsCommon, 0., 1., 50)));
  varMap.insert(std::pair<varType, PIDVariable>(CALO_MuSys, PIDVariable("CaloMuSys", "#mu system deposit", "GeV", nBinsCommon, binsCommon, 0., 10., 50)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_Chi2, PIDVariable("CluShapeChi2", "Cluster shape #chi^{2}", "", nBinsCommon, binsCommon, -2., 20., 44)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_DiscrL, PIDVariable("DiscrepancyL", "Shower max / EM shower max", "", nBinsCommon, binsCommon, -30., 80., 55)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_DiscrT, PIDVariable("DiscrepancyT", "Absorption length", "R_{m}", nBinsCommon, binsCommon, 0., 1., 50)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_xl20, PIDVariable("Xl20", "xl20", "?", nBinsCommon, binsCommon, 0., 50., 50)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2electron, PIDVariable("dEdxChi2electron", "#chi^{2}_{dE/dx} (electron)", "", 0, NULL, -50., 50., 100)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2muon, PIDVariable("dEdxChi2muon", "#chi^{2}_{dE/dx} (#mu)", "", 0, NULL, -50., 50., 100)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2pion, PIDVariable("dEdxChi2pion", "#chi^{2}_{dE/dx} (#pi)", "", 0, NULL, -50., 50., 100)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2kaon, PIDVariable("dEdxChi2kaon", "#chi^{2}_{dE/dx} (K)", "", 0, NULL, -50., 50., 100)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2proton, PIDVariable("dEdxChi2proton", "#chi^{2}_{dE/dx} (proton)", "", 0, NULL, -50., 50., 100)));
}


void PIDVariables::Update(EVENT::ReconstructedParticle* _particle) {
  EVENT::ClusterVec cluvec=_particle->getClusters();
  EVENT::TrackVec trk = _particle->getTracks();
  TVector3 p3(_particle->getMomentum());

  Update(cluvec, trk, p3);
}


void PIDVariables::Update(const EVENT::ClusterVec cluvec, const EVENT::TrackVec trax, const TVector3 p3){

  p = p3.Mag();

  //get deposit energy and shapes
  EVENT::FloatVec shapes;
  double ecal=0., hcal=0., mucal=0.;

  if(cluvec.size()>0){
    for(unsigned int i=0; i<cluvec.size(); i++){
      FloatVec sde = cluvec[i]->getSubdetectorEnergies();
      ecal += sde[0];
      hcal += sde[1];
      mucal+= sde[2];
    }
    // FIXME: DO we really want to use only the zeroth cluster?
    shapes=cluvec[0]->getShape();
  }
  else {
    ecal = hcal = mucal = 0.;
  }

  // IMPROVE HERE: use directly MCTruthTrackRelation and choose track with larger weight for dE/dx
/*    dEdx = 0;
    for (unsigned int itrack = 0; itrack < trax.size(); itrack++) {
      dEdx += trax.at(itrack)->getdEdx();
    }
*/
  if(trax.size() > 0) { dEdx = trax.at(0)->getdEdx(); }
  else { dEdx = -dEdx_MIP; } // Keep an eye on get_dEdxChi2()...

  // Normalise dEdx to MIP
  dEdx /= dEdx_MIP;

  /*************************************/
  /***      Calculate variables      ***/
  /*************************************/

  for (VarMap::iterator it=varMap.begin(); it!=varMap.end(); it++) it->second.SetValue(0);

  // Basic variables
  if(p3.Perp() > ptCut) {
    varMap.at(CALO_Total).SetValue( (ecal +hcal) / p );
  }
  else { varMap.at(CALO_Total).SetValue( -1. ); }

  if(ecal+hcal > caloCut) {
    // Avoid ECAL fraction == 1.0 (1.0 goes into the overflow bin)
    varMap.at(CALO_EFrac).SetValue( ecal/(ecal+hcal) - 1.e-6);
  }
  else { varMap.at(CALO_EFrac).SetValue( -1. ); }

  // Introducing a very small additive spread to avoid problems with TMVA for electrons
  varMap.at(CALO_MuSys).SetValue(mucal+_rand->Gaus(0.,1.e-6));

  // Shower shapes
  if(shapes.size()!=0){
    varMap.at(CLUSHAPE_Chi2).SetValue(shapes[0]);
    varMap.at(CLUSHAPE_DiscrL).SetValue(TMath::Sign(float(TMath::Log(fabs(shapes[5])+FLT_MIN)), shapes[5]));
//    varMap.at(CLUSHAPE_DiscrL).SetValue(shapes[5]);
//    varMap.at(CLUSHAPE_DiscrT).SetValue(fabs(shapes[3]/(shapes[6])));
    if(fabs(shapes[3]) < FLT_MAX)
    {  varMap.at(CLUSHAPE_DiscrT).SetValue(TMath::Sign(float(TMath::Log(shapes[3]+FLT_MIN)), shapes[3])); }
    else { varMap.at(CLUSHAPE_DiscrT).SetValue(_rand->Gaus(-1.,1.e-6)); }
//    varMap.at(CLUSHAPE_DiscrT).SetValue(shapes[3]);
    varMap.at(CLUSHAPE_xl20).SetValue(shapes[15]/(2.0*3.50));
  }
  else
  {  // If shapes empty, push the value out of bounds. Then these variables
     // give the same likelihood for all particle types
    varMap.at(CLUSHAPE_Chi2).SetValue(_rand->Gaus(-1.,1.e-6));
    varMap.at(CLUSHAPE_DiscrL).SetValue(_rand->Gaus(-100.,1.e-6));
    varMap.at(CLUSHAPE_DiscrT).SetValue(_rand->Gaus(-1.,1.e-6));
    varMap.at(CLUSHAPE_xl20).SetValue(_rand->Gaus(-1.,1.e-6));
  }

  // dE/dx
  /*
  varMap.at(DEDX_Chi2electron).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::electron))) );
  varMap.at(DEDX_Chi2muon).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::muon))) );
  varMap.at(DEDX_Chi2pion).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::pion))) );
  varMap.at(DEDX_Chi2kaon).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::kaon))) );
  varMap.at(DEDX_Chi2proton).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::proton))) );
*/
  varMap.at(DEDX_Chi2electron).SetValue( get_dEdxSignedLogChi2(&(particlePars->at(PIDParticles::electron))) );
  varMap.at(DEDX_Chi2muon).SetValue( get_dEdxSignedLogChi2(&(particlePars->at(PIDParticles::muon))) );
  varMap.at(DEDX_Chi2pion).SetValue( get_dEdxSignedLogChi2(&(particlePars->at(PIDParticles::pion))) );
  varMap.at(DEDX_Chi2kaon).SetValue( get_dEdxSignedLogChi2(&(particlePars->at(PIDParticles::kaon))) );
  varMap.at(DEDX_Chi2proton).SetValue( get_dEdxSignedLogChi2(&(particlePars->at(PIDParticles::proton))) );

  /* The following anyway does not correspond to histograms in the standard ILDConfiguration.
  for(VarMap::iterator it=varMap.find(dEdx_first); it!=varMap.find(N_VarTypes); it++)
    it->second.SetValue(-0.5*fabs(it->second.Value()));
    */
  // This might have been the formula used for the standard histos.
  //+TMath::Log(sqrt(2.0*TMath::Pi()*0.05*trk->getdEdx()));

}


void PIDVariables::SetOutOfRange() {
  for (VarMap::iterator it=varMap.begin(); it!=varMap.end(); it++)
    { it->second.SetValue(-FLT_MAX); }
}


float PIDVariables::get_dEdxChi2(PIDParticles::PIDParticle_base* hypothesis) const {

  //get expected dE/dx
  float ExpdEdx=BetheBloch(hypothesis);

  float chi2 = -100.;
  if(dEdx>FLT_MIN) {
    //get chi2!!(so far 5% error assumed. conservative)
    chi2=TMath::Power((dEdx-ExpdEdx)/(0.05*dEdx),2.0);
    if(dEdx-ExpdEdx<0.0) chi2=-chi2;    //get signed chi2
  }

  return chi2;
}

float PIDVariables::get_dEdxSignedLogChi2(PIDParticles::PIDParticle_base* hypothesis) const {

  //get expected dE/dx
  float ExpdEdx=BetheBloch(hypothesis);

  float result = _rand->Gaus(-100., 1e-6);
  if(dEdx>FLT_MIN) {
    float normdev = (dEdx-ExpdEdx)/(0.05*dEdx);
    result = TMath::Sign(float(2.*TMath::Log(fabs(normdev))+FLT_MIN), normdev);
  }

  return result;
}

double PIDVariables::BetheBloch(PIDParticles::PIDParticle_base* hypothesis) const {

  Double_t bg=p/hypothesis->mass;
  Double_t b=sqrt(bg*bg/(1.0+bg*bg));
  //Double_t g=bg/b;
  Double_t tmax=hypothesis->GetBBpars()[2]*TMath::Power(bg,2.0);   ///(1.0+pars[3]*g+pars[4]);

  return (0.5*hypothesis->GetBBpars()[0]*TMath::Log(hypothesis->GetBBpars()[1]*TMath::Power(bg,2.0)*tmax)
          - hypothesis->GetBBpars()[3]*b*b - hypothesis->GetBBpars()[4]*bg/2.0)/(b*b);
}
