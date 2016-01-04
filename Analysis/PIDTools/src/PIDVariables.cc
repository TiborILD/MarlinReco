#include "PIDVariables.hh"
//#include "PIDParticles.hh"
#include "LikelihoodPID.hh"

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
const double PIDVariables::caloCut = 0.05; // Minimum ECAL+HCAL to calculate ECAL/(ECAL+HCAL)
const double PIDVariables::muSysCut = 0.01;
const double PIDVariables::muSysPCut = 5.;

const double PIDVariables::dEdx_MIP = 1.35e-7;

PIDVariables::PIDVariables() : dEdx(0), p(0) {
  // Create map of particle properties.
  // Priors are irrelevant for this class. Set them to 0.
  std::vector<float> dummypriors;
  for(int i=0; i<PIDParticles::nParticleTypes; i++) dummypriors.push_back(0.);
  particlePars = PIDParticles::CreateMap(dummypriors);

  PopulateMap();
}

PIDVariables::PIDVariables(EVENT::ReconstructedParticle* _particle) {
  // Create map of particle properties. This must be done first
  // Priors are irrelevant for this class. Set them to 0.
  std::vector<float> dummypriors;
  for(int i=0; i<PIDParticles::nParticleTypes; i++) dummypriors.push_back(0.);
  particlePars = PIDParticles::CreateMap(dummypriors);

  PopulateMap();
  Update(_particle);
}

void PIDVariables::PopulateMap() {
  varMap.insert(std::pair<varType, PIDVariable>(CALO_Total, PIDVariable("CaloTotal", "(ECAL+HCAL)/p", 0., 5.)));
  varMap.insert(std::pair<varType, PIDVariable>(CALO_EFrac, PIDVariable("CaloEFrac", "ECAL/(ECAL+HCAL)", 0., 1.)));
  varMap.insert(std::pair<varType, PIDVariable>(CALO_MuSys, PIDVariable("CaloMuSys", "#mu system deposit (GeV)", 0., 10.)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_Chi2, PIDVariable("CluShapeChi2", "Cluster shape #chi^{2}", -2., 20.)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_DiscrL, PIDVariable("DiscrepancyL", "Shower max / EM shower max", -30., 80.)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_DiscrT, PIDVariable("DiscrepancyT", "Absorption length / R_{m}", 0., 1.)));
  varMap.insert(std::pair<varType, PIDVariable>(CLUSHAPE_xl20, PIDVariable("Xl20", "xl20", 0., 50.)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2electron, PIDVariable("dEdxChi2electron", "#chi^{2}_{dE/dx} (electron)", -50., 0.)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2muon, PIDVariable("dEdxChi2muon", "#chi^{2}_{dE/dx} (#mu)", -50., 0.)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2pion, PIDVariable("dEdxChi2pion", "#chi^{2}_{dE/dx} (#pi)", -50., 0.)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2kaon, PIDVariable("dEdxChi2kaon", "#chi^{2}_{dE/dx} (K)", -50., 0.)));
  varMap.insert(std::pair<varType, PIDVariable>(DEDX_Chi2proton, PIDVariable("dEdxChi2proton", "#chi^{2}_{dE/dx} (p)", -50., 0.)));
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
  if(cluvec.size()!=0){
    for(unsigned int i=0; i<cluvec.size(); i++){
      FloatVec sde = cluvec[i]->getSubdetectorEnergies();
      ecal += cluvec[i]->getSubdetectorEnergies()[0];
      hcal += cluvec[i]->getSubdetectorEnergies()[1];
      mucal+=cluvec[i]->getSubdetectorEnergies()[2];
    }
    shapes=cluvec[0]->getShape();
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
  // FIXME: There is a strange bug affecting CALO_Total, producing sharp peaks between 0.9 and 1.
  varMap.at(CALO_Total).SetValue( (ecal+hcal) / p );
  if(ecal+hcal > caloCut) { varMap.at(CALO_EFrac).SetValue( ecal/(ecal+hcal) ); }
  else { varMap.at(CALO_EFrac).SetValue( -1. ); }
  varMap.at(CALO_MuSys).SetValue(mucal);

  // Shower shapes
  if(shapes.size()!=0){
    varMap.at(CLUSHAPE_Chi2).SetValue(shapes[0]);
    varMap.at(CLUSHAPE_DiscrL).SetValue(shapes[5]);
    varMap.at(CLUSHAPE_DiscrT).SetValue(fabs(shapes[3])/(shapes[6]));
    varMap.at(CLUSHAPE_xl20).SetValue(shapes[15]/(2.0*3.50));
  }
  else
  {  // If shapes empty, push the value out of bounds. Then these variables
     // give the same likelihood for all particle types
    varMap.at(CLUSHAPE_Chi2).SetValue(-DBL_MAX);
    varMap.at(CLUSHAPE_DiscrL).SetValue(-DBL_MAX);
    varMap.at(CLUSHAPE_DiscrT).SetValue(-DBL_MAX);
    varMap.at(CLUSHAPE_xl20).SetValue(-DBL_MAX);
  }

  // dE/dx
  varMap.at(DEDX_Chi2electron).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::electron))) );
  varMap.at(DEDX_Chi2muon).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::muon))) );
  varMap.at(DEDX_Chi2pion).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::pion))) );
  varMap.at(DEDX_Chi2kaon).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::kaon))) );
  varMap.at(DEDX_Chi2proton).SetValue( get_dEdxChi2(&(particlePars->at(PIDParticles::proton))) );

  /* The following anyway does not correspond to histograms in the standard ILDConfiguration.
  for(VarMap::iterator it=varMap.find(dEdx_first); it!=varMap.find(N_VarTypes); it++)
    it->second.SetValue(-0.5*fabs(it->second.Value()));
    */
  // This might have been the formula used for the standard histos.
  //+TMath::Log(sqrt(2.0*TMath::Pi()*0.05*trk->getdEdx()));

}



double PIDVariables::get_dEdxChi2(PIDParticle* hypothesis) const {

  //get expected dE/dx
  double ExpdEdx=BetheBloch(hypothesis);

  //get chi2!!(so far 5% error assumed. conservative)
  Double_t chi2=TMath::Power((dEdx-ExpdEdx)/(0.05*dEdx),2.0);
  if(dEdx-ExpdEdx<0.0) chi2=-chi2;    //get signed chi2
  return chi2;
}

double PIDVariables::BetheBloch(PIDParticle* hypothesis) const {

  Double_t bg=p/hypothesis->mass;
  Double_t b=sqrt(bg*bg/(1.0+bg*bg));
  //Double_t g=bg/b;
  Double_t tmax=hypothesis->GetBBpars()[2]*TMath::Power(bg,2.0);   ///(1.0+pars[3]*g+pars[4]);

  return (0.5*hypothesis->GetBBpars()[0]*TMath::Log(hypothesis->GetBBpars()[1]*TMath::Power(bg,2.0)*tmax)
          - hypothesis->GetBBpars()[3]*b*b - hypothesis->GetBBpars()[4]*bg/2.0)/(b*b);
}
