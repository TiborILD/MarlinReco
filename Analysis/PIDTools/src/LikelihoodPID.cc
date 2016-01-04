/*
--------------- inplement Particle Identofication using Global likelihood -------------

This code includes class for Particle ID. The method is Naive Bayes Classfier
1. First, check whether a particle is electron or not. 
   If so, the particle is identified as electron and return electron status code.
2. And then, check whether a particle a particle is muonn or not. 
   If so, the particle is identified as muon and return muon status code.
3. If a particle is neither electron nor muon, 
   try to classify a particle to each Hadron type(Pion, Kaon, or Proton)

TODO:
Bethe-Bloch parameters should be moved to steer file
Threshold parameters should be moved to steer file

make flags to choose Particle ID method 

risk-minimization and MAP
 */

//#include <string>
//#include <sstream>
//#include <TFile.h>
//#include <TH1F.h>
#include "LikelihoodPID.hh"
#include "TString.h"
#include "TMath.h"
//#include "PIDVariables.hh"

using namespace std;

std::string itos(int i)  
{
  std::stringstream s;
  s << i;
  return s.str();
}





/*******************************************************
 *
 *   Implementation of the LikelihoodPID class
 *
 ******************************************************/

const short LikelihoodPID::MASK_Basic  = 1 ;
const short LikelihoodPID::MASK_dEdx   = 1 << 1;
const short LikelihoodPID::MASK_Shapes = 1 << 2;

const double LikelihoodPID::lowestLogL = -DBL_MAX/100;

LikelihoodPID::LikelihoodPID(string fname, std::vector<float> priors) :
_particle(NULL),
_algoFlags(MASK_Basic | MASK_dEdx | MASK_Shapes)
{
  // Prepare particle map and set priors
  double tot=0.;
  for (unsigned int i=0; i<priors.size(); i++) tot += priors.at(i);
  for (unsigned int i=0; i<priors.size(); i++) priors.at(i) /= tot;

  particlePars = PIDParticles::CreateMap(priors);
  // Initialize bestParticle to an invalid particle
  bestParticle = particlePars->end();

  fpdf=new TFile(fname.c_str());

  string hname,hname2;
  //for(Int_t i=0;i<6;i++) pdf[i] =new TH1F()[14];
  for(Int_t i=0;i<PIDParticles::nParticleTypes;i++){
    //ep
    hname="hep" + itos(i+1);
    hname2="hep" + itos(i+1) + "_2";
    pdf[i][0]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //ehad
    hname="hehad" + itos(i+1);
    hname2="hehad" + itos(i+1) + "_2";
    pdf[i][1]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //mucal
    hname="hmucal" + itos(i+1);
    hname2="hmucal" + itos(i+1) + "_2";
    pdf[i][2]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //chi2
    hname="hchi2" + itos(i+1);
    hname2="hchi2" + itos(i+1) + "_2";
    pdf[i][3]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //showermax/exp.showermax
    hname="hldiscrepancy" + itos(i+1);
    hname2="hldiscrepancy" + itos(i+1) + "_2";
    pdf[i][4]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //absorption length
    hname="htdiscrepancy" + itos(i+1);
    hname2="htdiscrepancy" + itos(i+1) + "_2";
    pdf[i][5]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //xl20
    hname="hxl20" + itos(i+1);
    hname2="hxl20" + itos(i+1) + "_2";
    pdf[i][6]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //likeliele
    hname="hlikeliele" + itos(i+1);
    hname2="hlikeliele" + itos(i+1) + "_2";
    pdf[i][7]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //likelimuo
    hname="hlikelimuo" + itos(i+1);
    hname2="hlikelimuo" + itos(i+1) + "_2";
    pdf[i][8]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //likelipi
    hname="hlikelipi" + itos(i+1);
    hname2="hlikelipi" + itos(i+1) + "_2";
    pdf[i][9]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //likelik
    hname="hlikelik" + itos(i+1);
    hname2="hlikelik" + itos(i+1) + "_2";
    pdf[i][10]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //likelip
    hname="hlikelip" + itos(i+1);
    hname2="hlikelip" + itos(i+1) + "_2";
    pdf[i][11]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //deltax
    hname="hdeltax" + itos(i+1);
    hname2="hdeltax" + itos(i+1) + "_2";
    pdf[i][12]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
    //deltaz
    hname="hdeltaz" + itos(i+1);
    hname2="hdeltaz" + itos(i+1) + "_2";
    pdf[i][13]=(TH1F*)fpdf->Get(hname.c_str())/*->Clone(hname2.c_str())*/;
  }
  
  //normalize histograms
  Double_t weight=1.0;
  for(Int_t i=0;i<PIDParticles::nParticleTypes;i++){
    for(Int_t j=0;j<14;j++){
      //normalize histograms
      weight=pdf[i][j]->Integral(0,pdf[i][j]->GetNbinsX()+1,"");
      pdf[i][j]->Scale(1.0/weight);
    }
  }

  //set threshold
  //this is original
  // FIXME: This is not used
  particlePars->at(PIDParticles::electron).SetThreshold(TMath::Exp(-0.55));
  particlePars->at(PIDParticles::muon).SetThreshold(TMath::Exp(-0.7));
  particlePars->at(PIDParticles::pion).SetThreshold(0.0);
  particlePars->at(PIDParticles::kaon).SetThreshold(TMath::Exp(-0.6));
  particlePars->at(PIDParticles::proton).SetThreshold(TMath::Exp(-0.6));
  
  return;
}

LikelihoodPID::~LikelihoodPID(){

  fpdf->Close();
  delete particlePars;
}

//public
Int_t LikelihoodPID::Classification(IMPL::ReconstructedParticleImpl* particle){
  _particle = particle;

  return Classification();
}

Int_t LikelihoodPID::Classification(){

  CalcPosteriors();

  double maxpost = -DBL_MAX;
  bestParticle = particlePars->end(); // Initialize to unset particle
  for(particle_c_iterator it=particlePars->begin(); it!=particlePars->end(); it++) {
//    std::cout << "Posterior for " << it->second.Name() << ": "
  //      << it->second.Posterior() << std::endl;
    if(it->second.Posterior() > maxpost) {
      maxpost = it->second.Posterior();
      bestParticle = it;
    }
  }

  return bestParticle->second.pdg;
}

int LikelihoodPID::GetBestType() const {
  if(!_particle) return 0;
  if(bestParticle == particlePars->end()) return 0; // Unset particle

  return bestParticle->first;
}

int LikelihoodPID::GetBestPDG() const {
  if(!_particle) return 0;
  if(bestParticle == particlePars->end()) return 0; // Unset particle

  return bestParticle->second.pdg;
}

double LikelihoodPID::GetBestLikelihood() const {
  if(bestParticle == particlePars->end()) return -DBL_MAX;
  return TMath::Exp(bestParticle->second.LogL());
}

double LikelihoodPID::GetBestProbability() const {
  if(bestParticle == particlePars->end()) return -1.;
  return bestParticle->second.Posterior();
}


Double_t LikelihoodPID::getCorrEnergy(TLorentzVector pp, parType parttype){
  Double_t tmpmass=particlePars->at(parttype).mass;
  return sqrt(pp.P()*pp.P()+tmpmass*tmpmass);
}


void LikelihoodPID::CalcPosteriors(){
  // Calculation of likelihoods and Bayesian probabilities

  //  Prepare sensitive variables  ***/
  variables.Update(_particle);

  // Initialize likelihoods
  for(particle_iterator it=particlePars->begin(); it!=particlePars->end(); it++) {
    it->second.ResetLogL();
  }


  /*********************************************/
  /***   Loop over the sensitive variables   ***/
  /*********************************************/

  // Basic variables
  if( (_algoFlags&MASK_Basic) ) {

    for (variable_c_iterator varIt = variables.GetMap()->find(PIDVariables::basic_first);
            varIt!=variables.GetMap()->find(PIDVariables::calo_beyond); varIt++)
    {
      if (variables.GetVariable(PIDVariables::CALO_Total) < PIDVariables::caloCut) continue;
      // Addup logL for all particles
      for(particle_iterator partIt=particlePars->begin(); partIt!=particlePars->end(); partIt++) {
        // If deposit in the muon-system below the threshold, muon calo deposit distributions are different
        // FIXME: Should this not be decided by the tracker pT rather than the mu-sys deposit?
        if(variables.GetVariable(PIDVariables::CALO_MuSys) < PIDVariables::muSysCut && partIt->first==PIDParticles::muon) {
          partIt->second.AddLogL(LogL(PIDParticles::lowEmuon, varIt->first, varIt->second.Value() ));
        }
        else { partIt->second.AddLogL(LogL(partIt->first, varIt->first, varIt->second.Value())); }
      }

    }

    // If momentum sufficient, use Muon system deposit (NEW)
    if(variables.GetP() > PIDVariables::muSysPCut) {

      for(particle_iterator partIt=particlePars->begin(); partIt!=particlePars->end(); partIt++) {
        partIt->second.AddLogL(LogL(partIt->first, PIDVariables::CALO_MuSys,
                                 variables.GetVariable(PIDVariables::CALO_MuSys))); }

    }
  }

  // Cluster shape vars
  if( (_algoFlags&MASK_Shapes) && variables.GetVariable(PIDVariables::CALO_Total)>PIDVariables::caloCut) {

    for (variable_c_iterator varIt = variables.GetMap()->find(PIDVariables::clushape_first);
          varIt!=variables.GetMap()->find(PIDVariables::clushape_beyond); varIt++)
    {
      // Addup logL for all particles
      for(particle_iterator partIt=particlePars->begin(); partIt!=particlePars->end(); partIt++) {
        // If deposit in the muon-system below the threshold, muon calo deposit distributions are different
        // FIXME: Should this not be decided by the tracker pT rather than the mu-sys deposit?
        if(variables.GetVariable(PIDVariables::CALO_MuSys) < PIDVariables::muSysCut
            && partIt->second.pdg==13) {
          partIt->second.AddLogL(LogL(PIDParticles::lowEmuon, varIt->first, varIt->second.Value()));
        }
        else { partIt->second.AddLogL(LogL(partIt->first, varIt->first, varIt->second.Value())); }
      }
    }
  }

  // dEdx vars
  if(_algoFlags&MASK_dEdx) {

    for (variable_c_iterator varIt = variables.GetMap()->find(PIDVariables::dEdx_first);
        varIt!=variables.GetMap()->find(PIDVariables::dEdx_beyond); varIt++)
    {
      // Addup logL for all particles
      for(particle_iterator partIt=particlePars->begin(); partIt!=particlePars->end(); partIt++) {
        // If no trace in the muon system, muons get their L from the tracker, etc
        if(variables.GetVariable(PIDVariables::CALO_MuSys) < PIDVariables::muSysCut
            && partIt->second.pdg==13) {
          partIt->second.AddLogL(LogL(PIDParticles::lowEmuon, varIt->first, varIt->second.Value()));
        }
        else { partIt->second.AddLogL(LogL(partIt->first, varIt->first, varIt->second.Value())); }
      }
    }
  }

  // Normalisation denominator for the posteriors
  double posteriorWt = 0.;
  for(particle_c_iterator partIt=particlePars->begin(); partIt!=particlePars->end(); partIt++) {
    posteriorWt += partIt->second.prior * TMath::Exp(partIt->second.LogL());
  }
    // Calculate posteriors
  for(particle_iterator partIt=particlePars->begin(); partIt!=particlePars->end(); partIt++) {
    partIt->second.SetPosterior(TMath::Exp(partIt->second.LogL() + TMath::Log(partIt->second.prior)
                            - TMath::Log(posteriorWt)));
  }

  return;
}

const Double_t LikelihoodPID::LogL(parType type, varType valtype, Double_t value){

  TH1F *ph = pdf[type][valtype];
  TAxis *pa = ph->GetXaxis();

  Int_t bin = pa->FindFixBin(value);
  // Value out of range condition occurs at the same time for all hypotheses
  // (At least it should be so by design of PDF histograms
  // Thus in this case we consider this test neutral and return 0.
  if (bin==0 || bin>pa->GetNbins()) return 0;

  //get likelihood
  Double_t val=ph->GetBinContent(bin);

  if(val == 0) return lowestLogL;
  return TMath::Log(val);
}

Double_t LikelihoodPID::getPenalty(Int_t ptype, Int_t hypothesis, Double_t p){
  Double_t par[3]={0.0,0.0,0.0};
  //set parameters
  switch(ptype){
  case 0:  //electron
    if(hypothesis==0){
      par[0]=0.0;
      par[1]=1.0;
      par[2]=1.00178;

    }else if(hypothesis==1){
      par[0]=0.101945;
      par[1]=2.18719;
      par[2]=0.990000;

    }else if(hypothesis==2){
      par[0]=0.104521;
      par[1]=1.98490;
      par[2]=0.990000;

    }else if(hypothesis==3){
      par[0]=0.234299;
      par[1]=0.276835;
      par[2]=0.965973;

    }else{
      par[0]=0.414613;
      par[1]=-0.132530;
      par[2]=0.949679;

    }

    break;
  case 1:   //muon
    if(hypothesis==0){
      par[0]=-0.00512190;
      par[1]=-0.211835;
      par[2]=1.00024;

    }else if(hypothesis==1){
      par[0]=0.0;
      par[1]=1.00;
      par[2]=0.999684;

    }else if(hypothesis==2){
      par[0]=0.00833283;
      par[1]=9.99995;
      par[2]=0.999027;

    }else if(hypothesis==3){
      par[0]=0.0964021;
      par[1]=-0.214469;
      par[2]=0.989688;

    }else{
      par[0]=0.318674;
      par[1]=-0.197755;
      par[2]=0.968436;
    }

    break;
  case 2:   //pion
    if(hypothesis==0){
      par[0]=-0.0123577;
      par[1]=-0.141521;
      par[2]=1.00273;

    }else if(hypothesis==1){
      par[0]=-0.00558462;
      par[1]=-0.136941;
      par[2]=1.00135;

    }else if(hypothesis==2){
      par[0]=0.0;
      par[1]=1.0;
      par[2]=1.00001;

    }else if(hypothesis==3){
      par[0]=0.122083;
      par[1]=-0.1333923;
      par[2]=0.976863;

    }else{
      par[0]=0.401111;
      par[1]=-0.116807;
      par[2]=0.930906;
    }

    break;
  case 3:    //kaon
    if(hypothesis==0){
      par[0]=-0.102300;
      par[1]=-0.139570;
      par[2]=1.01362;

    }else if(hypothesis==1){
      par[0]=-0.0973257;
      par[1]=-0.138932;
      par[2]=1.01293;

    }else if(hypothesis==2){
      par[0]=-0.0936450;
      par[1]=-0.138469;
      par[2]=1.01242;

    }else if(hypothesis==3){
      par[0]=0.0;
      par[1]=1.0;
      par[2]=0.999865;

    }else{
      par[0]=0.223317;
      par[1]=-0.101273;
      par[2]=0.973484;
    }

    break;
  case 4:   //proton
    if(hypothesis==0){
      par[0]=-0.260150;
      par[1]=-0.0990612;
      par[2]=1.02854;

    }else if(hypothesis==1){
      par[0]=-0.256503;
      par[1]=-0.0976728;
      par[2]=1.02811;

    }else if(hypothesis==2){
      par[0]=-0.253788;
      par[1]=-0.0966732;
      par[2]=1.02779;

    }else if(hypothesis==3){
      par[0]=-0.183031;
      par[1]=-0.0742123;
      par[2]=1.01965;

    }else{
      par[0]=0.0;
      par[1]=1.0;
      par[2]=0.999791;
    }

    break;
  }

  return par[0]/sqrt(p*p+par[1])+par[2]; 
}

