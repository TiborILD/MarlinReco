/*
 * PIDParticles.cc
 *
 *  Created on: Dec 23, 2015
 *      Author: strahinja
 */

#include "PIDParticles.hh"

PIDParticles::ParameterMap* PIDParticles::CreateMap(std::vector<float> priors) {

  ParameterMap *parameterMap = new  ParameterMap;

  double BBparsElectron[5] = {-2.40638e-03, 7.10337e-01, 2.87718e-01, -1.71591e+00, 0.0};
  parameterMap->insert(std::pair<particleType, PIDParticle>
      (electron, PIDParticle("electron", 11, .000510998, priors.at(electron), BBparsElectron)));

  double BBparsMuon[5] = {8.11408e-02, 9.92207e-01, 7.58509e+05, -1.70167e-01, 4.63670e-04};
  parameterMap->insert(std::pair<particleType, PIDParticle>
      (muon, PIDParticle("muon", 13, .105658, priors.at(muon), BBparsMuon)));

  double BBparsPion[5] = {8.10756e-02, -1.45051e+06, -3.09843e+04, 2.84056e-01, 3.38131e-04};
  parameterMap->insert(std::pair<particleType, PIDParticle>
      (pion, PIDParticle("pion", 211, .139570, priors.at(pion), BBparsPion)));

  double BBparsKaon[5] = {7.96117e-02, 4.13335e+03, 1.13577e+06, 1.80555e-01, -3.15083e-04};
  parameterMap->insert(std::pair<particleType, PIDParticle>
      (kaon, PIDParticle("kaon", 321, .493677, priors.at(kaon), BBparsKaon)));

  double BBparsProton[5] = {7.78772e-02, 4.49300e+04, 9.13778e+04, 1.50088e-01, -6.64184e-04};
  parameterMap->insert(std::pair<particleType, PIDParticle>
      (proton, PIDParticle("proton", 2212, .938272, priors.at(proton), BBparsProton)));

  return parameterMap;
}

