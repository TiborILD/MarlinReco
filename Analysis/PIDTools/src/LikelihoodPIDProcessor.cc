#include <vector>
#include <string>

#include <EVENT/LCCollection.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

#include "TLorentzVector.h"

#include "LikelihoodPIDProcessor.hh"
#include "LikelihoodPID.hh"
#include "LowMomentumMuPiSeparationPID_BDTG.hh"

LikelihoodPIDProcessor aLikelihoodPIDProcessor ;

LikelihoodPIDProcessor::LikelihoodPIDProcessor()
  : Processor("LikelihoodPIDProcessor"),
    _myPID(NULL), _pfoCol(NULL), _mupiPID(NULL),
    _basicFlg(true), _dEdxFlg(true), _showerShapesFlg(true),
    _nEvt(0)
    {
  
  // Processor description
  _description = "Particle ID code using Bayesian Classifier" ;
 
  std::vector< std::string > xmlfiles;
  xmlfiles.push_back( "TMVAClassification_BDTG_02GeVP_clusterinfo.weights.xml" );  
  xmlfiles.push_back( "TMVAClassification_BDTG_03GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_04GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_05GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_06GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_07GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_08GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_09GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_10GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_11GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_12GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_13GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_14GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_15GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_16GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_17GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_18GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_19GeVP_clusterinfo.weights.xml" );
  xmlfiles.push_back( "TMVAClassification_BDTG_20GeVP_clusterinfo.weights.xml" );
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
           "RecoParticleCollection" ,
           "Input collection of Reconstrcuted Particle",
           _inputPFOsCollection,
           std::string("PandoraPFOs"));
  
  registerProcessorParameter( "EnergyBoundaries" ,
			      "Boundaries for energy",
			      _energyBoundary,
			      EVENT::FloatVec(0,1.0e+07));
  
  registerProcessorParameter( "FilePDFName" ,
			      "rootfile of PDF",
			      _PDFName,
			      std::string("pdf_ParticleID_ok.root") );
  
  registerProcessorParameter( "FileWeightFormupiSeparationName" ,
            "weight file for low momentum mu pi separation",
            _weightFileName,
            xmlfiles );

  registerProcessorParameter( "UseBasicPID" ,
            "Selection flag to use basic PID",
            _basicFlg,
            bool(true) );

  registerProcessorParameter( "Use_dEdx" ,
            "Selection flag to use dEdx PID",
            _dEdxFlg,
            bool(true) );

  registerProcessorParameter( "UseShowerShapes" ,
            "Selection flag to use shower shapes",
            _showerShapesFlg,
            bool(true) );

  registerProcessorParameter( "AlgoName" ,
            "Algorithm name to be written to the slcio file",
            _algoName,
            std::string("LikelihoodPID") );

  std::vector<float> defaultPriors;
  defaultPriors.push_back(0.2);
  defaultPriors.push_back(0.2);
  defaultPriors.push_back(0.2);
  defaultPriors.push_back(0.2);
  defaultPriors.push_back(0.2);
  registerProcessorParameter( "ParticlePriors" ,
            "Prior probabilities for particle types. E.g. global frequencies in the experiment",
            _particlePriors,
            defaultPriors );
}  
  


void LikelihoodPIDProcessor::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  // Processor output parameters
  _parNames.push_back("electronLikelihood");
  _parNames.push_back("muonLikelihood");
  _parNames.push_back("pionLikelihood");
  _parNames.push_back("kaonLikelihood");
  _parNames.push_back("protonLikelihood");
  _parNames.push_back("MVAOutput_mupiSeparation");
  _parNames.push_back("electronProbability");
  _parNames.push_back("muonProbability");
  _parNames.push_back("pionProbability");
  _parNames.push_back("kaonProbability");
  _parNames.push_back("protonProbability");
 
  _myPID = new LikelihoodPID(_PDFName, _particlePriors);
  _myPID->setBasicFlg(_basicFlg);
  _myPID->setdEdxFlg(_dEdxFlg);
  _myPID->setShowerShapesFlg(_showerShapesFlg);

  streamlog_out(MESSAGE) << "Priors normalized: ";
  for (particle_c_iterator it = _myPID->GetParticlePars()->begin();
                                               it != _myPID->GetParticlePars()->end(); it++)
    streamlog_out(MESSAGE) << it->second.prior << " ";
  streamlog_out(MESSAGE) << std::endl;

  //mupi separation class
  _mupiPID = new LowMomentumMuPiSeparationPID_BDTG(_weightFileName);

  _nEvt = 0;

  printParameters();
  
}

void LikelihoodPIDProcessor::processRunHeader( LCRunHeader* run) { 
} 

void LikelihoodPIDProcessor::processEvent( LCEvent * evt ) { 
  try { _pfoCol = evt->getCollection( _inputPFOsCollection ) ; }
  catch ( lcio::DataNotAvailableException &e )
    {
        streamlog_out(WARNING) <<  _inputPFOsCollection
            << " collection not available. Skipping event." << std::endl;
        _pfoCol = NULL;
        return;
    }
  
  _nEvt++;
  if(_nEvt%1000 == 0)
    streamlog_out(MESSAGE) << "Processing event " << _nEvt
           << " (# on file " << evt->getEventNumber() << ")" << std::endl;

  PIDHandler pidh(_pfoCol);   //BasicPID
//  int algoID = pidh.addAlgorithm("BayesianPID", _parNames);
//  int algoID = pidh.addAlgorithm("NewBayesianPID", _parNames);
  int algoID = pidh.getAlgorithmID(_algoName.c_str());
//  streamlog_out(DEBUG) << "Added algorithm " << algoid << std::endl;
 
  int npfo = _pfoCol->getNumberOfElements();
  for (int i = 0; i < npfo; i++ ) {
    ReconstructedParticleImpl* part =
        dynamic_cast<ReconstructedParticleImpl*>( _pfoCol->getElementAt(i) );
    
    if(part->getCharge()==0) continue;  //avoid neutral particle ...
    
    
//////////////////////////////////////////////////////////////////////////////////
    streamlog_out(DEBUG) << "Calling Classification." << std::endl;
    int pdg = _myPID->Classification(part);
    streamlog_out(DEBUG) << "Done Classification." << std::endl;

    //mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
    // It would be good to also make MuPISeparation that takes (ReconstructedParticleImpl*)
    // as argument.
    // TODO: Move to using the new helper classes
    Float_t MVAoutput = -1.0;
    TLorentzVector pp(TVector3(part->getMomentum()), part->getEnergy());
    if((abs(pdg)==13 || abs(pdg)==211) && pp.P()<2.0){
      streamlog_out(DEBUG) << "Checking mu/pi." << std::endl;
      EVENT::ClusterVec clu=part->getClusters();
      lcio::Track* trk = part->getTracks()[0];
      int parttype=_mupiPID->MuPiSeparation(pp, trk, clu);
      _myPID->setBestParticle((parttype==1) ? PIDParticles::muon : PIDParticles::pion);
      MVAoutput = _mupiPID->getMVAOutput();
    }

    //create PIDHandler
    streamlog_out(DEBUG) << "Creating PID handler." << std::endl;
    createParticleIDClass(part, pidh, algoID, MVAoutput);

  }
  
}

void LikelihoodPIDProcessor::check( LCEvent * evt ) { 
}

void LikelihoodPIDProcessor::end() { 
  delete _myPID;
  delete _mupiPID;
}

void LikelihoodPIDProcessor::createParticleIDClass(ReconstructedParticle *part,
                                 PIDHandler &pidh, int algoID, float MVAoutput) {

  const ParticleMap *particlePars = _myPID->GetParticlePars();
  std::vector<float> pidPars;
  for(particle_c_iterator it=particlePars->begin(); it!=particlePars->end(); it++) {
    pidPars.push_back(it->second.LogL());
  }
  pidPars.push_back(MVAoutput);
  for(particle_c_iterator it=particlePars->begin(); it!=particlePars->end(); it++) {
    pidPars.push_back(it->second.Posterior());
  }
  for(particle_c_iterator it=particlePars->begin(); it!=particlePars->end(); it++) {
    pidPars.push_back(it->second.DEdxDist());
  }

  
  //set pid results
  // Perhaps user id should be used!
  pidh.setParticleID(part, int(_myPID->GetBestType()), _myPID->GetBestPDG(), _myPID->GetBestLikelihood(), algoID, pidPars);
  
  return;
}
