/*
 * MvaPidProcessor.cc
 *
 *
 *
 *  Created on: Feb 4, 2016
 *      Author: Strahinja Lukic
 */


#include "MvaPidProcessor.hh"

MvaPidProcessor::MvaPidProcessor() :
  Processor("MvaPidProcessor"),
  _reader(NULL), _variables(NULL), _hypotheses(NULL),
  _description("Particle ID using MVA"),
  _pfoCol(NULL), _mupiPID(NULL), _nEvt(0)
{
  // Defaults for lowE mu-pi separation
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

  registerProcessorParameter( "WeightFileName" ,
            "Weight file for general PID",
            _weightFileName,
            std::string("MvaPidBDT.weights.xml") );

  registerProcessorParameter( "MVAMethod" ,
            "MVA method name",
            _mvaMethod,
            std::string("BDT") );

  registerProcessorParameter( "FileWeightFormupiSeparationName" ,
            "Weight files for low momentum mu pi separation",
            _muPiWeightFileNames,
            xmlfiles );

  registerProcessorParameter( "EnergyBoundaries" ,
            "Boundaries for energy",
            _energyBoundary,
            EVENT::FloatVec(0,1.0e+07));

}


void MvaPidProcessor::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  _hypotheses = PIDParticles::CreateMVAPIDMap();

  // Processor output parameters: MVA output and Q-statistic for each hypothesis
  for (particle_c_iterator it=_hypotheses->begin(); it!=_hypotheses->end(); it++) {
    std::string mvaname(it->second.Name()); mvaname.append("MVAout");
    _pidPars.insert(std::pair<std::string, float>(mvaname, 0.));
    std::string qname(it->second.Name()); qname.append("Q");
    _pidPars.insert(std::pair<std::string, float>(qname, 0.));
  }

  _variables = new PIDVariables;
  _reader = new TMVA::Reader("Silent");

  //mupi separation class
  _mupiPID = new LowMomentumMuPiSeparationPID_BDTG(_muPiWeightFileNames);

  _nEvt = 0;

  printParameters();

}


void MvaPidProcessor::check( LCEvent * evt ) {
}

void MvaPidProcessor::end() {
  delete _mupiPID;
  delete _reader;
  delete _variables;
}
