/*
 * MvaPidProcessor.cc
 *
 *
 *
 *  Created on: Feb 4, 2016
 *      Author: Strahinja Lukic
 */


#include "MvaPidProcessor.hh"
#include "LowMomentumMuPiSeparationPID_BDTG.hh"



const char *MvaPidProcessor::algoName = "MvaPid";

MvaPidProcessor::MvaPidProcessor() :
  Processor("MvaPidProcessor"),
  _variables(NULL), _hypotheses(NULL),
  _description("Particle ID using MVA"),
  _pfoCol(NULL), _pidh(NULL), _mupiPID(NULL), _nEvt(0)
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


/**********************************************************
***********************************************************/
void MvaPidProcessor::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  _hypotheses = PIDParticles::CreateMVAPIDMap();
  _variables = new PIDVariables;

  // Prepare vectors of processor output parameters for the PIDHandler:
  // MVA output and Q-statistic for each hypothesis
  // Prepare one TMVA::Reader for each hypothesis
  for (hypotheses_c_iterator it=_hypotheses->begin(); it!=_hypotheses->end(); it++) {
    std::string mvaname(it->second.Name()); mvaname.append("MVAout");
    _pidPars.push_back(0.);
    _pidParNames.push_back(mvaname);
    std::string qname(it->second.Name()); qname.append("Q");
    _pidPars.push_back(0.);
    _pidParNames.push_back(qname);

    _readerMap.insert(std::pair<particleType, TMVA::Reader*>
                                (it->first, new TMVA::Reader("Silent")));
  }


  // Prepare maps of mva variables
  // (necessity because TMVA::Reader::AddVariable(...) does not take const pointers)
  // Add all mva variables to the MVA readers
  for (variable_c_iterator it=_variables->GetMap()->begin();
                           it!=_variables->GetMap()->end();
                           it++)
  { // Hope this works...
    _mvaVars.insert(std::pair<const char*, float>(it->second.Name(), 0.));
    for(ReaderMap::const_iterator itr = _readerMap.begin(); itr != _readerMap.end(); itr++)
      itr->second->AddVariable(it->second.Name(), &(_mvaVars.at(it->second.Name())));
  }
  _mvaVars.insert(std::pair<const char*, float>("p", 0.));
  for(ReaderMap::const_iterator itr = _readerMap.begin(); itr != _readerMap.end(); itr++)
    itr->second->AddVariable("p", &(_mvaVars.at("p")));

  //mupi separation class
  _mupiPID = new LowMomentumMuPiSeparationPID_BDTG(_muPiWeightFileNames);

  _nEvt = 0;

  printParameters();

}


/**********************************************************
***********************************************************/
void MvaPidProcessor::processRunHeader( LCRunHeader* run) {
}


/**********************************************************
***********************************************************/
void MvaPidProcessor::processEvent( LCEvent * evt ) {

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

  _pidh = new PIDHandler(_pfoCol);
  _pidh->addAlgorithm(algoName, _pidParNames);

  // Loop over PFO objects
  int npfo = _pfoCol->getNumberOfElements();
  for (int i = 0; i < npfo; i++ ) {
    ReconstructedParticleImpl* part =
        dynamic_cast<ReconstructedParticleImpl*>( _pfoCol->getElementAt(i) );

    // if( part->getCharge()==0) continue;  //avoid neutral particle ...
    if( part->getTracks().size() == 0 ) continue;  // skip particles without tracks ...


    //////////////////////////////////////////////////////////////////////////////////
    streamlog_out(DEBUG) << "Starting identification." << std::endl;
    _bestHypothesis=_hypotheses->end();
    Identify(part);
    streamlog_out(DEBUG) << "Done identification." << std::endl;

    //mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
    // It would be good to also make MuPISeparation that takes (ReconstructedParticleImpl*)
    // as argument.
    // TODO: Move this part to Identify()
    Float_t MVAoutput = -1.0;
    if((abs(_bestHypothesis->second.pdg)==13 || abs(_bestHypothesis->second.pdg)==211)
        && _variables->GetP()<2.0)
    {
      streamlog_out(DEBUG) << "Checking mu/pi." << std::endl;
      TLorentzVector pp(TVector3(part->getMomentum()), part->getEnergy());
      EVENT::ClusterVec clu=part->getClusters();
      lcio::Track* trk = part->getTracks()[0];
      int parttype=_mupiPID->MuPiSeparation(pp, trk, clu);
      // this will work in Identify: pdg = (parttype==1) ? 13 : 211;
      MVAoutput = _mupiPID->getMVAOutput();
    }

    _pidh->setParticleID(part, _bestHypothesis->first, _bestHypothesis->second.pdg,
                   GetQ(_bestHypothesis), _pidh->getAlgorithmID(algoName), _pidPars);

  } // Loop over PFO objects


  delete _pidh;
  _pidh = NULL;

}


/**********************************************************
***********************************************************/
void MvaPidProcessor::check( LCEvent * evt ) {
}

/**********************************************************
***********************************************************/
void MvaPidProcessor::end() {
  delete _mupiPID; _mupiPID = NULL;
  for(ReaderMap::const_iterator itr = _readerMap.begin(); itr != _readerMap.end(); itr++)
  { delete itr->second; }
  _readerMap.clear();
  delete _variables; _variables = NULL;
  delete _hypotheses; _hypotheses = NULL;
}


/**********************************************************
 *
 * Private utility functions
 *
***********************************************************/

// Updates _variables, sets _mvaVars, evaluates MVA, selects best hypothesis,
// Checks lowE mu-pi separation, fills _pidPars
void MvaPidProcessor::Identify(ReconstructedParticle* particle) {
  _variables->Update(particle);

  for(hypotheses_iterator hit=_hypotheses->begin(); hit != _hypotheses->end(); hit++) {
    for(variable_c_iterator it=_variables->GetMap()->begin();
                            it != _variables->GetMap()->end();
                            it++)
    {  _mvaVars.at(it->second.Name()) = it->second.Value();  }

    hit->second.SetMVAout(_readerMap.at(hit->first)->EvaluateMVA(_mvaMethod));
  }
}


