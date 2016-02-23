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

#include "TFile.h"


const char *MvaPidProcessor::algoName = "MvaPid";

MvaPidProcessor aMvaPidProcessor ;

MvaPidProcessor::MvaPidProcessor() :
  Processor("MvaPidProcessor"),
  _variables(NULL), _hypotheses(NULL),
  _description("Particle ID using MVA"),
  _pfoCol(NULL), _pidh(NULL), _mupiPID(NULL),
  _nEvt(0), _nPFO(0), _nUnidentified(0), _nDecisionQ(0),
_nEmptyClusters(0), _nEmptyTracks(0), _nEmptyShapes(0),
_nZerodEdx(0)
{
  // Defaults for lowE mu-pi separation
  std::vector< std::string > mupi_weightfiles;
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_02GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_03GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_04GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_05GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_06GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_07GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_08GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_09GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_10GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_11GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_12GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_13GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_14GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_15GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_16GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_17GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_18GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_19GeVP_clusterinfo.weights.xml" );
  mupi_weightfiles.push_back( "TMVAClassification_BDTG_20GeVP_clusterinfo.weights.xml" );

  std::vector< std::string > weightfiles;
  PIDParticles::ParticleMap* defaultmap = PIDParticles::CreateParticleMap();
  for(PIDParticles::ParticleMap::iterator pit=defaultmap->begin();
                                          pit!=defaultmap->end(); pit++)
  { weightfiles.push_back(std::string("MvaPid_") + pit->second.Name()); }


  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
           "RecoParticleCollection" ,
           "Input collection of Reconstructed Particle",
           _inputPFOsCollection,
           std::string("PandoraPFOs"));

  registerProcessorParameter( "WeightFileName" ,
            "Weight file for general PID",
            _weightFileNames,
            weightfiles );

  registerProcessorParameter( "MVAMethod" ,
            "MVA method name",
            _mvaMethod,
            std::string("BDT") );

  registerProcessorParameter( "FileWeightFormupiSeparationName" ,
            "Weight files for low momentum mu pi separation",
            _muPiWeightFileNames,
            mupi_weightfiles );

}


/**********************************************************
***********************************************************/
void MvaPidProcessor::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  _hypotheses = PIDParticles::CreateMVAPIDMap();
  _variables = new PIDVariables;

  // Prepare vectors of processor output parameters for the PIDHandler:
  // MVA output and Q-statistic for each hypothesis
  for (hypotheses_c_iterator it=_hypotheses->begin(); it!=_hypotheses->end(); it++) {
    std::string mvaname(it->second.Name()); mvaname.append("MVAout");
    _pidPars.push_back(0.);
    _pidParNames.push_back(mvaname);
    std::string qname(it->second.Name()); qname.append("Q");
    _pidPars.push_back(0.);
    _pidParNames.push_back(qname);
  }


  // Prepare separate maps of MVA variables
  // (necessity because TMVA::Reader::AddVariable(...) does not take const pointers)
  // Add MVA variables to the MVA readers
  for (variable_c_iterator it=_variables->GetMap()->begin();
                           it!=_variables->GetMap()->end();
                           it++)
  {
    _mvaVars.insert(std::pair<const char*, float>(it->second.Name(), 0.));
  }
  _mvaVars.insert(std::pair<const char*, float>("p", 0.));

  for(hypotheses_iterator ith = _hypotheses->begin(); ith != _hypotheses->end(); ith++)
  {
    for (variable_c_iterator it=_variables->GetMap()->begin();
                             it!=_variables->GetMap()->end();
                             it++)
    {
      ith->second.AddMVAVariable(it->second.Name(), &(_mvaVars.at(it->second.Name())));
    }
    ith->second.AddMVAVariable("seenP", &(_mvaVars.at("p")));
  }

  for(hypotheses_iterator ith = _hypotheses->begin(); ith != _hypotheses->end(); ith++) {
    // Dirty? Using particleType as vector index for _weightFileNames
    std::string wFileFullName(_weightFileNames.at(ith->first));
    wFileFullName += ".weights.xml";
    ith->second.BookMVA(_mvaMethod, wFileFullName);
    streamlog_out(DEBUG) << "Booked " << _mvaMethod << " with " << wFileFullName.c_str() << std::endl;
    std::string qFileName(_weightFileNames.at(ith->first));
    qFileName += ".Q.root";
    TFile qfile(qFileName.c_str());
    if (!qfile.IsOpen()) {
      streamlog_out(ERROR) << "Cannot open Q file " << qFileName.c_str() << std::endl;
      exit(0);
    }

    TH1F *histoQ;
    qfile.GetObject("histoQ", histoQ);
    if (!histoQ) {
      streamlog_out(ERROR) << "Cannot read histoQ from file " << qFileName.c_str() << std::endl;
      exit(0);
    }
    ith->second.SetHistoQ(histoQ);

    TH1F *histoSig;
    qfile.GetObject("sigMVA", histoSig);
    if (!histoSig) {
      streamlog_out(ERROR) << "Cannot read sigMVA from file " << qFileName.c_str() << std::endl;
      exit(0);
    }
    ith->second.SetHistoSig(histoSig);

    TH1F *histoBkg;
    qfile.GetObject("bkgMVA", histoBkg);
    if (!histoBkg) {
      streamlog_out(ERROR) << "Cannot read bkgMVA from file " << qFileName.c_str() << std::endl;
      exit(0);
    }
    ith->second.SetHistoBkg(histoBkg);

    TObjString *mvaCutString;
    qfile.GetObject("MVACut", mvaCutString);
    if (!mvaCutString) {
      streamlog_out(ERROR) << "Cannot read MVACut string object from file " << qFileName.c_str() << std::endl;
      exit(0);
    }
    ith->second.SetMVACut(mvaCutString->GetString().Atof());

    _mapNDecisionQ.insert(std::pair<particleType, unsigned int>(ith->first, 0));
    _mapNDecisionTot.insert(std::pair<particleType, unsigned int>(ith->first, 0));
    _mapNDecisionSigAbove.insert(std::pair<particleType, unsigned int>(ith->first, 0));
  }


  //mupi separation class
  _mupiPID = new LowMomentumMuPiSeparationPID_BDTG(_muPiWeightFileNames);

  _nEvt = _nPFO = _nUnidentified = _nDecisionQ = 0;
  _nEmptyClusters = _nEmptyShapes = _nEmptyTracks = _nZerodEdx = 0;

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
  if(_nEvt%100 == 0)
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

    _pidPars.clear();
    for (hypotheses_c_iterator ith=_hypotheses->begin(); ith!=_hypotheses->end(); ith++) {
      _pidPars.push_back(ith->second.GetMVAout());
      _pidPars.push_back(ith->second.GetQ());
    }


    if(_bestHypothesis != _hypotheses->end()) {
      _pidh->setParticleID(part, _bestHypothesis->first, _bestHypothesis->second.pdg,
                   GetQ(_bestHypothesis), _pidh->getAlgorithmID(algoName), _pidPars);
    }
    else {
      _pidh->setParticleID(part, -1, 0, -1., _pidh->getAlgorithmID(algoName), _pidPars);
    }
    _nPFO++;

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

  streamlog_out(MESSAGE) << "\n===========================================================\n";
  streamlog_out(MESSAGE) << "\nInfo on PID in this sample:\n";
  streamlog_out(MESSAGE) << "Analysed total " << _nPFO << " PFOs.\n";
  streamlog_out(MESSAGE) << "Total " << _nDecisionQ << " (" << TMath::Floor((1000.*_nDecisionQ)/_nPFO + .5)*.1
      << "%) decided by Q-statistic.\nOf that,\n";
  for(hypotheses_c_iterator ith=_hypotheses->begin(); ith!=_hypotheses->end(); ith++) {
    streamlog_out(MESSAGE) << ith->second.Name() << " " << _mapNDecisionQ.at(ith->first)
        << " times by Q and " << _mapNDecisionTot.at(ith->first) << " total.\n";
  }
  streamlog_out(MESSAGE) << "Total " << _nDecisionSigAbove << " (" << TMath::Floor((1000.*_nDecisionSigAbove)/_nPFO + .5)*.1
      << "%) decided by minimising S>.\nOf that,\n";
  for(hypotheses_c_iterator ith=_hypotheses->begin(); ith!=_hypotheses->end(); ith++) {
    streamlog_out(MESSAGE) << ith->second.Name() << " " << _mapNDecisionSigAbove.at(ith->first)
        << " times by minimising S> and " << _mapNDecisionTot.at(ith->first) << " total.\n";
  }
  streamlog_out(MESSAGE) << "Total " << _nUnidentified << " (" << TMath::Floor((1000.*_nUnidentified)/_nPFO + .5)*.1
      << "%) remained unidentified.\n";
  streamlog_out(MESSAGE) << "Empty clusters " << _nEmptyClusters << " times ("
                         << TMath::Floor(1000.0*_nEmptyClusters/_nPFO+.5)*.1 << "% of all PFOs).\n";
  streamlog_out(MESSAGE) << "Empty tracks " << _nEmptyTracks<< " times ("
                         << TMath::Floor(1000.0*_nEmptyTracks/_nPFO+.5)*.1 << "% of all PFOs).\n";
  streamlog_out(MESSAGE) << "Empty shapes " << _nEmptyShapes<< " times ("
                         << TMath::Floor(1000.0*_nEmptyShapes/_nPFO+.5)*.1 << "% of all PFOs).\n";
  streamlog_out(MESSAGE) << "dEdx < 1e-10 " << _nZerodEdx << " times ("
                         << TMath::Floor(1000.0*_nZerodEdx/_nPFO+.5)*.1 << "% of all PFOs).\n";

  delete _mupiPID; _mupiPID = NULL;
  delete _variables; _variables = NULL;
  _hypotheses->clear(); delete _hypotheses; _hypotheses = NULL;
}


/**********************************************************
 *
 * Private utility functions
 *
***********************************************************/

// Updates _variables, sets _mvaVars, evaluates MVA, selects best hypothesis,
// Checks lowE mu-pi separation, fills _pidPars
void MvaPidProcessor::Identify(ReconstructedParticle* particle) {

  short updateres = _variables->Update(particle);
  if (updateres & PIDVariables::MASK_EmptyClusters) _nEmptyClusters++;
  if (updateres & PIDVariables::MASK_EmptyTracks) _nEmptyTracks++;
  if (updateres & PIDVariables::MASK_EmptyShapes) _nEmptyShapes++;
  if (updateres & PIDVariables::MASK_ZerodEdx) _nZerodEdx++;

  if (updateres & (PIDVariables::MASK_EmptyClusters | PIDVariables::MASK_ZerodEdx) ) {
    _bestHypothesis = _hypotheses->end();
    return;
  }

  for(hypotheses_iterator ith=_hypotheses->begin(); ith != _hypotheses->end(); ith++) {
    // Not sure whether it is really necessary to repeatedly fill _mvaVars with the same
    // values within the hypothesis loop. However, TMVA::Reader does not promise to not
    // change the sensitive variables so I am just being cautious.
    for(variable_c_iterator it=_variables->GetMap()->begin();
                            it != _variables->GetMap()->end();
                            it++)
    {  _mvaVars.at(it->second.Name()) = it->second.Value();  }

//    ith->second.SetMVAout(_readerMap.at(ith->first)->EvaluateMVA(_mvaMethod));
    ith->second.Evaluate(TString(_mvaMethod.c_str()));
  }

  particleType bestH = PIDParticles::nParticleTypes;

  // Check which hypotheses pass the cut...
  std::vector<particleType> passH;
  for(hypotheses_iterator ith=_hypotheses->begin(); ith != _hypotheses->end(); ith++) {
    if(ith->second.PassesCut()) passH.push_back(ith->first);
  }

  // If neither hypothesis passes the cut, decide by minimising Q=-ln(S</B>)
  if(passH.size() == 0) {
    float minQ = FLT_MAX;
    for(hypotheses_iterator ith=_hypotheses->begin(); ith != _hypotheses->end(); ith++) {
      if( ith->second.GetQ() < minQ ) {
        minQ = ith->second.GetQ();
        bestH = ith->first;
      }
    }
    _nDecisionQ++;
    _mapNDecisionQ.at(bestH)++;
    _mapNDecisionTot.at(bestH)++;
  }

  // If only one passes the cut, select it
  if(passH.size() == 1) {
    bestH = passH.at(0);
    _mapNDecisionTot.at(bestH)++;
  }

  // If several pass the cut, select by minimising S>
  if(passH.size() > 1) {
    float minSigAbove = FLT_MAX;
    for(hypotheses_iterator ith=_hypotheses->begin(); ith != _hypotheses->end(); ith++) {
      if( ith->second.GetSigAbove() < minSigAbove ) {
        minSigAbove = ith->second.GetSigAbove();
        bestH = ith->first;
      }
    }
    _nDecisionSigAbove++;
    _mapNDecisionSigAbove.at(bestH)++;
    _mapNDecisionTot.at(bestH)++;
  }


  //mu-pi Separation for very low momentum tracks (from 0.2 GeV until 2 GeV)
  // It would be good to also make MuPISeparation that takes (ReconstructedParticleImpl*)
  // as argument.
  Float_t MVAoutput = -1.0;
  if((bestH == PIDParticles::muon || bestH == PIDParticles::pion)
      && _variables->GetP()<2.0)
  {
    streamlog_out(DEBUG) << "Checking mu/pi." << std::endl;
    TLorentzVector pp(TVector3(particle->getMomentum()), particle->getEnergy());
    EVENT::ClusterVec clu=particle->getClusters();
    lcio::Track* trk = particle->getTracks()[0];
    int parttype=_mupiPID->MuPiSeparation(pp, trk, clu);
    bestH = (parttype==1) ? PIDParticles::muon : PIDParticles::pion;
    MVAoutput = _mupiPID->getMVAOutput();
  }

  _bestHypothesis = _hypotheses->find(bestH);

}


