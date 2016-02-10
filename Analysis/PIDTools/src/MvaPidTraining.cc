/*
 * MvaPidTraining.cc
 *
 *  Created on: Feb 8, 2016
 *      Author: S. Lukic
 */

#include "TCut.h"

#include "MvaPidTraining.hh"

MvaPidTraining::MvaPidTraining() :
  Processor("LikelihoodPIDProcessor"),
  _description("Training of particle ID using MVA"),
  _factory(NULL), _rootfile(NULL), _tree(NULL),
  _seenP(0.), _truePDG(0), _isReconstructed(false),
  _signalPDG(0),
  _nEvt(0), _nMCPtot(0), _nRec(0), _nTrkCaloMismatch(0)
{
  registerProcessorParameter( "SignalPDG" ,
            "PDG of the signal hypothesis for this training",
            _signalPDG,
            11 );

  registerProcessorParameter( "RootFileName" ,
            "Root file with data for training",
            _rootFileName,
            std::string("MvaPidTraining.root") );

  registerProcessorParameter( "MVAMethod" ,
            "MVA method name",
            _mvaMethod,
            std::string("BDT") );

  registerProcessorParameter( "WeightFileName" ,
            "File name to write weights",
            _weightFileName,
            std::string("MvaPid_electron.xml") );
}


void MvaPidTraining::init() {

  streamlog_out(DEBUG) << "   init called  "
           << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nEvt = 0;
  _nMCPtot = _nRec = _nTrkCaloMismatch = 0;

  _tree = new TTree("varTree","varTree");

  for(variable_c_iterator it = _variables.GetMap()->begin();
      it != _variables.GetMap()->end();
      it++)
  {
//    _tree->Branch(it->second.Name(), &sensitiveVars[it->first]);
    _tree->Branch(it->second.Name(), it->second.Address());
  }
  _tree->Branch("seenP",&_seenP) ;
  _tree->Branch("truePDG",&_truePDG) ;
  _tree->Branch("isReconstructed",&_isReconstructed) ;

}


void MvaPidTraining::processRunHeader( LCRunHeader* run) {

}

void MvaPidTraining::processEvent( LCEvent * evt ) {

  _truePDG = 0;
  _seenP = 0.;


  streamlog_out(DEBUG) << " iterator and navigator " << std::endl;
  LCIterator<MCParticle> mcpIt( evt, _mcParticleCollectionName ) ;
  streamlog_out(DEBUG) << " got mcpIt " << mcpIt.size() << std::endl;
  LCRelationNavigator mc2recoNav(evt->getCollection( _trueToReco ));
  streamlog_out(DEBUG) << " got mc2recoNav from " << mc2recoNav.getFromType() << " to " << mc2recoNav.getToType() << std::endl;

  //----------------------------------------------------------------------------------------------------------------------------
  // loop over MCParticles
  //----------------------------------------------------------------------------------------------------------------------------

  while( MCParticle* mcp = mcpIt.next()  ) {

    if (!mcp) continue;

    streamlog_out(DEBUG) << "MCparticle id = " << mcp->getPDG() << ", genstat = " << mcp->getGeneratorStatus() << std::endl;

    if (mcp->getGeneratorStatus() != 1) continue;

    _nMCPtot++;
    _truePDG = mcp->getPDG();

    streamlog_out(DEBUG) << " get reco particle " << std::endl;
    const EVENT::LCObjectVec& recovec = mc2recoNav.getRelatedToObjects(mcp);
    streamlog_out(DEBUG) << " recovec has length " << recovec.size() << std::endl;
    const EVENT::FloatVec& recoweightvec = mc2recoNav.getRelatedToWeights(mcp);
    streamlog_out(DEBUG) << " recoweightvec has length " << recoweightvec.size() << std::endl;
    double maxtrckweight = 0;
    double maxcaloweight = 0;
    double maxweight = 0;
    int imaxtrckweight = -1;
    int imaxcaloweight = -1;
    int imaxweight = -1;
    for (unsigned int irel = 0; irel < recovec.size(); irel++) {
      double tmptrackw = double((int(recoweightvec.at(irel))%10000))/1000;
      double tmpcalow = double((int(recoweightvec.at(irel))/10000))/1000;

      streamlog_out(DEBUG) << " irel " << irel << ", recoweight = " << int(recoweightvec.at(irel))
                           << ", tmptrackw = " << tmptrackw
                           << ", tmpcalow = " << tmpcalow << std::endl;
      if (tmptrackw > maxtrckweight) {
        imaxtrckweight = irel;
        maxtrckweight = tmptrackw;
      }
      if (tmpcalow > maxcaloweight) {
        imaxcaloweight = irel;
        maxcaloweight = tmpcalow;
      }
    }

    if (maxtrckweight > 0) {
      imaxweight = imaxtrckweight;
      maxweight = maxtrckweight;
    }
    if (imaxcaloweight != imaxtrckweight) _nTrkCaloMismatch++;

    streamlog_out(DEBUG) << "Found reco particle for mcp at imaxweight = " << imaxweight << " with weight = " << maxweight << std::endl ;

    _isReconstructed = (imaxweight>=0);

    if (_isReconstructed) {

      _nRec++;

      ReconstructedParticle* rcp =  (ReconstructedParticle*) recovec.at(imaxweight);
      EVENT::TrackVec trax = rcp->getTracks();
      TVector3 p3(rcp->getMomentum());

      //  PID sensitive variables  ***/
      _variables.Update(rcp);
      _seenP = p3.Mag();

    } // if reco part
    // TODO: IMPROVE HERE: CHECK if track or cluster exists even if no PF!
    else {
      streamlog_out(DEBUG) << "No ReconstructedParticle found for this mcp!" << std::endl ;

      _variables.SetOutOfRange();
      _seenP = 0.;

    }
  }  // loop over MCPs


  _tree->Fill();
  _nEvt++;

  if (_nEvt%1000 == 0)
    streamlog_out(MESSAGE) <<  "=================== event " << _nEvt << std::endl ;

}




void MvaPidTraining::end() {

  // FIXME: Can be local?
  _factory = new TMVA::Factory( "TMVAClassification", _weightFileName,
      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");


  for(variable_c_iterator vit=_variables.GetMap()->begin(); vit!=_variables.GetMap()->end(); vit++)
  {
    _factory->AddVariable(vit->second.Name(), vit->second.Description(), vit->second.Unit(), 'F');
  }
  _factory->AddVariable("seenP", "Measured momentum", "GeV", 'F');
//  _factory->AddSpectator("truePDG", "True PDG", "", 'I');
//  _factory->AddSpectator("isReconstructed", "Boolean true if MC Particle reconstructed", "", 'B');

  TCut signalCut("signalCut", Form("truePDG==%d&&isReconstructed", _signalPDG));
  TCut backgroundCut("backgroundCut", Form("truePDG!=%d&&isReconstructed", _signalPDG));
  _factory->SetInputTrees(_tree, signalCut, backgroundCut);

  _factory->BookMethod(_mvaMethod, _mvaMethod); // Fixme: Fix this

  _factory->TrainAllMethods();

  // Test
  // Evaluate
  // Store info on "optimal" cuts, impact of variables etc.

  // Output messages mismatched track/calo etc.

}

