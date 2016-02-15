/*
 * MvaPidTraining.cc
 *
 *  Created on: Feb 8, 2016
 *      Author: S. Lukic
 */

#include "MvaPidTraining.hh"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCIterator.h>

#include "TFile.h"
#include "TCut.h"
#include "TMVA/Factory.h"
//#include "TObjArray.h"
#include "TCanvas.h"


MvaPidTraining aMvaPidTraining;

MvaPidTraining::MvaPidTraining() :
  Processor("MvaPidTraining"),
  _description("Training of particle ID using MVA"),
  //_rootfile(NULL),
  _tree(NULL),
  _seenP(0.), _truePDG(0), _isReconstructed(false),
  _signalPDG(0),
  _nEvt(0), _nMCPtot(0), _nRec(0), _nTrkCaloMismatch(0)
{
  registerInputCollection( LCIO::LCRELATION,
         "MCTruth2RecoLinkCollectionName" ,
         "true - reco relation collection"  ,
         _trueToReco,
         std::string("MCTruthRecoLink") ) ;


  registerInputCollection( LCIO::LCRELATION,
         "Reco2MCTruthLinkCollectionName" ,
         "reco - true relation collection"  ,
         _recoToTrue,
         std::string("RecoMCTruthLink") ) ;


  registerInputCollection( LCIO::MCPARTICLE,
         "MCParticleCollection" ,
         "Name of the MCParticle input collection"  ,
         _mcParticleCollectionName ,
         std::string("MCParticle") ) ;


  registerProcessorParameter( "SignalPDG" ,
            "PDG of the signal hypothesis for this training",
            _signalPDG,
            11 );


  registerProcessorParameter( "RootFileName" ,
            "Root file with MVA response data",
            _rootFileName,
            std::string("MvaPidTraining.response.root") );
/**/
  registerProcessorParameter( "MVAMethod" ,
            "MVA method name",
            _mvaMethod,
            std::string("BDT") );

  registerProcessorParameter( "MVAMethodOptions" ,
            "MVA method options",
            _mvaMethodOptions,
            std::string("") );

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
    _trainingVars.insert( std::pair<variableType, float>(it->first, 0.) );
//    _tree->Branch(it->second.Name(), it->second.Address());
    _tree->Branch(it->second.Name(), &(_trainingVars.at(it->first)));
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

    imaxweight = imaxtrckweight;
    maxweight = maxtrckweight;

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

    for(std::map<variableType, float>::iterator vit=_trainingVars.begin(); vit!=_trainingVars.end(); vit++)
    { vit->second = _variables.GetValue(vit->first); }
    _tree->Fill();

  }  // loop over MCPs


  _nEvt++;

  if (_nEvt%1000 == 0)
    streamlog_out(MESSAGE) <<  "=================== event " << _nEvt << std::endl ;

}


void MvaPidTraining::check( LCEvent * evt ) {

}



void MvaPidTraining::end() {

  TFile* outputFile = TFile::Open( _rootFileName.c_str(), "RECREATE" );
  TMVA::Factory * factory = new TMVA::Factory(_weightFileName, outputFile,
      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");



  for(variable_c_iterator vit=_variables.GetMap()->begin(); vit!=_variables.GetMap()->end(); vit++)
  {
    factory->AddVariable(vit->second.Name(), vit->second.Description(), vit->second.Unit(), 'F');
  }
  factory->AddVariable("seenP", "Measured momentum", "GeV", 'F');
  // Do we need these spectators for the cuts?
//  _factory->AddSpectator("truePDG", "True PDG", "", 'I');
//  _factory->AddSpectator("isReconstructed", "Boolean true if MC Particle reconstructed", "", 'B');

  TCut signalCut("signalCut", Form("truePDG==%d&&isReconstructed", _signalPDG));
  TCut backgroundCut("backgroundCut", Form("(truePDG!=%d)&&isReconstructed", _signalPDG));

/**/  TObjArray* listb = _tree->GetListOfBranches();

  std::cout << "Branches in the training tree:\n";
//  for(TObjArray::Iterator_t bit=listb->begin(); bit!=listb->end(); bit++) {
  for(int i=0; i<listb->GetEntries(); i++) {
    std::cout << ((TBranch*)(listb->At(i)))->GetName() << std::endl;
  }

  TCanvas *c = new TCanvas;
  c->SetLogy();
  for(variable_c_iterator vit=_variables.GetMap()->begin(); vit!=_variables.GetMap()->end(); vit++) {
    const char *var = vit->second.Name();
    _tree->Draw(var, signalCut);
    c->Print(Form("signal_%s.pdf", var));
    _tree->Draw(var, backgroundCut);
    c->Print(Form("background_%s.pdf", var));
    _tree->Draw(var, "truePDG==13&&isReconstructed");
    c->Print(Form("muons_%s.pdf", var));
  }
  delete c; c=NULL;

  factory->SetInputTrees(_tree, signalCut, backgroundCut);

  // Fixme: Is this ok or needs a map of available methods with options
  // Or should options be steerable?
  factory->BookMethod(_mvaMethod, _mvaMethod, _mvaMethodOptions);
  factory->PrintHelpMessage(_mvaMethod);

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  // TODO:
  // Output optimal cuts (need a definition of "optimal"), impact of variables etc.
  // Store optimal cuts

  // Output messages on mismatched track/calo etc.

  outputFile->Close();

}

