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

#include "TCanvas.h"
#include "TH1F.h"


MvaPidTraining aMvaPidTraining;


MvaPidTraining::MvaPidTraining() :
  Processor("MvaPidTraining"),
  _description("Training of particle ID using MVA"),
  _treeBkgTraining(NULL), _treeSigTraining(NULL),
  _treeBkgTest(NULL), _treeSigTest(NULL),
  _variables(NULL),
  _seenP(0.), _truePDG(0), _isReconstructed(false),
  _hasClusters(false), _hasShapes(false), _hasdEdx(false), _hasMomentum(false),
  _signalPDG(0), _pMin(0.), _pMax(0.), _usedVars(0),
  _nEvt(0), _nMCPtot(0), _nRec(0), _nTrkCaloMismatch(0),
  _nInvalidMomentum(0), _nEmptyClusters(0), _nEmptyTracks(0), _nEmptyShapes(0),
  _nZerodEdx(0)
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

  std::vector< std::string > usedVars;

  registerProcessorParameter( "UsedVariables" ,
            "List of used variables",
            _usedVars,
            usedVars );

  registerProcessorParameter( "PMin" ,
            "Minimum measured momentum for this training",
            _pMin,
            float(5.0) );

  registerProcessorParameter( "PMax" ,
            "Maximum measured momentum for this training",
            _pMax,
            float(20.) );


  registerProcessorParameter( "MVAResponseFileName" ,
            "Root file with MVA response data",
            _mvaResponseFileName,
            std::string("MvaPidTraining.response") );


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
            std::string("MvaPid_electron") );
}


void MvaPidTraining::init() {

  streamlog_out(DEBUG) << "   init called  "
           << std::endl ;

  // usually a good idea to
  printParameters() ;

  _variables = new PIDVariables_MvaPid(ParticleTypeByPDG(_signalPDG));

  _nEvt = 0;
  _nMCPtot = _nRec = _nTrkCaloMismatch = 0;
  _nEmptyClusters = _nEmptyShapes = _nEmptyTracks = _nZerodEdx = 0;

  _treeBkgTraining = new TTree("treeBkgTraining","treeBkgTraining");
  _treeSigTraining = new TTree("treeSigTraining","treeSigTraining");
  _treeBkgTest = new TTree("treeBkgTest","treeBkgTest");
  _treeSigTest = new TTree("treeSigTest","treeSigTest");

  for(unsigned int i = 0; i < _variables->GetVariables()->size(); i++)
  {
    _treeBkgTraining->Branch(_variables->GetVariables()->at(i)->Name(), &(_variables->GetMvaVariables()->operator [](i)), _variables->GetVariables()->at(i)->Name());
    _treeSigTraining->Branch(_variables->GetVariables()->at(i)->Name(), &(_variables->GetMvaVariables()->operator [](i)), _variables->GetVariables()->at(i)->Name());
    _treeBkgTest->Branch(_variables->GetVariables()->at(i)->Name(), &(_variables->GetMvaVariables()->operator [](i)), _variables->GetVariables()->at(i)->Name());
    _treeSigTest->Branch(_variables->GetVariables()->at(i)->Name(), &(_variables->GetMvaVariables()->operator [](i)), _variables->GetVariables()->at(i)->Name());
  }

  _treeBkgTraining->Branch("seenP",&_seenP, "seenP/F") ;
  _treeBkgTraining->Branch("truePDG",&_truePDG, "truePDG/I") ;
  _treeBkgTraining->Branch("isReconstructed",&_isReconstructed, "isReconstructed/O") ;

  _treeSigTraining->Branch("seenP",&_seenP, "seenP/F") ;
  _treeSigTraining->Branch("truePDG",&_truePDG, "truePDG/I") ;
  _treeSigTraining->Branch("isReconstructed",&_isReconstructed, "isReconstructed/O") ;

  _treeBkgTest->Branch("seenP",&_seenP, "seenP/F") ;
  _treeBkgTest->Branch("truePDG",&_truePDG, "truePDG/I") ;
  _treeBkgTest->Branch("isReconstructed",&_isReconstructed, "isReconstructed/O") ;

  _treeSigTest->Branch("seenP",&_seenP, "seenP/F") ;
  _treeSigTest->Branch("truePDG",&_truePDG, "truePDG/I") ;
  _treeSigTest->Branch("isReconstructed",&_isReconstructed, "isReconstructed/O") ;
  // The following are not needed as branches. Events with any of these false are not written
//  _treeBkgTraining->Branch("hasClusters",&_hasClusters, "hasClusters/O") ;
//  _treeBkgTraining->Branch("hasShapes",&_hasShapes, "hasShapes/O") ;
//  _treeBkgTraining->Branch("hasdEdx",&_hasdEdx, "hasdEdx/O") ;
//  _treeBkgTraining->Branch("hasMomentum",&_hasMomentum, "hasMomentum/O") ;

  // Not working:
 // _treeSigTraining->CopyAddresses(_treeBkgTraining);
 // _treeBkgTest->CopyAddresses(_treeBkgTraining);
 // _treeSigTest->CopyAddresses(_treeBkgTraining);

  // Randomiser for smearing discrete values of variables
  PIDVariable_base::varRand = new TRandom3;

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
    // TODO: IMPROVE HERE: CHECK if track or cluster exists even if no PF!

    /*** Reject if no tracks found ***/
    if (!_isReconstructed) {
      streamlog_out(DEBUG) << "No ReconstructedParticle found for this mcp!" << std::endl ;
      continue;
    }

    _nRec++;

    ReconstructedParticle* rcp =  (ReconstructedParticle*) recovec.at(imaxweight);
    EVENT::TrackVec trax = rcp->getTracks();
    TVector3 p3(rcp->getMomentum());
    _seenP = p3.Mag();

    /*** Reject if outside of the momentum range ***/
    if (_seenP<_pMin||_seenP>_pMax) continue;

    //  PID sensitive variables  ***/
    short updateres = _variables->Update(rcp);
    if (updateres & PIDVariable_base::MASK_InvalidMomentum) { _nInvalidMomentum++; _hasMomentum=false; }
    else { _hasMomentum = true; }
    if (updateres & PIDVariable_base::MASK_EmptyClusters) { _nEmptyClusters++; _hasClusters=false; }
    else { _hasClusters = true; }
    if (updateres & PIDVariable_base::MASK_EmptyTracks) _nEmptyTracks++;
    if (updateres & PIDVariable_base::MASK_EmptyShapes) { _nEmptyShapes++; _hasShapes=false; }
    else { _hasShapes = true; }
    if (updateres & PIDVariable_base::MASK_ZerodEdx) { _nZerodEdx++; _hasdEdx = false; }
    else { _hasdEdx = true; }

    /*** Reject if clusters or dEdx missing ***/
    if (!(_hasClusters&&_hasdEdx)) continue;

    /*** Select and fill the correct tree ***/
    if(TMath::Abs(_truePDG)==_signalPDG) {
      if(PIDVariable_base::varRand->Uniform() < .5) {
        _treeSigTraining->Fill();
      }
      else {
        _treeSigTest->Fill();
      }
    }
    else { // background PDG
      if(PIDVariable_base::varRand->Uniform() < .5) {
        _treeBkgTraining->Fill();
      }
      else {
        _treeBkgTest->Fill();
      }
    }

  }  // loop over MCPs


  _nEvt++;

  if (_nEvt%1000 == 0)
    streamlog_out(MESSAGE) <<  "=================== event " << _nEvt << std::endl ;

}


void MvaPidTraining::check( LCEvent * evt ) {

}



void MvaPidTraining::end() {

//  streamlog_out(DEBUG) << "Clearing.\n";
//  _variables->ClearVars();
//  streamlog_out(DEBUG) << "Cleared.\n";

/*  _treeSigTraining->Write();
  _treeBkgTraining->Write();
  _treeSigTest->Write();
  _treeBkgTest->Write();
  exit(0);
*/

  _weightFileName += TString::Format("_%.1f-%.1fGeV", _pMin, _pMax);
  _mvaResponseFileName += TString::Format("_%.1f-%.1fGeV.root", _pMin, _pMax);

  TFile* outputFile = TFile::Open( _mvaResponseFileName.c_str(), "RECREATE" );
//  TMVA::Factory * factory = new TMVA::Factory(_weightFileName, (TFile*)(gDirectory),
  TMVA::Factory * factory = new TMVA::Factory(_weightFileName, outputFile,
      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  factory->AddSignalTree(_treeSigTraining, 1., TMVA::Types::kTraining);
  factory->AddBackgroundTree(_treeBkgTraining, 1., TMVA::Types::kTraining);
  factory->AddSignalTree(_treeSigTest, 1., TMVA::Types::kTesting);
  factory->AddBackgroundTree(_treeBkgTest, 1., TMVA::Types::kTesting);

  // Add sensitive variables - they are known from the map
  if ( _usedVars.size() > 0 ) {
    for(unsigned int i=0; i<_usedVars.size(); i++)
    {
      PIDVariable_base* pvar = _variables->FindVariable(_usedVars.at(i));
      if(pvar) {
        factory->AddVariable(pvar->Name(), pvar->Description(), pvar->Unit(), 'F');
        streamlog_out(MESSAGE) << "Adding variable " << pvar->Name() << " to training.\n";
      }
      else {
        streamlog_out(WARNING) << "Requested addition of unknown variable "
            << _usedVars.at(i) << ". Skipping." << std::endl;
      }
    }
  }
  else {
    for(unsigned int i=0; i<_variables->GetVariables()->size(); i++)
    {
      PIDVariable_base* pvar = _variables->GetVariables()->at(i);
      factory->AddVariable(pvar->Name(), pvar->Description(), pvar->Unit(), 'F');
      streamlog_out(MESSAGE) << "Adding variable " << pvar->Name() << " to training.\n";
    }
  }
  factory->AddSpectator("seenP", "Measured momentum", "GeV", 'F');

//  factory->PrepareTrainingAndTestTree( "", "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  factory->BookMethod(_mvaMethod, _mvaMethod, _mvaMethodOptions);
//  factory->PrintHelpMessage(_mvaMethod);


  factory->TrainAllMethods();
  factory->TestAllMethods();
//  outputFile->Write();
  factory->EvaluateAllMethods();
//  outputFile->Write();


  /*** Create, fill and store the histogram of the Q-statistic ***/
  streamlog_out(MESSAGE) << "\n===========================================================\n";
  streamlog_out(MESSAGE) << "\nGenerating and storing the Q-statistic:\n\n";
  TTree *testTree = NULL;
  gDirectory->GetObject("TestTree", testTree);
  if(!testTree) {
    streamlog_out(MESSAGE) << "No TestTree in gDirectory. Sorry, abandoning.\n";
    return;
  }

  TH1F histoQ("histoQ", "Q statistic; MVA response; Q", nChanQ, -1., 1.);
  histoQ.SetDirectory(0);

  streamlog_out(MESSAGE) << "Projecting signal MVA from the test tree.\n";
  TH1F sigMVA(histoQ);
  sigMVA.SetName("sigMVA"); sigMVA.SetTitle("Signal MVA response; MVA; Count");
  testTree->Project("sigMVA", _mvaMethod.c_str(), "classID==0");
  sigMVA.Scale(1./sigMVA.Integral());

  streamlog_out(MESSAGE) << "Projecting background MVA from the test tree.\n";
  TH1F bkgMVA(histoQ);
  bkgMVA.SetName("bkgMVA"); bkgMVA.SetTitle("Background MVA response; MVA; Count");
  testTree->Project("bkgMVA", _mvaMethod.c_str(), "classID==1");
  bkgMVA.Scale(1./bkgMVA.Integral());

  sigMVA.SetDirectory(0);
  bkgMVA.SetDirectory(0);

  double nSig = sigMVA.Integral();
  double nBkg = bkgMVA.Integral();
  streamlog_out(MESSAGE) << "Done projecting. nSig = " << nSig << "; nBkg = " << nBkg << ".\n";

  double nSigAbove = nSig;
  double nBkgBelow = 0;
  double nSigBelow = 1.e-10; // We will never have 1e10 events in training
  double nBkgAbove = nBkg;
  float mvaCut = -1.;

  for (int ibin=1; ibin<histoQ.GetNbinsX(); ibin++) {
    // Above and below refer to the upper edge of the bin
    nSigAbove -= sigMVA.GetBinContent(ibin);
    if (nSigAbove < 1.e6*FLT_MIN) { nSigAbove = 1.e6*FLT_MIN; }
    nBkgBelow += bkgMVA.GetBinContent(ibin);
    nSigBelow += sigMVA.GetBinContent(ibin);
    nBkgAbove -= bkgMVA.GetBinContent(ibin);
    if (nBkgAbove < 1.e6*FLT_MIN) { nBkgAbove = 1.e6*FLT_MIN; }
    float effSig = nSigAbove/nSig;
    float effBkg = nBkgAbove/nSig;
//    if(effSig>.99) mvaCut = sigMVA.GetBinLowEdge(ibin+1);
    if(effBkg>.01) mvaCut = sigMVA.GetBinLowEdge(ibin);
    float q;
// Outer Q
//    if (nBkgBelow > 0) { q = effSig*nBkg/nBkgBelow; }
//    else { q = FLT_MAX; }
// Inner Q
//    q = - ( TMath::Log(nSigBelow) - TMath::Log(nBkgAbove) );
    q = - (TMath::Log(effSig) + TMath::Log(nSigAbove) - TMath::Log(nSigAbove+nBkgAbove));
    if (TMath::IsNaN(q)) {
      streamlog_out(DEBUG) << "Q is NaN. ibin = " << ibin << "; effSig = " << effSig << "; effBkg = "
          << effBkg << "; nSigAbove = " << nSigAbove << "; nBkgAbove = " << nBkgAbove << std::endl;
    }

    streamlog_out(DEBUG) << "Setting content in bin " << ibin << " to " << q << std::endl;
    histoQ.SetBinContent(ibin, q);
  }


  streamlog_out(MESSAGE) << "MVA cut for signal efficiency > 99% is " << mvaCut << ".\n";
  TObjString mvaCutString(TString::Format("%7.4f", mvaCut));

  TFile qfile(TString::Format("weights/%s_%s.Q.root",
      _weightFileName.c_str(), _mvaMethod.c_str()), "RECREATE");
  histoQ.Write();
  sigMVA.Write();
  bkgMVA.Write();
  mvaCutString.Write("MVACut");

  streamlog_out(MESSAGE) << "\n===========================================================\n";
  streamlog_out(MESSAGE) << "\nInfo on the training sample:\n";
  streamlog_out(MESSAGE) << "Analysed total " << _nMCPtot << " final MC particles.\n";
  streamlog_out(MESSAGE) << "Found total " << _nRec << " reconstructed PFO object linked to the analysed MC particles.\n";
  streamlog_out(MESSAGE) << "Reconstructed percentage: " << TMath::Floor(1000.0*_nRec/_nMCPtot +.5)*.1 << "%.\n";
  streamlog_out(MESSAGE) << "Best cluster and best track found in different PFOs " << _nTrkCaloMismatch << " times ("
                         << TMath::Floor(1000.0*_nTrkCaloMismatch/_nRec+.5)*.1 << "% of all PFOs).\n";
  streamlog_out(MESSAGE) << "Empty clusters " << _nEmptyClusters << " times ("
                         << TMath::Floor(1000.0*_nEmptyClusters/_nRec+.5)*.1 << "% of all PFOs).\n";
  streamlog_out(MESSAGE) << "Empty tracks " << _nEmptyTracks<< " times ("
                         << TMath::Floor(1000.0*_nEmptyTracks/_nRec+.5)*.1 << "% of all PFOs).\n";
  streamlog_out(MESSAGE) << "Empty shapes " << _nEmptyShapes<< " times ("
                         << TMath::Floor(1000.0*_nEmptyShapes/_nRec+.5)*.1 << "% of all PFOs).\n";
  streamlog_out(MESSAGE) << "dEdx < 1e-15 MIP " << _nZerodEdx << " times ("
                         << TMath::Floor(1000.0*_nZerodEdx/_nRec+.5)*.1 << "% of all PFOs).\n";

}



const PIDParticle_base* MvaPidTraining::ParticleTypeByPDG(int pdg) {
  switch(abs(pdg)) {
  case 11:
    return &PIDParticles::electronProperties;
  case 13:
    return &PIDParticles::muonProperties;
  case 211:
    return &PIDParticles::pionProperties;
  case 321:
    return &PIDParticles::kaonProperties;
  case 2212:
    return &PIDParticles::protonProperties;
  default:
    streamlog_out(ERROR) << "Unknown hypothesis with PDG " << pdg << ". Aborting.\n";
    exit(0);
  }
}
