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
  _tree(NULL),
  _seenP(0.), _truePDG(0), _isReconstructed(false),
  _hasClusters(false), _hasShapes(false), _hasdEdx(false),
  _signalPDG(0),
  _nEvt(0), _nMCPtot(0), _nRec(0), _nTrkCaloMismatch(0),
  _nEmptyClusters(0), _nEmptyTracks(0), _nEmptyShapes(0),
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


  registerProcessorParameter( "MVAResponseFileName" ,
            "Root file with MVA response data",
            _mvaResponseFileName,
            std::string("MvaPidTraining.response.root") );


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
  _nEmptyClusters = _nEmptyShapes = _nEmptyTracks = _nZerodEdx = 0;

  _tree = new TTree("varTree","varTree");

  for(unsigned int i = 0; i < _variables.GetVariables()->size(); i++)
  {
    _tree->Branch(_variables.GetVariables()->at(i)->Name(), &(_variables.GetMvaVariables()->operator [](i)));
  }

  _tree->Branch("seenP",&_seenP) ;
  _tree->Branch("truePDG",&_truePDG) ;
  _tree->Branch("isReconstructed",&_isReconstructed) ;
  _tree->Branch("hasClusters",&_hasClusters) ;
  _tree->Branch("hasShapes",&_hasShapes) ;
  _tree->Branch("hasdEdx",&_hasdEdx) ;

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
      short updateres = _variables.Update(rcp);
      if (updateres & PIDVariable_base::MASK_EmptyClusters) { _nEmptyClusters++; _hasClusters=false; }
      else { _hasClusters = true; }
      if (updateres & PIDVariable_base::MASK_EmptyTracks) _nEmptyTracks++;
      if (updateres & PIDVariable_base::MASK_EmptyShapes) { _nEmptyShapes++; _hasShapes=false; }
      else { _hasShapes = true; }
      if (updateres & PIDVariable_base::MASK_ZerodEdx) { _nZerodEdx++; _hasdEdx = false; }
      else { _hasdEdx = true; }
      _seenP = p3.Mag();

    } // if reco part
    // TODO: IMPROVE HERE: CHECK if track or cluster exists even if no PF!
    else {
      streamlog_out(DEBUG) << "No ReconstructedParticle found for this mcp!" << std::endl ;

      _variables.SetOutOfRange();
      _seenP = 0.;

    }

    _tree->Fill();

  }  // loop over MCPs


  _nEvt++;

  if (_nEvt%1000 == 0)
    streamlog_out(MESSAGE) <<  "=================== event " << _nEvt << std::endl ;

}


void MvaPidTraining::check( LCEvent * evt ) {

}



void MvaPidTraining::end() {


//  _tree->Write();
//  exit(0);

//  TFile* outputFile = TFile::Open( _mvaResponseFileName.c_str(), "RECREATE" );
  TMVA::Factory * factory = new TMVA::Factory(_weightFileName, (TFile*)(gDirectory),
      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  // Add sensitive variables - they are known from the map
  for(variable_c_iterator vit=_variables.GetVariables()->begin();
      vit!=_variables.GetVariables()->end(); vit++)
  {
//    std::cout << "Adding variable:\nExpression: " << (*vit)->Name() << "; title: "
//              << (*vit)->Description() << "; Unit: " << (*vit)->Unit() << std::endl;
    factory->AddVariable((*vit)->Name(), (*vit)->Description(), (*vit)->Unit(), 'F');
  }
//  exit(0);
  factory->AddSpectator("seenP", "Measured momentum", "GeV", 'F');

  // Recognise signal and background from the truePDG - only reconstructed particles.
  const char * basicCut = "isReconstructed&&hasClusters&&hasShapes&&hasdEdx";
  TCut signalCut("signalCut", Form("(abs(truePDG)==%d)&&%s", _signalPDG, basicCut));
  TCut backgroundCut("backgroundCut", Form("(abs(truePDG)!=%d)&&%s", _signalPDG, basicCut));

  factory->SetInputTrees(_tree, signalCut, backgroundCut);

  factory->BookMethod(_mvaMethod, _mvaMethod, _mvaMethodOptions);
//  factory->PrintHelpMessage(_mvaMethod);

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();


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
    nBkgBelow += bkgMVA.GetBinContent(ibin);
    nSigBelow += sigMVA.GetBinContent(ibin);
    nBkgAbove -= bkgMVA.GetBinContent(ibin);
    float effSig = nSigAbove/nSig;
    float effBkg = nBkgAbove/nSig;
//    if(effSig>.99) mvaCut = sigMVA.GetBinLowEdge(ibin+1);
    if(effBkg>.01) mvaCut = sigMVA.GetBinLowEdge(ibin);
    float q;
// Outer Q
//    if (nBkgBelow > 0) { q = effSig*nBkg/nBkgBelow; }
//    else { q = FLT_MAX; }
// Inner Q
    if (nBkgAbove < 1.e6*FLT_MIN) { nBkgAbove = 1.e6*FLT_MIN; }
//    q = - ( TMath::Log(nSigBelow) - TMath::Log(nBkgAbove) );
    if (nBkgAbove < 1.e6*FLT_MIN) { nBkgAbove = 1.e6*FLT_MIN; }
    q = - (TMath::Log(effSig) + TMath::Log(nSigAbove) - TMath::Log(nSigAbove+nBkgAbove));

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
  streamlog_out(MESSAGE) << "dEdx < 1e-10 " << _nZerodEdx << " times ("
                         << TMath::Floor(1000.0*_nZerodEdx/_nRec+.5)*.1 << "% of all PFOs).\n";

}

