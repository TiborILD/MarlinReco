#include <CluShapeTreeProcessor.hh>

#include <UTIL/LCRelationNavigator.h>
#include "UTIL/LCIterator.h"
#include <UTIL/PIDHandler.h>

#include "TROOT.h"
#include "TString.h"
#include <TH2F.h>



CluShapeTreeProcessor aCluShapeTreeProcessor ;


CluShapeTreeProcessor::CluShapeTreeProcessor() : Processor("CluShapeTreeProcessor"),
    nEvt(0),
    varTree(NULL),
    nMCParticles(0), nMCPtot(0), nRec(0), nTrkCaloMismatch(0)
{
  
  // modify processor description
  _description = "Creates ROOT PDF histograms of variables used by the "
      "LikelihoodPIDProcessor. Also fills a tree "
      "with the same variables for detailed analysis of the PDF.";

  
  // register steering parameters: name, description, class-variable, default value


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


  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "PFOs" , 
                            "particle flow objects"  ,
                            _pandoraPFOs ,
                            std::string("PandoraPFOs") ) ;
 
  registerInputCollection( LCIO::TRACK,
			   "StudiedTracks" , 
			   "Name of the FullLDC track collection"  ,
			   _trackColName ,
			   std::string("MarlinTrkTracks") ) ;  

}


void CluShapeTreeProcessor::init() {

  streamlog_out(DEBUG) << "   init called  "
		       << std::endl ;

  // usually a good idea to
  printParameters() ;
  nEvt = 0;
  nMCParticles = 0;
  nMCPtot = nRec = nTrkCaloMismatch = 0;

  gROOT->ProcessLine("#include <vector>");

  varTree = new TTree("varTree","varTree");
  varTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I") ;

  // Cluster shape variables
  for(unsigned int i=0; i<_variables.Size(); i++) {
    _treeVars.push_back(new vector<double>);
    varTree->Branch(_variables.Name(i), _treeVars.back());
  }

  // Other variables for study
  varTree->Branch("trueP",&trueP) ;
  varTree->Branch("truePt",&truePt) ;
  varTree->Branch("trueTheta",&trueTheta) ;
  varTree->Branch("truePhi",&truePhi) ;
  varTree->Branch("trueCharge",&trueCharge) ;
  varTree->Branch("truePDG",&truePDG) ;

  varTree->Branch("isReconstructed",&isReconstructed) ;
  varTree->Branch("isSeen",&isSeen) ;
  varTree->Branch("seenP",&seenP) ;
  varTree->Branch("seenPt",&seenPt) ;
  varTree->Branch("seenTheta",&seenTheta) ;
  varTree->Branch("seenPhi",&seenPhi) ;
  varTree->Branch("seenCharge",&seenCharge) ;

}

void CluShapeTreeProcessor::processRunHeader( LCRunHeader* run) {
    
} 

void CluShapeTreeProcessor::processEvent( LCEvent * evt ) {

  //streamlog_out(MESSAGE) << " start processing event # " << evt->getEventNumber() << std::endl;

  
  streamlog_out(DEBUG) << " clearing vectors " << std::endl;

  nMCParticles = 0;

  for(unsigned int i=0; i<_treeVars.size(); i++) {
    _treeVars.at(i)->clear();
  }

  trueP.clear();
  truePt.clear();
  trueTheta.clear();
  truePhi.clear();
  trueCharge.clear();
  truePDG.clear();
  isReconstructed.clear();
  isSeen.clear();
  seenP.clear();
  seenPt.clear();
  seenTheta.clear();
  seenPhi.clear();
  seenCharge.clear();


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
    
    streamlog_out(DEBUG) << " mcparticle id = " << mcp->getPDG() << ", genstat = " << mcp->getGeneratorStatus() << std::endl;
    
    if (mcp->getGeneratorStatus() != 1) continue;
    
    nMCParticles++;
    nMCPtot++;

    TVector3 v( mcp->getVertex() );
    TVector3 e( mcp->getEndpoint() );
    TVector3 p( mcp->getMomentum() );

    streamlog_out(DEBUG) << " start push_back " << std::endl;
    trueP.push_back(p.Mag());
    truePt.push_back(p.Perp());
    trueTheta.push_back(p.Theta());
    truePhi.push_back(p.Phi());
    trueCharge.push_back(mcp->getCharge());
    truePDG.push_back(mcp->getPDG());



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
    
    streamlog_out(DEBUG) << " found reco particle for mcp at imaxtrackweight = " << imaxtrckweight << " with weight = " << maxtrckweight << std::endl ;
    streamlog_out(DEBUG) << " found reco particle for mcp at imaxcaloweight = " << imaxcaloweight << " with weight = " << maxcaloweight << std::endl ;

//    imaxweight = imaxcaloweight;
//    maxweight = maxcaloweight;
    if (maxtrckweight > 0) {
      imaxweight = imaxtrckweight;
      maxweight = maxtrckweight;
    }
    if (imaxcaloweight != imaxtrckweight) nTrkCaloMismatch++;
 //      streamlog_out(WARNING) << " imaxcaloweight != imaxtrckweight, choosing recoparticle with larger fraction" << std::endl ;


    streamlog_out(DEBUG) << " found reco particle for mcp at imaxweight = " << imaxweight << " with weight = " << maxweight << std::endl ;


    isSeen.push_back(maxweight);
    isReconstructed.push_back(imaxweight>=0);

    if (imaxweight >= 0) {

      nRec++;

      ReconstructedParticle* rcp =  (ReconstructedParticle*) recovec.at(imaxweight);
      EVENT::TrackVec trax = rcp->getTracks();
      TVector3 p3(rcp->getMomentum());

      //  PID sensitive variables  ***/
      _variables.Update(rcp);
      for(unsigned int i=0; i<_variables.Size(); i++)
      {
        double value = _variables.GetValue(i);
        _treeVars.at(i)->push_back(value);
      }


      // Other variables
      seenP.push_back(p3.Mag());
      seenPt.push_back(p3.Perp());
      seenTheta.push_back(p3.Theta());
      seenPhi.push_back(p3.Phi());
      seenCharge.push_back(rcp->getCharge());

    } // if reco part
    // IMPROVE HERE: CHECK if track or cluster exists!
    else {
      streamlog_out(DEBUG) << "no ReconstructedParticle found for this mcp!" << std::endl ;

      for(unsigned int i=0; i<_variables.Size(); i++)
      {
        double value = _variables.GetValue(i);
        _treeVars.at(i)->push_back(-FLT_MAX);
      }

      seenP.push_back(0);
      seenPt.push_back(0);
      seenTheta.push_back(0);
      seenPhi.push_back(0);
      seenCharge.push_back(0);

    }
  }  // loop over MCPs
  
  
  varTree->Fill();
  nEvt++;
  
  if (nEvt%1000 == 0)
    streamlog_out(MESSAGE) <<  "=================== event " << nEvt << std::endl ;
  
}



void CluShapeTreeProcessor::check( LCEvent * evt ) {

}


void CluShapeTreeProcessor::end(){
  streamlog_out(MESSAGE) << "Found " << nMCPtot << " MC particles.\n";
  streamlog_out(MESSAGE) << "Matched " << nRec << " reconstructed particles.\n";
  streamlog_out(MESSAGE) << "Found " << nTrkCaloMismatch << " cases of track/calo mismatch.\n";
}



