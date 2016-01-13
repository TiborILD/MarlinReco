#include <PIDvarPDF.hh>

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include "TROOT.h"
#include "TString.h"



PIDvarPDF aPIDvarPDF ;


PIDvarPDF::PIDvarPDF() : Processor("PIDvarPDF"),
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

PIDvarPDF::~PIDvarPDF() {
  // Is this necessary, or does marlin magically take care of this?
  for(variable_c_iterator it = pidVars.GetMap()->begin(); it != pidVars.GetMap()->end(); it++) {
    for (particle_c_iterator jt = pidVars.GetParticleMap()->begin(); jt != pidVars.GetParticleMap()->end(); jt++)
    delete sensVarHistos[jt->first][it->first];
  }

}

void PIDvarPDF::init() {

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

  // Sensitive variables used by LikelihoodPID
  for(variable_c_iterator it = pidVars.GetMap()->begin(); it != pidVars.GetMap()->end(); it++) {
    varTree->Branch(it->second.Name(), &sensitiveVars[it->first]);

    for(particle_c_iterator jt  = pidVars.GetParticleMap()->begin();
                            jt != pidVars.GetParticleMap()->end(); jt++)
    {
      sensVarHistos[jt->first][it->first] =
          new TH1F(Form("%s.%s", jt->second.Name(), it->second.Name()),
                   Form("%s.%s;%s;", jt->second.Name(), it->second.Name(), it->second.AxisTitle()),
                          200, it->second.LoLim(), it->second.HiLim());
    }
  }


  static const char *bg = "(x/[5])";
  static const char *bg2 = Form("(%s**2)", bg);
  static const char *b2 = Form("(%s/(1+%s))", bg2, bg2);
  static const char *tmax = Form("[2]*%s", bg2);
  static const char *bbfun = Form("(0.5*[0]*log([1]*%s*%s)-[3]*%s-[4]*%s/2)/%s", bg2, tmax, b2, bg, b2);

  for(particle_c_iterator jt  = pidVars.GetParticleMap()->begin();
                          jt != pidVars.GetParticleMap()->end(); jt++)
  {
    bbFunction[jt->first] = new TF1(Form("%s_BetheBlochFun", jt->second.Name()), bbfun, 0., 300.);
    for (int i=0; i<5; i++)
      bbFunction[jt->first]->FixParameter(i, jt->second.GetBBpars()[i]);

    bbFunction[jt->first]->FixParameter(5, jt->second.mass);

    bbHistos[jt->first] = new TH1F(Form("%s_BetheBloch", jt->second.Name()),
                   Form("%s Bethe-Bloch; p (GeV); dE/dx (MIP)", jt->second.Name()),
                          100, 0., 300.);
    bbHistos[jt->first]->GetListOfFunctions()->Add(bbFunction[jt->first]);
    bbHistos[jt->first]->SetMinimum(0.);
    bbHistos[jt->first]->SetMaximum(10.);

  }
  // Other variables for study
  varTree->Branch("trueP",&trueP) ;
  varTree->Branch("truePt",&truePt) ;
  varTree->Branch("trueTheta",&trueTheta) ;
  varTree->Branch("truePhi",&truePhi) ;
  varTree->Branch("trueCharge",&trueCharge) ;
  varTree->Branch("truePDG",&truePDG) ;
  varTree->Branch("trueMother",&trueMother) ;

  varTree->Branch("isReconstructed",&isReconstructed) ;
  varTree->Branch("isSeen",&isSeen) ;
  varTree->Branch("seenP",&seenP) ;
  varTree->Branch("seenPt",&seenPt) ;
  varTree->Branch("seenTheta",&seenTheta) ;
  varTree->Branch("seenPhi",&seenPhi) ;
  varTree->Branch("seenDEdx",&seenDEdx) ;
  varTree->Branch("seenCharge",&seenCharge) ;

}

void PIDvarPDF::processRunHeader( LCRunHeader* run) {
    
} 

void PIDvarPDF::processEvent( LCEvent * evt ) {

  //streamlog_out(MESSAGE) << " start processing event # " << evt->getEventNumber() << std::endl;

  
  streamlog_out(DEBUG) << " clearing vectors " << std::endl;

  nMCParticles = 0;

  for(variable_c_iterator it = pidVars.GetMap()->begin(); it != pidVars.GetMap()->end(); it++)
    { sensitiveVars[it->first].clear(); }

  trueP.clear();
  truePt.clear();
  trueTheta.clear();
  truePhi.clear();
  trueCharge.clear();
  truePDG.clear();
  trueMother.clear();
  isReconstructed.clear();
  isSeen.clear();
  seenP.clear();
  seenPt.clear();
  seenTheta.clear();
  seenPhi.clear();
  seenDEdx.clear();
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


    if (mcp->getParents().size() > 0) {
  streamlog_out(DEBUG) << " Pushing back true mother " << std::endl;
      trueMother.push_back(mcp->getParents()[0]->getPDG());
    }
    else {
  streamlog_out(DEBUG) << " Pushing back -1 " << std::endl;
      trueMother.push_back(-1);
    }

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

      //  PID sensitive variables  ***/
      pidVars.Update(rcp);
      for(variable_c_iterator it = pidVars.GetMap()->begin(); it != pidVars.GetMap()->end(); it++)
      {
        double value = pidVars.GetVariable(it->first);
        sensitiveVars[it->first].push_back(value);
        // Find which histo to fill, if any, by looking at PDG
        for (particle_c_iterator jt = pidVars.GetParticleMap()->begin(); jt != pidVars.GetParticleMap()->end(); jt++) {
          if(jt->second.pdg == TMath::Abs(mcp->getPDG())) {
            sensVarHistos[jt->first][it->first]->Fill(value);
            break;
          }
        }
      }


      // Other variables
      EVENT::TrackVec trax = rcp->getTracks();
      TVector3 p3(rcp->getMomentum());

      seenP.push_back(p3.Mag());
      seenPt.push_back(p3.Perp());
      seenTheta.push_back(p3.Theta());
      seenPhi.push_back(p3.Phi());
      seenCharge.push_back(rcp->getCharge());

      seenDEdx.push_back(pidVars.GetDEdx());

    } // if reco part
    // IMPROVE HERE: CHECK if track or cluster exists!
    else {
      streamlog_out(DEBUG) << "no ReconstructedParticle found for this mcp!" << std::endl ;

      for(variable_c_iterator it = pidVars.GetMap()->begin(); it != pidVars.GetMap()->end(); it++)
      { sensitiveVars[it->first].push_back(-DBL_MAX); }

      seenP.push_back(0);
      seenPt.push_back(0);
      seenTheta.push_back(0);
      seenPhi.push_back(0);
      seenCharge.push_back(0);
      seenDEdx.push_back(0);

    }
  }  // loop over MCPs
  
  
  varTree->Fill();
  nEvt++;
  
  if (nEvt%1000 == 0)
    streamlog_out(MESSAGE) <<  "=================== event " << nEvt << std::endl ;
  
}



void PIDvarPDF::check( LCEvent * evt ) {

}


void PIDvarPDF::end(){
  streamlog_out(MESSAGE) << "Found " << nMCPtot << " MC particles.\n";
  streamlog_out(MESSAGE) << "Matched " << nRec << " reconstructed particles.\n";
  streamlog_out(MESSAGE) << "Found " << nTrkCaloMismatch << " cases of track/calo mismatch.\n";
}



