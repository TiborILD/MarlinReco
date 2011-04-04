
#include "MaterialDB_F77.hh"
#include "material_structsF77.h"

#include <stdexcept>
#include <sstream>
#include <vector>

#include <marlin/Global.h>

#include "streamlog/streamlog.h"

#include "lcio.h"

#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/BField.h>


using marlin::Global ;

namespace marlin_delphiF77{

  /** helper function to get double values from GEAR parameters - 
      reads values from DoubleVec - if thius does not exist try to read from single double
      value and add n times to vector.
      Needed for backward compatibility of SIT/SET/ETD Mokka drivers ... */
  
  void getDoubleValues( EVENT::DoubleVec& v, const gear::GearParameters& p, 
			const std::string& name, unsigned n){
    
    // first try to read vector:
    try{
      
      const EVENT::DoubleVec& dv = p.getDoubleVals( name )  ;
      
      
      if( n != dv.size() ){ 
	
	std::stringstream em ;
	
	em << " wrong size of array " << name << " - expected: " << n 
	   << " - got: " << dv.size() << std::endl ;
	
	throw gear::Exception( em.str() );
      }
      
      std::copy( dv.begin() , dv.end() , std::back_inserter( v ) ) ; 
      
      
      return ;
      
    } catch(gear::UnknownParameterException& e){}
    
    
    // try to read single double value
    
    double d = p.getDoubleVal( name ) ;
    
    // no catch here - if parameter does not exist we want the exception to be thrown
    
    v.resize( n ) ;
    
    for( unsigned i=0 ; i < n ; ++i ) { v[i] = d ; } 
    
  }
  
// Global static pointer used to ensure a single instance of the class.
  MaterialDB_F77* MaterialDB_F77::_pInstance = NULL; 

  MaterialDB_F77* MaterialDB_F77::Instance() {

    if( ! _pInstance ) {
      _pInstance = new MaterialDB_F77 ;
      _pInstance->checkCommonBlocks() ;
      _pInstance->initialise() ;
    } 

    _pInstance->checkCommonBlocks() ;

    return _pInstance ;
  }

  MaterialDB_F77::~MaterialDB_F77(){
  }

  void MaterialDB_F77::initialise(){
    
    this->buildBeamPipe() ;
    this->buildVXD();
    this->buildSIT();
    this->buildTPC();

    streamlog_out(MESSAGE) << "MaterialDB_F77: Detector building finnished" << std::endl ; 

    switchONMaterial();
    
    // extrapolation surfaces
    fkexts_.nexs = _Nexs;

    // set the magnetic field constants 
    float bField = float(Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z());
    fkfild_.consb = 2.997924e-3*bField;
    coildims_.bfield = bField * 10.0 ; // F77 code uses kGauss
    
    streamlog_out(MESSAGE) << "MaterialDB_F77 " 
			   << " Npmat " << _Npmat << ":" <<  fkddes_.npmat
			   << " Ncmat " << _Ncmat << ":" << fkddes_.ncmat
			   << " Nconmat " << _Nconmat << ":" << fkddes1_.nconmat
			   << " Nplmat " << _Nplmat << ":" << fkddes2_.nplmat
			   << " Nexs " << _Nexs << ":" << fkexts_.nexs 
			   << std::endl;
    
  }


  void MaterialDB_F77::checkCommonBlocks(){
    
    MaterialDB_F77exception exp;  
    // check that nothing has changed in the common blocks since we last checked	

    if( _Nexs != fkexts_.nexs )  {
      streamlog_out(ERROR) << "Miss-match in common blocks: " 
			   << " Nexs = " << _Nexs
			   << " fkexts_.nexs = " << fkexts_.nexs 
			   << std::endl ;
      throw exp ; 
    }

    if( _useMaterials ) {
      if( 
	 fkddes_.ncmat != _Ncmat
	 ||
	 fkddes_.npmat != _Npmat  
	 ||
	 fkddes1_.nconmat != _Nconmat
	 ||
	 fkddes2_.nplmat != _Nplmat  
	  )
	{ 
	  streamlog_out(ERROR) << "Miss-match in common blocks: " 
			       << " Ncmat = " << _Ncmat
			       << " fkddes_.ncmat = " << fkddes_.ncmat
			       << " Npmat = " << _Npmat
			       << " fkddes_.npmat = " << fkddes_.npmat 
			       << " Nconmat = " << _Nconmat
			       << " fkddes1_.nconmat = " << fkddes1_.nconmat 
			       << " Nplmat = " << _Nplmat
			       << " fkddes2_.nplmat = " << fkddes2_.nplmat 
			       << std::endl ;
	  throw exp ; 
	}      
    }
    else{
      if( 
	 fkddes_.ncmat != 0
	 ||
	 fkddes_.npmat != 0 
	 ||
	fkddes1_.nconmat != 0
	 ||
	 fkddes2_.nplmat != 0
	  )      
	{ 
	  throw exp ; 
	}      
    }
  }

  void MaterialDB_F77::switchOFFMaterial(){
    // set the number of material planes to 0
    fkddes_.ncmat = 0 ;
    fkddes_.npmat = 0 ;
    fkddes1_.nconmat = 0 ;
    fkddes2_.nplmat = 0 ;

    _useMaterials = false ;

  }

  void MaterialDB_F77::switchONMaterial(){

    // setting numbers of planar, cyllinder, conical and realistic ladder surfaces
    fkddes_.npmat = _Npmat;
    fkddes_.ncmat = _Ncmat;
    fkddes1_.nconmat = _Nconmat;
    fkddes2_.nplmat = _Nplmat;

    _useMaterials = true ;
    
  }


  bool MaterialDB_F77::buildBeamPipe(){

    // **************************************** //
    // **  Build Database for Beam Pipe   ** //
    // **************************************** //
    streamlog_out(DEBUG) << "build beampipe  ..." << std::endl;
    
    try{
      
      const gear::GearParameters& pBeamPipe = Global::GEAR->getGearParameters("BeamPipe");
      // all values need to be converted from mm to cm
      float beamPipeRadius = 0.1*float(pBeamPipe.getDoubleVal("BeamPipeRadius"));
      float beamPipeHalfZ  = 0.1*float(pBeamPipe.getDoubleVal("BeamPipeHalfZ"));
      float beamPipe_thickness = 0.1*float(pBeamPipe.getDoubleVal("BeamPipeThickness"));
      float beamPipe_radLength = 0.1*float(pBeamPipe.getDoubleVal("BeamPipeProperties_RadLen"));
      //SJA: not sure why this is 10.0*float I guess this is 1/mm
      float beamPipe_dedx = 10.0*float(pBeamPipe.getDoubleVal("BeamPipeProperties_dEdx"));
    
      fkddes_.rcmat[_Ncmat]  =  beamPipeRadius;
      fkddes_.zcmin[_Ncmat]  = -beamPipeHalfZ;
      fkddes_.zcmax[_Ncmat]  =  beamPipeHalfZ;
      fkddes_.xrlc[_Ncmat]   =  beamPipe_thickness / beamPipe_radLength;
      fkddes_.xelosc[_Ncmat] =  beamPipe_thickness * beamPipe_dedx;
    
      _Ncmat++;
    
    
      fkexts_.itexts[_Nexs] = 0;  
      fkexts_.rzsurf[_Nexs] = 0.05*beamPipeRadius ; // SJA?
      fkexts_.zrmin[_Nexs] = -10000.;
      fkexts_.zrmax[_Nexs] = 10000.;
      
      _Nexs++;

    }
    catch(gear::UnknownParameterException){
      streamlog_out(MESSAGE) << "Beampipe not present in Gear file -- Beampipe Material not created" << std::endl;
      return false ;
    }
    return true;
  }

  bool MaterialDB_F77::buildTPC(){
    // **************************************** //
    // ** Build Database for TPC Detector ** //
    // **************************************** //
    streamlog_out(DEBUG) << "build TPC  ..." << std::endl;

    try{
    
      const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;

      // all dimentions need to be converted from mm to cm
        
      float RTPCINN = 0.1*float(gearTPC.getDoubleVal("tpcInnerRadius"));
      float RTPCOUT = 0.1*float(gearTPC.getDoubleVal("tpcOuterRadius"));
      float TPCTHBI = 0.1*float(gearTPC.getDoubleVal("tpcInnerWallThickness"));
      float TPCTHBO = 0.1*float(gearTPC.getDoubleVal("tpcOuterWallThickness"));

      float TPCHLFZ = 0.1 * float( gearTPC.getMaxDriftLength() );

      float xralu     = 0.1*float(gearTPC.getDoubleVal("TPCWallProperties_RadLen"));
      float dedxalu   = 10.*float(gearTPC.getDoubleVal("TPCWallProperties_dEdx"));
      float xrargon   = 0.1*float(gearTPC.getDoubleVal("TPCGasProperties_RadLen"));
      float dedxargon = 10.*float(gearTPC.getDoubleVal("TPCGasProperties_dEdx"));

      // inner tube    
      fkddes_.rcmat[_Ncmat]  =  RTPCINN + TPCTHBI/2.;
      fkddes_.zcmin[_Ncmat]  = -TPCHLFZ;
      fkddes_.zcmax[_Ncmat]  =  TPCHLFZ;
      fkddes_.xrlc[_Ncmat]   =  TPCTHBI/xralu;
      fkddes_.xelosc[_Ncmat] =  TPCTHBI*dedxalu;
    
      fkexts_.itexts[_Nexs]  = 0;
      fkexts_.rzsurf[_Nexs]  = fkddes_.rcmat[_Ncmat]+0.1; //SJA?
      fkexts_.zrmin[_Nexs]   = fkddes_.zcmin[_Ncmat];
      fkexts_.zrmax[_Nexs]   = fkddes_.zcmax[_Ncmat];
    
      _Nexs++;
      _Ncmat++;
    
      int ncyl = 50; // SJA: why 50?
    
      // Gas volume in TPC
      float xstep = (RTPCOUT-RTPCINN-TPCTHBO-TPCTHBI)/float(ncyl);
      for (int icyl=0;icyl<ncyl;++icyl) {
	fkddes_.rcmat[_Ncmat]  =   RTPCINN + TPCTHBI + (float(icyl)+0.5)*xstep;
	fkddes_.zcmin[_Ncmat]  =  -TPCHLFZ;
	fkddes_.zcmax[_Ncmat]  =   TPCHLFZ;
	fkddes_.xrlc[_Ncmat]   =   xstep/xrargon;
	fkddes_.xelosc[_Ncmat] =   xstep*dedxargon;
	_Ncmat++;
      }
    
      // outer tube
      fkddes_.rcmat[_Ncmat]  =   RTPCOUT - TPCTHBO/2.;
      fkddes_.zcmin[_Ncmat]  =  -TPCHLFZ;
      fkddes_.zcmax[_Ncmat]  =   TPCHLFZ;
      fkddes_.xrlc[_Ncmat]   =   TPCTHBO/xralu;
      fkddes_.xelosc[_Ncmat] =   TPCTHBO*dedxalu;
    
      fkexts_.itexts[_Nexs] = 0;
      fkexts_.rzsurf[_Nexs] = fkddes_.rcmat[_Ncmat]-0.1; //SJA ?
      fkexts_.zrmin[_Nexs]  = fkddes_.zcmin[_Ncmat];
      fkexts_.zrmax[_Nexs]  = fkddes_.zcmax[_Ncmat];
  
      _Nexs++;
      _Ncmat++;
    
    }
  
    catch(gear::UnknownParameterException){
      streamlog_out(MESSAGE) << "TPC not present in Gear file -- TPC Material not created" << std::endl;
      return false ;
    }
    return true ;
  }

  bool MaterialDB_F77::buildVXD(){ //SJA:: Extension Surfaces are missing

    streamlog_out(DEBUG) << "building VXD ... " << std::endl;

    bool build_VXD = true;

  try{ 
    const gear::VXDParameters& pVXDDetMain = Global::GEAR->getVXDParameters();
    pVXDDetMain.getVXDLayerLayout();
    Global::GEAR->getGearParameters("VXDInfra");
  }
  catch(gear::UnknownParameterException){
    streamlog_out(MESSAGE) << "VXD not present in Gear file -- VXD Material not created" << std::endl;
    build_VXD = false;
  }

  if( build_VXD ) {

    //--Get GEAR Parameters--
    const gear::VXDParameters& pVXDDetMain = Global::GEAR->getVXDParameters();
    const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();
    const gear::GearParameters& pVXDDet = Global::GEAR->getGearParameters("VXDInfra");
    
    //  Cyllinders and cones in VXD00
    
    //--Get additional parameters from GEAR--
    int nLayersVTX = pVXDLayerLayout.getNLayers();
    int nLadderGaps = int(pVXDDet.getDoubleVals("LadderGaps").size());
    int nStripLines = int(pVXDDet.getDoubleVals("StripLineFinalZ").size());

    float dedx_si = 10.0*float(pVXDDet.getDoubleVal("ActiveLayerProperties_dEdx"));
    float dedx_ber = 10.0*float(pVXDDet.getDoubleVal("SupportLayerProperties_dEdx"));
    
    // Check values

    if (nLadderGaps != nLayersVTX) {
      std::stringstream errorMsg ;
      errorMsg << "MaterialDB_F77 : vector size of LadderGaps vector ("
	       << nLadderGaps << ")  not equal to number of VXD Layers (" << nLayersVTX << ")" << std::endl;
      throw gear::Exception(errorMsg.str());
    }
    if (nStripLines != nLayersVTX) {
      std::stringstream errorMsg ; 
      errorMsg << "MaterialDB_F77 : vector size of StripLines vector ("
	       << nStripLines << ")  not equal to number of VXD Layers (" << nLayersVTX << ")" << std::endl;
      throw gear::Exception(errorMsg.str());
    }
    if (pVXDDetMain.getShellOuterRadius()<pVXDDetMain.getShellInnerRadius()) {
      std::stringstream errorMsg ;
      errorMsg << "Outer Shell Radius (" << pVXDDetMain.getShellOuterRadius() 
	       << ") is smaller than inner Shell radius (" << pVXDDetMain.getShellInnerRadius()
	       << ")" << std::endl;
      throw gear::Exception(errorMsg.str());
    }

    float shell_thickness = float(pVXDDetMain.getShellOuterRadius()-pVXDDetMain.getShellInnerRadius());
    float VTXShell_thickness = shell_thickness;
    float VTXShell_Radius = float(pVXDDetMain.getShellInnerRadius());
    float VTXShell_HalfZ = float(pVXDDetMain.getShellHalfLength());
    float VTXEndPlate_innerRadius = float(pVXDDet.getDoubleVal("VXDEndPlateInnerRadius"));
    float VTXShell_radLength = float(pVXDDetMain.getShellRadLength());
    
    //--The Ladder structure (realistic ladder)--
    int nLadders;
    float Pi = acos(-1);
    //  float deg2rad = Pi / 180.0;
    
    for (int i=0; i<nLayersVTX; ++i) {
      nLadders = pVXDLayerLayout.getNLadders(i);
      
      float ladder_phi0 = float(pVXDLayerLayout.getPhi0(i));
      float ladder_distance = float(pVXDLayerLayout.getLadderDistance(i));
      float ladder_thickness = float(pVXDLayerLayout.getLadderThickness(i));
      float ladder_width = float(pVXDLayerLayout.getLadderWidth(i));
      float ladder_length = float (pVXDLayerLayout.getLadderLength(i));
      float ladder_offset = float (pVXDLayerLayout.getLadderOffset(i));
      float ladder_radLength = 0.1*float(pVXDLayerLayout.getLadderRadLength(i));
      
      float sensitive_distance = float(pVXDLayerLayout.getSensitiveDistance(i));
      float sensitive_thickness = float(pVXDLayerLayout.getSensitiveThickness(i));
      float sensitive_width = float(pVXDLayerLayout.getSensitiveWidth(i));
      float sensitive_length = float(pVXDLayerLayout.getSensitiveLength(i));
      float sensitive_offset = float (pVXDLayerLayout.getSensitiveOffset(i));
      float sensitive_radLength = 0.1*float(pVXDLayerLayout.getSensitiveRadLength(i));

      float halfLadderGaps = float(pVXDDet.getDoubleVals("LadderGaps")[i]);

      //--Realistic ladder structure--
      
      float currPhi;
      float angleLadders = 2*Pi / nLadders;
      float cosphi, sinphi;
      
      ladder_distance += 0.5* ladder_thickness;
      sensitive_distance +=0.5* sensitive_thickness;
	
      for (int j=0; j<nLadders; ++j) {
	
	currPhi = ladder_phi0 + (angleLadders * j);
	cosphi = cos(currPhi);
	sinphi = sin(currPhi);
	
	//Beryllium support
	fkddes2_.xplmat[_Nplmat] = 0.1*(ladder_distance*cosphi - ladder_offset*sinphi);
	fkddes2_.yplmat[_Nplmat] = 0.1*(ladder_distance*sinphi + ladder_offset*cosphi);
	fkddes2_.zplmat[_Nplmat] = 0.0; //Ladder is centered around z=0
	
	fkddes2_.widplmat[_Nplmat] = 0.1*(ladder_width); 
	fkddes2_.lenplmat[_Nplmat] = 0.1*(2.*(ladder_length+halfLadderGaps));
	fkddes2_.phiplmat[_Nplmat] = currPhi;

	fkddes2_.xrlpl[_Nplmat] = 0.1*ladder_thickness/ladder_radLength;
	fkddes2_.xelospl[_Nplmat] = 0.1*ladder_thickness*dedx_ber;
	
	_Nplmat++;
      
	
	//Sensitive Si part
	fkddes2_.xplmat[_Nplmat] = 0.1*(sensitive_distance*cosphi - sensitive_offset*sinphi);
	fkddes2_.yplmat[_Nplmat] = 0.1*(sensitive_distance*sinphi + sensitive_offset*cosphi);
	fkddes2_.zplmat[_Nplmat] = 0.0; //Ladder is centered around z=0
	
	fkddes2_.widplmat[_Nplmat] = 0.1*(sensitive_width); 
	fkddes2_.lenplmat[_Nplmat] = 0.1*(2.*(sensitive_length+halfLadderGaps));
	fkddes2_.phiplmat[_Nplmat] = currPhi;
	
	fkddes2_.xrlpl[_Nplmat] = 0.1*sensitive_thickness/sensitive_radLength;
	fkddes2_.xelospl[_Nplmat] = 0.1*sensitive_thickness*dedx_si;
	
	_Nplmat++;
      }
      
      //Extension Surface (used by the Kalman filter); Radius is set to the distance from IP to the center of the support
      fkexts_.itexts[_Nexs] = 0;
      fkexts_.rzsurf[_Nexs] = 0.1*ladder_distance;
      fkexts_.zrmin[_Nexs] = - 0.1*ladder_length;
      fkexts_.zrmax[_Nexs] = 0.1*ladder_length;
      _Nexs++;	
    
    }

    //--Cryostat--
    float AlRadius = float(pVXDDet.getDoubleVal("CryostatAlRadius"));
    float AlHalfLength = float(pVXDDet.getDoubleVal("CryostatAlHalfZ"));
    float AlThickness = float(pVXDDet.getDoubleVal("CryostatAlThickness"));
    float AlZEndCap = float(pVXDDet.getDoubleVal("CryostatAlZEndCap"));
    float AlRinEndCap = float(pVXDDet.getDoubleVal("CryostatAlInnerR"));
    float xrad_cryo = 0.1*float(pVXDDet.getDoubleVal("Cryostat_RadLen"));
    float dedx_cryo = 10.0*float(pVXDDet.getDoubleVal("Cryostat_dEdx"));
    
    // Al cryostat barrel
    fkddes_.rcmat[_Ncmat] = 0.1*(AlRadius+0.5*AlThickness);
    fkddes_.zcmin[_Ncmat] = -0.1*AlHalfLength;
    fkddes_.zcmax[_Ncmat] = 0.1*AlHalfLength;
    fkddes_.xrlc[_Ncmat] = 0.1*AlThickness/xrad_cryo;
    fkddes_.xelosc[_Ncmat] = 0.1*AlThickness*dedx_cryo;
    _Ncmat++;
  

    // Al cryostat endcaps
    fkddes_.zpmat[_Npmat] = -0.1*(AlZEndCap+0.5*AlThickness);
    fkddes_.rpmin[_Npmat] = 0.1*AlRinEndCap;
    fkddes_.rpmax[_Npmat] = 0.1*(AlRadius+AlThickness);;
    fkddes_.xrlp[_Npmat] = 0.1*AlThickness/xrad_cryo;
    fkddes_.xelosp[_Npmat] = 0.1*AlThickness*dedx_cryo;
    _Npmat++;

    fkddes_.zpmat[_Npmat] = 0.1*(AlZEndCap+0.5*AlThickness);
    fkddes_.rpmin[_Npmat] = 0.1*AlRinEndCap;
    fkddes_.rpmax[_Npmat] = 0.1*(AlRadius+AlThickness);;
    fkddes_.xrlp[_Npmat] = 0.1*AlThickness/xrad_cryo;
    fkddes_.xelosp[_Npmat] = 0.1*AlThickness*dedx_cryo;
    _Npmat++;
    

    //  Outer support cyllinder for VTX
    fkddes_.rcmat[_Ncmat] = 0.1*(VTXShell_Radius+0.5*VTXShell_thickness);
    fkddes_.zcmin[_Ncmat] = -0.1*VTXShell_HalfZ;
    fkddes_.zcmax[_Ncmat] = 0.1*VTXShell_HalfZ;
    fkddes_.xrlc[_Ncmat] = 0.1*VTXShell_thickness/VTXShell_radLength;
    fkddes_.xelosc[_Ncmat] = 0.1*VTXShell_thickness*dedx_ber;
    _Ncmat++;


    //  EndPlate support disk for VTX ; left part
    fkddes_.zpmat[_Npmat] = -0.1*(VTXShell_HalfZ+0.5*VTXShell_thickness);
    fkddes_.rpmin[_Npmat] = 0.1*VTXEndPlate_innerRadius;
    fkddes_.rpmax[_Npmat] = 0.1*(VTXShell_Radius+VTXShell_thickness);
    fkddes_.xrlp[_Npmat] = 0.1*VTXShell_thickness/VTXShell_radLength;
    fkddes_.xelosp[_Npmat] = 0.1*VTXShell_thickness*dedx_ber;  
    _Npmat++;    


    //  EndPlate support disk for VTX ; right part
    fkddes_.zpmat[_Npmat] = 0.1*(VTXShell_HalfZ+0.5*VTXShell_thickness);
    fkddes_.rpmin[_Npmat] = 0.1*VTXEndPlate_innerRadius;
    fkddes_.rpmax[_Npmat] = 0.1*(VTXShell_Radius+VTXShell_thickness);
    fkddes_.xrlp[_Npmat] = 0.1*VTXShell_thickness/VTXShell_radLength;
    fkddes_.xelosp[_Npmat] = 0.1*VTXShell_thickness*dedx_ber;  
    _Npmat++;    

    return build_VXD ;

  }

  // *********************************************** //
  // ** Build Database for VXD_Simple Detector ** //
  // *********************************************** //
  streamlog_out(DEBUG) << "build VXD_Simple ..." << std::endl;
  
  bool build_VXD_Simple = true;

  try{ 
    Global::GEAR->getGearParameters("VXD_Simple");
  }
  catch(gear::UnknownParameterException){
    streamlog_out(MESSAGE) << "VXD_Simple not present in Gear file -- VXD_Simple Material not created" << std::endl;
    build_VXD_Simple = false;
  }
  
  if( build_VXD_Simple ) 
    {
      
      const gear::GearParameters& theVXD_Simple = Global::GEAR->getGearParameters("VXD_Simple");
      
      int nLayers = int(theVXD_Simple.getDoubleVals("VXDSensitiveLayerInnerRadius").size());
      
      float sensitiveThickness = float(theVXD_Simple.getDoubleVal("VXDSensitiveLayerThickness")); 
      float supportThickness   = float(theVXD_Simple.getDoubleVal("VXDSupportLayerThickness")); 
      
      EVENT::DoubleVec SensitiveLayerHalfLengths ;
      getDoubleValues(SensitiveLayerHalfLengths, theVXD_Simple, "VXDLayerHalfLength",  nLayers);
      
      EVENT::DoubleVec SensitiveLayerInnerRadii ;
      getDoubleValues(SensitiveLayerInnerRadii, theVXD_Simple, "VXDSensitiveLayerInnerRadius",  nLayers);
      
      float sensitiveLayerRadLen    = float(theVXD_Simple.getDoubleVal("VXDSensitiveLayer_RadLen"));
      float supportLayerRadLen      = float(theVXD_Simple.getDoubleVal("VXDSupportLayer_RadLen"));
      float sensitiveLayerdEdx      = float(theVXD_Simple.getDoubleVal("VXDSensitiveLayer_dEdx"));
      float supportLayerdEdx        = float(theVXD_Simple.getDoubleVal("VXDSupportLayer_dEdx"));
      
      for (int iL=0;iL<nLayers;++iL) {
	// sensitive layers
	float radius_sen = float(SensitiveLayerInnerRadii[iL]) + 0.5 * sensitiveThickness;
	float halfz      = float(SensitiveLayerHalfLengths[iL]);
	fkddes_.rcmat[_Ncmat] = 0.1*radius_sen; // convert to cm
	fkddes_.zcmin[_Ncmat] =  -0.1*halfz; // convert to cm
	fkddes_.zcmax[_Ncmat] = 0.1*halfz; // convert to cm
	fkddes_.xrlc[_Ncmat] = sensitiveThickness/sensitiveLayerRadLen;
	fkddes_.xelosc[_Ncmat] = sensitiveThickness*sensitiveLayerdEdx;      
	_Ncmat++;
	
	// support layers
	float radius_sup = radius_sen + 0.5 * sensitiveThickness + 0.5 * supportThickness ;
	
	fkddes_.rcmat[_Ncmat] = 0.1*radius_sup; // convert to cm
	fkddes_.zcmin[_Ncmat] = -0.1*halfz; // convert to cm
	fkddes_.zcmax[_Ncmat] = 0.1*halfz;  // convert to cm
	fkddes_.xrlc[_Ncmat] = supportThickness/supportLayerRadLen;
	fkddes_.xelosc[_Ncmat] = supportThickness*supportLayerdEdx;
	_Ncmat++;

	fkexts_.itexts[_Nexs] =  0;
	fkexts_.rzsurf[_Nexs] =  0.1*radius_sen-0.1; //SJA ?
	fkexts_.zrmin[_Nexs]  = -0.1*halfz;
	fkexts_.zrmax[_Nexs]  =  0.1*halfz;
  
	_Nexs++;

	
      }
    }
  return build_VXD_Simple;
  }



  bool MaterialDB_F77::buildSIT(){ //SJA:: Extension Surfaces are missing

    streamlog_out(DEBUG) << "build SIT  ..." << std::endl;
    bool build_SIT = true;

    try{
      Global::GEAR->getGearParameters("SIT");      
    }

    catch(gear::UnknownParameterException){
      streamlog_out(MESSAGE) << "SIT not present in Gear file -- SIT Material not created" << std::endl;
      build_SIT = false;
    }

    if(build_SIT){
      // SIT layers
      const gear::GearParameters& pSITDet = Global::GEAR->getGearParameters("SIT");

      std::vector<float> rSIT;
      std::vector<float> halfZSIT;
      std::vector<float> rSITSupport;
      std::vector<float> halfZSITSupport;

      std::vector<double> SITLayer_thickness ;
      std::vector<double> SITLayerSupport_thickness ;

      int nSITR = int(pSITDet.getDoubleVals("SITLayerRadius").size());
      int nSITHL = int(pSITDet.getDoubleVals("SITLayerHalfLength").size());
      //  int SITModel = int(pSITDet.getIntVal("SITModel"));
      int nLayersSIT = 0;
      
      if (nSITR == nSITHL) {
	nLayersSIT = nSITR;
	rSIT.resize(nLayersSIT);
	halfZSIT.resize(nLayersSIT);
	rSITSupport.resize(nLayersSIT);
	halfZSITSupport.resize(nLayersSIT);
      }
      else {
	std::stringstream errorMsg;
	errorMsg << "Size of SITLayerRadius vector (" << nSITR 
		  << ") is not equal to the size of SITHalfLength vector ("
		  << nSITHL << ")" << std::endl;
	throw gear::Exception(errorMsg.str());
      }
      
      
      getDoubleValues(SITLayer_thickness,        pSITDet, "SITLayerThickness", nLayersSIT );
      getDoubleValues(SITLayerSupport_thickness, pSITDet, "SITSupportLayerThickness", nLayersSIT );
      
      for (int iL=0;iL<nLayersSIT;++iL) {
	rSIT[iL] = float(pSITDet.getDoubleVals("SITLayerRadius")[iL]);
	halfZSIT[iL] = float(pSITDet.getDoubleVals("SITLayerHalfLength")[iL]);
	rSITSupport[iL] = float(pSITDet.getDoubleVals("SITSupportLayerRadius")[iL]);
	halfZSITSupport[iL] = float(pSITDet.getDoubleVals("SITSupportLayerHalfLength")[iL]);
      }
      
      float radlen_si = 0.1*float(pSITDet.getDoubleVal("SITLayer_RadLen"));
      float dedx_si = 10.*float(pSITDet.getDoubleVal("SITLayer_dEdx"));
      
      
      float radlen_ber = 0.1*float(pSITDet.getDoubleVal("SITSupportLayer_RadLen"));
      float dedx_ber = 10.*float(pSITDet.getDoubleVal("SITSupportLayer_dEdx"));
      
      
      for (int iL = 0; iL < nLayersSIT; ++iL) {
	
		
	fkddes_.rcmat[_Ncmat] = 0.1*rSIT[iL];
	fkddes_.zcmin[_Ncmat] = -0.1*halfZSIT[iL];
	fkddes_.zcmax[_Ncmat] = 0.1*halfZSIT[iL];
	fkddes_.xrlc[_Ncmat] = 0.1*SITLayer_thickness[iL]/radlen_si;
	fkddes_.xelosc[_Ncmat] = 0.1*SITLayer_thickness[iL]*dedx_si;     
	++_Ncmat;
	
	fkddes_.rcmat[_Ncmat] = 0.1*rSITSupport[iL];
	fkddes_.zcmin[_Ncmat] = -0.1*halfZSITSupport[iL];
	fkddes_.zcmax[_Ncmat] = 0.1*halfZSITSupport[iL];
	fkddes_.xrlc[_Ncmat] = 0.1*SITLayerSupport_thickness[iL]/radlen_ber;
	fkddes_.xelosc[_Ncmat] = 0.1*SITLayerSupport_thickness[iL]*dedx_ber;
	++_Ncmat;
	
    
	fkexts_.itexts[_Nexs] = 0;
	fkexts_.rzsurf[_Nexs] = 0.1*rSIT[iL];
	fkexts_.zrmin[_Nexs]  = -halfZSIT[iL];
	fkexts_.zrmax[_Nexs]  = halfZSIT[iL];
	
	++_Nexs;	      	      
      } 
      return build_SIT;
    }    
    
  


    // *********************************************** //
    // ** Build Database for SIT_Simple Detector ** //
    // *********************************************** //
    streamlog_out(DEBUG) << "build SIT_Simple ..." << std::endl;
    
    bool build_SIT_Simple = true;
    try{ 
      Global::GEAR->getGearParameters("SIT_Simple");
    }
    catch(gear::UnknownParameterException){
      streamlog_out(MESSAGE) << "SIT_Simple not present in Gear file -- SIT_Simple Material not created" << std::endl;
      build_SIT_Simple = false;
    }
    
    if( build_SIT_Simple ) 
      {
	
	const gear::GearParameters& theSIT_Simple = Global::GEAR->getGearParameters("SIT_Simple");
	
	int nLayers = int(theSIT_Simple.getDoubleVals("SITSensitiveLayerInnerRadius").size());
	
	float sensitiveThickness = float(theSIT_Simple.getDoubleVal("SITSensitiveLayerThickness")); 
	float supportThickness   = float(theSIT_Simple.getDoubleVal("SITSupportLayerThickness")); 
	
	EVENT::DoubleVec SensitiveLayerHalfLengths ;
	getDoubleValues(SensitiveLayerHalfLengths, theSIT_Simple, "SITLayerHalfLength",  nLayers);
	
	EVENT::DoubleVec SensitiveLayerInnerRadii ;
	getDoubleValues(SensitiveLayerInnerRadii, theSIT_Simple, "SITSensitiveLayerInnerRadius",  nLayers);
	
	float sensitiveLayerRadLen    = float(theSIT_Simple.getDoubleVal("SITSensitiveLayer_RadLen"));
	float supportLayerRadLen      = float(theSIT_Simple.getDoubleVal("SITSupportLayer_RadLen"));
	float sensitiveLayerdEdx      = float(theSIT_Simple.getDoubleVal("SITSensitiveLayer_dEdx"));
	float supportLayerdEdx        = float(theSIT_Simple.getDoubleVal("SITSupportLayer_dEdx"));
	
	for (int iL=0;iL<nLayers;++iL) {
	  // sensitive layers
	  float radius_sen = float(SensitiveLayerInnerRadii[iL]) + 0.5 * sensitiveThickness;
	  float halfz      = float(SensitiveLayerHalfLengths[iL]);
	  fkddes_.rcmat[_Ncmat] = 0.1*radius_sen; // convert to cm
	  fkddes_.zcmin[_Ncmat] =  -0.1*halfz; // convert to cm
	  fkddes_.zcmax[_Ncmat] = 0.1*halfz; // convert to cm
	  fkddes_.xrlc[_Ncmat] = sensitiveThickness/sensitiveLayerRadLen;
	  fkddes_.xelosc[_Ncmat] = sensitiveThickness*sensitiveLayerdEdx;      
	  _Ncmat++;
	  
	  // support layers
	  float radius_sup = radius_sen + 0.5 * sensitiveThickness + 0.5 * supportThickness ;
	  fkddes_.rcmat[_Ncmat] = 0.1*radius_sup; // convert to cm
	  fkddes_.zcmin[_Ncmat] = -0.1*halfz; // convert to cm
	  fkddes_.zcmax[_Ncmat] = 0.1*halfz;  // convert to cm
	  fkddes_.xrlc[_Ncmat] = supportThickness/supportLayerRadLen;
	  fkddes_.xelosc[_Ncmat] = supportThickness*supportLayerdEdx;
	  _Ncmat++;

	  fkexts_.itexts[_Nexs] =  0;
	  fkexts_.rzsurf[_Nexs] =  0.1*radius_sen-0.1; //SJA ?
	  fkexts_.zrmin[_Nexs]  = -0.1*halfz;
	  fkexts_.zrmax[_Nexs]  =  0.1*halfz;
	  
	  _Nexs++;
	
	}	
      }
    return build_SIT_Simple;
  }


} // end of marlin_delphiF77 namespace 
