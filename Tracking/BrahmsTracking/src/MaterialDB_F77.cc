
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
    this->buildSET();
    this->buildTPC();
    this->buildFTD();

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
      fkexts_.rzsurf[_Nexs] = beamPipeRadius - 0.1 ; // place an extrapolation surface just inside the beampipe
      fkexts_.zrmin[_Nexs] = -10000.;
      fkexts_.zrmax[_Nexs] = 10000.;
      
      _Nexs++;
        
      fkexts_.itexts[_Nexs] = 0;  
      fkexts_.rzsurf[_Nexs] = 0.05*beamPipeRadius ; // place an extrapolation surface close to the IP
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
      const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
      const gear::Vector2D padCoord = padLayout.getPadCenter(1) ;

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
    
//      fkexts_.itexts[_Nexs]  = 0;
//      fkexts_.rzsurf[_Nexs]  = fkddes_.rcmat[_Ncmat]+0.1; // place extrapolation surface just on the outside of the inner field cage material 
//      fkexts_.zrmin[_Nexs]   = fkddes_.zcmin[_Ncmat];
//      fkexts_.zrmax[_Nexs]   = fkddes_.zcmax[_Ncmat];
      _Ncmat++;

      fkexts_.itexts[_Nexs]  = 0;
      fkexts_.rzsurf[_Nexs]  =  RTPCINN ; // place extrapolation surface one inner side (wrt IP) of the inner field cage material
      fkexts_.zrmin[_Nexs]   = -TPCHLFZ;
      fkexts_.zrmax[_Nexs]   =  TPCHLFZ;
      _Nexs++;

      // put another extrapolation surface on the gas side of the field cage
      fkexts_.itexts[_Nexs]  =  0;
      fkexts_.rzsurf[_Nexs]  =  RTPCINN + TPCTHBI ; // place extrapolation surface one outer side (wrt IP) of the inner field cage material
      fkexts_.zrmin[_Nexs]   = -TPCHLFZ;
      fkexts_.zrmax[_Nexs]   =  TPCHLFZ;
      _Nexs++;
    
      // put another extrapolation surface at the center of the first pad
      fkexts_.itexts[_Nexs]  =  0;
      fkexts_.rzsurf[_Nexs]  =  0.1 * padLayout.getPlaneExtent()[0];
      fkexts_.zrmin[_Nexs]   = -TPCHLFZ;
      fkexts_.zrmax[_Nexs]   =  TPCHLFZ;
      _Nexs++;

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

      _Ncmat++;
  
//      fkexts_.itexts[_Nexs] = 0;
//      fkexts_.rzsurf[_Nexs] = fkddes_.rcmat[_Ncmat]-0.1; // place extrapolation surface just on the inside of the outer field cage material 
//      fkexts_.zrmin[_Nexs]  = fkddes_.zcmin[_Ncmat];
//      fkexts_.zrmax[_Nexs]  = fkddes_.zcmax[_Ncmat];

      // put another extrapolation surface at the inner radius of the pad plane
      fkexts_.itexts[_Nexs]  =  0;
      fkexts_.rzsurf[_Nexs]  =  0.1 * padLayout.getPlaneExtent()[1];
      fkexts_.zrmin[_Nexs]   = -TPCHLFZ;
      fkexts_.zrmax[_Nexs]   =  TPCHLFZ;
      _Nexs++;
  
      fkexts_.itexts[_Nexs]  = 0;
      fkexts_.rzsurf[_Nexs]  =  RTPCOUT - TPCTHBO; // place extrapolation surface one inner side (wrt IP) of the outer field cage material
      fkexts_.zrmin[_Nexs]   = -TPCHLFZ;
      fkexts_.zrmax[_Nexs]   =  TPCHLFZ;
      _Nexs++;
      
      // put another extrapolation surface on the gas side of the field cage
      fkexts_.itexts[_Nexs]  =  0;
      fkexts_.rzsurf[_Nexs]  =  RTPCOUT; // place extrapolation surface one outer side (wrt IP) of the outer field cage material
      fkexts_.zrmin[_Nexs]   = -TPCHLFZ;
      fkexts_.zrmax[_Nexs]   =  TPCHLFZ;
      _Nexs++;
    
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
	fkexts_.rzsurf[_Nexs] =  0.1*radius_sen-0.1; // place an extrapolation surface just inside the inner edge sensitive material 
	fkexts_.zrmin[_Nexs]  = -0.1*halfz;
	fkexts_.zrmax[_Nexs]  =  0.1*halfz;
  
	_Nexs++;

	
      }
      return build_VXD_Simple;
    }
  return false ; // if this point is reached nothing was sucessfully built 
  }

  bool MaterialDB_F77::buildFTD(){

    // ************************************* //
    // ** Build Database for FTD_Detector ** //
    // ************************************* //
    streamlog_out(DEBUG) << "build FTD ..." << std::endl;
    
    bool build_FTD = true;

    try{ 
      Global::GEAR->getGearParameters("FTD");
    }
    catch(gear::UnknownParameterException){
      streamlog_out(MESSAGE) << "FTD not present in Gear file -- FTD Material not created" << std::endl;
      build_FTD = false;
    }
  
    if( build_FTD ) 
      {
      
	const gear::GearParameters& theFTD = Global::GEAR->getGearParameters("FTD");

	std::vector<float> zFTD ;
	std::vector<float> rInFTD ;
	std::vector<float> rOutFTD ;
	std::vector<float> dzSiFTD ;
	std::vector<float> dzSupportFTD ;
	
	// Check which version of the FTD this is:
	// i)  Support Rings
	// ii) Support Disks
      
	// Planar detectors in FTD

	bool has_si872 (false);
	bool has_SupportDisks (false);

	try{
	  int nFTDSupport_dZ = int(theFTD.getDoubleVals("FTDDiskSupportThickness").size());
	  nFTDSupport_dZ = 0; // avoid compiler warning 
	  has_SupportDisks = true;
	}
	catch(gear::UnknownParameterException &e){}

	try{
	  float dedx_si872 = float(theFTD.getDoubleVal("Silicon872_dEdx"));
	  dedx_si872 = 0.; // avoid compiler warning 
	  has_si872 = true;
	}
	catch(gear::UnknownParameterException &e){}


      
	// Planar detectors in FTD

	if(has_SupportDisks==true && has_si872==false) {
 
	  streamlog_out(DEBUG) << "build FTD according to the SupportDisk structure" << std::endl;

	  int nLayersFTD = 0;
	  int nFTDZ = int(theFTD.getDoubleVals("FTDZCoordinate").size());
	  int nFTDRin = int(theFTD.getDoubleVals("FTDInnerRadius").size());
	  int nFTDRout = int(theFTD.getDoubleVals("FTDOuterRadius").size());
	  int nFTDSi_dZ = int(theFTD.getDoubleVals("FTDDiskSiThickness").size());
	  int nFTDSupport_dZ = int(theFTD.getDoubleVals("FTDDiskSupportThickness").size());
  

	  if (nFTDZ == nFTDRin && nFTDRin == nFTDRout && nFTDRout == nFTDSi_dZ && nFTDSi_dZ == nFTDSupport_dZ) {
	    nLayersFTD = nFTDZ;
	    zFTD.resize(nLayersFTD);
	    rInFTD.resize(nLayersFTD);
	    rOutFTD.resize(nLayersFTD);
	    dzSiFTD.resize(nLayersFTD);
	    dzSupportFTD.resize(nLayersFTD);
	  }
	  else {
	    std::stringstream errorMsg ;
	    errorMsg << "Size of vectors FTDZCoordinate, FTDInnerRadius, FTDInnerRadius, FTDDiskSiThickness and FTDDiskSupportThickness are not equal --->" 
		     << " # FTDZCoordinate : " << nFTDZ
		     << " # FTDInnerRadius : " << nFTDRin
		     << " # FTDOuterRadius : " << nFTDRout 
		     << " # FTDDiskSiThickness : " << nFTDSi_dZ
		     << " # FTDDiskSupportThickness : " << nFTDSupport_dZ
		     << std::endl ;
	    throw gear::Exception(errorMsg.str());
	  }
	  for (int i=0;i<nLayersFTD;++i) {
	    zFTD[i] = float(theFTD.getDoubleVals("FTDZCoordinate")[i]);
	    rInFTD[i] = float(theFTD.getDoubleVals("FTDInnerRadius")[i]);
	    rOutFTD[i] = float(theFTD.getDoubleVals("FTDOuterRadius")[i]);
	    dzSiFTD[i] = float(theFTD.getDoubleVals("FTDDiskSiThickness")[i]);
	    dzSupportFTD[i] = float(theFTD.getDoubleVals("FTDDiskSupportThickness")[i]);
	  } 
	  
	  float zFTDOuterCyllinderStart = float(theFTD.getDoubleVal("zFTDOuterCylinderStart"));
	  float zFTDOuterCyllinderEnd = float(theFTD.getDoubleVal("zFTDOuterCylinderEnd"));
//	  float zFTDInnerConeStart = float(theFTD.getDoubleVal("zFTDInnerConeStart"));
//	  float zFTDInnerConeEnd = float(theFTD.getDoubleVal("zFTDInnerConeEnd"));
	  float FTD_copper_thickness = float(theFTD.getDoubleVal("FTDCopperThickness"));
	  float FTD_kaptonCyl_thickness = float(theFTD.getDoubleVal("FTDOuterCylinderThickness"));
	
	  float dedx_si = 10.0*float(theFTD.getDoubleVal("Silicon_dEdx"));
	  float dedx_kapton = 10.0*float(theFTD.getDoubleVal("Kapton_dEdx"));
	  float dedx_copper = 10.0*float(theFTD.getDoubleVal("Copper_dEdx"));
	  float radlen_si = 0.1*float(theFTD.getDoubleVal("Silicon_RadLen"));  
	  float radlen_kapton = 0.1*float(theFTD.getDoubleVal("Kapton_RadLen"));
	  float radlen_copper = 0.1*float(theFTD.getDoubleVal("Copper_RadLen"));
	  
	  for (int i=0;i<nLayersFTD;++i) {
	    // FTD Si Disks
	    float dedx = dedx_si;
	    float radlen  = radlen_si;
	    
	    // right-hand part
	    fkddes_.zpmat[_Npmat] = 0.1*zFTD[i];
	    fkddes_.rpmin[_Npmat] = 0.1*rInFTD[i];
	    fkddes_.rpmax[_Npmat] = 0.1*rOutFTD[i];
	    fkddes_.xrlp[_Npmat] = 0.1*dzSiFTD[i]/radlen;
	    fkddes_.xelosp[_Npmat] = 0.1*dzSiFTD[i]*dedx;
	    
	    fkexts_.itexts[_Nexs] = 1;
	    fkexts_.rzsurf[_Nexs] = fkddes_.zpmat[_Npmat];
	    fkexts_.zrmin[_Nexs] = fkddes_.rpmin[_Npmat];
	    fkexts_.zrmax[_Nexs] = fkddes_.rpmax[_Npmat];
	    
	    _Nexs++;
	    _Npmat++;
	    
	    // left-hand part
	    fkddes_.zpmat[_Npmat] = -0.1*zFTD[i];
	    fkddes_.rpmin[_Npmat] = 0.1*rInFTD[i];
	    fkddes_.rpmax[_Npmat] = 0.1*rOutFTD[i];
	    fkddes_.xrlp[_Npmat] = 0.1*dzSiFTD[i]/radlen;
	    fkddes_.xelosp[_Npmat] = 0.1*dzSiFTD[i]*dedx;
	    
	    fkexts_.itexts[_Nexs] = 1;
	    fkexts_.rzsurf[_Nexs] = fkddes_.zpmat[_Npmat];
	    fkexts_.zrmin[_Nexs] = fkddes_.rpmin[_Npmat];
	    fkexts_.zrmax[_Nexs] = fkddes_.rpmax[_Npmat];
	    
	    _Nexs++;
	    _Npmat++;
	    
	    // Support Disks
	    // right-hand part
	    fkddes_.zpmat[_Npmat] = 0.1*(zFTD[i]+(dzSiFTD[i]/2.0)+(dzSupportFTD[i]/2.0));
	    fkddes_.rpmin[_Npmat] = 0.1*rInFTD[i];
	    fkddes_.rpmax[_Npmat] = 0.1*rOutFTD[i];
	    fkddes_.xrlp[_Npmat] = 0.1*(dzSupportFTD[i]/radlen_kapton);
	    fkddes_.xelosp[_Npmat] = 0.1*(dzSupportFTD[i]*dedx_kapton);
	    _Npmat++;
	
	    // left-hand part
	    fkddes_.zpmat[_Npmat] = -0.1*(zFTD[i]+(dzSiFTD[i]/2.0)+(dzSupportFTD[i]/2.0));
	    fkddes_.rpmin[_Npmat] = 0.1*rInFTD[i];
	    fkddes_.rpmax[_Npmat] = 0.1*rOutFTD[i];
	    fkddes_.xrlp[_Npmat] = 0.1*dzSupportFTD[i]/radlen_kapton;
	    fkddes_.xelosp[_Npmat] = 0.1*dzSupportFTD[i]*dedx_kapton;
	    _Npmat++;
	    
	  }
	  // Outer support cyllinders (FTD)
	  
	    // copper cables; left part
	  fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+0.5+
				       0.5*FTD_copper_thickness);
	  fkddes_.zcmin[_Ncmat] = -0.1*zFTDOuterCyllinderEnd;
	  fkddes_.zcmax[_Ncmat] = -0.1*zFTDOuterCyllinderStart;
	  fkddes_.xrlc[_Ncmat] = 0.1*FTD_copper_thickness/radlen_copper;
	  fkddes_.xelosc[_Ncmat] = 0.1*FTD_copper_thickness*dedx_copper;
	  _Ncmat++;
	  // copper cables; right part
	  fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+0.5+ 
				       0.5*FTD_copper_thickness);
	  fkddes_.zcmin[_Ncmat] = 0.1*zFTDOuterCyllinderStart;
	  fkddes_.zcmax[_Ncmat] = 0.1*zFTDOuterCyllinderEnd;
	  fkddes_.xrlc[_Ncmat] = 0.1*FTD_copper_thickness/radlen_copper;
	  fkddes_.xelosc[_Ncmat] = 0.1*FTD_copper_thickness*dedx_copper;
	  _Ncmat++;  
	  // kapton cyllinder; left part
	  fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+0.5+
				       FTD_copper_thickness+
				       0.5*FTD_kaptonCyl_thickness);
	  fkddes_.zcmin[_Ncmat] = -0.1*zFTDOuterCyllinderEnd;
	  fkddes_.zcmax[_Ncmat] = -0.1*zFTDOuterCyllinderStart;
	  fkddes_.xrlc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness/radlen_kapton;
	  fkddes_.xelosc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness*dedx_kapton;
	  _Ncmat++;
	  // kapton cyllinder; right part
	  fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+0.5+
				       FTD_copper_thickness+
				       0.5*FTD_kaptonCyl_thickness);
	  fkddes_.zcmin[_Ncmat] = 0.1*zFTDOuterCyllinderStart;
	  fkddes_.zcmax[_Ncmat] = 0.1*zFTDOuterCyllinderEnd;
	  fkddes_.xrlc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness/radlen_kapton;
	  fkddes_.xelosc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness*dedx_kapton;
	  _Ncmat++;
	  
	}
	
	else if(has_si872==true && has_SupportDisks==false) {
	    
	    streamlog_out(DEBUG) << "build FTD according to the SupportRing structure" << std::endl;
	    
	    int nLayersFTD = 0;
	    int nFTDZ = int(theFTD.getDoubleVals("FTDZCoordinate").size());
	    int nFTDRin = int(theFTD.getDoubleVals("FTDInnerRadius").size());
	    int nFTDRout = int(theFTD.getDoubleVals("FTDOuterRadius").size());
	    
	    
	    if (nFTDZ == nFTDRin && nFTDRin == nFTDRout) {
	      nLayersFTD = nFTDZ;
	      zFTD.resize(nLayersFTD);
	      rInFTD.resize(nLayersFTD);
	      rOutFTD.resize(nLayersFTD);
	    }
	    else {
	      streamlog_out(DEBUG) << "Size of vectors FTDZCoordinate, FTDInnerRadius and  FTDInnerRadius are not equal --->" << std::endl;
	      streamlog_out(DEBUG) << "# FTDZCoordinate : " << nFTDZ << std::endl;
	      streamlog_out(DEBUG) << "# FTDInnerRadius : " << nFTDRin << std::endl;
	      streamlog_out(DEBUG) << "# FTDOuterRadius : " << nFTDRout << std::endl;
	      exit(1);
	    }
	    for (int i=0;i<nLayersFTD;++i) {
	      zFTD[i] = float(theFTD.getDoubleVals("FTDZCoordinate")[i]);
	      rInFTD[i] = float(theFTD.getDoubleVals("FTDInnerRadius")[i]);
	      rOutFTD[i] = float(theFTD.getDoubleVals("FTDOuterRadius")[i]);
	    } 
	    float FTDdisk_thickness = float(theFTD.getDoubleVal("FTDDiskThickness"));
	    float FTD_innerSupport_dR = float(theFTD.getDoubleVal("FTDInnerSupportdR"));
	    float FTD_outerSupport_dR = float(theFTD.getDoubleVal("FTDOuterSupportdR"));
	    float FTD_innerSupport_thickness = float(theFTD.getDoubleVal("FTDInnerSupportThickness"));
	    float FTD_outerSupport_thickness = float(theFTD.getDoubleVal("FTDOuterSupportThickness"));
	    float zFTDOuterCyllinderStart = float(theFTD.getDoubleVal("zFTDOuterCylinderStart"));
	    float zFTDOuterCyllinderEnd = float(theFTD.getDoubleVal("zFTDOuterCylinderEnd"));
//	    float zFTDInnerConeStart = float(theFTD.getDoubleVal("zFTDInnerConeStart"));
//	    float zFTDInnerConeEnd = float(theFTD.getDoubleVal("zFTDInnerConeEnd"));
	    float FTD_copper_thickness = float(theFTD.getDoubleVal("FTDCopperThickness"));
	    float FTD_kaptonCyl_thickness = float(theFTD.getDoubleVal("FTDOuterCylinderThickness"));
	    int iLast = theFTD.getIntVal("LastHeavyLayer");
	    float dedx_si = 10.0*float(theFTD.getDoubleVal("Silicon_dEdx"));
	    float dedx_si872(0);
	    float radlen_si872(0);
	    if (iLast>0) {
	      dedx_si872 = 10.0*float(theFTD.getDoubleVal("Silicon872_dEdx"));
	      radlen_si872 = 0.1*float(theFTD.getDoubleVal("Silicon872_RadLen"));
	    }
	    float dedx_kapton = 10.0*float(theFTD.getDoubleVal("Kapton_dEdx"));
	    float dedx_copper = 10.0*float(theFTD.getDoubleVal("Copper_dEdx"));
	    float radlen_si = 0.1*float(theFTD.getDoubleVal("Silicon_RadLen"));  
	    float radlen_kapton = 0.1*float(theFTD.getDoubleVal("Kapton_RadLen"));
	    float radlen_copper = 0.1*float(theFTD.getDoubleVal("Copper_RadLen"));
	    
	    for (int i=0;i<nLayersFTD;++i) {
	      // FTD Si Disks
	      float dedx = dedx_si;
	      float radlen  = radlen_si;
	      if (i<iLast) {
		dedx = dedx_si872;
		radlen = radlen_si872;
	      }
	      // right-hand part
	      fkddes_.zpmat[_Npmat] = 0.1*zFTD[i];
	      fkddes_.rpmin[_Npmat] = 0.1*rInFTD[i];
	      fkddes_.rpmax[_Npmat] = 0.1*rOutFTD[i];
	      fkddes_.xrlp[_Npmat] = 0.1*FTDdisk_thickness/radlen;
	      fkddes_.xelosp[_Npmat] = 0.1*FTDdisk_thickness*dedx;
	      
	      fkexts_.itexts[_Nexs] = 1;
	      fkexts_.rzsurf[_Nexs] = fkddes_.zpmat[_Npmat];
	      fkexts_.zrmin[_Nexs] = fkddes_.rpmin[_Npmat];
	      fkexts_.zrmax[_Nexs] = fkddes_.rpmax[_Npmat];
	      
	      _Nexs++;
	      _Npmat++;
	      
	      // left-hand part
	      fkddes_.zpmat[_Npmat] = -0.1*zFTD[i];
	      fkddes_.rpmin[_Npmat] = 0.1*rInFTD[i];
	      fkddes_.rpmax[_Npmat] = 0.1*rOutFTD[i];
	      fkddes_.xrlp[_Npmat] = 0.1*FTDdisk_thickness/radlen;
	      fkddes_.xelosp[_Npmat] = 0.1*FTDdisk_thickness*dedx;
	      
	      fkexts_.itexts[_Nexs] = 1;
	      fkexts_.rzsurf[_Nexs] = fkddes_.zpmat[_Npmat];
	      fkexts_.zrmin[_Nexs] = fkddes_.rpmin[_Npmat];
	      fkexts_.zrmax[_Nexs] = fkddes_.rpmax[_Npmat];
	      
	      _Nexs++;
	      _Npmat++;
	      
	      // Inner Support Rings
	      // right-hand part
	      fkddes_.zpmat[_Npmat] = 0.1*zFTD[i];
	      fkddes_.rpmin[_Npmat] = 0.1*(rInFTD[i]-FTD_innerSupport_dR);
	      fkddes_.rpmax[_Npmat] = 0.1*rInFTD[i];
	      fkddes_.xrlp[_Npmat] = 0.1*FTD_innerSupport_thickness/radlen_kapton;
	      fkddes_.xelosp[_Npmat] = 0.1*FTD_innerSupport_thickness*dedx_kapton;
	      _Npmat++;
	      
	      // left-hand part
	      fkddes_.zpmat[_Npmat] = -0.1*zFTD[i];
	      fkddes_.rpmin[_Npmat] = 0.1*(rInFTD[i]-FTD_innerSupport_dR);
	      fkddes_.rpmax[_Npmat] = 0.1*rInFTD[i];
	      fkddes_.xrlp[_Npmat] = 0.1*FTD_innerSupport_thickness/radlen_kapton;
	      fkddes_.xelosp[_Npmat] = 0.1*FTD_innerSupport_thickness*dedx_kapton;
	      _Npmat++;
	      
	      // Outer Support
	      // right-hand part
	      fkddes_.zpmat[_Npmat] = 0.1*zFTD[i];
	      fkddes_.rpmin[_Npmat] = 0.1*rOutFTD[i];
	      fkddes_.rpmax[_Npmat] = 0.1*(rOutFTD[i]+FTD_outerSupport_dR);
	      fkddes_.xrlp[_Npmat] = 0.1*FTD_outerSupport_thickness/radlen_kapton;
	      fkddes_.xelosp[_Npmat] = 0.1*FTD_outerSupport_thickness*dedx_kapton;
	      _Npmat++;
	      
	      // left-hand part
	      fkddes_.zpmat[_Npmat] = -0.1*zFTD[i];
	      fkddes_.rpmin[_Npmat] = 0.1*rOutFTD[i];
	      fkddes_.rpmax[_Npmat] = 0.1*(rOutFTD[i]+FTD_outerSupport_dR);
	      fkddes_.xrlp[_Npmat]  = 0.1*FTD_outerSupport_thickness/radlen_kapton;
	      fkddes_.xelosp[_Npmat]= 0.1*FTD_outerSupport_thickness*dedx_kapton;
	      _Npmat++;
	    }
	    
	    
	    
	    // copper cables; left part
	    fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+
					FTD_outerSupport_dR+0.5+
					0.5*FTD_copper_thickness);
	    fkddes_.zcmin[_Ncmat] = -0.1*zFTDOuterCyllinderEnd;
	    fkddes_.zcmax[_Ncmat] = -0.1*zFTDOuterCyllinderStart;
	    fkddes_.xrlc[_Ncmat] = 0.1*FTD_copper_thickness/radlen_copper;
	    fkddes_.xelosc[_Ncmat] = 0.1*FTD_copper_thickness*dedx_copper;
	    _Ncmat++;
	    
	    // copper cables; right part
	    fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+
					FTD_outerSupport_dR+0.5+ 
					0.5*FTD_copper_thickness);
	    fkddes_.zcmin[_Ncmat] = 0.1*zFTDOuterCyllinderStart;
	    fkddes_.zcmax[_Ncmat] = 0.1*zFTDOuterCyllinderEnd;
	    fkddes_.xrlc[_Ncmat] = 0.1*FTD_copper_thickness/radlen_copper;
	    fkddes_.xelosc[_Ncmat] = 0.1*FTD_copper_thickness*dedx_copper;
	    _Ncmat++;  
	    
	    // kapton cyllinder; left part
	    fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+
					FTD_outerSupport_dR+0.5+
					FTD_copper_thickness+
					0.5*FTD_kaptonCyl_thickness);
	    fkddes_.zcmin[_Ncmat] = -0.1*zFTDOuterCyllinderEnd;
	    fkddes_.zcmax[_Ncmat] = -0.1*zFTDOuterCyllinderStart;
	    fkddes_.xrlc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness/radlen_kapton;
	    fkddes_.xelosc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness*dedx_kapton;
	    _Ncmat++;
	    
	    // kapton cyllinder; right part
	    fkddes_.rcmat[_Ncmat] = 0.1*(rOutFTD[nLayersFTD-1]+
					FTD_outerSupport_dR+0.5+
					FTD_copper_thickness+
					0.5*FTD_kaptonCyl_thickness);
	    fkddes_.zcmin[_Ncmat] = 0.1*zFTDOuterCyllinderStart;
	    fkddes_.zcmax[_Ncmat] = 0.1*zFTDOuterCyllinderEnd;
	    fkddes_.xrlc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness/radlen_kapton;
	    fkddes_.xelosc[_Ncmat] = 0.1*FTD_kaptonCyl_thickness*dedx_kapton;
	    _Ncmat++;
	  }
	  
	  else {
	    std::stringstream errorMsg;
	    errorMsg << "MaterialDB Processor : FTD Geometery not correctly described. \n"
		      << " It is neither SupportRing or SupportDisk based."
		      <<  std::endl;
	    throw gear::Exception(errorMsg.str());
	  }
	  
      }
    return build_FTD;	 
  }
   

  bool MaterialDB_F77::buildSET(){ 

    streamlog_out(DEBUG) << "build SET  ..." << std::endl;
    bool build_SET = true;

    try{
      Global::GEAR->getGearParameters("SET");      
    }

    catch(gear::UnknownParameterException){
      streamlog_out(MESSAGE) << "SET not present in Gear file -- SET Material not created" << std::endl;
      build_SET = false;
    }

    if(build_SET){
          // SET layers
      const gear::GearParameters& pSETDet = Global::GEAR->getGearParameters("SET");

      std::vector<float> rSET;
      std::vector<float> halfZSET;
      std::vector<float> rSETSupport;
      std::vector<float> halfZSETSupport;

      std::vector<double> SETLayer_thickness ;
      std::vector<double> SETLayerSupport_thickness ;

      int nSETR = int(pSETDet.getDoubleVals("SETLayerRadius").size());
      int nSETHL = int(pSETDet.getDoubleVals("SETLayerHalfLength").size());
      //  int SETModel = int(pSETDet.getIntVal("SETModel"));
      int nLayersSET = 0;
      
      if (nSETR == nSETHL) {
	nLayersSET = nSETR;
	rSET.resize(nLayersSET);
	halfZSET.resize(nLayersSET);
	rSETSupport.resize(nLayersSET);
	halfZSETSupport.resize(nLayersSET);
      }
      else {
	std::stringstream errorMsg;
	errorMsg << "Size of SETLayerRadius vector (" << nSETR 
		  << ") is not equal to the size of SETHalfLength vector ("
		  << nSETHL << ")" << std::endl;
	throw gear::Exception(errorMsg.str());
      }
      
      
      getDoubleValues(SETLayer_thickness,        pSETDet, "SETLayerThickness", nLayersSET );
      getDoubleValues(SETLayerSupport_thickness, pSETDet, "SETSupportLayerThickness", nLayersSET );
      
      for (int iL=0;iL<nLayersSET;++iL) {
	rSET[iL] = float(pSETDet.getDoubleVals("SETLayerRadius")[iL]);
	halfZSET[iL] = float(pSETDet.getDoubleVals("SETLayerHalfLength")[iL]);
	rSETSupport[iL] = float(pSETDet.getDoubleVals("SETSupportLayerRadius")[iL]);
	halfZSETSupport[iL] = float(pSETDet.getDoubleVals("SETSupportLayerHalfLength")[iL]);
      }
      
      float radlen_si = 0.1*float(pSETDet.getDoubleVal("SETLayer_RadLen"));
      float dedx_si = 10.*float(pSETDet.getDoubleVal("SETLayer_dEdx"));
      
      
      float radlen_ber = 0.1*float(pSETDet.getDoubleVal("SETSupportLayer_RadLen"));
      float dedx_ber = 10.*float(pSETDet.getDoubleVal("SETSupportLayer_dEdx"));
      
      
      for (int iL = 0; iL < nLayersSET; ++iL) {
	
		
	fkddes_.rcmat[_Ncmat] = 0.1*rSET[iL];
	fkddes_.zcmin[_Ncmat] = -0.1*halfZSET[iL];
	fkddes_.zcmax[_Ncmat] = 0.1*halfZSET[iL];
	fkddes_.xrlc[_Ncmat] = 0.1*SETLayer_thickness[iL]/radlen_si;
	fkddes_.xelosc[_Ncmat] = 0.1*SETLayer_thickness[iL]*dedx_si;     
	++_Ncmat;
	
	fkddes_.rcmat[_Ncmat] = 0.1*rSETSupport[iL];
	fkddes_.zcmin[_Ncmat] = -0.1*halfZSETSupport[iL];
	fkddes_.zcmax[_Ncmat] = 0.1*halfZSETSupport[iL];
	fkddes_.xrlc[_Ncmat] = 0.1*SETLayerSupport_thickness[iL]/radlen_ber;
	fkddes_.xelosc[_Ncmat] = 0.1*SETLayerSupport_thickness[iL]*dedx_ber;
	++_Ncmat;
	
    
	fkexts_.itexts[_Nexs] = 0;
	fkexts_.rzsurf[_Nexs] = 0.1*rSET[iL];
	fkexts_.zrmin[_Nexs]  = -halfZSET[iL];
	fkexts_.zrmax[_Nexs]  = halfZSET[iL];
	
	++_Nexs;	      	      
      } 
    }
    return build_SET;
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
	  fkexts_.rzsurf[_Nexs] =  0.1*radius_sen-0.1; // place an extrapolation surface just inside the inner edge sensitive material 
	  fkexts_.zrmin[_Nexs]  = -0.1*halfz;
	  fkexts_.zrmax[_Nexs]  =  0.1*halfz;
	  
	  _Nexs++;
	
	}	
	return build_SIT_Simple;
      }
    return false ; // if this point is reached nothing was sucessfully built
  }

} // end of marlin_delphiF77 namespace 
