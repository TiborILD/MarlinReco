#ifndef SimDigital_HHH
#define SimDigital_HHH

#include "marlin/Processor.h"
#include "lcio.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <stdlib.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <marlin/Global.h>
#include <gear/GearParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <TROOT.h>
#include <TFile.h>
#include <TF2.h>
#include <TH1.h>
#include "TH1F.h"

#include "CalorimeterHitType.h" //in MarlinUtil
#include "marlinutil/LCGeometryTypes.h"

class TH1F;
class TF1;
class TTree;

namespace AIDA {
  class ITuple;
}

using namespace lcio ;
using namespace marlin ;

// Gerald Grenier :
// start removing the ECAL part whih has no use for SDHCAL
// this can be fully activated when there is an ECAL digitizer 
// that don't do HCAL digitization
#define SIMDIGITAL_WITHECAL

#ifdef SIMDIGITAL_WITHECAL
const int MAX_LAYERS = 200;
const int MAX_STAVES =  16;
#endif 

/** Digitization for the SDHcal - based on NewLDCCaloDigi. 
 *
 *  @author  G.Grenier, INPL
 *  @author  R.Han, INPL
 *  @version $Id$
 */



class SimDigitalGeomRPCFrame;

struct StepAndCharge
{
  StepAndCharge() : step(), charge(0) {}
  StepAndCharge(LCVector3D vec) : step(vec), charge(0) {}
  LCVector3D step;
  float charge;
};


//helper class to manage cellId and local geometry
class SimDigitalGeomCellId 
{
 public:
  static void bookTuples(const marlin::Processor* proc);
  SimDigitalGeomCellId(LCCollection *inputCol, LCCollectionVec *outputCol);
  ~SimDigitalGeomCellId();
  //return the list of step positions in coordinates corresponding to 'I' ,'J' and 'layer'
  std::vector<StepAndCharge> decode(SimCalorimeterHit *hit);
  void encode(CalorimeterHitImpl *hit,int delta_I, int delta_J);
  void setLayerLayout( CHT::Layout layout);
  float getCellSize();
  const LCVector3D& normalToRPCPlane() {return _normal;}
  const LCVector3D& Iaxis() {return _Iaxis;}
  const LCVector3D& Jaxis() {return _Jaxis;}

  int I() {return _Iy;}
  int J() {return _Jz;}
  int K() {return _trueLayer;}
  int stave() {return _stave;}
  int module() {return _module;}

 private:
  enum HCAL_GEOM {VIDEAU,TESLA};
  HCAL_GEOM _geom;
  int _trueLayer;
  int _stave;
  int _module;
  int _Iy;
  int _Jz;
  const float* _hitPosition; 
  CellIDDecoder<SimCalorimeterHit> _decoder;
  CellIDEncoder<CalorimeterHitImpl> _encoder;
  const gear::LayerLayout* _layerLayout;
  SimDigitalGeomRPCFrame* _normal_I_J_setter;
  CHT::Layout _currentHCALCollectionCaloLayout;
  LCVector3D _normal;
  LCVector3D _Iaxis;
  LCVector3D _Jaxis;
  static AIDA::ITuple* _tupleHit;
  enum {TH_DETECTOR,TH_CHTLAYOUT,TH_MODULE,TH_STAVE,TH_LAYER,TH_I,TH_J,
	TH_X,TH_Y,TH_Z,
	TH_NORMALX,TH_NORMALY,TH_NORMALZ,
	TH_IX,TH_IY,TH_IZ,
	TH_JX,TH_JY,TH_JZ};
  static AIDA::ITuple* _tupleStep;
  enum {TS_DETECTOR,TS_CHTLAYOUT,TS_HITCELLID,TS_NSTEP,
	TS_HITX,TS_HITY,TS_HITZ,
	TS_STEPX,TS_STEPY,TS_STEPZ,
	TS_DELTAI,TS_DELTAJ,TS_DELTALAYER};
  friend class SimDigitalGeomRPCFrame;
};


//hierarchy of classes to determine the RPC reference frame 
class SimDigitalGeomRPCFrame 
{
 public:
  SimDigitalGeomRPCFrame(SimDigitalGeomCellId& h) : _layerInfo(h) {;}
  virtual ~SimDigitalGeomRPCFrame() {}  
  virtual void setRPCFrame()=0;
 private:
  SimDigitalGeomCellId& _layerInfo;
 protected:
  int stave() {return _layerInfo._stave;}
  int module() {return _layerInfo._module;}
  LCVector3D& normal() {return _layerInfo._normal;}
  LCVector3D& Iaxis() {return _layerInfo._Iaxis;}
  LCVector3D& Jaxis() {return _layerInfo._Jaxis;}
};

class SimDigitalGeomRPCFrame_TESLA_BARREL : public SimDigitalGeomRPCFrame
{
 public:
  SimDigitalGeomRPCFrame_TESLA_BARREL(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
  void setRPCFrame();
};
class SimDigitalGeomRPCFrame_VIDEAU_BARREL : public SimDigitalGeomRPCFrame
{
 public:
  SimDigitalGeomRPCFrame_VIDEAU_BARREL(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
  void setRPCFrame();
};
class SimDigitalGeomRPCFrame_TESLA_ENDCAP : public SimDigitalGeomRPCFrame
{
 public:
  SimDigitalGeomRPCFrame_TESLA_ENDCAP(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
  void setRPCFrame();
};
class SimDigitalGeomRPCFrame_VIDEAU_ENDCAP : public SimDigitalGeomRPCFrame
{
 public:
  SimDigitalGeomRPCFrame_VIDEAU_ENDCAP(SimDigitalGeomCellId& h) : SimDigitalGeomRPCFrame(h) {}
  void setRPCFrame();
};



//helper class to manage multiplicity
class multiplicityChargeSplitterBase
{
 public:
  multiplicityChargeSplitterBase();
  virtual ~multiplicityChargeSplitterBase() {}
  typedef std::pair<int,int> I_J_Coordinates;
  virtual void addCharge(float charge, float pos_I, float pos_J)=0;
  void newHit(float cell_Size) {_chargeMap.clear();_cellSize=cell_Size;}
  const std::map<I_J_Coordinates,float>& chargeMap() {return _chargeMap;}
 protected:
  float _cellSize;
  std::map<I_J_Coordinates,float> _chargeMap;
};


class multiplicityChargeSplitterUniform : public multiplicityChargeSplitterBase
{
 public:
  multiplicityChargeSplitterUniform();
  void addCharge(float charge, float pos_I, float pos_J);
 private:
  float _edgeSize;
  friend class SimDigital;
};


class multiplicityChargeSplitterFunction : public multiplicityChargeSplitterBase
{
 public:
  multiplicityChargeSplitterFunction();
  virtual ~multiplicityChargeSplitterFunction();
  void init();
  void addCharge(float charge, float pos_I, float pos_J);
 private:
  TF2* _f2;
  float _range;
  std::string _function_description;
  std::vector<float> _functionParameters;
  float _normalisation;
  float _RPC_PadSeparation; //distance between cells
  friend class SimDigital;
};

class multiplicityChargeSplitterErfFunction : public multiplicityChargeSplitterBase
{
 public:
  multiplicityChargeSplitterErfFunction();
  virtual ~multiplicityChargeSplitterErfFunction();
  void init();
  void addCharge(float charge, float pos_I, float pos_J);
 private:
  float _range;
  std::vector<float> _erfWidth;
  std::vector<float> _erfWeigth;
  float _normalisation;
  float _RPC_PadSeparation; //distance between cells
  friend class SimDigital;
};


//helper class to manage efficiency maps
class effMapBase
{
 public:
  virtual float getEfficiency(int I, int J, int K, int stave, int module)=0;
};

class effMapConstant : public effMapBase
{
 public:
 effMapConstant(float val=1.0) : value(val) {}
  virtual float getEfficiency(int I, int J, int K, int stave, int module) { return value;}
 private:
  float value;
};


class effMapProtoByAsic : public effMapBase
{
 public:
  effMapProtoByAsic(std::string fileName);
  virtual float getEfficiency(int I, int J, int K, int stave, int module) 
  {
    std::map<int,float>::iterator it=_effMap.find( (I-1)/8+((J-1)/8)*12+K*1000 );
    return (it != _effMap.end() ? it->second : 1.0);
  }
 private:
  std::map<int,float> _effMap;
};


class SimDigital : public Processor {
 public:
  virtual Processor*  newProcessor() { return new SimDigital;}
  SimDigital();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

/** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ;

  virtual void check( LCEvent * evt ) ;

  /** Called after data processing for clean up.
   */
  virtual void end() ;

#ifdef SIMDIGITAL_WITHECAL
  virtual void fillECALGaps() ;
#endif 
 
 private:
#ifdef SIMDIGITAL_WITHECAL
  std::vector<std::string> _ecalCollections;
  std::string _outputEcalCollection0;
  std::string _outputEcalCollection1;
  std::string _outputEcalCollection2;
  std::vector<std::string> _outputEcalCollections;
  std::vector<float> _calibrCoeffEcal;
  std::vector<int> _ecalLayers;
  int _digitalEcal;
  float _thresholdEcal;
  int _ecalGapCorrection;
  float _ecalGapCorrectionFactor;
  float _ecalModuleGapCorrectionFactor;
  float _ecalEndcapCorrectionFactor;
  std::vector<CalorimeterHitImpl*> _calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int> _calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  float _zOfEcalEndcap;
  float _barrelPixelSizeT[MAX_LAYERS];
  float _barrelPixelSizeZ[MAX_LAYERS];
  float _endcapPixelSizeX[MAX_LAYERS];
  float _endcapPixelSizeY[MAX_LAYERS];
  float _barrelStaveDir[MAX_STAVES][2];
#endif 
  //std::vector<int> _hcalLayers;
  //std::vector<std::string> _calorimeterHitCollections;
  std::vector<std::string> _hcalCollections;
  std::vector<std::string> _outputHcalCollections;
  std::map<std::string, int> _counters;
  std::vector<float> _thresholdHcal;
  std::vector<float>_calibrCoeffHcal;
  std::string _outputRelCollection;
  bool _printSimDigital;
  TF1 * _QPolya;
  double _polyaAverageCharge, _polyaFunctionWidthParameter;
  
  LCCollectionVec* _relcol;
#ifdef SIMDIGITAL_WITHECAL
  void registerECALparameters();
  void setECALgeometry();
  void processECAL(LCEvent* evt, LCFlagImpl& flag);
#endif 
  void processHCAL(LCEvent* evt, LCFlagImpl& flag);

  static bool sortStepWithCharge(StepAndCharge s1, StepAndCharge s2){return s1.charge>s2.charge;}

  //intermediate storage class
  struct hitMemory
  {
    hitMemory() : ahit(0),relatedHits(), maxEnergydueToHit(-1), rawHit(-1) {} 
    CalorimeterHitImpl *ahit;
    std::set<int> relatedHits;
    float maxEnergydueToHit;
    int rawHit;
  };
  float depositedEnergyInRPC;
  multiplicityChargeSplitterUniform _chargeSplitterUniform;
  multiplicityChargeSplitterFunction _chargeSplitterFunction;
  multiplicityChargeSplitterErfFunction _chargeSplitterErfFunction;
  multiplicityChargeSplitterBase* _theChosenSplitter;
  bool _doThresholds;
  std::string  _chargeSplitterOption;
  std::string _effMapFileName;
  float _constEffMapValue;
  std::string _effMapOption;
  effMapBase *_effMap;
  float _absZstepFilter;
  bool _keepAtLeastOneStep;
  float _minXYdistanceBetweenStep;
  AIDA::ITuple* _debugTupleStepFilter;
  AIDA::ITuple* _tupleStepFilter;
  AIDA::ITuple* _tupleCollection;
  multiplicityChargeSplitterBase& getSplitter() { return *_theChosenSplitter; }

  //predicate class to remove potential hit below threshold
  class ThresholdIsBelow
  {
    float value;
  public:
    ThresholdIsBelow(float f) : value(f) {;}
    bool operator()(std::pair<int,hitMemory> f) { return f.second.ahit->getEnergy()<value;}
  };

  //predicate class to remove steps
  class absZGreaterThan
  {
  public:
    absZGreaterThan(float val) : _value(val) {}
      bool operator()(StepAndCharge& v) { return fabs( v.step.z() ) >_value;}
  private:
    float _value;
  };
  class randomGreater
  {
  public:
  randomGreater(float val) : _value(val) {}
    bool operator()(StepAndCharge& v) { return double(rand())/RAND_MAX>_value;}
  private:
    float _value;
  };
  //helper function to remove steps too close in I,J
  void remove_adjacent_step(std::vector<StepAndCharge>& vec);
  void fillTupleStep(std::vector<StepAndCharge>& vec,int level);

  LCCollectionVec * processHCALCollection(LCCollection * col ,CHT::Layout layout, LCFlagImpl& flag);
  void createPotentialOutputHits(std::map<int,hitMemory>& myHitMap, LCCollection *col, SimDigitalGeomCellId& aGeomCellId);
  void removeHitsBelowThreshold(std::map<int,hitMemory>& myHitMap, float threshold);
  void applyThresholds(std::map<int,hitMemory>& myHitMap);

};
  
#endif
