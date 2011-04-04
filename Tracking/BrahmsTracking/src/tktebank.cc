#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tktebank.h"

using namespace std;

Tk_Te_Bank::~Tk_Te_Bank()
{
}

void Tk_Te_Bank::clear(){
  te_bank.clear();
  itedat_bank.clear();
}

void Tk_Te_Bank::add_te(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
{

  tk_te ate;
  ate.subdetector_ID = subid;
  ate.submodule = submod;
  ate.unused = unused;
  ate.measurement_code = MesrCode;
  ate.pointer_to_end_of_TE = PnteTE;
  ate.charge = Q;
  ate.ndf = ndf;
  ate.chi2 = chi2;
  ate.length = L;
  ate.coord1_of_ref_point = cord1;
  ate.coord2_of_ref_point = cord2;
  ate.coord3_of_ref_point = cord3;
  ate.theta = theta;
  ate.phi = phi;
  ate.invp = invp;
  ate.de_dx = dedx;
  ate.covmatrix1 = cov1;
  ate.covmatrix2 = cov2;
  ate.covmatrix3 = cov3;
  ate.covmatrix4 = cov4;
  ate.covmatrix5 = cov5;
  ate.covmatrix6 = cov6;
  ate.covmatrix7 = cov7;
  ate.covmatrix8 = cov8;
  ate.covmatrix9 = cov9;
  ate.covmatrix10 = cov10;
  ate.covmatrix11 = cov11;
  ate.covmatrix12 = cov12;
  ate.covmatrix13 = cov13;
  ate.covmatrix14 = cov14;
  ate.covmatrix15 = cov15;


  te_bank.push_back(ate);

  unsigned int size = itedat_bank.size() ;
  itedat_bank.resize(size+1) ;

}

// due to the fact that the bank will be accessed by integer value of the te number this is probaly dangerous
// so leave it out for now

//void Tk_Te_Bank::remove_te(int te)

//{
//  te_bank.erase(te_bank.begin()+te);
//}

//int tkmktecpp(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
//{
// 
//  Tk_Te_Bank::Instance().add_te(subid, submod, unused, MesrCode, PnteTE, Q, ndf, chi2, L, cord1, cord2, cord3, theta, phi, invp, dedx, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
// 
//  return 0;
//}

int tkmktecpp(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float* cov)
{

  float cov1 = cov[0]; 
  float cov2 = cov[1];
  float cov3 = cov[2];
  float cov4 = cov[3];
  float cov5 = cov[4];
  float cov6 = cov[5];
  float cov7 = cov[6];
  float cov8 = cov[7];
  float cov9 = cov[8];
  float cov10 = cov[9];
  float cov11 = cov[10];
  float cov12 = cov[11];
  float cov13 = cov[12];
  float cov14 = cov[13];
  float cov15 = cov[14];

//   cout << "cov1 = " << cov1 << endl; 
//   cout << "cov2 = " << cov2 << endl; 
//   cout << "cov3 = " << cov3 << endl; 
//   cout << "cov4 = " << cov4 << endl; 
//   cout << "cov5 = " << cov5 << endl; 
//   cout << "cov6 = " << cov6 << endl; 
//   cout << "cov15 = " << cov15 << endl; 

  Tk_Te_Bank::Instance().add_te(subid, submod, unused, MesrCode, PnteTE, Q, ndf, chi2, L, cord1, cord2, cord3, theta, phi, invp, dedx, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
  return 0;
}
FCALLSCFUN17(INT,tkmktecpp,TKMKTECPP,tkmktecpp, INT , INT ,INT ,INT ,INT ,INT ,INT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOATV)

float rreadtktecpp(int attribute, int te)
{

  te = te - 1;

  if(te>Tk_Te_Bank::Instance().size()) return 0.;

  switch (attribute) {
  case 1: 
    return Tk_Te_Bank::Instance().getSubdetector_ID(te);
    break;
  case 2: 
    return Tk_Te_Bank::Instance().getSubmodule(te);
     break;
  case 3: 
    return Tk_Te_Bank::Instance().getUnused(te);
    break;
  case 4: 
    return Tk_Te_Bank::Instance().getMeasurement_code(te);
    break;
  case 5: 
    return Tk_Te_Bank::Instance().getPointer_to_end_of_TE(te);
    break;
  case 6: 
    return Tk_Te_Bank::Instance().getCharge(te);
    break;
  case 7: 
    return Tk_Te_Bank::Instance().getNdf(te);
    break;
  case 8: 
    return Tk_Te_Bank::Instance().getChi2(te);
    break;
  case 9: 
    return Tk_Te_Bank::Instance().getLength(te);
    break;
  case 10: 
    return Tk_Te_Bank::Instance().getCoord1_of_ref_point(te);
    break;
  case 11: 
    return Tk_Te_Bank::Instance().getCoord2_of_ref_point(te);
    break;
  case 12: 
    return Tk_Te_Bank::Instance().getCoord3_of_ref_point(te);
    break;
  case 13: 
    return Tk_Te_Bank::Instance().getTheta(te);
    break;
  case 14: 
    return Tk_Te_Bank::Instance().getPhi(te);
    break;
  case 15: 
    return Tk_Te_Bank::Instance().getInvp(te);
    break;
  case 16: 
    return Tk_Te_Bank::Instance().getDe_dx(te);
    break;
  case 17: 
    return Tk_Te_Bank::Instance().getCovmatrix1(te);
    break;
  case 18: 
    return Tk_Te_Bank::Instance().getCovmatrix2(te);
    break;
  case 19: 
    return Tk_Te_Bank::Instance().getCovmatrix3(te);
    break;
  case 20: 
    return Tk_Te_Bank::Instance().getCovmatrix4(te);
    break;
  case 21: 
    return Tk_Te_Bank::Instance().getCovmatrix5(te);
    break;
  case 22: 
    return Tk_Te_Bank::Instance().getCovmatrix6(te);
    break;
  case 23: 
    return Tk_Te_Bank::Instance().getCovmatrix7(te);
    break;
  case 24: 
    return Tk_Te_Bank::Instance().getCovmatrix8(te);
    break;
  case 25: 
    return Tk_Te_Bank::Instance().getCovmatrix9(te);
    break;
  case 26: 
    return Tk_Te_Bank::Instance().getCovmatrix10(te);
    break;
  case 27: 
    return Tk_Te_Bank::Instance().getCovmatrix11(te);
    break;
  case 28: 
    return Tk_Te_Bank::Instance().getCovmatrix12(te);
    break;
  case 29: 
    return Tk_Te_Bank::Instance().getCovmatrix13(te);
    break;
  case 30: 
    return Tk_Te_Bank::Instance().getCovmatrix14(te);
    break;
  case 31: 
    return Tk_Te_Bank::Instance().getCovmatrix15(te);
    break;
  default: 
    throw runtime_error("te attribute not valid");
  } 
}
FCALLSCFUN2(FLOAT,rreadtktecpp,RREADTKTECPP,rreadtktecpp, INT, INT)

int ireadtktecpp(int attribute, int te)
{

  te = te - 1;

  if(te>Tk_Te_Bank::Instance().size()) return 0;

  switch (attribute) {
  case 1: 
    return (int)Tk_Te_Bank::Instance().getSubdetector_ID(te);
    break;
  case 2: 
    return (int)Tk_Te_Bank::Instance().getSubmodule(te);
     break;
  case 3: 
    return (int)Tk_Te_Bank::Instance().getUnused(te);
    break;
  case 4: 
    return (int)Tk_Te_Bank::Instance().getMeasurement_code(te);
    break;
  case 5: 
    return (int)Tk_Te_Bank::Instance().getPointer_to_end_of_TE(te);
    break;
  case 6: 
    return (int)Tk_Te_Bank::Instance().getCharge(te);
    break;
  case 7: 
    return (int)Tk_Te_Bank::Instance().getNdf(te);
    break;
  case 8: 
    return (int)Tk_Te_Bank::Instance().getChi2(te);
    break;
  case 9: 
    return (int)Tk_Te_Bank::Instance().getLength(te);
    break;
  case 10: 
    return (int)Tk_Te_Bank::Instance().getCoord1_of_ref_point(te);
    break;
  case 11: 
    return (int)Tk_Te_Bank::Instance().getCoord2_of_ref_point(te);
    break;
  case 12: 
    return (int)Tk_Te_Bank::Instance().getCoord3_of_ref_point(te);
    break;
  case 13: 
    return (int)Tk_Te_Bank::Instance().getTheta(te);
    break;
  case 14: 
    return (int)Tk_Te_Bank::Instance().getPhi(te);
    break;
  case 15: 
    return (int)Tk_Te_Bank::Instance().getInvp(te);
    break;
  case 16: 
    return (int)Tk_Te_Bank::Instance().getDe_dx(te);
    break;
  case 17: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix1(te);
    break;
  case 18: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix2(te);
    break;
  case 19: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix3(te);
    break;
  case 20: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix4(te);
    break;
  case 21: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix5(te);
    break;
  case 22: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix6(te);
    break;
  case 23: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix7(te);
    break;
  case 24: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix8(te);
    break;
  case 25: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix9(te);
    break;
  case 26: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix10(te);
    break;
  case 27: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix11(te);
    break;
  case 28: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix12(te);
    break;
  case 29: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix13(te);
    break;
  case 30: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix14(te);
    break;
  case 31: 
    return (int)Tk_Te_Bank::Instance().getCovmatrix15(te);
    break;
  default: 
    throw runtime_error("te attribute not valid");
  } 
}
FCALLSCFUN2(INT,ireadtktecpp,IREADTKTECPP,ireadtktecpp, INT, INT)


int writetktecpp(float value, int attribute, int te){

  if(te>Tk_Te_Bank::Instance().size()) return 1;

  te = te - 1;

  switch (attribute) {
  case 1: 
    Tk_Te_Bank::Instance().setSubdetector_ID((int)value,te);
    return 0;
    break;
  case 2: 
    Tk_Te_Bank::Instance().setSubmodule((int)value,te);
    return 0;
    break;
  case 3: 
    Tk_Te_Bank::Instance().setUnused((int)value,te);
    return 0;
    break;
  case 4: 
    Tk_Te_Bank::Instance().setMeasurement_code((int)value,te);
    return 0;
    break;
  case 5: 
    Tk_Te_Bank::Instance().setPointer_to_end_of_TE((int)value,te);
    return 0;
    break;
  case 6: 
    Tk_Te_Bank::Instance().setCharge((int)value,te);
    return 0;
    break;
  case 7: 
    Tk_Te_Bank::Instance().setNdf((int)value,te);
    return 0;
    break;
  case 8: 
    Tk_Te_Bank::Instance().setChi2(value,te);
    return 0;
    break;
  case 9: 
    Tk_Te_Bank::Instance().setLength(value,te);
    return 0;
    break;
  case 10: 
    Tk_Te_Bank::Instance().setCoord1_of_ref_point(value,te);
    return 0;
    break;
  case 11: 
    Tk_Te_Bank::Instance().setCoord2_of_ref_point(value,te);
    return 0;
    break;
  case 12: 
    Tk_Te_Bank::Instance().setCoord3_of_ref_point(value,te);
    return 0;
    break;
  case 13: 
    Tk_Te_Bank::Instance().setTheta(value,te);
    return 0;
    break;
  case 14: 
    Tk_Te_Bank::Instance().setPhi(value,te);
    return 0;
    break;
  case 15: 
    Tk_Te_Bank::Instance().setInvp(value,te);
    return 0;
    break;
  case 16: 
    Tk_Te_Bank::Instance().setDe_dx(value,te);
    return 0;
    break;
  case 17: 
    Tk_Te_Bank::Instance().setCovmatrix1(value,te);
    return 0;
    break;
  case 18: 
    Tk_Te_Bank::Instance().setCovmatrix2(value,te);
    return 0;
    break;
  case 19: 
    Tk_Te_Bank::Instance().setCovmatrix3(value,te);
    return 0;
    break; 
  case 20: 
    Tk_Te_Bank::Instance().setCovmatrix4(value,te);
    return 0;
    break;
  case 21: 
    Tk_Te_Bank::Instance().setCovmatrix5(value,te);
    return 0;
    break;
  case 22: 
    Tk_Te_Bank::Instance().setCovmatrix6(value,te);
    return 0;
    break;
  case 23: 
    Tk_Te_Bank::Instance().setCovmatrix7(value,te);
    return 0;
    break;
  case 24: 
    Tk_Te_Bank::Instance().setCovmatrix8(value,te);
    return 0;
    break;
  case 25: 
    Tk_Te_Bank::Instance().setCovmatrix9(value,te);
    return 0;
    break;
  case 26: 
    Tk_Te_Bank::Instance().setCovmatrix10(value,te);
    return 0;
    break;
  case 27: 
    Tk_Te_Bank::Instance().setCovmatrix11(value,te);
    return 0;
    break;
  case 28: 
    Tk_Te_Bank::Instance().setCovmatrix12(value,te);
    return 0;
    break;
  case 29: 
    Tk_Te_Bank::Instance().setCovmatrix13(value,te);
    return 0;
    break;
  case 30: 
    Tk_Te_Bank::Instance().setCovmatrix14(value,te);
    return 0;
    break;
  case 31: 
    Tk_Te_Bank::Instance().setCovmatrix15(value,te);
    return 0;
    break;
  default: 
    cout << "attribute = " << attribute << endl;  
    throw runtime_error("te attribute not valid");
  } 
  
}
FCALLSCFUN3(INT,writetktecpp,WRITETKTECPP,writetktecpp, FLOAT, INT, INT)

int addhittktecpp(int hit, int te) 
{
  te = te - 1;
  hit = hit - 1;

  Tk_Te_Bank::Instance().addHit(hit,te);
  return 0;
}
FCALLSCFUN2(INT,addhittktecpp,ADDHITTKTECPP,addhittktecpp, INT, INT)

int writetkitedatcpp(int value, int attribute, int te){
  te = te - 1;
  
  switch (attribute) {
  case 1: 
    Tk_Te_Bank::Instance().setPosOfFirstHitInHitList(value,te);
    return 0;
    break;
  case 2: 
    Tk_Te_Bank::Instance().setNumOfHits(value,te);
    return 0;
    break;
  case 3: 
    Tk_Te_Bank::Instance().setPointrToFirstExclusion(value,te);
    return 0;
    break;
  case 4: 
    Tk_Te_Bank::Instance().setNumOfExclusions(value,te);
    return 0;
    break;
  case 5: 
    Tk_Te_Bank::Instance().setTrackNo(value,te);
    return 0;
    break;
  default: 
    std::cout << "attribute = " << attribute <<  std::endl;
    throw runtime_error("writetkteitedatcpp: te attribute not valid");
  } 

}
FCALLSCFUN3(INT,writetkitedatcpp,WRITETKITEDATCPP,writetkitedatcpp, INT, INT, INT)


int readtkitedatcpp(int attribute, int te)
{

  te = te - 1;

  if(te>Tk_Te_Bank::Instance().size()) return 0;

  switch (attribute) {
  case 1: 
    return Tk_Te_Bank::Instance().getPosOfFirstHitInHitList(te);
    break;
  case 2: 
    return Tk_Te_Bank::Instance().getNumOfHits(te);
    break;
  case 3: 
    return Tk_Te_Bank::Instance().getPointrToFirstExclusion(te);
    break;
  case 4: 
    return Tk_Te_Bank::Instance().getNumOfExclusions(te);
    break;
  case 5: 
    return Tk_Te_Bank::Instance().gettrackNo(te);
    break;
  default:
    std::cout << "attribute = " << attribute << std::endl ;
    std::cout << "te = " << te << std::endl ;
    throw runtime_error("te attribute not valid");
  } 

}
FCALLSCFUN2(INT,readtkitedatcpp,READTKITEDATCPP,readtkitedatcpp, INT, INT)
