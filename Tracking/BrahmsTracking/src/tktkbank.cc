#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tktkbank.h"

using namespace std;

// Global pointer to the Tk_Tk_Bank structure
//Tk_Tk_Bank * TkTkBank = NULL;


Tk_Tk_Bank::~Tk_Tk_Bank()
{
}

void Tk_Tk_Bank::clear(){
  tk_bank.clear();
  itkdat_bank.clear();
}

void Tk_Tk_Bank::add_tk(int modid,int subdetbits,int MesrCode,int tracktype,int numtes,int Charge,int unused,int ndf,float chi2,float L,float xstart,float ystart,float zstart,float xend,float yend, float zend,float cord1,float cord2,float cord3,float theta,float phi,float invp,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
{

  tk_tk atk;
  atk.modid = modid;
  atk.subdetbits = subdetbits;
  atk.measurement_code = MesrCode;
  atk.tracktype = tracktype;
  atk.numtes = numtes;
  atk.charge = Charge;
  atk.unused = unused;
  atk.ndf = ndf;
  atk.chi2 = chi2;
  atk.length = L;
  atk.xstart = xstart ;
  atk.ystart = ystart ;
  atk.zstart = zstart ;
  atk.xend = xend ;
  atk.yend = yend ;
  atk.zend = zend ;
  atk.coord1_of_ref_point = cord1;
  atk.coord2_of_ref_point = cord2;
  atk.coord3_of_ref_point = cord3;
  atk.theta = theta;
  atk.phi = phi;
  atk.invp = invp;
  atk.covmatrix1 = cov1;
  atk.covmatrix2 = cov2;
  atk.covmatrix3 = cov3;
  atk.covmatrix4 = cov4;
  atk.covmatrix5 = cov5;
  atk.covmatrix6 = cov6;
  atk.covmatrix7 = cov7;
  atk.covmatrix8 = cov8;
  atk.covmatrix9 = cov9;
  atk.covmatrix10 = cov10;
  atk.covmatrix11 = cov11;
  atk.covmatrix12 = cov12;
  atk.covmatrix13 = cov13;
  atk.covmatrix14 = cov14;
  atk.covmatrix15 = cov15;


  tk_bank.push_back(atk);

  unsigned int size = itkdat_bank.size() ;
  itkdat_bank.resize(size+1) ;

}

// due to the fact that the bank will be accessed by integer value of the tk number this is probaly dangerous
// so leave it out for now

//void Tk_Tk_Bank::remove_tk(int tk)

//{
//  tk_bank.erase(tk_bank.begin()+tk);
//}

//int tkmktkcpp(int subid,int submod,int unused,int MesrCode,int PntkTK,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
//{
// 
//  TkTkBank->add_tk(subid, submod, unused, MesrCode, PnteTK, Q, ndf, chi2, L, cord1, cord2, cord3, theta, phi, invp, dedx, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
// 
//  return 0;
//}


// The following are global functions 

FCALLSCFUN23(INT,tkmktkcpp,TKMKTKCPP,tkmktkcpp, INT , INT ,INT ,INT ,INT ,INT ,INT ,INT, FLOAT, FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT ,FLOAT, FLOAT, FLOAT, FLOATV)

int tkmktkcpp(int modid,int subdetbits,int MesrCode,int tracktype, int numtes,int Charge,int unused,int ndf,float chi2,float L,float xstart, float ystart, float zstart, float xend, float yend, float zend, float cord1,float cord2,float cord3,float theta,float phi,float invp,float* cov)
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

  Tk_Tk_Bank::Instance().add_tk(modid, subdetbits, MesrCode, tracktype, numtes, Charge, unused, ndf, chi2, L, xstart, ystart, zstart, xend, yend, zend, cord1, cord2, cord3, theta, phi, invp, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
  return 0;
}

FCALLSCFUN2(FLOAT,rreadtktkcpp,RREADTKTKCPP,rreadtktkcpp, INT, INT)
float rreadtktkcpp(int attribute, int tk)
{

  tk = tk - 1;

  if(tk>Tk_Tk_Bank::Instance().size()) return 0.;

  switch (attribute) {
  case 1: 
    return Tk_Tk_Bank::Instance().getMod_ID(tk);
    break;
  case 2: 
    return Tk_Tk_Bank::Instance().getSubdetbits(tk);
     break;
  case 3: 
    return Tk_Tk_Bank::Instance().getMeasurement_code(tk);
    break;
  case 4: 
    return Tk_Tk_Bank::Instance().getTrackType(tk);
    break;
  case 5: 
    return Tk_Tk_Bank::Instance().getNumTEs(tk);
    break;
  case 6: 
    return Tk_Tk_Bank::Instance().getCharge(tk);
    break;
  case 7: 
    return Tk_Tk_Bank::Instance().getUnused(tk);
    break;
  case 8: 
    return Tk_Tk_Bank::Instance().getNdf(tk);
    break;
  case 9: 
    return Tk_Tk_Bank::Instance().getChi2(tk);
    break;
  case 10: 
    return Tk_Tk_Bank::Instance().getLength(tk);
    break;
  case 11:
    return Tk_Tk_Bank::Instance().getStartX(tk); 
    break;
  case 12:
    return Tk_Tk_Bank::Instance().getStartY(tk); 
    break;
  case 13:
    return Tk_Tk_Bank::Instance().getStartZ(tk); 
    break;
  case 14:
    return Tk_Tk_Bank::Instance().getEndX(tk); 
    break;
  case 15:
    return Tk_Tk_Bank::Instance().getEndY(tk); 
    break;
  case 16:
    return Tk_Tk_Bank::Instance().getEndZ(tk); 
    break;
  case 17: 
    return Tk_Tk_Bank::Instance().getCoord1_of_ref_point(tk);
    break;
  case 18: 
    return Tk_Tk_Bank::Instance().getCoord2_of_ref_point(tk);
    break;
  case 19: 
    return Tk_Tk_Bank::Instance().getCoord3_of_ref_point(tk);
    break;
  case 20: 
    return Tk_Tk_Bank::Instance().getTheta(tk);
    break;
  case 21: 
    return Tk_Tk_Bank::Instance().getPhi(tk);
    break;
  case 22: 
    return Tk_Tk_Bank::Instance().getInvp(tk);
    break;
  case 23: 
    return Tk_Tk_Bank::Instance().getCovmatrix1(tk);
    break;
  case 24: 
    return Tk_Tk_Bank::Instance().getCovmatrix2(tk);
    break;
  case 25: 
    return Tk_Tk_Bank::Instance().getCovmatrix3(tk);
    break;
  case 26: 
    return Tk_Tk_Bank::Instance().getCovmatrix4(tk);
    break;
  case 27: 
    return Tk_Tk_Bank::Instance().getCovmatrix5(tk);
    break;
  case 28: 
    return Tk_Tk_Bank::Instance().getCovmatrix6(tk);
    break;
  case 29: 
    return Tk_Tk_Bank::Instance().getCovmatrix7(tk);
    break;
  case 30: 
    return Tk_Tk_Bank::Instance().getCovmatrix8(tk);
    break;
  case 31: 
    return Tk_Tk_Bank::Instance().getCovmatrix9(tk);
    break;
  case 32: 
    return Tk_Tk_Bank::Instance().getCovmatrix10(tk);
    break;
  case 33: 
    return Tk_Tk_Bank::Instance().getCovmatrix11(tk);
    break;
  case 34: 
    return Tk_Tk_Bank::Instance().getCovmatrix12(tk);
    break;
  case 35: 
    return Tk_Tk_Bank::Instance().getCovmatrix13(tk);
    break;
  case 36: 
    return Tk_Tk_Bank::Instance().getCovmatrix14(tk);
    break;
  case 37: 
    return Tk_Tk_Bank::Instance().getCovmatrix15(tk);
    break;
  default: 
    throw runtime_error("tk attribute not valid");
  } 
}

FCALLSCFUN2(INT,ireadtktkcpp,IREADTKTKCPP,ireadtktkcpp, INT, INT)
int ireadtktkcpp(int attribute, int tk)
{

  tk = tk - 1;

  if(tk>Tk_Tk_Bank::Instance().size()) return 0;

  switch (attribute) {
  case 1: 
    return (int)Tk_Tk_Bank::Instance().getMod_ID(tk);
    break;
  case 2: 
    return (int)Tk_Tk_Bank::Instance().getSubdetbits(tk);
     break;
  case 3: 
    return (int)Tk_Tk_Bank::Instance().getMeasurement_code(tk);
    break;
  case 4: 
    return (int)Tk_Tk_Bank::Instance().getTrackType(tk);
    break;
  case 5: 
    return (int)Tk_Tk_Bank::Instance().getNumTEs(tk);
    break;
  case 6:
    return (int)Tk_Tk_Bank::Instance().getCharge(tk);
    break;
  case 7: 
    return (int)Tk_Tk_Bank::Instance().getUnused(tk);
    break;
  case 8: 
    return (int)Tk_Tk_Bank::Instance().getNdf(tk);
    break;
  case 9: 
    return (int)Tk_Tk_Bank::Instance().getChi2(tk);
    break;
  case 10:
    return (int)Tk_Tk_Bank::Instance().getLength(tk);
    break;
  case 11: 
    return (int)Tk_Tk_Bank::Instance().getStartX(tk);
    break;
  case 12: 
    return (int)Tk_Tk_Bank::Instance().getStartY(tk);
    break;
  case 13: 
    return (int)Tk_Tk_Bank::Instance().getStartZ(tk);
    break;
  case 14: 
    return (int)Tk_Tk_Bank::Instance().getEndX(tk);
    break;
  case 15: 
    return (int)Tk_Tk_Bank::Instance().getEndY(tk);
    break;
  case 16: 
    return (int)Tk_Tk_Bank::Instance().getEndZ(tk);
    break;
  case 17: 
    return (int)Tk_Tk_Bank::Instance().getCoord1_of_ref_point(tk);
    break;
  case 18: 
    return (int)Tk_Tk_Bank::Instance().getCoord2_of_ref_point(tk);
    break;
  case 19: 
    return (int)Tk_Tk_Bank::Instance().getCoord3_of_ref_point(tk);
    break;
  case 20: 
    return (int)Tk_Tk_Bank::Instance().getTheta(tk);
    break;
  case 21: 
    return (int)Tk_Tk_Bank::Instance().getPhi(tk);
    break;
  case 22: 
    return (int)Tk_Tk_Bank::Instance().getInvp(tk);
    break;
  case 23: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix1(tk);
    break;
  case 24: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix2(tk);
    break;
  case 25: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix3(tk);
    break;
  case 26: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix4(tk);
    break;
  case 27: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix5(tk);
    break;
  case 28: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix6(tk);
    break;
  case 29: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix7(tk);
    break;
  case 30: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix8(tk);
    break;
  case 31: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix9(tk);
    break;
  case 32: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix10(tk);
    break;
  case 33: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix11(tk);
    break;
  case 34: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix12(tk);
    break;
  case 35: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix13(tk);
    break;
  case 36: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix14(tk);
    break;
  case 37: 
    return (int)Tk_Tk_Bank::Instance().getCovmatrix15(tk);
    break;
  default: 
    throw runtime_error("tk attribute not valid");
  } 
}

FCALLSCFUN3(INT,writetktkcpp,WRITETKTKCPP,writetktkcpp, FLOAT, INT, INT)
int writetktkcpp(float value, int attribute, int tk){

  if(tk>Tk_Tk_Bank::Instance().size()) return 1;

  tk = tk - 1;

  switch (attribute) {
  case 1: 
    Tk_Tk_Bank::Instance().setMod_ID((int)value,tk);
    return 0;
    break;
  case 2: 
    Tk_Tk_Bank::Instance().setSubdetbits((int)value,tk);
    return 0;
    break;
  case 3: 
    Tk_Tk_Bank::Instance().setMeasurement_code((int)value,tk);
    return 0;
    break;
  case 4: 
    Tk_Tk_Bank::Instance().setTrackType((int)value,tk);
    return 0;
    break;
  case 5: 
    Tk_Tk_Bank::Instance().setNumTEs((int)value,tk);
    return 0;
    break;
  case 6: 
    Tk_Tk_Bank::Instance().setCharge((int)value,tk);
    return 0;
    break;
  case 7: 
    Tk_Tk_Bank::Instance().setUnused((int)value,tk);
    return 0;
    break;
  case 8: 
    Tk_Tk_Bank::Instance().setNdf((int)value,tk);
    return 0;
    break;
  case 9: 
    Tk_Tk_Bank::Instance().setChi2(value,tk);
    return 0;
    break;
  case 10: 
    Tk_Tk_Bank::Instance().setLength(value,tk);
    return 0;
    break;
  case 11:
    Tk_Tk_Bank::Instance().setStartX(value,tk);
    return 0;
    break;
  case 12:
    Tk_Tk_Bank::Instance().setStartY(value,tk);
    return 0;
    break;
  case 13:
    Tk_Tk_Bank::Instance().setStartZ(value,tk);
    return 0;
    break;
  case 14:
    Tk_Tk_Bank::Instance().setEndX(value,tk);
    return 0;
    break;
  case 15:
    Tk_Tk_Bank::Instance().setEndY(value,tk);
    return 0;
    break;
  case 16:
    Tk_Tk_Bank::Instance().setEndZ(value,tk);
    return 0;
    break;
  case 17: 
    Tk_Tk_Bank::Instance().setCoord1_of_ref_point(value,tk);
    return 0;
    break;
  case 18: 
    Tk_Tk_Bank::Instance().setCoord2_of_ref_point(value,tk);
    return 0;
    break;
  case 19: 
    Tk_Tk_Bank::Instance().setCoord3_of_ref_point(value,tk);
    return 0;
    break;
  case 20: 
    Tk_Tk_Bank::Instance().setTheta(value,tk);
    return 0;
    break;
  case 21: 
    Tk_Tk_Bank::Instance().setPhi(value,tk);
    return 0;
    break;
  case 22: 
    Tk_Tk_Bank::Instance().setInvp(value,tk);
    return 0;
    break;
  case 23: 
    Tk_Tk_Bank::Instance().setCovmatrix1(value,tk);
    return 0;
    break;
  case 24: 
    Tk_Tk_Bank::Instance().setCovmatrix2(value,tk);
    return 0;
    break;
  case 25: 
    Tk_Tk_Bank::Instance().setCovmatrix3(value,tk);
    return 0;
    break; 
  case 26: 
    Tk_Tk_Bank::Instance().setCovmatrix4(value,tk);
    return 0;
    break;
  case 27: 
    Tk_Tk_Bank::Instance().setCovmatrix5(value,tk);
    return 0;
    break;
  case 28: 
    Tk_Tk_Bank::Instance().setCovmatrix6(value,tk);
    return 0;
    break;
  case 29: 
    Tk_Tk_Bank::Instance().setCovmatrix7(value,tk);
    return 0;
    break;
  case 30: 
    Tk_Tk_Bank::Instance().setCovmatrix8(value,tk);
    return 0;
    break;
  case 31: 
    Tk_Tk_Bank::Instance().setCovmatrix9(value,tk);
    return 0;
    break;
  case 32: 
    Tk_Tk_Bank::Instance().setCovmatrix10(value,tk);
    return 0;
    break;
  case 33: 
    Tk_Tk_Bank::Instance().setCovmatrix11(value,tk);
    return 0;
    break;
  case 34: 
    Tk_Tk_Bank::Instance().setCovmatrix12(value,tk);
    return 0;
    break;
  case 35: 
    Tk_Tk_Bank::Instance().setCovmatrix13(value,tk);
    return 0;
    break;
  case 36: 
    Tk_Tk_Bank::Instance().setCovmatrix14(value,tk);
    return 0;
    break;
  case 37: 
    Tk_Tk_Bank::Instance().setCovmatrix15(value,tk);
    return 0;
    break;
  default: 
    cout << "attribute = " << attribute << endl;  
    throw runtime_error("tk attribute not valid");
  } 
  
}

FCALLSCFUN2(INT,addtetktkcpp,ADDTETKTKCPP,addtetktkcpp, INT, INT)
int addtetktkcpp(int te , int tk) 
{
  tk = tk - 1;
  te = te - 1;

  Tk_Tk_Bank::Instance().addTE(te ,tk);
  return 0;
}

FCALLSCFUN3(INT,writetkitkdatcpp,WRITETKITKDATCPP,writetkitkdatcpp, INT, INT, INT)
int writetkitkdatcpp(int value, int attribute, int tk){
  tk = tk - 1;
  
  switch (attribute) {
  case 1: 
    Tk_Tk_Bank::Instance().setPosOfFirstTEInTEList(value,tk);
    return 0;
    break;
  case 2: 
    Tk_Tk_Bank::Instance().setNumOfTEs(value,tk);
    return 0;
    break;
  case 3: 
    Tk_Tk_Bank::Instance().setTrackNo(value,tk);
    return 0;
    break;
  default: 
    std::cout << "attribute = " << attribute <<  std::endl;
    throw runtime_error("writetktkitkdatcpp: tk attribute not valid");
  } 

}

FCALLSCFUN2(INT,readtkitkdatcpp,READTKITKDATCPP,readtkitkdatcpp, INT, INT)
int readtkitkdatcpp(int attribute, int tk)
{

  tk = tk - 1;

  if(tk>Tk_Tk_Bank::Instance().size()) return 0;

  switch (attribute) {
  case 1: 
    return Tk_Tk_Bank::Instance().getPosOfFirstTEInTEList(tk);
    break;
  case 2: 
    return Tk_Tk_Bank::Instance().getNumOfTEs(tk);
    break;
  case 3: 
    return Tk_Tk_Bank::Instance().getTrackNo(tk);
    break;
  default:
    std::cout << "attribute = " << attribute << std::endl ;
    std::cout << "tk = " << tk << std::endl ;
    throw runtime_error("tk attribute not valid");
  } 

}
