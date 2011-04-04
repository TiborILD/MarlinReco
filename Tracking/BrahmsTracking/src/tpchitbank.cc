#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tpchitbank.h"

using namespace std;


// Declaration of Common Block used to store the number of hits in the hit bank defined in TPChitbank.h
//CNTPC_DEF CNTPC;


TPC_Hit_Bank::~TPC_Hit_Bank()
{
}

void TPC_Hit_Bank::clear()
{
  hit_bank.clear();
}

//void TPC_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID, int TrkID, int PntToEx, int NEx, int ResC, float Res1, float Res2)
void TPC_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID, float Res1, float Res2, int TrkID)
{

  TPC_hit ahit;

  ahit.x = X;
  ahit.y = Y;
  ahit.z = Z;
  ahit.energy = E;
  ahit.subdetector_ID = SubID;
  //  ahit.pointer_to_first_exclusion = PntToEx;
  //  ahit.number_of_exclusions = NEx;
  //  ahit.resolution_code = ResC;
  ahit.resolution_1 = Res1;
  ahit.resolution_2 = Res2;
  ahit.track_ID = TrkID;

  hit_bank.push_back(ahit);

}

// due to the fact that the bank will be accessed by integer value of the hit number this is probaly dangerous
// so leave it out for now

//void TPC_Hit_Bank::remove_hit(int hit)
//{
//  hit_bank.erase(hit_bank.begin()+hit);
//}

float readtpchitscpp(int attribute, int hit)
{

  if(hit>TPC_Hit_Bank::Instance().size()) return 0.;

  hit = hit - 1;

  switch (attribute) {
  case 1: 
    return TPC_Hit_Bank::Instance().getX(hit);
    break;
  case 2: 
    return TPC_Hit_Bank::Instance().getY(hit);
     break;
  case 3: 
    return TPC_Hit_Bank::Instance().getZ(hit);
    break;
  case 4: 
    return TPC_Hit_Bank::Instance().getEnergy(hit);
    break;
  case 5: 
    return TPC_Hit_Bank::Instance().getSubdetectorID(hit);
    break;
  case 6: 
    return TPC_Hit_Bank::Instance().getResolution1(hit);
    break;
  case 7: 
    return TPC_Hit_Bank::Instance().getResolution2(hit);
    break;
  case 8: 
    return TPC_Hit_Bank::Instance().getTrackID(hit);
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 
}
FCALLSCFUN2(FLOAT,readtpchitscpp,READTPCHITSCPP,readtpchitscpp, INT, INT)

int writetpccpp(float value, int attribute, int hit)
{

  hit = hit - 1;

  if(hit>TPC_Hit_Bank::Instance().size()) return 1;

  switch (attribute) {
  case 1: 
    TPC_Hit_Bank::Instance().setX(value, hit);
    return 0;
    break;
  case 2: 
    TPC_Hit_Bank::Instance().setY(value, hit);
    return 0;
    break;
  case 3: 
    TPC_Hit_Bank::Instance().setZ(value, hit);
    return 0;
    break;
  case 4: 
    TPC_Hit_Bank::Instance().setEnergy(value, hit);
    return 0;
    break;
  case 5: 
    TPC_Hit_Bank::Instance().setSubdetectorID(((int)value),hit);
    return 0;
    break;
  case 6: 
    TPC_Hit_Bank::Instance().setResolution1(value,hit);
    return 0;
    break;
  case 7: 
    TPC_Hit_Bank::Instance().setResolution2(value,hit);
    return 0; 
    break;
  case 8: 
    TPC_Hit_Bank::Instance().setTrackID(((int)value),hit);
    return 0;
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 

}
FCALLSCFUN3(INT,writetpccpp,WRITETPCCPP,writetpccpp, FLOAT, INT, INT)
