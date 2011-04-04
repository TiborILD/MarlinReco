#include <iostream>
#include <vector>
#include <string>
//stl exception handler
#include<stdexcept>

#include "tkhitbank.h"

using namespace std;

// Global pointer to the Tk_Hit_Bank structure
Tk_Hit_Bank * TkHitBank;

// global functions needed for the F77/C++ interface

int subdetfirsthitindex(string subdet);

int numofsubdethits(string subdet);

float rreadtkhitscpp(int a, int b); 

int ireadtkhitscpp(int a, int b); 

int writetkhitcpp(float c, int a, int b); 


Tk_Hit_Bank::~Tk_Hit_Bank()
{
}

void Tk_Hit_Bank::clear()
{
  hit_bank.clear();
  firstHitEntry.clear();
  lastHitEntry.clear();
}

void Tk_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID,  int TrkID, int PntToEx, int NEx, int ResC, float Res1, float Res2)
  //void Tk_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID, float Res1, float Res2, int TrkID)
{

  tk_hit ahit;

  ahit.x = X;
  ahit.y = Y;
  ahit.z = Z;
  ahit.energy = E;
  ahit.subdetector_ID = SubID;
  ahit.track_ID = TrkID;
  ahit.pointer_to_first_exclusion = PntToEx;
  ahit.number_of_exclusions = NEx;
  ahit.resolution_code = ResC;
  ahit.resolution_1 = Res1;
  ahit.resolution_2 = Res2;


  hit_bank.push_back(ahit);

}

// due to the fact that the bank will be accessed by integer value of the hit number this is probaly dangerous
// so leave it out for now

//void Tk_Hit_Bank::remove_hit(int hit)
//{
//  hit_bank.erase(hit_bank.begin()+hit);
//}





int subdetfirsthitindex(string subdet)
{
  return Tk_Hit_Bank::Instance().getFirstHitIndex( subdet )+1;
}
FCALLSCFUN1(INT,subdetfirsthitindex,SUBDETFIRSTHITINDEX,subdetfirsthitindex, STRING)

int numofsubdethits(string subdet)
{
  return Tk_Hit_Bank::Instance().getNumOfSubDetHits( subdet );
}
FCALLSCFUN1(INT,numofsubdethits,NUMOFSUBDETHITS,numofsubdethits, STRING)

float rreadtkhitscpp(int attribute, int hit)
{

  if(hit>Tk_Hit_Bank::Instance().size()) return 0.;

  hit = hit - 1;

  switch (attribute) {
  case 1: 
    //    std::cout << "getX = " << Tk_Hit_Bank::Instance().getX(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getX(hit);
    break;
  case 2: 
    //    std::cout << "getY = " << Tk_Hit_Bank::Instance().getY(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getY(hit);
     break;
  case 3: 
    //    std::cout << "getZ = " << Tk_Hit_Bank::Instance().getZ(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getZ(hit);
    break;
  case 4: 
    //    std::cout << "getEnergy = " << Tk_Hit_Bank::Instance().getEnergy(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getEnergy(hit);
    break;
  case 5: 
    //    std::cout << "getSubdetectorID = " << Tk_Hit_Bank::Instance().getSubdetectorID(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getSubdetectorID(hit)+0.5;
    break;
  case 6: 
    //    std::cout << "getTrackID = " << Tk_Hit_Bank::Instance().getTrackID(hit)<< std::endl;
    return Tk_Hit_Bank::Instance().getTrackID(hit)+0.5;
    break;
  case 7: 
    //    std::cout << "getPntToFirstExclusion = " << Tk_Hit_Bank::Instance().getPntToFirstExclusion(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getPntToFirstExclusion(hit)+0.5;
    break;
  case 8: 
    //    std::cout << "Tk_Hit_Bank::Instance().getNExclusion = " << Tk_Hit_Bank::Instance().getNExclusion(hit)<< std::endl;
    return Tk_Hit_Bank::Instance().getNExclusion(hit)+0.5;
    break;
  case 9: 
    //    std::cout << "getResolutionCode = " << Tk_Hit_Bank::Instance().getResolutionCode(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getResolutionCode(hit)+0.5;
    break;
  case 10: 
    //    std::cout << "getResolution1 = " << Tk_Hit_Bank::Instance().getResolution1(hit)<< std::endl;
    return Tk_Hit_Bank::Instance().getResolution1(hit);
    break;
  case 11: 
    //    std::cout << "getResolution2 = " << Tk_Hit_Bank::Instance().getResolution2(hit) << std::endl;
    return Tk_Hit_Bank::Instance().getResolution2(hit);
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 
}
FCALLSCFUN2(FLOAT,rreadtkhitscpp,RREADTKHITSCPP,rreadtkhitscpp, INT, INT)

int ireadtkhitscpp(int attribute, int hit)
{

  if(hit>Tk_Hit_Bank::Instance().size()) return 0;

  hit = hit - 1;

  switch (attribute) {
  case 1: 
    //    std::cout << "getX = " << Tk_Hit_Bank::Instance().getX(hit) << std::endl;
    return (int)Tk_Hit_Bank::Instance().getX(hit);
    break;
  case 2: 
    //    std::cout << "getY = " << Tk_Hit_Bank::Instance().getY(hit) << std::endl;
    return (int)Tk_Hit_Bank::Instance().getY(hit);
     break;
  case 3: 
    //    std::cout << "getZ = " << Tk_Hit_Bank::Instance().getZ(hit) << std::endl;
    return (int)Tk_Hit_Bank::Instance().getZ(hit);
    break;
  case 4: 
    //    std::cout << "getEnergy = " << Tk_Hit_Bank::Instance().getEnergy(hit) << std::endl;
    return (int)Tk_Hit_Bank::Instance().getEnergy(hit);
    break;
  case 5: 
    //    std::cout << "getSubdetectorID = " << Tk_Hit_Bank::Instance().getSubdetectorID(hit) << std::endl;
    return (int)(Tk_Hit_Bank::Instance().getSubdetectorID(hit));
    break;
  case 6: 
    //    std::cout << "getTrackID = " << Tk_Hit_Bank::Instance().getTrackID(hit)<< std::endl;
    return (int)Tk_Hit_Bank::Instance().getTrackID(hit);
    break;
  case 7: 
    //    std::cout << "getPntToFirstExclusion = " << Tk_Hit_Bank::Instance().getPntToFirstExclusion(hit) << std::endl;
    return (int)Tk_Hit_Bank::Instance().getPntToFirstExclusion(hit);
    break;
  case 8: 
    //    std::cout << "Tk_Hit_Bank::Instance().getNExclusion = " << Tk_Hit_Bank::Instance().getNExclusion(hit)<< std::endl;
    return (int)Tk_Hit_Bank::Instance().getNExclusion(hit);
    break;
  case 9: 
    //    std::cout << "getResolutionCode = " << Tk_Hit_Bank::Instance().getResolutionCode(hit) << std::endl;
    return (int)Tk_Hit_Bank::Instance().getResolutionCode(hit);
    break;
  case 10: 
    //    std::cout << "getResolution1 = " << Tk_Hit_Bank::Instance().getResolution1(hit)<< std::endl;
    return (int)Tk_Hit_Bank::Instance().getResolution1(hit);
    break;
  case 11: 
    //    std::cout << "getResolution2 = " << Tk_Hit_Bank::Instance().getResolution2(hit) << std::endl;
    return (int)Tk_Hit_Bank::Instance().getResolution2(hit);
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 
}
FCALLSCFUN2(INT,ireadtkhitscpp,IREADTKHITSCPP,ireadtkhitscpp, INT, INT)

int writetkhitcpp(float value, int attribute, int hit)
{

  hit = hit - 1;

  if(hit>Tk_Hit_Bank::Instance().size()) return 1;

  switch (attribute) {
  case 1: 
    Tk_Hit_Bank::Instance().setX(value, hit);
    return 0;
    break;
  case 2: 
    Tk_Hit_Bank::Instance().setY(value, hit);
    return 0;
    break;
  case 3: 
    Tk_Hit_Bank::Instance().setZ(value, hit);
    return 0;
    break;
  case 4: 
    Tk_Hit_Bank::Instance().setEnergy(value, hit);
    return 0;
    break;
  case 5: 
    Tk_Hit_Bank::Instance().setSubdetectorID(((int)value),hit);
    return 0;
    break;
  case 6: 
    Tk_Hit_Bank::Instance().setTrackID(((int)value),hit);
    return 0;
    break;
  case 7: 
    Tk_Hit_Bank::Instance().setPntToFirstExclusion(((int)value),hit);
    return 0;
    break;
  case 8: 
    Tk_Hit_Bank::Instance().setNExclusion(((int)value),hit);
    return 0;
    break;
  case 9: 
    Tk_Hit_Bank::Instance().setResolutionCode(((int)value),hit);
    return 0;
    break;
  case 10: 
    Tk_Hit_Bank::Instance().setResolution1(value,hit);
    return 0;
    break;
  case 11: 
    Tk_Hit_Bank::Instance().setResolution2(value,hit);
    return 0;
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 

}
FCALLSCFUN3(INT,writetkhitcpp,WRITETKHITCPP,writetkhitcpp, FLOAT, INT, INT)

