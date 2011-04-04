#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tkmcbank.h"

using namespace std;

Tk_MC_Bank::~Tk_MC_Bank()
{
}

void Tk_MC_Bank::clear(){
  mc_bank.clear();
}

void Tk_MC_Bank::add_mc(float PX, float PY, float PZ, float E, float X, float Y, float Z, int GEANT_PID, int MCTRACK, int NHITS, int HEPEVT_NUM)

{

  tk_mc amc;

  amc.px = PX;
  amc.py = PY;
  amc.pz = PZ;
  amc.energy = E;
  amc.x = X;
  amc.y = Y;
  amc.z = Z;
  amc.geant_pid = GEANT_PID; 
  amc.mctrack = MCTRACK;
  amc.nhits = NHITS;
  amc.hepevt_num = HEPEVT_NUM;

  mc_bank.push_back(amc);

}

// due to the fact that the bank will be accessed by integer value of the mc number this is probaly dangerous
// so leave it out for now

//void Tk_MC_Bank::remove_mc(int mc)
//{
//  mc_bank.erase(mc_bank.begin()+mc);
//}

float readtkmccpp(int attribute, int mc)
{

  if(mc>Tk_MC_Bank::Instance().size()) return 0.;

  mc = mc - 1;

  switch (attribute) {
  case 1: 
    //    std::cout << "getX = " << Tk_MC_Bank::Instance().getPX(mc) << std::endl;
    return Tk_MC_Bank::Instance().getPX(mc);
    break;
  case 2: 
    //    std::cout << "getY = " << Tk_MC_Bank::Instance().getPY(mc) << std::endl;
    return Tk_MC_Bank::Instance().getPY(mc);
     break;
  case 3: 
    //    std::cout << "getPZ = " << Tk_MC_Bank::Instance().getPZ(mc) << std::endl;
    return Tk_MC_Bank::Instance().getPZ(mc);
    break;
  case 4: 
    //    std::cout << "getEnergy = " << Tk_MC_Bank::Instance().getEnergy(mc) << std::endl;
    return Tk_MC_Bank::Instance().getEnergy(mc);
    break;
  case 5: 
    //    std::cout << "getX = " << Tk_MC_Bank::Instance().getPX(mc) << std::endl;
    return Tk_MC_Bank::Instance().getPX(mc);
    break;
  case 6: 
    //    std::cout << "getY = " << Tk_MC_Bank::Instance().getPY(mc) << std::endl;
    return Tk_MC_Bank::Instance().getPY(mc);
     break;
  case 7: 
    //    std::cout << "getPZ = " << Tk_MC_Bank::Instance().getPZ(mc) << std::endl;
    return Tk_MC_Bank::Instance().getPZ(mc);
    break;
  case 8: 
    //    std::cout << "getGEANT_PID = " << Tk_MC_Bank::Instance().getGEANT_PID(mc) << std::endl;
    return Tk_MC_Bank::Instance().getGEANT_PID(mc)+0.5;
    break;
  case 9: 
    //    std::cout << "getGEANT_TRACK_NUMBER = " << Tk_MC_Bank::Instance().getGEANT_TRACK_NUMBER(mc)<< std::endl;
    return Tk_MC_Bank::Instance().getMCTRACK(mc)+0.5;
    break;
  case 10: 
    //    std::cout << "getNHITS = " << Tk_MC_Bank::Instance().getNHITS(mc)<< std::endl;
    return Tk_MC_Bank::Instance().getNHITS(mc);
    break;
  case 11: 
    //    std::cout << "getHEPEVT_NUM = " << Tk_MC_Bank::Instance().getHEPEVT_NUM(mc) << std::endl;
    return Tk_MC_Bank::Instance().getHEPEVT_NUM(mc);
    break;
  default: 
    throw runtime_error("mc attribute not valid");
  } 
}

int writetkmccpp(float value, int attribute, int mc)
{

  mc = mc - 1;

  if(mc>Tk_MC_Bank::Instance().size()) return 1;

  switch (attribute) {
  case 1: 
    Tk_MC_Bank::Instance().setPX(value, mc);
    return 0;
    break;
  case 2: 
    Tk_MC_Bank::Instance().setPY(value, mc);
    return 0;
    break;
  case 3: 
    Tk_MC_Bank::Instance().setPZ(value, mc);
    return 0;
    break;
  case 4: 
    Tk_MC_Bank::Instance().setEnergy(value, mc);
    return 0;
    break;
  case 5: 
    Tk_MC_Bank::Instance().setX(value, mc);
    return 0;
    break;
  case 6: 
    Tk_MC_Bank::Instance().setY(value, mc);
    return 0;
    break;
  case 7: 
    Tk_MC_Bank::Instance().setZ(value, mc);
    return 0;
    break;
  case 8: 
    Tk_MC_Bank::Instance().setGEANT_PID(((int)value),mc);
    return 0;
    break;
  case 9: 
    Tk_MC_Bank::Instance().setMCTRACK(((int)value),mc);
    return 0;
    break;
  case 10: 
    Tk_MC_Bank::Instance().setNHITS(((int)value),mc);
    return 0;
    break;
  case 11: 
    Tk_MC_Bank::Instance().setHEPEVT_NUM(((int)value),mc);
    return 0;
    break;
  default: 
    throw runtime_error("mc attribute not valid");
  } 

}
