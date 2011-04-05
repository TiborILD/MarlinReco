#ifndef MaterialDB_F77_h
#define MaterialDB_F77_h

#include <string>
#include <map>
#include <exception>

namespace marlin_delphiF77{

  class MaterialDB_F77exception: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "MaterialDB_F77exception occurred";
    }
  } ;
  
  class MaterialDB_F77 
  {

  public:

    static MaterialDB_F77* Instance() ;
    
    // Other non-static member functions
  
  public:
  
    ~MaterialDB_F77();   

    void switchOFFMaterial() ;
    void switchONMaterial() ;

  private:

    MaterialDB_F77():            // Private constructor 
    _useMaterials(false), _Ncmat(0), _Npmat(0), _Nconmat(0), _Nexs(0), _Nplmat(0)
    { } ;

    MaterialDB_F77(const MaterialDB_F77&) ;                 // Prevent copy-construction
    MaterialDB_F77& operator=(const MaterialDB_F77&) ;      // Prevent assignment

    void initialise() ;  

    void checkCommonBlocks() ;

    bool buildBeamPipe() ;
    bool buildVXD();
    bool buildFTD();
    bool buildSIT();
    bool buildSET();
    bool buildTPC();

    // private memeber variables
    static MaterialDB_F77* _pInstance ;  

    bool _useMaterials;

    int _Ncmat ;
    int _Npmat ;
    int _Nconmat ;
    int _Nexs ;
    int _Nplmat ;
  
  };

} // end of marlin_delphiF77 namespace 

#endif
