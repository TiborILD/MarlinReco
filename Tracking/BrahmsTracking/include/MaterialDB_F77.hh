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
    static MaterialDB_F77& Instance()
    {
      static MaterialDB_F77 singleton;
      return singleton;
    }
    
    // Other non-static member functions
  
  public:
  
    ~MaterialDB_F77();   

    void initialise(bool withMSOn) ;
    bool isInitialise() ;

  private:
    MaterialDB_F77(): 
      _isInitialised(false), _Ncmat(0), _Npmat(0), _Nconmat(0), _Nexs(0), _Nplmat(0)
    { } ;           // Private constructor
    MaterialDB_F77(const MaterialDB_F77&) ;                 // Prevent copy-construction
    MaterialDB_F77& operator=(const MaterialDB_F77&) ;      // Prevent assignment
  
    void finaliseCommonBlocks() ;
    bool buildBeamPipe() ;
    bool buildVXD();
    bool buildSIT();
    bool buildTPC();

    // private memeber variables
  
    bool _isInitialised;
    bool _useMaterials;

    int _Ncmat ;
    int _Npmat ;
    int _Nconmat ;
    int _Nexs ;
    int _Nplmat ;
  
  };

} // end of marlin_kaltest namespace 

#endif
