#ifndef MaterialDB_h
#define MaterialDB_h 1

extern "C" {
  extern struct { // size defined by NCMAMX in fkparm.inc
    int ncmat;
    float rcmat[100];
    float zcmax[100];
    float xrlc[100];
    int npmat;
    float zpmat[60];
    float rpmin[60];
    float rpmax[60];
    float xrlp[60];
    float xelosc[100];
    float xelosp[60];   
    float zcmin[100];
  } fkddes_; 
}

extern "C" {
  extern struct {  // size defined by NCMAMX in fkparm.inc
    int nconmat;
    float z1conmat[100];
    float z2conmat[100];
    float r1conmat[100];
    float r2conmat[100];
    float xrl1con[100];
    float xrl2con[100];
    float xel1con[100];
    float xel2con[100];
  } fkddes1_; 
}

extern "C" {
  extern struct {  // size defined by NPLMAMX in fkparm.inc
    int nplmat;          //total number of ladders
    float xplmat[200];   //X coordinate of the center of the ladder
    float yplmat[200];   //Y coordinate of the center of the ladder
    float zplmat[200];   //Z coordinate of the center of the ladder
    float widplmat[200]; //width of the ladder
    float lenplmat[200]; //length of the ladder (ladder parallel to z axis)
    float phiplmat[200]; //angle between normal vector of ladder and x axis
    float xrlpl[200];    //radiation length
    float xelospl[200];  //energy loss
  } fkddes2_; 
}

extern "C" {
  extern struct { // size defined by NEXSMX in fkparm.inc
    int nexs;
    float rzsurf[60];
    float zrmin[60];
    float zrmax[60];
    int itexts[60];
    int nexhpc;
  } fkexts_;

}

//extern "C" {
//  extern struct { // size defined by NEXSMX in fkparm.inc
//    char typx[60];
//  } fkexty_;
//
//}




extern "C" {
  extern struct {
    float consb;
  } fkfild_;
}

extern "C" { // only bfield is used anywhere in the F77 code 
  extern struct {
    float bfield;
    float bfldrtn;
    float rclin;
    float rclout;
    float bfieldq;
    float clthsc;
  } coildims_;
}

#endif
