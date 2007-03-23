#ifndef ROOT_EdbTestAl
#define ROOT_EdbTestAl

#include "TVector3.h"
#include "TH2F.h"
#include "EdbPattern.h"
#include "EdbAffine.h"

class EdbTestAl : public TObject
{
public:
  //EdbPattern *eP1;
  //EdbPattern *eP2;

  Int_t eITMAX;   // angular step (def=50)
  Int_t eOCMAX;   // occupancy (def=100)
  
  float eOffset;
  float eBinSize; //microns

  TObjArray *eS1, *eS2;   // pointers to segments selected by HDistance

  TTree *eT;
  TFile *eFile;
  TH2F  *HD;
  TH2F  *HDF;
  TH2F  *HDF2;

  Int_t    eN[4];                          // nx,ny,nz,nphi - number of divisions
  Float_t  eDmin[4], eDmax[4], eD0[4];     // limits, and found value of the peak for dx,dy,dz
  //TVector3 eVmin,  eVmax,  eVbin,  eV0;     
  //Float_t  ePHImin,ePHImax,ePHIbin,ePHI0;   // the same for the rotation

public:
  EdbTestAl();
  virtual ~EdbTestAl();

  void  HDistance(EdbPattern &p1, EdbPattern &p2);
  int   FillTree(float gdz=0);
  int   MakeTrans(EdbAffine2D &aff, float dz=0, const char *ccut="abs(dy)<400&&abs(dx)<400");
  int   CheckMaxBin(float dz, float phi, float &meanbin);
  int   CheckMaxBin();

  ClassDef(EdbTestAl,1)  // alignment test class
};

#endif /* ROOT_EdbTestAl */
