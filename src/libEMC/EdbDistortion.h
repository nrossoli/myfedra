#ifndef ROOT_EdbDistortion
#define ROOT_EdbDistortion
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDistortion - use overlapped frames as "cormtx"                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TObject.h"
#include "TFile.h"
#include "TEnv.h"
#include "TNtuple.h"
#include "EdbCell2.h"
#include "EdbRun.h"

//______________________________________________________________________________
class EdbClDist : public TObject {
public:
  Float_t  eX,eY,eZ;  // cluster coords in view RS
  Float_t  eXv,eYv;   // view coords
  Int_t    eView;
  Int_t    eFrame;
  Short_t  eIsCenter;   
  
  EdbClDist(){};
  EdbClDist(float x,float y,float z, float xv,float yv, int view, int frame, short i, short j)
  { eX=x; eY=y; eZ=z; eXv=xv; eYv=yv; eView=view; eFrame=frame;}
  virtual ~EdbClDist(){}
  
  ClassDef(EdbClDist,1)  // service structure for views correction
};

//______________________________________________________________________________
class EdbDistortionMap : public TObject
{
private:
  float eXpix,  eYpix;     // camera pixel size in microns
  int   eNXpix, eNYpix;    // camera matrix size in pixels
  EdbCell2   eMap;     // view divided into cells [microns space]
  
public:
  EdbDistortionMap(){}
  virtual ~EdbDistortionMap(){}
  
  void InitMap(int nxpix,int nypix, float xpix, float ypix, float stepX, float stepY);
  void SetDX( int j, double dx );
  void SetDY( int j, double dy );
  void AddDX( int j, double dx );
  void AddDY( int j, double dy );
  void Fill(float x, float y, float dx, float dy);
  void Norm();
  double DX(int j) const;
  double DY(int j) const;
  int Bin(int j) const {return eMap.Bin(j);}
  int Jcell(float x, float y) const {return eMap.Jcell(x,y);}
  
  float StepX() {return eMap.Xbin();}
  float StepY() {return eMap.Ybin();}
  
  void Add(const EdbDistortionMap &map, float k=1.);
  void Substract(const EdbDistortionMap &map) { Add(map,-1.); }
  void Scale(const float k);
  
  TH2F *GetH2dX(const char *name="hdx");
  TH2F *GetH2dY(const char *name="hdy");
  void PutH2dX( const TH2F &h2 );
  void PutH2dY( const TH2F &h2 );
  void Smooth(int n=1, Option_t *opt="k5a");
  
  void DrawCorrMap(TFile *file=0, const char *name=0);
  void ReadMatrix2Map(const char *file);
  void GenerateCorrectionMatrix(const char *file);
  void ApplyCorr(EdbClDist &c);
  void Save(TFile *f, const char *suffix="");
  
  ClassDef(EdbDistortionMap,1)  // distortion map for the views correction  
};

//______________________________________________________________________________
class EdbDistortionCorr : public TObject {
private:
  int        eNClMin;      // minimal number of clusters inside grain to be used for corrections
  float      eR2CenterMax; // the maximal distance to the matrix center for the reference cluster
  int        eNCenterMin;  // the minimal number of cluster appearance in the central region
  float      eRmax;        // acceptance for clusters matching
  
  TClonesArray  eCl;       // array of EdbClDist objects
  TClonesArray  eGr;       // array of EdbSegment objects (cluster "clouds" in this case)
  EdbCell2   eGMap;        // map of grains (EdbSegments)
  EdbDistortionMap  eCorrMap;     // map of corrections (differential)
  EdbDistortionMap  eCorrMap0;    // input map of corrections
  EdbDistortionMap  eCorrMapTot;    // eCorrMap+eCorrMap0
  
  float eXpix,  eYpix;     // pixel size in microns
  int   eNXpix, eNYpix;    // camera matrix size in pixels
  int   eDirX, eDirY;      // microscope steps directions {-1, 0, 1}
  
  float  eCorrectionMatrixStepX;
  float  eCorrectionMatrixStepY;
  
public:
  TFile *eOutputFile;
  bool   eDumpGr;                     // if(1) dump grains tree
  float  eAreaMin, eAreaMax;          // clusters Area limits
  float  eVolumeMin, eVolumeMax;      // clusters Volume limits
  
public:
  EdbDistortionCorr();
  virtual ~EdbDistortionCorr(){}
  
  void InitGMap();
  void InitCorrMap();
  
  void ReadClusters( const char *fname, TClonesArray &arr );
  void AddCluster( EdbClDist *c );
  
  void CalculateGrRef();
  void CalculateGrRefMean();
  void CalculateCorr();
  void SaveCorrMap();
  void GenerateCorrectionMatrix(bool do_add);
  void SetPar(TEnv &env);
  void SetPixelSize(float xpix, float ypix) { eXpix=xpix; eYpix=ypix; }
  void CutGrRef();
  
  TNtuple *DumpGr(const char *name);
  
  void MakeDistortionMap(const char *fname, TEnv &env, const char *usefile=0, const char *addfile=0);
  int MakeViewDirectionList( EdbRun &run, int dirx, int diry, TArrayI &entries );
  
  void Print();
  
  ClassDef(EdbDistortionCorr,1)  // 
};

#endif /* ROOT_EdbDistortion */
