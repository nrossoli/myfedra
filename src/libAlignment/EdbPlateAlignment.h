#ifndef ROOT_EdbPlateAlignment
#define ROOT_EdbPlateAlignment

#include "EdbAlignmentV.h"

//-------------------------------------------------------------------------------------------------
class EdbPlateAlignment : public EdbAlignmentV
{
 public:
  
  Float_t   eSigma[2];     // sigma of the bt useful for the fine alignment ie:(10,0.01)
  Float_t   eOffsetMax;    // the maximal offset to be looked for

  Float_t   eDZ;           // the range +- dz   will be scanned by coarce align
  Float_t   eDPHI;         // the range +- dphi will be scanned by coarce align

  Bool_t   eDoTestAl, eTestAlOK;
  Bool_t   eDoCoarse, eCoarseOK;
  Bool_t   eDoFine,   eFineOK;
  Bool_t   eSaveCouples;          // save couples tree with the report file
  Bool_t   eStatus;               // overall alignment result (true - OK)
  Int_t    eNcoins;               // final number of coinsidence used for affine calc

  Int_t    eFineMin;              // minimum coinsidences to accept alignment
  Int_t    eCoarseMin;            // minimum coinsidences to accept alignment

  EdbPeak2  eH_zphi_coarse;   // the results of the coarse alignment
  EdbPeak2  eH_xy_coarse;
  EdbPeak2  eH_xy_final;      // the final alignment peak

 public:
  EdbPlateAlignment();
  virtual ~EdbPlateAlignment();

  void Align(EdbPattern &p1, EdbPattern &p2, float dz );
  void TestAl(EdbPattern &p1, EdbPattern &p2);
  void CoarseAl(EdbPattern &p1, EdbPattern &p2);
  void FineAl(EdbPattern &p1, EdbPattern &p2);
  void FineAlAff(EdbPattern &p1, EdbPattern &p2, EdbLayer &la1);
  void DoubletsFilterOut(EdbPattern &p1, EdbPattern &p2);

  void SetParTestAl( float zcorr, float dz=500, float dphi=0.03 );
  void SetParCoarseAl( float zcorr, float dpos=300, float dang=0.015, float dz=122, float dphi=0.01 );
  void SetParFineAl();
  void ProduceReport();
  void SaveCouplesTree();
  void SetSigma(float spos, float sang) { eSigma[0]=spos; eSigma[1]=sang; }
  
  void SlowAlignXY(EdbPattern &p1, EdbPattern &p2, EdbH2 &hxy, EdbH1 &hphi, const char *name="slowal" );

  ClassDef(EdbPlateAlignment,1)  // plate-to-plate alignment
};
#endif /* ROOT_EdbAlignmentV */
