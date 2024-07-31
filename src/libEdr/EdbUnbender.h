#ifndef ROOT_EdbUnbender
#define ROOT_EdbUnbender
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbUnbender                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <TEnv.h>
#include "EdbPVRec.h"
#include "EdbScanSet.h"
#include "EdbCouplesTree.h"

//______________________________________________________________________________
class EdbUnbender : public TObject {
  private:
    EdbPVRec   *eV;   // volume to be unbended (assumed to be with the selected tracks only)
    EdbScanSet *eSS;  // scanset with layers to be modified
    TEnv       *eEnv;

    EdbCouplesTree eCPT; // debug output
  public:
    bool do_save_hist;
    bool do_save_tree;
    
  public:
   EdbUnbender(){ Set0(); }
   virtual ~EdbUnbender(){}

   void Set0();
   void Unbend3a(EdbPVRec &pvr, EdbScanSet &ss, TEnv &env);
   void Unbend3g(EdbPVRec &pvr, EdbScanSet &ss, TEnv &env);
   void Unbend5g(EdbPVRec &pvr, EdbScanSet &ss, TEnv &env);
   void CalcMeanPath( int cycle );
   void GlobalCorr3( int offset, int step, int cycle );

   void GlobalCorrN( int length, int position, int offset, int step, int flag );
   void CalculateCorrections(const TObjArray &acorr, int length, EdbPattern &p, EdbPlateP &plate, int flag);
   void CheckResolutions(TObjArray &a1, TObjArray &a2, int plate, int flag);
   
  ClassDef(EdbUnbender,1)  // unbend tracks group
};

#endif /* ROOT_EdbUnbender */
