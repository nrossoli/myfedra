#ifndef ROOT_EdbMosaicIO
#define ROOT_EdbMosaicIO

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbMosaicIO                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TFile.h>
#include "EdbID.h"
#include "EdbPattern.h"
#include "EdbLayer.h"

//______________________________________________________________________________
class EdbMosaicIO : public TObject {
  public:

    TFile *eFile;
    
  public:
   EdbMosaicIO(){ eFile=0; }
   virtual ~EdbMosaicIO(){ if(eFile) eFile->Close(); }
    
    void Init(const char *file, Option_t* option = "");
    void SaveFragment(EdbPattern &p);
    EdbPattern *GetFragment( int plate, int side, int id, bool do_corr );
    
    void SaveFragmentObj(TObject *ob, int plate, int side, int id, const char *pref);
    void SaveSideObj(TObject *ob, int plate, int side, const char *pref);

    void SaveCorrMap( int plate, int side, EdbLayer &l, const char *file);
    void SaveCorrMap( int plate, int side, EdbLayer &l );
    EdbLayer *GetCorrMap( int plate, int side );
    
    char *FileName(int brick, int plate, int major, int minor, const char *pref="", const char *suff="");

    void DrawFragment(EdbPattern &p);

    ClassDef(EdbMosaicIO,1)  //Mosaic IO
};

#endif /* ROOT_EdbMosaicIO */
