#ifndef ROOT_EdbMosaicPath
#define ROOT_EdbMosaicPath

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbMosaicPath                                                            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TFile.h>
#include <TEnv.h>
#include "EdbID.h"
#include "EdbMosaicIO.h"
#include "EdbPattern.h"
#include "EdbLayer.h"
#include "EdbViewMap.h"
#include "EdbRunAccess.h"

//______________________________________________________________________________
class EdbMosaicPath : public TObject {
  private:
   Int_t     eN0;
   TObjArray eHarr;         // full array of a view headers (length eN0)
   TArrayI   eOK;           // alignment success 0,1        (length eN0)
   TArrayF   eDist;         // distances from headers to (X0,Y0) (eN elements useful)
   TArrayI   eInd;          // order (eN elements useful)
   Int_t     eN;            // number of selected elements
  public:
   Float_t   eX0,eY0;       // fragment center
   Float_t   eR0;           // visinity range
   
  public:
   EdbMosaicPath( int n ): eHarr(n),eOK(n),eDist(n),eInd(n) {eN0=n;}
   virtual ~EdbMosaicPath(){}

   EdbViewHeader *GetHeader(int i)  const { return (EdbViewHeader *)(eHarr.At(i)); }
   EdbViewHeader *FindNearest( const TObjArray &harr, const float x0, const float y0 );
   void    InitArea(const TObjArray &harr, const float x0, const float y0, int nsegmin);
   float   InitLineX(const TObjArray &harr, const float y0, const float dy);
   float   InitLineY(const TObjArray &harr, const float x0, const float dx);
   void    SetOK( const int i) { eOK[i]++; }
   int     GetAlignedNeighbours(const int i0, TArrayI &list) const;
   int     I(const int i)    const { return eInd[i]; }
   int     N0()        const { return eN0; }
   int     N()         const { return eN;  }
   int     OK(const int i)   const { return eOK[i]; }
   Float_t Dist(const int i) const { return eDist[i]; }
  
  ClassDef(EdbMosaicPath,1)  //Mosaic path manager 
};

#endif /* ROOT_EdbMosaicPath */
