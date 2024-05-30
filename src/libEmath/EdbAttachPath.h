#ifndef ROOT_EdbAttachPath
#define ROOT_EdbAttachPath
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbAttachPath                                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <TArrayF.h>
#include <TArrayI.h>
//______________________________________________________________________________
class EdbAttachPath : public TObject {
  private:
   Int_t     eN0;
   TArrayF   eX;         // input points [eN0]
   TArrayF   eY;         //
   TArrayI   eID;        // external id's of input elements [eN0]
   TArrayI   eOK;        // alignment success levels 0,1,2 etc       [eN0]
 
   Float_t   eX0,eY0;    // starting point
   Float_t   eR;         // visinity range for neighbours search
  
   TArrayF   eDist;      // distances from headers to (X0,Y0) (eN elements useful)
   TArrayI   eInd;       // ordered index (eN elements useful)
   Int_t     eN;         // number of selected elements

    
  public:
   EdbAttachPath( int n ): eX(n),eY(n),eID(n),eOK(n),eDist(n),eInd(n) {eN0=n; eN=0;}
   virtual ~EdbAttachPath(){}

   void SetStartingPosition(const int i) {eX0=eX[i]; eY0=eY[i];}
   void SetStartingPosition(const float x, const float y);
   void SetStartingAtCenter();
   void SetRange(const float r) {eR=r;}
   void SetPoint(const int i, const float x, const float y, const int id);
   void OrderPointsRadial();                           // here define eN, fill eInd and eDist
   void RiseOK( const int i) { eOK[i]++; }
 
   int     I(const int i)    const { return eInd[i];  }
   int     N0()              const { return eN0;      }
   int     N()               const { return eN;       }
   int     OK(const int i)   const { return eOK[i];   }
   int     ID(const int i)   const { return eID[i];   }
   Float_t Dist(const int i) const { return eDist[i]; }
   int     GetAlignedNeighbours(const int i0, TArrayI &list) const;
   float   GetMin(TArrayF &a);
   float   GetMax(TArrayF &a);
   float   Xmin(){return GetMin( eX ); }
   float   Xmax(){return GetMax( eX ); }
   float   Ymin(){return GetMin( eY ); }
   float   Ymax(){return GetMax( eY ); }

  ClassDef(EdbAttachPath,1)  //Mosaic path manager 
};

#endif /* ROOT_EdbAttachPath */
