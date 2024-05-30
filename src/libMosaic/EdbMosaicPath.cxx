//-- Author :  Valeri Tioukov   15/02/2024

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbMosaic                                                          //
//                                                                      //
//                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TCut.h>
#include "EdbMosaicIO.h"
#include "EdbMosaic.h"
#include "EdbPlateAlignment.h"
#include "EdbScanProc.h"
#include "EdbLog.h"

using namespace TMath;
ClassImp(EdbMosaicPath);

//-----------------------------------------------------------------------
int EdbMosaicPath::GetAlignedNeighbours(const int i0, TArrayI &list) const
{
  // for a given view return list of already aligned neighbours
  EdbViewHeader *h0 = (EdbViewHeader *)(eHarr.At(i0));
  int n=0;
  for(int i=0; i<eN0; i++)
  {
    if(eOK[i]) {
      EdbViewHeader *h=(EdbViewHeader *)(eHarr.At(i));
      float  dx = h->GetXview() - h0->GetXview();
      float  dy = h->GetYview() - h0->GetYview();
      float  d  = Sqrt(dx*dx+dy*dy);
      if( d<eR0 && d>0.1 ) { list[n] = i; n++; }
    }
  }
  return n;
}

//-----------------------------------------------------------------------
float EdbMosaicPath::InitLineX( const TObjArray &harr, const float y0, const float dy0 )
{
  eN0 = harr.GetEntries();
  int cnt=0;
  for(int i=0; i<eN0; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    eHarr.Add(h);
    if( Abs( h->GetYview()-y0 ) < dy0 ) {
      eDist[i] = h->GetXview();
      cnt++;
    } else    eDist[i] = kMaxLong;
  }
  eN  = cnt;
  TMath::Sort(eN0,eDist.GetArray(),eInd.GetArray(),0);
  eY0 = y0;
  eX0 = ((EdbViewHeader *)(eHarr.At(I(0))))->GetXview();
  return ((EdbViewHeader *)(eHarr.At(I(eN-1))))->GetXview() - eX0;
}

//-----------------------------------------------------------------------
float EdbMosaicPath::InitLineY( const TObjArray &harr, const float x0, const float dx0 )
{
  eN0 = harr.GetEntries();
  int cnt=0;
  for(int i=0; i<eN0; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    eHarr.Add(h);
    if( Abs( h->GetXview()-x0 ) < dx0 ) {
      eDist[i] = h->GetYview();
      cnt++;
    } else    eDist[i] = kMaxLong;    
  }
  eN  = cnt;
  TMath::Sort(eN0,eDist.GetArray(),eInd.GetArray(),0);
  eX0 = x0;
  eY0 = ((EdbViewHeader *)(eHarr.At(I(0))))->GetYview();
  return ((EdbViewHeader *)(eHarr.At(I(eN-1))))->GetYview() - eY0;
}

//-----------------------------------------------------------------------
EdbViewHeader *EdbMosaicPath::FindNearest( const TObjArray &harr, const float x0, const float y0 )
{
  int n = harr.GetEntries();
  float r2min=kMaxLong;
  EdbViewHeader *h0=0;
  for(int i=0; i<n; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    float dx = h->GetXview()-x0;
    float dy = h->GetYview()-y0;
    float r2 = dx*dx+dy*dy;
    if(r2<r2min)
    {
      r2min=r2;
      h0=h;
    }
  }
  return h0;
}

//-----------------------------------------------------------------------
void EdbMosaicPath::InitArea( const TObjArray &harr, const float x0, const float y0 )
{
  eX0 = x0;
  eY0 = y0;
  int n = harr.GetEntries();
  if(n>eN0) {Log(1,"EdbMosaicPath::InitArea","ERROR: array length %d > %d",n,eN0); return;}
  eN0=n;
  for(int i=0; i<eN0; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    eHarr.Add(h);
    float dx = h->GetXview()-eX0;
    float dy = h->GetYview()-eY0;
    eDist[i] = Sqrt(dx*dx+dy*dy);
  }
  TMath::Sort(eN0,eDist.GetArray(),eInd.GetArray(),0);
}
