//-- Author :  Valeri Tioukov   19.05.2002
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbPattern                                                           //
//                                                                      //
// Segments pattern                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TIndexCell.h"
#include "EdbAffine.h"
#include "EdbPattern.h"

ClassImp(EdbSegP)
ClassImp(EdbSegmentsBox)
ClassImp(EdbTrackP)
ClassImp(EdbPattern)
ClassImp(EdbPatternsVolume)

//______________________________________________________________________________
EdbSegP::EdbSegP() 
{
  ePID=0;
  eID=0;
  eVid[0]=0;
  eVid[1]=0;
  eFlag=0;
  eProb=0;
  eW=0;
}

//______________________________________________________________________________
EdbSegP& EdbSegP::operator += (EdbSegP const& s) 
{
 
  float dz = s.Z() - Z();

  float tx,ty;
  float w1,w2;
  float x1,y1,x2,y2;
  float xa,ya,sxa,sya,xb,yb,sxb,syb;

  // project to one side (a):

  x1 = X();
  w1 = 1./(SX()*SX());
  x2 = s.X() - s.TX()*dz;
  w2 = 1./(s.SX()*s.SX() + s.STX()*dz*s.STX()*dz); 
  xa = (x1*w1+x2*w2)/(w1+w2);
  sxa = TMath::Sqrt(1./(w1+w2));

  y1 = Y();
  w1 = 1./(SY()*SY());
  y2 = s.Y() - s.TY()*dz;
  w2 = 1./(s.SY()*s.SY() + s.STY()*dz*s.STY()*dz); 
  ya = (y1*w1+y2*w2)/(w1+w2);
  sya = TMath::Sqrt(1./(w1+w2));

  // project to another side (b):

  x1 = s.X();
  w1 = 1./(s.SX()*s.SX());
  x2 = X() + TX()*dz;
  w2 = 1./(SX()*SX() + STX()*dz*STX()*dz); 
  xb = (x1*w1+x2*w2)/(w1+w2);
  sxb = TMath::Sqrt(1./(w1+w2));
  
  y1 = s.Y();
  w1 = 1./(s.SY()*s.SY());
  y2 = Y() + TY()*dz;
  w2 = 1./(SY()*SY() + STY()*dz*STY()*dz); 
  yb = (y1*w1+y2*w2)/(w1+w2);
  syb = TMath::Sqrt(1./(w1+w2));

  float z  = (Z()+s.Z())/2.;
  tx = (xb-xa)/dz;
  ty = (yb-ya)/dz;
  float x  = (xb+xa)/2.;
  float y  = (yb+ya)/2.;

  Set(s.ID(),x,y,tx,ty,W()+s.W(),s.Flag());
  SetZ(z);
  //SetErrors(sx,sy,sz,stx,sty);

  return *this;
}

///______________________________________________________________________________
void EdbSegP::LinkMT(const EdbSegP* s1,const EdbSegP* s2, EdbSegP* s)
{
  /// Segments fit by Andrey Aleksandrov (Jul-2003)

  register Double_t dz = s2->Z() - s1->Z();
  Double_t dz2 = dz*dz;
 
  Double_t q1,q2,w1,w2;
  Double_t d1,d2,dxx11,dxx22;
  Double_t dtt01,dtt02,dtx01,dtx02;
  Double_t dxx01,dxx02,dxt01,dxt02;
  Double_t xm1,xm2,sx0,sy0,stx0,sty0;
 
  register Double_t q;

  if(dz==0.0) {
    s->SetZ(s1->Z());
    s->SetID(s1->ID());
      
    q1 = s1->SX()*s1->SX();
    q2 = s2->SX()*s2->SX();
    w1 = s1->STX()*s1->STX();
    w2 = s2->STX()*s2->STX();
    
    sx0 = q1*q2/(q1+q2);
    q = (s1->X()/q1+s2->X()/q2)*sx0;
    s->SetX(q);
    stx0 = w1*w2/(w1+w2);
    q = (s1->TX()/w1+s2->TX()/w2)*stx0;
    s->SetTX(q);
 
    q1 = s1->SY()*s1->SY();
    q2 = s2->SY()*s2->SY();
    w1 = s1->STY()*s1->STY();
    w2 = s2->STY()*s2->STY();
 
    sy0 = q1*q2/(q1+q2);
    q = (s1->Y()/q1+s2->Y()/q2)*sy0;
    s->SetY(q);
    sty0 = w1*w2/(w1+w2);
    q = (s1->TY()/w1+s2->TY()/w2)*sty0;
    s->SetTY(q);
 
    s->SetErrors(TMath::Sqrt(sx0),TMath::Sqrt(sy0),0.0,TMath::Sqrt(stx0),TMath::Sqrt(sty0));
    s->SetW( s1->W()+s2->W() );
    return;
  }

  q = 0.5*(s1->Z()+s2->Z());
  register Double_t dzr = 1.0/dz;
 
  s->SetZ(q);
  s->SetID(s1->ID());
 
  q1 = s1->SX()*s1->SX();
  q2 = s2->SX()*s2->SX();
  w1 = s1->STX()*s1->STX();
  w2 = s2->STX()*s2->STX();
 
  q = dz2*w2+q2;
  d1 = 1.0/(q+q1);
  xm1 = (q*s1->X()+(s2->X()-dz*s2->TX())*q1)*d1;
 
  q = dz2*w1+q1;
  d2 = 1.0/(q+q2);
  xm2 = (q*s2->X()+(s1->X()+dz*s1->TX())*q2)*d2;


  dtt01 = q2*d2;
  dtt02 = q1*d1;
  dxx11 = 1.0-dtt02;
  dxx22 = 1.0-dtt01;
  dxx01 = 0.5*(dxx11+dtt01);
  dxx02 = 0.5*(dxx22+dtt02);
  dxt01 = 0.5*dz*dtt01;
  dxt02 = -0.5*dz*dtt02;
  dtx01 = dzr*(dtt01-dxx11);
  dtx02 = dzr*(dxx22-dtt02);
 
  q = (xm1+xm2)*0.5;
  s->SetX(q);
  q = (xm2-xm1)*dzr;
  s->SetTX(q);
  sx0 = TMath::Sqrt(dxx01*dxx01*q1+dxx02*dxx02*q2+dxt01*dxt01*w1+dxt02*dxt02*w2);
  stx0 = TMath::Sqrt(dtx01*dtx01*q1+dtx02*dtx02*q2+dtt01*dtt01*w1+dtt02*dtt02*w2);
 
  q1 = s1->SY()*s1->SY();
  q2 = s2->SY()*s2->SY();
  w1 = s1->STY()*s1->STY();
  w2 = s2->STY()*s2->STY();
 
  q = dz2*w2+q2;
  d1 = 1.0/(q+q1);
  xm1 = (q*s1->Y()+(s2->Y()-dz*s2->TY())*q1)*d1;
 
  q = dz2*w1+q1;
  d2 = 1.0/(q+q2);
  xm2 = (q*s2->Y()+(s1->Y()+dz*s1->TY())*q2)*d2;

  dtt01 = q2*d2;
  dtt02 = q1*d1;
  dxx11 = 1.0-dtt02;
  dxx22 = 1.0-dtt01;
  dxx01 = 0.5*(dxx11+dtt01);
  dxx02 = 0.5*(dxx22+dtt02);
  dxt01 = 0.5*dz*dtt01;
  dxt02 = -0.5*dz*dtt02;
  dtx01 = dzr*(dtt01-dxx11);
  dtx02 = dzr*(dxx22-dtt02);
 
  q = (xm1+xm2)*0.5;
  s->SetY(q);
  q = (xm2-xm1)*dzr;
  s->SetTY(q);
  sy0 = TMath::Sqrt(dxx01*dxx01*q1+dxx02*dxx02*q2+dxt01*dxt01*w1+dxt02*dxt02*w2);
  sty0 = TMath::Sqrt(dtx01*dtx01*q1+dtx02*dtx02*q2+dtt01*dtt01*w1+dtt02*dtt02*w2);
 
  s->SetErrors(sx0,sy0,0.0,stx0,sty0);
  s->SetW( s1->W()+s2->W() );
}

//______________________________________________________________________________
double EdbSegP::Chi2( EdbSegP &s ) const
{
  //calculated at the s.Z();

  double dz = s.Z()-Z();
  double x1  = X() + dz * TX();
  double y1  = Y() + dz * TY();

  double sx = TMath::Sqrt( SX()*SX() + dz*STX()*dz*STX() );
  double sy = TMath::Sqrt( SY()*SY() + dz*STY()*dz*STY() );

  double stx = TMath::Sqrt( STX()*STX() + s.STX()*s.STX() );
  double sty = TMath::Sqrt( STY()*STY() + s.STY()*s.STY() );

  double dx  = (s.X()-x1)/sx;
  double dy  = (s.Y()-y1)/sy;
  double dtx = (s.TX()-TX())/stx;
  double dty = (s.TY()-TY())/sty;

  return TMath::Sqrt(dx*dx + dy*dy + dtx*dtx + dty*dty)/2.;
}

//______________________________________________________________________________
float EdbSegP::Chi2A( EdbSegP &s ) const
{
  // ignore the position errors here;

  float dz  = s.Z()-Z();
  float tx  = (s.X()-X())/dz;
  float ty  = (s.Y()-Y())/dz;

  float dtx1 = (TX()-tx)*(TX()-tx)/STX()/STX();
  float dty1 = (TY()-ty)*(TY()-ty)/STY()/STY();
  float dtx2 = (s.TX()-tx)*(s.TX()-tx)/s.STX()/s.STX();
  float dty2 = (s.TY()-ty)*(s.TY()-ty)/s.STY()/s.STY();

  return TMath::Sqrt(dtx1+dty1+dtx2+dty2)/2.;
}

//______________________________________________________________________________
float EdbSegP::Chi2Aprob( EdbSegP &s ) const
{
  // ignore the position errors here;
  // use segments probablility information

  float dz  = s.Z()-Z();
  float tx  = (s.X()-X())/dz;
  float ty  = (s.Y()-Y())/dz;

  float dtx1 = (TX()-tx)*(TX()-tx)/STX()/STX()/Prob();
  float dty1 = (TY()-ty)*(TY()-ty)/STY()/STY()/Prob();
  float dtx2 = (s.TX()-tx)*(s.TX()-tx)/s.STX()/s.STX()/s.Prob();
  float dty2 = (s.TY()-ty)*(s.TY()-ty)/s.STY()/s.STY()/s.Prob();

  return TMath::Sqrt(dtx1+dty1+dtx2+dty2)/2.;
}

//______________________________________________________________________________
void EdbSegP::PropagateTo( float z ) 
{
  float dz = z-Z();

  eX  = X() + TX()*dz;
  eY  = Y() + TY()*dz;
  eZ  = z;
  eSX = TMath::Sqrt( SX()*SX() + STX()*dz*STX()*dz );
  eSY = TMath::Sqrt( SY()*SY() + STY()*dz*STY()*dz );
}

//______________________________________________________________________________
float EdbSegP::ProbLink( EdbSegP &s1, EdbSegP &s2 )
{
  // return probability of the correct link in case Up/Down - specifics is 
  // that the position errors are neglected)

  double dz = s2.Z() - s1.Z();

  double tx = (s2.X() - s1.X()) / dz;
  double ty = (s2.Y() - s1.Y()) / dz;

  double dtx1 = (s1.TX()-tx)/s1.STX();
  double dty1 = (s1.TY()-ty)/s1.STY();
  double dtx2 = (s2.TX()-tx)/s2.STX();
  double dty2 = (s2.TY()-ty)/s2.STY();

  double chi2 = TMath::Sqrt(dtx1*dtx1 + dty1*dty1 + dtx2*dtx2 + dty2*dty2);
  double p3 = TMath::Prob(chi2,4);

  return s1.Prob()*s2.Prob()*p3;
}

//______________________________________________________________________________
void EdbSegP::MergeTo( EdbSegP &s ) 
{
  // create linked segment at Z of s2 
  // TODO - navesti nauku covariantnuiu 
  
  PropagateTo( s.Z() );

  float wx1,wx2, wy1,wy2;
  float wtx1,wtx2, wty1,wty2;

  wx1 = 1/SX()/SX();
  wx2 = 1/s.SX()/s.SX();
  wy1 = 1/SY()/SY();
  wy2 = 1/s.SY()/s.SY();
  wtx1 = 1/STX()/STX();
  wtx2 = 1/s.STX()/s.STX();
  wty1 = 1/STY()/STY();
  wty2 = 1/s.STY()/s.STY();

  eX = (X()*wx1 + s.X()*wx2)/(wx1+wx2);
  eY = (Y()*wy1 + s.Y()*wy2)/(wy1+wy2);
  eSX = TMath::Sqrt( 1./(wx1+wx2) );
  eSY = TMath::Sqrt( 1./(wy1+wy2) );

  eTX = (TX()*wtx1 + s.TX()*wtx2)/(wtx1+wtx2);
  eTY = (TY()*wty1 + s.TY()*wty2)/(wty1+wty2);
  eSTX = TMath::Sqrt( 1./(wtx1+wtx2) );
  eSTY = TMath::Sqrt( 1./(wty1+wty2) );

  eZ = s.Z();
  eSZ = TMath::Sqrt(( SZ()*SZ() + s.SZ()*s.SZ())/2);

  eW = W()+s.W();

  eID   = s.ID();
  ePID  = s.PID();
  eFlag = s.Flag();
}

//______________________________________________________________________________
void EdbSegP::Print(Option_t *opt) const
{
  printf("EdbSegP: %d  %f %f %f  %f %f  %f  %d\n", ID(),X(),Y(),Z(),TX(),TY(),W(),Flag() );
} 

//______________________________________________________________________________
EdbSegmentsBox::EdbSegmentsBox()
{
  eSegments = new TClonesArray("EdbSegP");
  Set0();

}
 
//______________________________________________________________________________
EdbSegmentsBox::EdbSegmentsBox(float x0, float y0, float z0)
{
  eSegments = new TClonesArray("EdbSegP");
  eX=x0;  eY=y0;  eZ=z0;
  eDZkeep=0;
}
 
//______________________________________________________________________________
EdbSegmentsBox::~EdbSegmentsBox( )
{
  if(eSegments) delete eSegments;
}
 
//______________________________________________________________________________
void EdbSegmentsBox::Set0()
{
  eX=0;  eY=0;  eZ=0;
  eDZkeep=0;
}
 
//______________________________________________________________________________
void EdbSegmentsBox::AddSegment(int id, float x, float y, float tx, float ty, 
			    float w, int flag)
{
  new((*eSegments)[N()])  EdbSegP( id,x,y,tx,ty,w,flag );
}
 
//______________________________________________________________________________
void EdbSegmentsBox::AddSegment( EdbSegP &s )
{
  new((*eSegments)[N()])  EdbSegP( s );
}
 
//______________________________________________________________________________
void EdbSegmentsBox::AddSegment( EdbSegP &s1, EdbSegP &s2 )
{
  EdbSegP *s = new((*eSegments)[N()])  EdbSegP( s1 );
  s->MergeTo(s2);
}
 
//______________________________________________________________________________
void EdbSegmentsBox::Reset()
{
  Set0();
  if(eSegments) eSegments->Clear();
}
 
//______________________________________________________________________________
float EdbSegmentsBox::DiffAff(EdbAffine2D *aff)
{
  EdbSegmentsBox p,p1;
  p.AddSegment(0,Xmin(),Ymin(),0.,0.);
  p.AddSegment(0,Xmax(),Ymin(),0.,0.);
  p.AddSegment(0,Xmin(),Ymax(),0.,0.);
  p.AddSegment(0,Xmax(),Ymax(),0.,0.);

  p1.AddSegment(0,Xmin(),Ymin(),0.,0.);
  p1.AddSegment(0,Xmax(),Ymin(),0.,0.);
  p1.AddSegment(0,Xmin(),Ymax(),0.,0.);
  p1.AddSegment(0,Xmax(),Ymax(),0.,0.);

  p1.Transform(aff);
  return p.Diff(p1);
}
 
//______________________________________________________________________________
float EdbSegmentsBox::Diff(EdbSegmentsBox &p)
{
  // return the mean difference beteween pattern elements
  int nseg = TMath::Min( N(), p.N() );
  if(nseg<1) return 0;

  EdbSegP *s1=0, *s2=0;

  float dx=0, dy=0;
  double sdx=0;
  for(int i=0; i<nseg; i++ ) {
    s1 =   GetSegment(i);
    s2 = p.GetSegment(i);
    dx = s2->X() - s1->X();
    dy = s2->Y() - s1->Y();
    sdx += TMath::Sqrt( dx*dx + dy*dy);
  }
  return sdx/nseg;
}
 
//______________________________________________________________________________
void EdbSegmentsBox::ProjectTo(const float dz)
{
  eZ += dz;  eDZkeep += dz;
 
  EdbSegP *p;
  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);
    p->SetX( p->X() + p->TX()*dz );
    p->SetY( p->Y() + p->TY()*dz );
  }
}
 
//______________________________________________________________________________
int EdbSegmentsBox::CalculateXY( EdbSegmentsBox *pat, EdbAffine2D *aff )
{
  int n=N();
  if( n>pat->N() )  n=pat->N();
  if(n<2) return 0;

  aff->Calculate(this,pat);
  return 1;
}

//______________________________________________________________________________
int EdbSegmentsBox::CalculateAXAY( EdbSegmentsBox *pat, EdbAffine2D *aff )
{
  int n=N();
  if( n>pat->N() )  n=pat->N();
  if(n<2) return 0;

  float *ax1 = new float[n];
  float *ay1 = new float[n];
  float *ax2 = new float[n];
  float *ay2 = new float[n];

  EdbSegP *p1,*p2;
  for(int i=0; i<n; i++ ) {
    p1 = GetSegment(i);
    p2 = pat->GetSegment(i);
    ax1[i] = p1->TX();
    ay1[i] = p1->TY();
    ax2[i] = p2->TX();
    ay2[i] = p2->TY();
  }
  aff->Calculate(n,ax1,ay1,ax2,ay2);

  delete ax1;
  delete ay1;
  delete ax2;
  delete ay2;

  return 1;
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformA( const EdbAffine2D *aff )
{
  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = aff->A11()*p->TX() + aff->A12()*p->TY() + aff->B1();
    ty = aff->A21()*p->TX() + aff->A22()*p->TY() + aff->B2();
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformARot( const EdbAffine2D *aff )
{
  // apply to the angles only rotation members of transformation

  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = aff->A11()*p->TX() + aff->A12()*p->TY();
    ty = aff->A21()*p->TX() + aff->A22()*p->TY();
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::Print(Option_t *opt) const
{
  int nseg=GetN();
  printf("EdbSegmentsBox: %d segments\n", nseg );
  for(int i=0; i<nseg; i++) GetSegment(i)->Print();
} 

//______________________________________________________________________________
//______________________________________________________________________________
EdbTrackP::EdbTrackP()
{
  eID = 0;
  eVid = 0;
}
 
//______________________________________________________________________________
EdbTrackP::~EdbTrackP()
{

}

//______________________________________________________________________________
void EdbTrackP::Copy(EdbTrackP &tr)
{
  Reset();
  eID = tr.ID();
  int nseg=tr.N();
  for(int i=0; i<nseg; i++)
    AddSegment(*tr.GetSegment(i));
}

//______________________________________________________________________________
//______________________________________________________________________________
EdbPattern::EdbPattern()
{
  eCell     = new TIndexCell();
  Set0();
}

//______________________________________________________________________________
EdbPattern::EdbPattern(float x0, float y0, float z0) : EdbSegmentsBox(x0,y0,z0) 
{
  eCell     = new TIndexCell();
  Set0();
}
 
//______________________________________________________________________________
EdbPattern::~EdbPattern( )
{
  if(eCell)     delete eCell;
}

//______________________________________________________________________________
void EdbPattern::Set0()
{
  eID = 0;
}

//______________________________________________________________________________
void EdbPattern::Reset()
{
  ((EdbSegmentsBox *)this)->Reset();
  Set0();
  if(eCell)     eCell->Drop();
}
 
////////////////////////////////////////////////////////////////////////////////
//
//   EdbPatternsVolume
//
////////////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
EdbPatternsVolume::EdbPatternsVolume()
{
  ePatterns   = new TObjArray();
  eTracksCell = 0;
  eTracks     = new TObjArray();
  eSV  = 0;
  Set0();
}

//______________________________________________________________________________
EdbPatternsVolume::EdbPatternsVolume(EdbPatternsVolume &pvol)
{
  ePatterns   = new TObjArray();
  eTracksCell = 0;
  eTracks     = new TObjArray();
  eSV  = 0;
  Set0();

  pvol.PassProperties(*this);
  int npat,nseg;
  npat = Npatterns();
  for(int j=0; j<npat; j++) {
    nseg = GetPattern(j)->N();
    for(int i=0; i<nseg; i++ ) {
      pvol.GetPattern(j)->AddSegment( *(GetPattern(j)->GetSegment(i)) );
    }
  }
}

//______________________________________________________________________________
EdbPatternsVolume::~EdbPatternsVolume()
{
  if(ePatterns) {
    ePatterns->Delete();
    delete ePatterns;
  }
}

//______________________________________________________________________________
void EdbPatternsVolume::Set0()
{

}
 
//______________________________________________________________________________
int EdbPatternsVolume::DropCouples()
{
  int count=0;
  int npat=Npatterns();
  for(int i=0; i<npat; i++ )
    count += GetPattern(i)->Cell()->DropCouples(4);
  if(count) printf("%d couples are dropped in volume cells\n",count);
  return count;
}
 
//______________________________________________________________________________
void EdbPatternsVolume::SetPatternsID()
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ )
    GetPattern(i)->SetID(i);
}
 
//______________________________________________________________________________
void EdbPatternsVolume::Transform( const EdbAffine2D *aff )
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ )  {
    GetPattern(i)->Transform(aff);
    GetPattern(i)->TransformARot(aff);
  }
}
 
//______________________________________________________________________________
void EdbPatternsVolume::Centralize()
{
  // find geometrical center (XY) of all patterns and set it as the center of 
  // coordinates  to simplify transformations
  // To be used before any operations on patterns

  float xc=0;
  float yc=0;
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    xc += GetPattern(i)->Xmax() + GetPattern(i)->Xmin();
    yc += GetPattern(i)->Ymax() + GetPattern(i)->Ymin();
  }
  xc = xc/Npatterns()/2;
  yc = yc/Npatterns()/2;

  eX = xc;  eY=yc;

  Shift(-xc,-yc);
  npat = Npatterns();
  for(int i=0; i<npat; i++ ) 
    GetPattern(i)->SetKeep(1,0,0,1,0,0);
}

//______________________________________________________________________________
void EdbPatternsVolume::PrintAff() const
{
  EdbAffine2D a;
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    GetPattern(i)->GetKeep(a);
    printf(" %d ",i); a.Print();
  }
}

//______________________________________________________________________________
void EdbPatternsVolume::PrintStat( Option_t *opt="") const
{
  int npat = Npatterns();
  printf("\nVolume statistics for %d patterns\n",npat);

  float dx,dy;
  EdbPattern *pat=0;
  printf("pat# \t segments \t dX \t\tdY \t meanDist \n");
  for(int i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    dx = pat->Xmax() - pat->Xmin();
    dy = pat->Ymax()- pat->Ymin();
    printf(" %d\t %d\t %10.2f \t %10.2f \t %10.4f \n", 
	   i, pat->GetN(),dx,dy, TMath::Sqrt(dx*dy/pat->GetN()) );
  }

  npat=Npatterns();
  for(int i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    pat->Cell()->PrintStat();
  }
}
 
//______________________________________________________________________________
void EdbPatternsVolume::PrintStat(EdbPattern &pat) const
{
} 

//______________________________________________________________________________
void EdbPatternsVolume::AddPattern( EdbPattern *pat )
{
  ePatterns->Add(pat);
}

//______________________________________________________________________________
EdbPattern *EdbPatternsVolume::GetPattern( int id ) const
{
  if(Npatterns()>id) return (EdbPattern*)ePatterns->At(id);
  else return 0;
}


//______________________________________________________________________________
void EdbPatternsVolume::PassProperties( EdbPatternsVolume &pvol )
{
  pvol.SetXYZ(eX,eY,eZ);
  EdbAffine2D a;
  EdbPattern *p=0;

  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    p = GetPattern(i);
    p->GetKeep(a);
    EdbPattern *psel = new EdbPattern( p->X(),p->Y(),p->Z() );
    psel->SetKeep( a.A11(), a.A12(), a.A21(), a.A22(), a.B1(),a.B2() );
    pvol.AddPattern(psel);
  }
  pvol.SetPatternsID();
}


//______________________________________________________________________________
void EdbPatternsVolume::Shift( float x, float y )
{
  EdbAffine2D aff;
  aff.ShiftX(x);
  aff.ShiftY(y);

  int npat = Npatterns();
  for(int i=0; i<npat; i++)
    GetPattern(i)->Transform(&aff);
}

//______________________________________________________________________________
void EdbPatternsVolume::DropCell()
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) GetPattern(i)->Cell()->Drop();
}
