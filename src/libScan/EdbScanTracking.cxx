  //-- Author :  Valeri Tioukov   19/02/2011
  //////////////////////////////////////////////////////////////////////////
  //                                                                      //
  // EdbScanTracking                                                      //
  //                                                                      //
  // To handle tracking in the scanset                                    //
  //                                                                      //
  //////////////////////////////////////////////////////////////////////////
#include "EdbLog.h"
#include "EdbDataSet.h"
#include "EdbScanTracking.h"
#include "EdbAlignmentV.h"
#include "EdbPlateAlignment.h"
#include "EdbTrackFitter.h"
  
using namespace TMath;
  
ClassImp(EdbTrackAssembler)
ClassImp(EdbScanTracking)
  
//--------------------------------------------------------------------------------------
EdbTrackAssembler::~EdbTrackAssembler()
{
}
  
  //--------------------------------------------------------------------------------------
  EdbTrackAssembler::EdbTrackAssembler()
{
  eMapMarg = 50.; // [microns]
  eZ = 0;
  eCellN=10;    //mean n/cell
  eDTmax=0.07;
  eDRmax=45.;
  eDZGapMax = 5000;
  eCollisionsRate=0;
    
    // for basetracks:
  eCond.SetDefault();
  eCond.SetSigma0( 4, 4, 0.005, 0.005 );
  eCond.SetPulsRamp0(14., 21.);
  eCond.SetPulsRamp04(14., 21.);
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::DoubletsFilterOut(EdbPattern &p)
{
  EdbAlignmentV adup;
  adup.eDVsame[0]=adup.eDVsame[1]=10.;
  adup.eDVsame[2]=adup.eDVsame[3]=0.08;
  adup.FillGuessCell(p,p,1.);
  adup.FillCombinations();
  adup.DoubletsFilterOut(0);   // assign flag -10 to the duplicated segments
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::CheckPatternAlignment(EdbPattern &p, int nsegmin)
{
  ExtrapolateTracksToZ( p.Z(), nsegmin);
  int ntr = eTrZ.GetEntries();
  EdbPattern ptr( 0, 0, p.Z(), ntr ); 
  for(int i=0; i<ntr; i++) ptr.AddSegment( *((EdbSegP*)(eTrZ.UncheckedAt(i))) );
    
  EdbPlateAlignment al;
  al.SetSigma( 10, 0.005 );
  al.eOffsetMax = 500.;
  al.eDZ        = 0;
  al.eDPHI      = 0.00;
    //al.eDoCoarse=1;
  al.Align(ptr,p,0);
  EdbAffine2D *aff = al.eCorrL[0].GetAffineXY();
  aff->Print();
  for(int i=0; i<ntr; i++) ((EdbSegP*)(eTrZ.UncheckedAt(i)))->Transform(aff);
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::AddPattern(EdbPattern &p)
{
    //int ntrBefore=eTracks.GetEntries();
    //if(ntrBefore>0) ExtrapolateTracksToZ(p.Z());
    
    //DoubletsFilterOut();
  int nseg = p.N();
  int attached=0;
  for(int j=0; j<nseg; j++) {
    EdbSegP *s = p.GetSegment(j);
    if(s->Flag()==-10)   continue;
    if( !AddSegment( *(p.GetSegment(j)) ) )
      AddSegmentAsTrack( *(p.GetSegment(j)) );
    else attached++;
  }
  int ntrAfter = eTracks.GetEntries();
    //int totSegTr=0;
    //for(int i=0; i<ntrAfter; i++) totSegTr += ((EdbTrackP*)(eTracks.At(i)))->N();
  Log(2,"EdbTrackAssembler::AddPattern","with z=%10.2f   %d/%d attached/tried; total collisions: %d;  tracks: %d",
      p.Z(), attached, nseg, eCollisionsRate, ntrAfter );
}
  
  //--------------------------------------------------------------------------------------
  EdbTrackP *EdbTrackAssembler::AddSegment(EdbSegP &s)
{
  TObjArray trsel;
  float v[2] = { s.X(), s.Y() };
  int nsel =  eTrZMap.SelectObjectsC( v, eDRmax+50 , trsel ); 
  if(!nsel) { 
      //printf("nsel=%d\n",nsel); s.PrintNice(); 
    return 0;  }
    //printf("nsel=%d\n",nsel);
    float prob, probmax = 0.00000001;
    EdbSegP *ssbest = 0;
    for(int i=0; i<nsel; i++) {
      EdbSegP *ss = (EdbSegP*)(trsel.At(i));
      //printf("\nbefore: ");    s.PrintNice();
      prob = ProbSeg( *ss, s );
      //printf("prob probmax:  %f %f\n",  prob, probmax);
      if( prob > probmax ) { ssbest = ss; probmax=prob; }
    }
    if(!ssbest)  return 0;
    //printf("probmax = %10.7f \n", probmax);
    s.SetProb(probmax);
    EdbTrackP *t = (EdbTrackP*)(ssbest);
    EdbSegP *sz = t->GetSegmentWithClosestZ( t->Z(), 45. );
    //if(sz) printf("%f %f \n", s.Prob(),sz->Prob());
    if(!sz) t->AddSegment( eSegments.AddSegment(s) );
    else {
      if( !SameSegment(s,*sz) ) {
        if( s.Prob() > sz->Prob() )  t->SubstituteSegment( sz ,  eSegments.AddSegment(s) );
      //printf("%f %f \n", s.Prob(),sz->Prob());
      //s.PrintNice();
        eCollisionsRate++;
      }
    }
    return t;
}
  
  //--------------------------------------------------------------------------------------
  bool EdbTrackAssembler::SameSegment( EdbSegP &s1, EdbSegP &s2 )
{
  if( Abs( s1.X() - s2.X() )   <0.0001 &&
      Abs( s1.Y() - s2.Y() )   <0.0001 &&
      Abs( s1.TX()- s2.TX() )  <0.0001 &&
      Abs( s1.TY()- s2.TY() )  <0.0001 &&
      Abs( s1.W() - s2.W() )   <0.0001  )    return true;
  return false;
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::RecalculateSegmentsProb( EdbTrackP &tr )
{
    // assumed that track is fitted: reset the segmets probabilities
  int n=tr.N();
  for(int i=0; i<n; i++) 
    tr.GetSegment(i)->SetProb( ProbSeg( *(tr.GetSegmentF(i)),  *(tr.GetSegment(i)) ) ); 
}
  
  //--------------------------------------------------------------------------------------
  float EdbTrackAssembler::ProbSeg( EdbSegP &s1, EdbSegP &s2 )
{
  // return the probability that the second segment can belong to track defined by s1
  float dtx = s1.TX() - s2.TX();
  if( Abs( dtx ) > eDTmax )    return 0;
  float dty = s1.TY() - s2.TY();
  if( Abs( dty ) > eDTmax )    return 0;
  double dt2 = dtx*dtx +  dty*dty;
  if(dt2>eDTmax*eDTmax)        return 0;
  
  float dz = s2.Z()-s1.Z();
  float dx = s2.X() - (s1.X() + dz*s1.TX());
  if( Abs( dx ) > eDRmax )     return 0;
  float dy = s2.Y() - (s1.Y() + dz*s1.TY());
  if( Abs( dy ) > eDRmax )     return 0;
  double dr2 = dx*dx +  dy*dy;
  if(dr2>eDRmax*eDRmax)        return 0;
  
  EdbSegP s;
  float chi2 = eFitter.Chi2SegM( s1, s2, s, eCond, eCond );
    
  float prob = (float)TMath::Prob( chi2*chi2, 4 );
    
  prob *= eCond.ProbSeg( s2.Theta(), s2.W() );            // the proability component depending on the grains number
  prob *= (float)TMath::Prob( s2.Chi2()*s2.Chi2(), 4 );   // the proability component depending on the segment strength
    
  return prob;
}
  
  //--------------------------------------------------------------------------------------
  EdbTrackP *EdbTrackAssembler::AddSegmentAsTrack(EdbSegP &s)
{
  if(s.W()<16    )  return 0;
  if(s.Chi2()>2.5)  return 0;
  EdbTrackP *t = new EdbTrackP( eSegments.AddSegment(s), 0.139);    // EdbTrackAssembler is owner of segments 
  eTracks.Add(t);
    //EdbSegP ss;
    //t->MakePredictionTo(eZ,ss);
    //eTrZ.AddSegment(ss);
  return t;
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::ExtrapolateTracksToZ(float z, int nsegmin)
{
  eTrZ.Clear();
  eTrZMap.CleanCells();
    
  eZ=z;
  int n=eTracks.GetEntries();
  for(int i=0; i<n; i++) {
    EdbTrackP *t = (EdbTrackP*)(eTracks.At(i));
      
    if( t->N() < nsegmin )      continue;
    if( !AcceptDZGap(*t, z) )   continue;
      
    t->MakePredictionTo(eZ,*t);                        // extrapolation of tracks itself
      //((EdbSegP *)t)->PrintNice();
    eTrZ.Add(t);
  }
    
    //FillTrZMap();
    //if(gEDBDEBUGLEVEL>2) eTrZMap.PrintStat();
}
  
  //--------------------------------------------------------------------------------------
  bool EdbTrackAssembler::AcceptDZGap(EdbTrackP &t, float z)
{
  float z1 = t.GetSegmentFirst()->Z();
  float z2 = t.GetSegmentLast()->Z();
  if(Min( Abs(z1-z), Abs(z2-z)) > eDZGapMax ) return false;
  return true;
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::FillTrZMap()
{
  int n=eTrZ.GetEntries();
  Log(2,"EdbTrackAssembler::FillTrZMap", "with %d tracks",n);
  for(int i=0; i<n; i++) {
    EdbSegP *s = (EdbSegP*)(eTrZ.At(i));
    eTrZMap.AddObject( s->X(), s->Y(), s );
  }
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::InitTrZMap( const char *str )
{
  int   nx=0, ny=0, ncell=0;
  float xmi,xma, ymi, yma;
  sscanf(str,"%d %f %f %d %f %f %d",&nx,&xmi,&xma,&ny,&ymi,&yma,&ncell);
  InitTrZMap( nx,xmi,xma,ny,ymi,yma,ncell );
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::InitTrZMap( int nx, float xmi, float xma, 
                                      int ny, float ymi, float yma,  int ncell)
{
  float  mi[2] = {  xmi, ymi };
  float  ma[2] = {  xma, yma };
  int     n[2] = { nx, ny };
  eTrZMap.InitCell(ncell, n, mi, ma);
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::InitTrZMap()
{
  /*  float  mi[2] = {  eTrZ.Xmin()-eMapMarg, eTrZ.Ymin()-eMapMarg };
  float  ma[2] = {  eTrZ.Xmax()+eMapMarg, eTrZ.Ymax()+eMapMarg };
  float  dens = eTrZ.N()/( (ma[0]-mi[0])*(ma[1]-mi[1]));
  float  step=10000;
  if(dens>0.00000001) step = Sqrt( eCellN/dens );
  int    n[2] = { int((ma[0]-mi[0])/step)+1, int((ma[1]-mi[1])/step)+1 };
  float stepX = (ma[0]-mi[0])/n[0];
  float stepY = (ma[1]-mi[1])/n[1];
  n[0] = int((ma[0]-mi[0]+1.)/stepX);
  n[1] = int((ma[1]-mi[1]+1.)/stepY);
  eTrZMap.InitCell(3*eCellN, n, mi, ma);*/
}

//--------------------------------------------------------------------------------------
void EdbTrackAssembler::SetSegmentsErrors()
{
  int ntr = eTracks.GetEntries();
    for( int i=0; i<ntr; i++ )     {
      EdbTrackP *t = (EdbTrackP*)(eTracks.At(i));
      if(t->Flag()==-10) continue;
      int nseg=t->N();
      if(nseg>0)  { 
        for(int j=0; j<nseg; j++) {
          EdbSegP *s = t->GetSegment(j);
          s->SetErrors();
          eCond.FillErrorsCov(s->TX(),s->TY(),s->COV());
        }
      }
    }
}

//--------------------------------------------------------------------------------------
void EdbTrackAssembler::FitTracks()
{
  EdbTrackFitter fit;

  int ntr = eTracks.GetEntries();
  for( int i=0; i<ntr; i++ )     {
    EdbTrackP *t = (EdbTrackP*)(eTracks.At(i));
    if(t->Flag()==-10) continue;
    int nseg=t->N();
    t->SetP(1.);
    t->FitTrackKFS(0,10000);
        //fit.FitTrackLine(*t);
    if(nseg>1) RecalculateSegmentsProb(*t);
  }
}

//______________________________________________________________________________
void EdbTrackAssembler::CombTracks( TObjArray &selected )
{
    // eliminate crossing&overlapping tracks with multiple segments usage
  
  int nsegMin=2;
  int nGapMax=50;
  
  int ntr = eTracks.GetEntries();
  Log(3,"EdbTrackAssembler::CombTracks","Comb %d tracks");
  
    // *** sort tracks by quality
  
  TIndexCell cn;  //"nseg:prob:entry"
  Long_t v[3];
  
  int nsegtot=0;
  EdbTrackP *tr=0;
  for(int i=0; i<ntr; i++) {
    tr = (EdbTrackP*)(eTracks.At(i));
    if( tr->Flag() == -10 ) continue;
    if( tr->N() < nsegMin ) continue;
    tr->SetID(i);
    tr->SetCounters();
    nsegtot += tr->SetSegmentsTrack(-1);
    v[0]= -(tr->N());
    v[1]= (Long_t)((1.-tr->Prob())*100);
    v[2]= i;
    cn.Add(3,v);
  }
  cn.Sort();
  
  Log(3,"EdbTrackAssembler::CombTracks","%d tracks with %d segments for processing...",ntr,nsegtot);
  
    // *** set track ID for segments attached to
  
  TIndexCell *cp=0, *c=0;
  int nn=cn.GetEntries();
  for(int i=nn-1; i>=0; i--) {
    cp = cn.At(i);                              // tracks with fixed npl
    int np = cp->GetEntries();
    for(int ip=np-1; ip>=0; ip--) {
      c = cp->At(ip);                           // tracks with fixed Npl & Prob
      int nt = c->GetEntries();
      for(int it=0; it<nt; it++) {
        tr = (EdbTrackP*)(eTracks.At( c->At(it)->Value() ) );
        tr->SetSegmentsTrack();
      }
    }
  }
  
  
  cp=0;   c=0;
  nn=cn.GetEntries();
  for(int i=0; i<nn; i++) {
    cp = cn.At(i);                              // tracks with fixed npl

    int np = cp->GetEntries();
    for(int ip=0; ip<np; ip++) {
      c = cp->At(ip);                           // tracks with fixed Npl & Prob

      int nt = c->GetEntries();
      for(int it=0; it<nt; it++) {
  
        tr = (EdbTrackP*)(eTracks.At( c->At(it)->Value() ) );

        if(tr->RemoveAliasSegments()>0){
          if(tr->N()<nsegMin)             tr->SetFlag(-10);
          if(tr->CheckMaxGap()>nGapMax)   tr->SetFlag(-10);
        }

        if( tr->Flag() != -10 ) selected.Add(tr);
      }
    }
  }

}

//=======================================================================================
  EdbScanTracking::EdbScanTracking()
{
  eNsegMin=2;
  eNgapMax=50;
}

//--------------------------------------------------------------------------------------
void EdbScanTracking::TrackSetBT(EdbID idset, TEnv &env)
{
  
  // read scanset object
  EdbScanSet *ss = eSproc->ReadScanSet(idset);
  if(!ss) { Log(1,"EdbScanTracking::TrackSetBT",
    "Error! set for %s do not found",idset.AsString()); return; }
  
    int npl = ss->eIDS.GetEntries();
    if(npl<2) { Log(1,"EdbScanTracking::TrackSetBT", "Warning! npl<2 : %d stop tracking!",npl); return; }
  
  // create and init tracking object
    EdbTrackAssembler etra;
  
    etra.eCond.SetSigma0(        env.GetValue("fedra.track.Sigma0"         , "3 3 0.005 0.005") );
    etra.eCond.SetPulsRamp0(     env.GetValue("fedra.track.PulsRamp0"      , "15 20") );
    etra.eCond.SetPulsRamp04(    env.GetValue("fedra.track.PulsRamp04"     , "15 20") );
    etra.eCond.SetDegrad(        env.GetValue("fedra.track.Degrad"          , 4) );
      
    etra.eDTmax                 = env.GetValue("fedra.track.DTmax"          ,     0.07 );
    etra.eDRmax                 = env.GetValue("fedra.track.DRmax"          ,    45.   );
    etra.eDZGapMax              = env.GetValue("fedra.track.DZGapMax"       ,  5000.   );
    bool        do_erase        = env.GetValue("fedra.track.erase"          ,  false   );
    const char  *cut            = env.GetValue("fedra.readCPcut"            ,     "1"  );
    bool        do_misalign     = env.GetValue("fedra.track.do_misalign"    ,      0   );
    int         npass           = env.GetValue("fedra.track.npass"          ,      1   );
    float       misalign_offset = env.GetValue("fedra.track.misalign_offset",    500.  );
    bool        do_local_corr   = env.GetValue("fedra.track.do_local_corr"  ,      1   );
    bool        eDoRealign      = env.GetValue("fedra.track.do_realign"     ,      0   );
    bool        do_comb         = env.GetValue("fedra.track.do_comb"        ,      0   );
    eNsegMin                    = env.GetValue("fedra.track.NsegMin"        ,      2   );
  
  //  etra.InitTrZMap(  2400, 0, 120000,   2000, 0, 100000,   30 );
    etra.InitTrZMap(  env.GetValue("fedra.track.TrZmap", "2400 0 120000   2000 0 100000   30" ) );
  
    EdbPattern p;
    
    EdbAffine2D misalign[60];
    if(do_misalign) {
        //           1 2 3 4 5 6  7  8  9
      int dx[9] = {0,0,1,1,1,0,-1,-1,-1};
      int dy[9] = {0,1,1,0,-1,-1,-1,0,1};
      for(int i=0; i<60; i++) {
        misalign[i].ShiftX( dx[i%9] * misalign_offset );
        misalign[i].ShiftY( dy[i%9] * misalign_offset );
        printf("%d |  %d  %d\n",i, dx[i%9], dy[i%9]);
      }
    }
  
  // read segments and use them for tracking
    for(int ipass=0; ipass<npass; ipass++) {
      printf("\n\n*************** ipass=%d ************\n",ipass);
      etra.eCollisionsRate=0;
      for(int i=0; i<npl; i++) {
        EdbID *id = ss->GetID(i);
      
        EdbPlateP *plate = ss->GetPlate(id->ePlate);
      
        EdbPattern p;
        eSproc->ReadPatCPnopar(p,*id, cut, do_erase);
        p.SetZ(plate->Z());
        p.SetSegmentsZ();
        p.SetID(i);
        p.SetPID(i);
        p.SetSegmentsPID();
      //plate->Print();
        p.Transform(    plate->GetAffineXY()   );
        p.TransformShr( plate->Shr() );
        p.TransformA(   plate->GetAffineTXTY() );
        p.SetSegmentsPlate(id->ePlate);
      //printf("pattern with z: %f\n", p.Z());
      
        if(do_local_corr) {
          int nseg = p.N();
          for(int j=0; j<nseg; j++) {
            EdbSegP *s = p.GetSegment(j);
            plate->CorrectSegLocal(*s);
          }
        }
      
      
        if(do_misalign) {
          p.Transform(&misalign[i]);
          Log(2,"EdbScanTracking::TrackSetBT","apply misalignment of %f",misalign_offset);
        //misalign[i].Print();
        }
  
        if(i>0) etra.ExtrapolateTracksToZ(p.Z());
        if( eDoRealign && i==1 ) etra.CheckPatternAlignment(p,1);
        if( eDoRealign && i>1  ) etra.CheckPatternAlignment(p,2);
        etra.FillTrZMap();
        etra.AddPattern(p);
      }
    }

    int ntr = etra.Tracks().GetEntries();
    
    for(int i=0; i<ntr; i++) {
      EdbTrackP *t = (EdbTrackP *)(etra.Tracks().At(i));
      if(t->N()<eNsegMin)  t->SetFlag(-10);
    }
    
    etra.SetSegmentsErrors();
    etra.FitTracks();

    TObjArray selectedTracks(ntr);
    if(do_comb) {
      etra.CombTracks(selectedTracks);
    } else {
      int cnt=0;
      for(int i=0; i<ntr; i++) {
        EdbTrackP *t = (EdbTrackP *)(etra.Tracks().At(i));
        if(t->Flag()!=-10)  {
          t->SetID(cnt++);
          t->SetCounters();
          t->SetSegmentsTrack();
          selectedTracks.Add(t);
        }
      }
    }
    EdbDataProc::MakeTracksTree( selectedTracks, 0., 0., Form("b%s.trk.root", idset.AsString()) );
    
    TFile f( Form("b%s.trk.root", idset.AsString()) ,"UPDATE");
    env.Write();
    f.Close();
}
