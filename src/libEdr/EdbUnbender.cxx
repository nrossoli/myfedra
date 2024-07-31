//-- Author :  Valeri Tioukov   17/24/2024
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbUnbender - unbend tracks group                                    //
//                                                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TMath.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TPolyMarker3D.h>
#include "EdbLog.h"
#include "EdbAlignmentV.h"
#include "EdbTrackFitter.h"
#include "EdbUnbender.h"

using namespace TMath;
ClassImp(EdbUnbender);

//-----------------------------------------------------------------------
void EdbUnbender::Set0()
{
  eV=0;
  eSS=0;
  eEnv=0;
  do_save_hist=0;
  do_save_tree=0;
}

//-----------------------------------------------------------------------
void EdbUnbender::Unbend3a(EdbPVRec &pvr, EdbScanSet &ss, TEnv &env)
{
 eV=&pvr; eSS=&ss; eEnv=&env;
  
  for( int iter=0; iter<5; iter++ )
  {
    CalcMeanPath(iter);   
    for(int i=0; i<3; i++) 
    {
      GlobalCorr3(0,2, 100*iter+i);
      GlobalCorr3(1,2, 100*iter+i);
    }
  }
}

//-----------------------------------------------------------------------
void EdbUnbender::Unbend3g(EdbPVRec &pvr, EdbScanSet &ss, TEnv &env)
{
  eV=&pvr; eSS=&ss; eEnv=&env;
  for( int iter=0; iter<1; iter++ )
  {
    CalcMeanPath(iter);   
    for(int i=0; i<1; i++) 
    {
      GlobalCorrN( 3,1,0,2, 100*iter+i); // length, position, offset, step, iteration
      GlobalCorrN( 3,1,1,2, 200*iter+i); // length, position, offset, step, iteration
    }
  }
}

//-----------------------------------------------------------------------
void EdbUnbender::Unbend5g(EdbPVRec &pvr, EdbScanSet &ss, TEnv &env)
{
  eV=&pvr; eSS=&ss; eEnv=&env;
  if(do_save_tree) eCPT.InitCouplesTree( "couples","unbend.root","RECREATE");
  
  int niter=1;
  for( int iter=0; iter<1; iter++ )
  {
    //CalcMeanPath(iter);
    for(int i=0; i<niter; i++) GlobalCorrN( 5,2,0, 1,  520100+i); // length, position, offset, step, iteration 
    //for(int i=0; i<niter; i++) GlobalCorrN( 5,1,0, 1,  510100+i); // length, position, offset, step, iteration 
    //for(int i=0; i<niter; i++) GlobalCorrN( 5,3,0, 1,  530100+i); // length, position, offset, step, iteration 
    //for(int i=0; i<niter; i++) GlobalCorrN( 5,0,0, 5,  500500+i); // length, position, offset, step, iteration 
    //for(int i=0; i<niter; i++) GlobalCorrN( 5,4,0,-5, -540500+i); // length, position, offset, step, iteration 
    //for(int i=0; i<2; i++)     GlobalCorrN( 5,2,0, 1,  520110+i); // length, position, offset, step, iteration
  }
}
//-----------------------------------------------------------------------
void EdbUnbender::CalcMeanPath( int cycle )
{
  EdbPVRec   &ali = *eV;
  EdbScanSet &ss  = *eSS;

  int npat = ali.Npatterns();
  
  double x0[npat],y0[npat], z0[npat], tx0[npat], ty0[npat], w0[npat];
  for( int ipat=0; ipat<npat; ipat++ )
  {
    EdbPattern *p = ali.GetPattern(ipat);
    if(p)
    {
      z0[ipat] = p->Z();
      x0[ipat]=y0[ipat]=tx0[ipat]=ty0[ipat]=w0[ipat]=0;
    }
  }

  EdbSegP s0;
  int ntr = ali.Ntracks();
  for( int i=0; i<ntr; i++ )
  {
    EdbTrackP *t = ali.GetTrack(i);
    for( int ipat=0; ipat<npat; ipat++ )
    {
      EdbPattern *p = ali.GetPattern(ipat);
      if(p)
      {
	t->EstimatePositionAt( z0[ipat], s0 );
	x0[ipat]  += s0.X();
	y0[ipat]  += s0.Y();
	tx0[ipat] += s0.TX();
	ty0[ipat] += s0.TY();
	w0[ipat]  += 1;
      }
    }
  }
  
  TH3F hmp("hmp","mean path", 100, -50, 50, 100, -50,50, npat, 0, z0[npat-1] );
  TPolyMarker3D mp(npat);
  EdbAffine2D aff;
  for( int ipat=0; ipat<npat; ipat++ )
  {
    EdbPattern *p = ali.GetPattern(ipat);
    if(p)
    {
      aff.Reset();
      x0[ipat]  /= w0[ipat];
      y0[ipat]  /= w0[ipat];
      tx0[ipat] /= w0[ipat];
      ty0[ipat] /= w0[ipat];
      mp.SetPoint(ipat, x0[ipat] - x0[0],y0[ipat]-y0[0],z0[ipat]-z0[0] );
      hmp.Fill(x0[ipat] - x0[0],y0[ipat]-y0[0],z0[ipat]-z0[0]);
      aff.ShiftX( -(x0[ipat]-x0[0]) );
      aff.ShiftY( -(y0[ipat]-y0[0]) );
      p->Transform(&aff);
      EdbPlateP *plate = ss.GetPlate( p->ScanID().ePlate );
      plate->GetAffineXY()->Transform(&aff);
    }
  }
  if(do_save_hist)
  {
    TFile f(Form("pm_%d.root",cycle),"RECREATE");
    mp.Write("pm");
    hmp.Write("hmp");
    f.Close();
  }
}

//----------------------------------------------------------------------------
void EdbUnbender::GlobalCorrN( int length, int position, int offset, int step, int flag )
{
  EdbPVRec   &ali = *eV;
  EdbScanSet &ss  = *eSS;
  int ntr  = ali.Ntracks();
  int npat = ali.Npatterns();
  Log(2,"EdbUnbender::GlobalCorrN","l_p_o_s_f: %d_%d_%d_%d_%d with %d tracks and %d patterns",length,position,offset,step,flag,  ntr, npat);
  if(npat<length) return;
  if(ntr<3)       return;

  int start,end,incr;
  if(step>0) {
    start=0; 
    end=npat-1;
    incr=1;
  }
  else if(step<0) {
    start=npat-1; 
    end=0;
    incr=-1;
  } else return;

  int cntg=start + incr*offset;
  printf("start end incr: %d %d %d\n",start,end,incr);

  do {
    EdbPattern *p[length];
    int nfill=0;
    int cnt=cntg;
    for(int i=0; i<length; i++)
    {
      printf("i cnt nfill: %d %d %d\n", i,cnt,nfill );
      p[i] = ali.GetPattern(cnt); 
      if(p[i]) {
	p[i]->SetPID(cnt);
	nfill++;
      }
      if(cnt==end) break;
      cnt+=incr;
    }
    if(nfill<(length)) break;
    
    TObjArray   acorr; acorr.SetOwner();    
    for(int itr=0; itr<ntr; itr++) {
      EdbTrackP *t = ali.GetTrack(itr);
      EdbTrackP *tcorr = new EdbTrackP();
      for(int j=0; j<t->N(); j++) {
	EdbSegP *s  = t->GetSegment(j);
	for(int i=0; i<length; i++)
	{
	  if(p[i]) if(s->PID()==p[i]->PID()) tcorr->AddSegment(s);
	}
      }
      acorr.Add(tcorr);
    }

    EdbPlateP *plate = ss.GetPlate( p[position]->ScanID().ePlate );
    CalculateCorrections(acorr, length, *(p[position]), *plate, flag );
    if(cnt==end+incr) break;
    else  cntg+=step;
  }
  while(cntg!=end);
}

//-----------------------------------------------------------------------
void EdbUnbender::CalculateCorrections(const TObjArray &acorr, int length, EdbPattern &p, EdbPlateP &plate, int flag)
{
  EdbAlignmentV al; al.eS[1].SetOwner();
  EdbTrackFitter ft;
  int ntr = acorr.GetEntries();
  Log(2,"EdbUnbender::CalculateCorrections","ntr length: %d %d  pid: %d  plate: %d/%d",
      ntr, length, p.PID(), p.ScanID().ePlate, plate.ID() );
  for(int i=0; i<ntr; i++)
  {
    EdbTrackP *t = (EdbTrackP*)(acorr.At(i));
    if(t->N()!=length) continue;

    EdbSegP   *s1=0;
    EdbTrackP *trest = new EdbTrackP();
    for(int j=0; j<length; j++) 
    {
      EdbSegP *s=t->GetSegment(j);
      if(s->ScanID().ePlate==p.ScanID().ePlate) s1=s;
      else trest->AddSegment(s);
    }
    if(s1)
    {
      ft.Fit3Pos(*trest);
      trest->PropagateTo(s1->Z());
      al.eS[0].Add(s1);
      al.eS[1].Add( new EdbSegP(*trest) );
    }
  }
  
  CheckResolutions(al.eS[0], al.eS[1], p.ScanID().ePlate, flag);
  
  EdbAffine2D affXYcum; 
  EdbAffine2D affTXTYcum;

  for(int iter=0; iter<3; iter++)
  {
    EdbAffine2D affXY;
    al.CalculateAffXY( al.eS[0], al.eS[1], affXY);
    p.Transform(&affXY);
    affXYcum.Transform(&affXY);
    
    EdbAffine2D affTXTY;
    al.CalculateAffTXTYTurn( al.eS[0], al.eS[1], affTXTY);
    p.TransformA(&affTXTY);
    affTXTYcum.Transform(&affTXTY);
    Log(2,"EdbUnbender::CalculateCorrections","%d \n coord: %s\n angle: %s",iter,affXY.AsString(),affTXTY.AsString());
  }
  
  plate.GetAffineTXTY()->Transform(&affTXTYcum);
  plate.GetAffineXY()->Transform(&affXYcum);  
}

//-----------------------------------------------------------------------
void EdbUnbender::CheckResolutions(TObjArray &a1, TObjArray &a2, int plate, int flag)
{ 
  //TH2F hxy( Form("hxy_%d",plate), "hxy", 100, -2.5,2.5, 100, -2.5,2.5 );
  //TH2F htxty( Form("htxty_%d",plate), "htxty", 200, -.01,0.01, 200, -0.01,0.01 );
  int n=a1.GetEntries();
  for(int i=0; i<n; i++)
  {
    EdbSegP *s1 = (EdbSegP*)(a1.At(i));
    EdbSegP *s2 = (EdbSegP*)(a2.At(i));
    //hxy.Fill(s2->X()-s1->X(),s2->Y()-s1->Y());
    //htxty.Fill(s2->TX()-s1->TX(),s2->TY()-s1->TY());
    if(do_save_tree) eCPT.Fill( s1, s2, 0,0, 0,0, 0, flag );

  }
  //hxy.Write();
  //htxty.Write();
  if(do_save_tree) eCPT.SaveTree();
}

//-----------------------------------------------------------------------
void EdbUnbender::GlobalCorr3( int offset, int step, int cycle )
{
  EdbPVRec   &ali = *eV;
  EdbScanSet &ss  = *eSS;
  
  bool doAng = true;
  
  int ntr  = ali.Ntracks();
  int npat = ali.Npatterns();
  Log(2,"EdbUnbender::GlobalCorr3","%d_%d_%d with %d tracks and %d patterns",offset,step,cycle, ntr, npat);
  if(npat<3) return;
  
  TH3F hdxyp("hdxyp","triplet position residuals", 150,-15,15, 150,-15,15,  57,0.5,57.5 );
  TH3F hdtxtyp("hdtxtyp","triplet angular residuals", 150,-0.015,0.015, 150,-0.015,0.015,  57,0.5,57.5 );
  TH2F h2theta13("h2theta13","theta plates 1 3", 57,0.5,57.5, 100, 0., 1.);
  TH2F h2theta2("h2theta2","theta plate 2", 57,0.5,57.5, 100, 0., 1.);
  
  for( int ipat=offset; ipat<npat; ipat+=step )
  {
    int id1=ipat;
    int id2=ipat+1;
    int id3=ipat+2;
    if(id3>=npat-1) break;
    
    EdbPattern *p1 = ali.GetPattern(id1);  p1->SetPID(id1);
    EdbPattern *p2 = ali.GetPattern(id2);  p2->SetPID(id2);
    EdbPattern *p3 = ali.GetPattern(id3);  p3->SetPID(id3);
    
    float lm[4] = { 
      Max(Max( p1->Xmin(), p2->Xmin()), p3->Xmin()), 
      Min(Min( p1->Xmax(), p2->Xmax()), p3->Xmax() ), 
      Max(Max( p1->Ymin(), p2->Ymin()), p3->Ymin() ), 
      Min(Min( p1->Ymax(), p2->Ymax()), p3->Ymax() )
    };
    
    TObjArray p1corr;
    TObjArray p2corr;
    TObjArray p3corr;
    int ncp=0;
    for(int itr=0; itr<ntr; itr++) {
      EdbTrackP *t = ali.GetTrack(itr);
      EdbSegP *s1=0, *s2=0, *s3=0;
      for(int j=0; j<t->N(); j++) {
	EdbSegP *s  = t->GetSegment(j);
	if(s->PID()==p1->PID()) s1=s;
	if(s->PID()==p2->PID()) s2=s;
	if(s->PID()==p3->PID()) s3=s;
      }
      if(s1&&s2&&s3) { p1corr.Add(s1); p2corr.Add(s2); p3corr.Add(s3); ncp++; }
      //FillEfficiency( p2->ScanID().ePlate, s1,s2,s3, lm, h2theta13, h2theta2 );
    }
    
    TObjArray p13corr; 
    p13corr.SetOwner();
    for(int i=0; i<ncp; i++)
    {
      EdbSegP *s1 = ((EdbSegP *)(p1corr.At(i)));
      EdbSegP *s2 = ((EdbSegP *)(p2corr.At(i)));
      EdbSegP *s3 = ((EdbSegP *)(p3corr.At(i)));
      EdbSegP *s  = new EdbSegP( *s2 );
      float r = (s2->Z()-s1->Z())/(s3->Z()-s1->Z());
      s->SetX( s1->X() + (s3->X()-s1->X())*r );
      s->SetY( s1->Y() + (s3->Y()-s1->Y())*r );
      s->SetTX( (s3->X()-s1->X())/(s3->Z()-s1->Z()) );
      s->SetTY( (s3->Y()-s1->Y())/(s3->Z()-s1->Z()) );
      p13corr.Add(s);
    }
    
    printf("ipat: %d  plate: %d pid: %d  z: %.1f    N3 = %d\n", ipat, p2->ScanID().ePlate, p2->PID(), p2->Z(), ncp );
    for(int i=0; i<ncp; i++)
    {
      EdbSegP *s2  = ((EdbSegP*)(p2corr.At(i)));
      EdbSegP *s13 = ((EdbSegP*)(p13corr.At(i)));
      hdxyp.Fill( s13->X()-s2->X(), s13->Y()-s2->Y(),  (float)(s2->ScanID().ePlate) );
      hdtxtyp.Fill( s13->TX()-s2->TX(), s13->TY()-s2->TY(),  (float)(s2->ScanID().ePlate) );
    }
    
    EdbAlignmentV al;
    for(int i=0; i<ncp; i++)
    {
      al.eS[0].Add(p2corr.At(i));
      al.eS[1].Add(p13corr.At(i));
    }
    
    EdbAffine2D aff2XYcum; 
    EdbAffine2D aff2TXTYcum;
    float dz02corr=0;
    float z02set=(p3->Z()+p1->Z())/2.;
    
    for(int iter=0; iter<3; iter++) {
      printf("***** iter %d\n",iter);
      
      EdbAffine2D aff1to2XY;
      al.CalculateAffXY( al.eS[0], al.eS[1], aff1to2XY);
      aff1to2XY.Print();
      p2->Transform(&aff1to2XY);
      aff2XYcum.Transform(&aff1to2XY);
      
      if(doAng) {
	EdbAffine2D aff1to2TXTY;
	al.CalculateAffTXTYTurn( al.eS[0], al.eS[1], aff1to2TXTY);
	aff1to2TXTY.Print();
	p2->TransformA(&aff1to2TXTY);
	aff2TXTYcum.Transform(&aff1to2TXTY);
      }
    }
    
    EdbPlateP *plate2 = ss.GetPlate( p2->ScanID().ePlate );
    if(doAng) plate2->GetAffineTXTY()->Transform(&aff2TXTYcum);
    plate2->GetAffineXY()->Transform(&aff2XYcum);
  }
  
  if(do_save_hist)
  {
    TFile fout(Form("glob3_%d_%d_%d.root",offset,step,cycle),"RECREATE");
    hdxyp.Write();
    hdtxtyp.Write();
    h2theta13.Write();
    h2theta2.Write();
    fout.Close();
  }
}