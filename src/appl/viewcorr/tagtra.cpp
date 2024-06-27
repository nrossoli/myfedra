//-- Author :  Valeri Tioukov   3/06/2024

#include <string.h>
#include <iostream>
#include <TEnv.h>
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbScanTracking.h"
//#include "EdbMosaic.h"
#include "EdbMosaicIO.h"
using namespace std;

void TrackSetTag(EdbID idset, TEnv &env, EdbScanProc &sproc);

void print_help_message()
{
  cout<< "\nUsage: \n\t  tagtra -set=ID [ -o=DATA_DIRECTORY -v=DEBUG] \n\n";
  cout<< "\t\t  ID              - id of the dataset formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout<< "\t\t  DEBUG           - verbosity level: 0-print nothing, 1-errors only, 2-normal, 3-print all messages\n";
  cout<< "\t\t  -o              - the data directory\n";
  cout<< "\nExample: \n";
  cout<< "\t  tagtra -id=4554.10.1.0 -o=/scratch/BRICKS \n";
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<<endl;
}

void set_default(TEnv &cenv)
{
  // default parameters for tracking
  cenv.SetValue("fedra.track.TrZmap", "4000 0 200000   4000 0 200000   30" );
  cenv.SetValue("fedra.track.npass"     ,   1 );
  cenv.SetValue("fedra.track.minPlate"  ,-999 );
  cenv.SetValue("fedra.track.maxPlate"  , 999 );
  cenv.SetValue("fedra.track.refPlate"  , 999 );
  cenv.SetValue("fedra.track.nsegmin"   ,   2 );
  cenv.SetValue("fedra.track.ngapmax"   ,   4 );
  cenv.SetValue("fedra.track.DZGapMax"  , 5000. );
  cenv.SetValue("fedra.track.DRmax"     , 100. );
  cenv.SetValue("fedra.track.DTmax"     , 0.07 );
  
  cenv.SetValue("fedra.track.Sigma0" , "50 50 0.005 0.005");
  cenv.SetValue("fedra.track.PulsRamp0"      , "0.1 0.1");
  cenv.SetValue("fedra.track.PulsRamp04"     , "0.1 0.1");
  cenv.SetValue("fedra.track.Degrad"         ,  4   );
  
  cenv.SetValue("fedra.track.probmin"   , 0.001 );
  cenv.SetValue("fedra.track.momentum"  , 2 );
  cenv.SetValue("fedra.track.mass"      , 0.14 );
  cenv.SetValue("fedra.track.do_use_mcs", 0 );
  cenv.SetValue("fedra.track.RadX0"     , 5810.);
  
  cenv.SetValue("fedra.track.do_shtag"       , false );
  cenv.SetValue("fedra.track.do_misalign"    , false );
  cenv.SetValue("fedra.track.misalign_offset", 5000.);
  cenv.SetValue("fedra.track.do_comb"        ,    0 );
  cenv.SetValue("fedra.track.NsegMin"        ,    2 );
  cenv.SetValue("tagtra.outdir"          , "..");
  cenv.SetValue("tagtra.env"             , "tagtra.rootrc");
  cenv.SetValue("tagtra.EdbDebugLevel"   , 1);
}

int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("trackenv");
  set_default(cenv);
  gEDBDEBUGLEVEL        = cenv.GetValue("tagtra.EdbDebugLevel" , 1);
  const char *env       = cenv.GetValue("tagtra.env"            , "tagtra.rootrc");
  const char *outdir    = cenv.GetValue("tagtra.outdir"         , "..");
  
  bool      do_set              = false;
  bool      do_pred             = false;
  Int_t     pred_plate  = 0, to_plate=0;
  Int_t     brick=0, plate=0, major=0, minor=0;
  
  for(int i=1; i<argc; i++ ) {
    char *key  = argv[i];
    
    if(!strncmp(key,"-set=",5))
    {
      if(strlen(key)>5)	sscanf(key+5,"%d.%d.%d.%d",&brick,&plate,&major,&minor);
      do_set=true;
    }
    else if(!strncmp(key,"-o=",3)) 
    {
      if(strlen(key)>3)	outdir=key+3;
    }
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
  }
  
  if(!do_set)   { print_help_message(); return 0; }
  
  cenv.SetValue("tagtra.env"                  , env);
  cenv.ReadFile( cenv.GetValue("tagtra.env"   , "tagtra.rootrc") ,kEnvLocal);
  cenv.SetValue("tagtra.outdir"               , outdir);
  cenv.WriteFile("tagtra.save.rootrc");
  
  if(do_set) {
    EdbScanProc sproc;
    sproc.eProcDirClient=outdir;
    printf("\n----------------------------------------------------------------------------\n");
    printf("tracking set %d.%d.%d.%d\n", brick,plate, major,minor);
    printf("----------------------------------------------------------------------------\n\n");
    
    EdbID id(brick,plate,major,minor);
    EdbScanSet *ss = sproc.ReadScanSet(id);
    ss->Brick().SetID(brick);
    TrackSetTag(id,cenv,sproc);
  }
  cenv.WriteFile("tagtra.save.rootrc");
  
  return 1;
}

//-----------------------------------------------------------------------------------
int GetPattern( EdbPattern &pat, EdbID id, int side )
{
  int n=0;
  EdbMosaicIO mio;
  TString file;
  file.Form("p%3.3d/%d.%d.%d.%d.tag.root",
	    id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor);
  TFile f( file.Data() );
  EdbPattern *p = (EdbPattern*)(f.Get( Form("gpat%d_%d_0",id.ePlate, side) ) );
  if(p) {
    n = p->N();
    //p->SetSegmentsW(1.);                        //TODO set correct weights peak weights in mostag
    pat.AddPattern(*p);
    delete p;
  }
  //  p->SetScanID(id); p->SetSide(side);
  printf( "read %d segments from %s\n",n, file.Data() );
  return n;
}


//--------------------------------------------------------------------------------------
void TrackSetTag(EdbID idset, TEnv &env, EdbScanProc &sproc)
{
  
  // read scanset object
  EdbScanSet *ss = sproc.ReadScanSet(idset);
  if(!ss) { Log(1,"EdbScanTracking::TrackSetBT", "Error! set for %s do not found",idset.AsString()); return; }
  
  int npl = ss->eIDS.GetSize();
  if(npl<2) { Log(1,"EdbScanTracking::TrackSetBT", "Warning! npl<2 : %d stop tracking!",npl); return; }
  
  // create and init tracking object
  EdbTrackAssembler etra;
  
  etra.eCond.SetSigma0(        env.GetValue("fedra.track.Sigma0"         , "50 50 0.005 0.005") );
  etra.eCond.SetPulsRamp0(     env.GetValue("fedra.track.PulsRamp0"      , "0 0") );
  etra.eCond.SetPulsRamp04(    env.GetValue("fedra.track.PulsRamp04"     , "0 0") );
  etra.eCond.SetDegrad(        env.GetValue("fedra.track.Degrad"         , 0) );
  etra.SetRadLength(           env.GetValue("fedra.track.RadX0"          , 5810.) );
  etra.eDoUseMCS              = env.GetValue("fedra.track.do_use_mcs"    , 0 );
  
  etra.eDTmax                 = env.GetValue("fedra.track.DTmax"          ,     0.07 );
  etra.eDRmax                 = env.GetValue("fedra.track.DRmax"          ,    100.   );
  etra.eDZGapMax              = env.GetValue("fedra.track.DZGapMax"       ,   5000.   );
  etra.eProbMin               = env.GetValue("fedra.track.probmin"        ,  0.001   );
  bool        do_erase        = env.GetValue("fedra.track.erase"          ,  false   );
  const char  *cut            = env.GetValue("fedra.readCPcut"            ,     "1"  );
  bool        do_shtag        = env.GetValue("fedra.track.do_shtag"       ,  false   );  // no tracking only shower tag info
  bool        do_misalign     = env.GetValue("fedra.track.do_misalign"    ,      0   );
  int         npass           = env.GetValue("fedra.track.npass"          ,      1   );
  float       misalign_offset = env.GetValue("fedra.track.misalign_offset",   2000.  );
  bool        do_local_corr   = env.GetValue("fedra.track.do_local_corr"  ,      0   );
  bool        eDoRealign      = env.GetValue("fedra.track.do_realign"     ,      0   );
  bool        do_comb         = env.GetValue("fedra.track.do_comb"        ,      0   );
  int         nsegmin         = env.GetValue("fedra.track.NsegMin"        ,      2   );
  float       momentum        = env.GetValue("fedra.track.momentum"       ,      2.  );  
  etra.SetMomentum (momentum);
  
  etra.InitTrZMap(  env.GetValue("fedra.track.TrZmap", "2400 0 120000   2000 0 100000   30" ) );
  
  EdbAffine2D misalign[60];
  if(do_misalign) {
    //           1 2 3 4 5 6  7  8  9
    int dx[9] = {0,0,1,1,1,0,-1,-1,-1};
    int dy[9] = {0,1,1,0,-1,-1,-1,0,1};
    for(int i=0; i<60; i++) {
      misalign[i].ShiftX( dx[i%9] * misalign_offset );
      misalign[i].ShiftY( dy[i%9] * misalign_offset );
      if(gEDBDEBUGLEVEL>1) printf("%d |  %d  %d\n",i, dx[i%9], dy[i%9]);
    }
  }
  
  // read segments and use them for tracking
  int all_volume_segments=0;
  for(int ipass=0; ipass<npass; ipass++) {
    if(gEDBDEBUGLEVEL>1) printf("\n\n*************** ipass=%d ************\n",ipass);
			   etra.eCollisionsRate=0;
    for(int i=0; i<npl; i++) {
      EdbID *id = ss->GetID(i);
      
      EdbPlateP *plate = ss->GetPlate(id->ePlate);
      
      EdbPattern p;
      p.SetScanID(*id);
      
      int side =1;
      all_volume_segments += GetPattern( p, *id, side );
      
      p.SetZ(plate->Z());
      p.SetSegmentsZ();
      p.SetID(i);
      p.SetSide(side);
      p.SetPID(i);
      p.SetSegmentsPID();
      p.Transform(    plate->GetAffineXY()   );
      p.TransformShr( plate->Shr() );
      p.TransformA(   plate->GetAffineTXTY() );
      p.SetSegmentsPlate(id->ePlate);
      
      
      if(do_misalign) {
	p.Transform(&misalign[i]);
	Log(2,"EdbScanTracking::TrackSetBT","apply misalignment of %f",misalign_offset);
      }
      
      if(i>0) etra.ExtrapolateTracksToZ(p.Z());

      etra.FillTrZMap();
      etra.AddPattern(p);

      //etra.FillXYseg(p);
    }
  }

  
  
  int ntr = etra.Tracks().GetEntriesFast();
  
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = (EdbTrackP *)(etra.Tracks().At(i));
    if(t->N()<nsegmin)  t->SetFlag(-10);
  }
  
  etra.SetSegmentsErrors();
  etra.FitTracks();
  
  int cnt_badtrk=0;
  int cnt_attached_segments=0;
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
	t->SetP(momentum);
	selectedTracks.Add(t);
	cnt_attached_segments += t->N();
      } else cnt_badtrk++;
    }
  }
  printf("\n segments all/attached: %d/%d     tracks selected/bad = %d/%d \n",
	 all_volume_segments, cnt_attached_segments, selectedTracks.GetEntriesFast(), cnt_badtrk );
  EdbDataProc::MakeTracksTree( selectedTracks, 0., 0., Form("b%s.trk.root", idset.AsString()) );
  TFile f( Form("b%s.trk.root", idset.AsString()) ,"UPDATE");
  env.Write();
  f.Close();
  
  //SaveHist(idset,etra);
}
