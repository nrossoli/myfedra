//-- Author :  Valeri Tioukov   20/04/2023

#include <string.h>
#include <iostream>
#include <TRint.h>
#include <TEnv.h>
#include <TChain.h>
#include <TList.h>
#include "EdbLog.h"
#include "EdbRunAccess.h"
#include "EdbLinking.h"
#include "EdbScanProc.h"
#include "EdbPlateAlignment.h"
#include "EdbMosaic.h"
#include "EdbMosaicIO.h"
#include "EdbAttachPath.h"

using namespace std;
using namespace TMath;
int  AlignSide( EdbID id, int side, TEnv &env );
bool Align2Fragments( EdbPattern &p1, EdbPattern &p2, EdbLayer &l1);

void print_help_message()
{
  cout<< "\nUsage: \n";
  cout<< "\t  mosalign  -id=ID  [-from=frag0 -nfrag=N  -merge -v=DEBUG] \n";

  cout<< "\t\t  ID    - id of the raw.root file formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout<< "\t\t  frag0 - the first fragment (default: 0) \n";
  cout<< "\t\t  N     - number of fragments to be processed (default: upto 1000000, stop at first empty) \n";
  cout<< "\t\t  merge - merge all fragments into one cp file \n";
  
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<< "\n If the parameters file (mosalign.rootrc) is not presented - the default \n";
  cout<< " parameters will be used. After the execution them are saved into mosalign.save.rootrc file\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default_link(TEnv &cenv)
{
  // default parameters for the new linking
  cenv.SetValue("fedra.link.AFID"                ,  1   );   // 1 is usually fine for scanned data; for the db-read data use 0!
  cenv.SetValue("fedra.link.DoImageCorr"           , 0  );
  cenv.SetValue("fedra.link.ImageCorrSide1"           , "1. 1. 0.");
  cenv.SetValue("fedra.link.ImageCorrSide2"           , "1. 1. 0.");
  
  cenv.SetValue("fedra.link.DoImageMatrixCorr"              , 0  );
  cenv.SetValue("fedra.link.ImageMatrixCorrSide1"           , "");
  cenv.SetValue("fedra.link.ImageMatrixCorrSide2"           , "");
  
  cenv.SetValue("fedra.link.CheckUpDownOffset"   ,  1   );   // check dXdY offsets between up and correspondent down views
  cenv.SetValue("fedra.link.BinOK"               , 6.   );
  cenv.SetValue("fedra.link.NcorrMin"            , 100  );
  cenv.SetValue("fedra.link.DoCorrectShrinkage"  , false );
  cenv.SetValue("fedra.link.read.InvertSides"    , 0  );
  cenv.SetValue("fedra.link.read.HeaderCut"      , "1"  );
  cenv.SetValue("fedra.link.read.UseDensityAsW"      , false  );
  cenv.SetValue("fedra.link.read.ICUT"           , "-1     -500. 500.   -500.   500.    -1.   1.      -1.   1.       0.  50.");
  cenv.SetValue("fedra.link.RemoveDoublets"      , "1    2. .01   1");  //yes/no   dr  dt  checkview(0,1,2)
  cenv.SetValue("fedra.link.DumpDoubletsTree"    , true );
  cenv.SetValue("fedra.link.shr.NsigmaEQ"        , 7.5  );
  cenv.SetValue("fedra.link.shr.Shr0"            , .990  );
  cenv.SetValue("fedra.link.shr.DShr"            , .3   );
  cenv.SetValue("fedra.link.shr.ThetaLimits"     ,"0.0  1.");
  cenv.SetValue("fedra.link.DoCorrectAngles"     , true );
  cenv.SetValue("fedra.link.ang.Chi2max"         , 1.5  );
  cenv.SetValue("fedra.link.DoFullLinking"       , true );
  cenv.SetValue("fedra.link.full.NsigmaEQ"       , 5.5  );
  cenv.SetValue("fedra.link.full.DR"             , 20.  );
  cenv.SetValue("fedra.link.full.DT"             , 0.1  );
  cenv.SetValue("fedra.link.full.CHI2Pmax"       , 3.   );
  cenv.SetValue("fedra.link.DoSaveCouples"       , true );
  cenv.SetValue("fedra.link.Sigma0"         , "1 1 0.007 0.007");
  cenv.SetValue("fedra.link.PulsRamp0"      , "6 9");
  cenv.SetValue("fedra.link.PulsRamp04"     , "6 9");
  cenv.SetValue("fedra.link.Degrad"         ,  5   );
  
  cenv.SetValue("fedra.link.LLfunction"     , "0.256336-0.16489*x+2.11098*x*x" );
  cenv.SetValue("fedra.link.CPRankingAlg"   , 0 );

  cenv.SetValue("emlink.reportfileformat"   , "pdf" );

  cenv.SetValue("emlink.outdir"          , "..");
  cenv.SetValue("emlink.env"             , "link.rootrc");
  cenv.SetValue("emlink.EdbDebugLevel"   , 1);
}

int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("mosalignenv");
  gEDBDEBUGLEVEL     = cenv.GetValue("mosalign.EdbDebugLevel" , 1);
  const char *env    = cenv.GetValue("mosalign.env"            , "mosalign.rootrc");
  const char *outdir = cenv.GetValue("mosalign.outdir"         , "..");
  
  bool      do_single   = false;
  bool      do_set      = false;
  bool      do_merge    = false;
  EdbID     id;
  int       from_fragment=0;
  int       n_fragments=1000000;

  for(int i=1; i<argc; i++ ) {
    char *key  = argv[i];

    if(!strncmp(key,"-set=",5))
    {
       if(strlen(key)>5) if(id.Set(key+5))   do_set=true;
    }
    else if(!strncmp(key,"-id=",4))
    {
      if(strlen(key)>4) if(id.Set(key+4))   do_single=true;
    }
    else if(!strncmp(key,"-from=",6))
    {
      if(strlen(key)>6) from_fragment = atoi(key+6);
    }
    else if(!strncmp(key,"-nfrag=",7))
    {
      if(strlen(key)>7) n_fragments = atoi(key+7);
    }
    else if(!strncmp(key,"-merge",6))
    {
      do_merge=true;
    }
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
  }

  if(!(do_single||do_set))   { print_help_message(); return 0; }
  if(  do_single&&do_set )   { print_help_message(); return 0; }

  set_default_link(cenv);
  cenv.SetValue("mosalign.env"            , env);
  cenv.ReadFile( cenv.GetValue("mosalign.env"   , "mosalign.rootrc") ,kEnvLocal);
  cenv.SetValue("mosalign.outdir"         , outdir);

  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("mosalign.outdir","..");
  cenv.WriteFile("mosalign.save.rootrc");

  printf("\n----------------------------------------------------------------------------\n");
  printf("mosalign  %s\n"      ,id.AsString()	   );
  printf(  "----------------------------------------------------------------------------\n\n");


  if(do_single)
  {
    AlignSide(id, 1, cenv);
    AlignSide(id, 2, cenv);
  }

  cenv.WriteFile("mosalign.save.rootrc");
  return 1;
}

//-----------------------------------------------------------------------------------
int AlignSide( EdbID id, int side, TEnv &cenv )
{
  EdbMosaicIO mio;
  TString file;
  file.Form("p%3.3d/%d.%d.%d.%d.mos.root",
	      id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor);
  mio.Init( file.Data() );
  
  EdbLayer *mapside = mio.GetCorrMap( id.ePlate, side );
  
  int nc=mapside->Map().Ncell();
  EdbAttachPath atp(nc);
  for( int i=0; i<nc; i++ )
  {
    EdbLayer *l = mapside->Map().GetLayer(i);
    atp.SetPoint(i,l->X(),l->Y(),l->ID());
   }
  //atp.SetStartingPosition(0);
  atp.SetStartingAtCenter();
  atp.OrderPointsRadial();
  atp.SetRange(11000);
  atp.RiseOK( atp.I(0) );
  
  int n = atp.N();
  if(n>1)
  {
    for( int i=1; i<n; i++ )
    {
      int id1 = atp.ID( atp.I(i) );
      EdbPattern *p1 = mio.GetFragment( id.ePlate, side, id1, false );
      p1->SetScanID(id);
      EdbLayer *l1 = mapside->Map().GetLayer( atp.I(i) );
      p1->Transform( l1->GetAffineXY() );
      printf("%d %f (%f %f), %f (%f %f) \n",l1->ID(),
	   l1->X(), p1->Xmin(), p1->Xmax(),l1->Y(), p1->Ymin(), p1->Ymax() );
     
      TArrayI list(10);
      int nn = atp.GetAlignedNeighbours( atp.I(i), list );
      printf( "nn = %d\n", nn);
      EdbPattern p2;
      p2.SetScanID(id);
      if(nn)
      {
	for( int j=0; j<nn; j++ )
	{
	  EdbPattern *p = mio.GetFragment( id.ePlate, side, list[j], false );
	  EdbLayer *l2 = mapside->Map().GetLayer( list[j] );
	  p2.Transform( l2->GetAffineXY() );
	  p2.AddPattern( *p );
	  delete p;
	}
      }
      if( Align2Fragments(*p1, p2, *l1) ) atp.RiseOK( atp.I(i) );
      delete p1;
    }
   }

  
   for( int i=0; i<nc; i++ )
  {
    EdbLayer *l = mapside->Map().GetLayer(i);
    l->GetAffineXY()->Print();
  }
  mio.SaveCorrMap( id.ePlate, side, *mapside,  file.Data() );
  
}

//-----------------------------------------------------------------------
bool Align2Fragments( EdbPattern &p1, EdbPattern &p2, EdbLayer &l1)
{
  bool success=false;
  int eMinPeak=100;
  bool do_transform = true; 
  
  EdbPlateAlignment av;
  av.eNoScale     = 1;         // calculate shift and rotation
  av.eNoScaleRot  = 0;         // calculate shift only
  av.eOffsetMax   = 10.;
  av.eDZ          = 0.;
  av.eDPHI        = 0.02;
  av.eDoFine      = 1;
  av.eSaveCouples = 1;
  av.SetSigma( 0.3, 0.007 );
  av.eDoublets[0] = av.eDoublets[1]=0.01;
  av.eDoublets[2] = av.eDoublets[3]=0.0001;

  av.eSaveCouples=1;
  av.InitOutputFile( Form( "p%.3d/%d_%d.alf.root", p1.ScanID().ePlate, p1.ID(), p2.ID() ) ); 
  av.Align( p1, p2, 0); //-190
  EdbAffine2D *affXY = av.eCorrL[0].GetAffineXY();
  EdbAffine2D *affTXTY = av.eCorrL[0].GetAffineTXTY();

  if(av.eNcoins > eMinPeak )
  {
    if(do_transform) p1.Transform( affXY );
    l1.GetAffineXY()->Transform( affXY );
    success=true;
  }
  
  av.CloseOutputFile();
  return success;
}

