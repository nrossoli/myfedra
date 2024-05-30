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
int  AlignToBeam( EdbID id, TEnv &env );
bool AlignFragmentToBeam0( EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2);

void print_help_message()
{
  cout<< "\nUsage: \n";
  cout<< "\t  mosalignebeam  -id=ID  [-from=frag0 -nfrag=N  -merge -v=DEBUG] \n";

  cout<< "\t\t  ID    - id of the raw.root file formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout<< "\t\t  frag0 - the first fragment (default: 0) \n";
  cout<< "\t\t  N     - number of fragments to be processed (default: upto 1000000, stop at first empty) \n";
  cout<< "\t\t  merge - merge all fragments into one cp file \n";
  
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<< "\n If the parameters file (mosalignebeam.rootrc) is not presented - the default \n";
  cout<< " parameters will be used. After the execution them are saved into mosalignebeam.save.rootrc file\n";
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
  
  TEnv cenv("mosalignebeamenv");
  gEDBDEBUGLEVEL     = cenv.GetValue("mosalignebeam.EdbDebugLevel" , 1);
  const char *env    = cenv.GetValue("mosalignebeam.env"            , "mosalignebeam.rootrc");
  const char *outdir = cenv.GetValue("mosalignebeam.outdir"         , "..");
  
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
  cenv.SetValue("mosalignebeam.env"            , env);
  cenv.ReadFile( cenv.GetValue("mosalignebeam.env"   , "mosalignebeam.rootrc") ,kEnvLocal);
  cenv.SetValue("mosalignebeam.outdir"         , outdir);

  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("mosalignebeam.outdir","..");
  cenv.WriteFile("mosalignebeam.save.rootrc");

  printf("\n----------------------------------------------------------------------------\n");
  printf("mosalignebeam  %s\n"      ,id.AsString()	   );
  printf(  "----------------------------------------------------------------------------\n\n");


  if(do_single)
  {
    AlignToBeam(id, cenv);
  }

  cenv.WriteFile("mosalignebeam.save.rootrc");
  return 1;
}

//-----------------------------------------------------------------------------------
int AlignToBeam( EdbID id, TEnv &cenv )
{
  EdbMosaicIO mio;
  TString file;
  file.Form("p%3.3d/%d.%d.%d.%d.mos.root",
	      id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor);
  mio.Init( file.Data() );
  
  bool use_saved_alignment=true;
  
  EdbLayer *mapside1 = mio.GetCorrMap( id.ePlate, 1 ); 
  EdbLayer *mapside2 = mio.GetCorrMap( id.ePlate, 2 ); // align side 2 to side 1
  
  int nc=mapside2->Map().Ncell();
  for( int i=0; i<nc; i++ )
  {
    EdbLayer   *l1 = mapside1->Map().GetLayer(i);
    EdbLayer   *l2 = mapside2->Map().GetLayer(i);
    if(!use_saved_alignment) 
    {
      l1->GetAffineXY()->Reset();
      l2->GetAffineXY()->Reset();
      l1->GetAffineTXTY()->Reset();
      l2->GetAffineTXTY()->Reset();
    }
    EdbPattern *p1 = mio.GetFragment( id.ePlate, 1, i, use_saved_alignment ); //get side 1
    EdbPattern *p2 = mio.GetFragment( id.ePlate, 2, i, use_saved_alignment ); //get side 2
    p1->SetScanID(id);
    p2->SetScanID(id);
    AlignFragmentToBeam0(*p2, *p1, *l2,*l1);  //align 2 to 1 using parallel beam tracks
    delete p1;
    delete p2;
   }

/*  
   for( int i=0; i<nc; i++ )
  {
    EdbLayer *l = mapside->Map().GetLayer(i);
    l->GetAffineXY()->Print();
  }
  */

  mio.SaveCorrMap( id.ePlate, 1, *mapside1,  file.Data() );
  mio.SaveCorrMap( id.ePlate, 2, *mapside2,  file.Data() );
 
}

//-----------------------------------------------------------------------
bool AlignFragmentToBeam0( EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2)
{
  // Assume 0 angle beam here
  //
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
  av.eDoCorrectAngle = false;  
  av.eSaveCouples=1;
  
  av.InitOutputFile( Form( "p%.3d/%d_%d.ab0.root", p1.ScanID().ePlate, p1.ID(), p2.ID() ) ); 
  av.Align( p1, p2, 0); //-190
  EdbAffine2D *affXY = av.eCorrL[0].GetAffineXY();
  EdbAffine2D *affTXTY = av.eCorrL[0].GetAffineTXTY();
  
  float dtx1 = av.CalcMeanDiff2Const(2,0,0); 
  float dty1 = av.CalcMeanDiff2Const(3,0,0);
  float dtx2 = av.CalcMeanDiff2Const(2,1,0); 
  float dty2 = av.CalcMeanDiff2Const(3,1,0);

  printf("\n angular offsets found: %f %f %f %f\n\n", dtx1,dty1,dtx2,dty2);
  
  if(av.eNcoins > eMinPeak )
  {
    if(do_transform) p1.Transform( affXY );
    l1.GetAffineXY()->Transform( affXY );
    if(1) {
      l1.GetAffineTXTY()->ShiftX(-dtx1);
      l1.GetAffineTXTY()->ShiftY(-dty1);
      l2.GetAffineTXTY()->ShiftX(-dtx2);
      l2.GetAffineTXTY()->ShiftY(-dty2);
    }
    l1.GetAffineTXTY()->Print();
    l2.GetAffineTXTY()->Print();
    success=true;
  }
  
  av.CloseOutputFile();
  return success;
}

