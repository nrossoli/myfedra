//-- Author :  Valeri Tioukov   30/05/2024

#include <string.h>
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TEnv.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TText.h>
#include <TEllipse.h>
#include "EdbLog.h"
#include "EdbAlignmentV.h"
#include "EdbRunAccess.h"
#include "EdbScanProc.h"
//#include "EdbPlateAlignment.h"
#include "EdbMosaic.h"
#include "EdbMosaicIO.h"

using namespace std;
using namespace TMath;

void DrawOut(const EdbPattern &ptag, EdbMosaicIO &omio);
void FindPeaks(EdbH2 &h2p, EdbPattern &ptag, EdbMosaicIO &omio);
void DrawEllipse(const EdbPattern &ptag, int col );
void DoubletsFilterOut(EdbPattern &p, EdbMosaicIO &omio);
void Pat2H2( const EdbPattern &p, EdbH2 &h2p);
void TagSide( EdbID id, int from, int nfrag, int side, TEnv &cenv, EdbMosaicIO &mio, EdbMosaicIO &omio );

EdbPattern gPat;  // global tags pattern
struct OutHist
{
  TH2F *dblxy;
  TH2F *dbltxty;
  TH1F *spectrOriginal; // original, no smooth
  TH1F *spectrSmooth;   // after smooth
  TH1F *spectrProc;     // after smooth and peaks erasing
  TH2F *xy;
};
OutHist gOH = {0,0,0,0};

bool do_save_gif=false;
bool do_save_canvas=false;

void print_help_message()
{
  cout<< "\nUsage: \n";
  cout<< "\t  mostag  -id=ID  [-from=frag0 -nfrag=N -v=DEBUG] \n";
  
  cout<< "\t\t  ID    - id of the raw.root file formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout<< "\t\t  frag0 - the first fragment (default: 0) \n";
  cout<< "\t\t  N     - number of fragments to be processed (default: upto 1000000, stop at first empty) \n";
  
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<< "\n If the parameters file (mostag.rootrc) is not presented - the default \n";
  cout<< " parameters will be used. After the execution them are saved into mostag.save.rootrc file\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default(TEnv &cenv)
{
  // default parameters for showers tagging
  cenv.SetValue("fedra.tag.save_gif"            ,  false   );   // save gif with plots
  cenv.SetValue("fedra.tag.save_canvas"          ,  false   );   // save canvas objects into root file

  cenv.SetValue("fedra.tag.RemoveDoublets"      , "1    2. .01   1");  //yes/no   dr  dt  checkview(0,1,2)
  cenv.SetValue("fedra.tag.DumpDoubletsTree"    , false );

  cenv.SetValue("mostag.outdir"          , "..");
  cenv.SetValue("mostag.env"             , "mostag.rootrc");
  cenv.SetValue("mostag.EdbDebugLevel"   , 1);
}

int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("mostagenv");
  gEDBDEBUGLEVEL     = cenv.GetValue("mostag.EdbDebugLevel" , 1);
  const char *env    = cenv.GetValue("mostag.env"            , "mostag.rootrc");
  const char *outdir = cenv.GetValue("mostag.outdir"         , "..");
  
  bool      do_single   = false;
  bool      do_set      = false;
  EdbID     id;
  int       from_fragment=0;
  int       n_fragments=0;
  
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
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
  }
  
  if(!(do_single||do_set))   { print_help_message(); return 0; }
  if(  do_single&&do_set )   { print_help_message(); return 0; }
  
  set_default(cenv);
  cenv.SetValue("mostag.env"            , env);
  cenv.ReadFile( cenv.GetValue("mostag.env"   , "mostag.rootrc") ,kEnvLocal);
  cenv.SetValue("mostag.outdir"         , outdir);
  
  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("mostag.outdir","..");
  cenv.WriteFile("mostag.save.rootrc");
  do_save_gif   = cenv.GetValue("fedra.tag.save_gif"            ,  false   );
  do_save_canvas = cenv.GetValue("fedra.tag.save_canvas"          ,  false   );
  
  printf("\n----------------------------------------------------------------------------\n");
  printf("mostag  %s\n"      ,id.AsString()	   );
  printf(  "----------------------------------------------------------------------------\n\n");
  
  
  if(do_single) 
  {
     EdbMosaicIO mio;
     mio.Init( Form("p%3.3d/%d.%d.%d.%d.mos.root",
		 id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor) );
  
     EdbMosaicIO omio;
     omio.Init( Form("p%3.3d/%d.%d.%d.%d.tag.root",
		  id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor),"RECREATE");
    TagSide(id, from_fragment, n_fragments, 1, cenv, mio,omio);
    TagSide(id, from_fragment, n_fragments, 2, cenv, mio,omio);
  }
  
  cenv.WriteFile("mostag.save.rootrc");
  return 1;
}


void DrawOut(const EdbPattern &p, EdbMosaicIO &omio)
{
  gROOT->SetBatch();
  TCanvas *c = new TCanvas(Form("peaks%d_%d_%d",p.Plate(), p.Side(), p.ID() ),
			   Form("peaks%d_%d_%d",p.Plate(), p.Side(), p.ID() ),
			   1100,600);
  c->Divide(2,1);
  c->cd(1)->SetGrid();
  gOH.xy->Draw("colz");
  gStyle->SetOptStat("n");
  DrawEllipse(p, 1);

  TVirtualPad *pad2 = c->cd(2);
  pad2->Divide(1,2);
  pad2->cd(1)->SetLogy();
  gOH.spectrOriginal->Draw();
  gOH.spectrSmooth->SetLineColor(2);
  gOH.spectrSmooth->Draw("same");
  gOH.spectrProc->SetLineColor(3);
  gOH.spectrProc->Draw("same");

  TVirtualPad *pdb = pad2->cd(2);
  pdb->Divide(2,1);
  pdb->cd(1);  gOH.dblxy->Draw("colz");
  gStyle->SetOptStat("n");
  pdb->cd(2);  gOH.dbltxty->Draw("colz");
  gStyle->SetOptStat("n");

  if(do_save_canvas) omio.SaveFragmentObj( c, p.Plate(), p.Side(), p.ID(), "peaks");
  if(do_save_gif) c->Print( omio.FileName( p.Brick(), p.Plate(), 0, 0, "", ".gif+30") );
  delete c;
}

void TagSide( EdbID id, int from, int nfrag, int side, TEnv &cenv, EdbMosaicIO &mio, EdbMosaicIO &omio )
{
  EdbScanProc sproc;
  sproc.eProcDirClient="..";
  
  EdbID id0=id; id0.ePlate=0;
  EdbScanSet *ss = sproc.ReadScanSet(id0);
  EdbPlateP *plate = ss->GetPlate(id.ePlate);
  
  EdbLayer *mapside = mio.GetCorrMap( id.ePlate, side );  
  int nc=mapside->Map().Ncell();
  if(from==0&&nfrag==0) nfrag=nc;
  
  for( int i=from; i<from+nfrag; i++ )
  {
    EdbPattern *p = mio.GetFragment( id.ePlate, side, i, true); // (plate,side,id)
    if(p) 
    {
      p->SetScanID(id);
      DoubletsFilterOut(*p, omio);
      
      EdbH2 h2p;
      Pat2H2(*p, h2p);
      //printf("Integral = %d \n",h2p.Integral());
      //gOH.xy = h2p.DrawH2("h","original mt");
      //omio.SaveFragmentObj( gOH.xy,  p->Plate(), p->Side(), p->ID(), "hxy");

      EdbPattern ptag;      ptag.SetScanID(id); ptag.SetSide(side); ptag.SetID(i);

      FindPeaks(h2p,ptag, omio);

      gPat.AddPattern(ptag);

      DrawOut(ptag,omio);
      SafeDelete(gOH.xy);
      SafeDelete(gOH.spectrOriginal);
      SafeDelete(gOH.spectrSmooth);
      SafeDelete(gOH.spectrProc);
      SafeDelete(gOH.dblxy);
      SafeDelete(gOH.dbltxty);
      delete p;      
    }
  }
  omio.SaveFragmentObj( &gPat,	id.ePlate, side, 0, "gpat");
}

void DoubletsFilterOut(EdbPattern &p, EdbMosaicIO &omio)
{
  int   checkview =1;
  int   fillhist  =1;
  float dr        =1;
  float dt        =0.008;
  TH2F *hxy=0, *htxty=0;
  if(fillhist)  {
    gOH.dblxy   = new TH2F("dblXY"  ,"Doublets DX DY"  ,50,-dr,dr,50,-dr,dr);
    gOH.dbltxty = new TH2F("dblTXTY","Doublets DTX DTY",50,-dt,dt,50,-dt,dt);
  }
  EdbAlignmentV adup;
  adup.eDVsame[0]=adup.eDVsame[1]= dr;
  adup.eDVsame[2]=adup.eDVsame[3]= dt;
  adup.FillGuessCell(p,p,1.);
  adup.FillCombinations();
  adup.DoubletsFilterOut(checkview, gOH.dblxy, gOH.dbltxty);   // assign flag -10 to the duplicated segments
  
  //omio.SaveFragmentObj(   hxy, p.Plate(), p.Side(), p.ID(), "dblxy");
  //omio.SaveFragmentObj( htxty, p.Plate(), p.Side(), p.ID(), "dbltxty");

  //if(eDoDumpDoubletsTree) DumpDoubletsTree(adup,"doublets");
  //SafeDelete(hxy);
  //SafeDelete(htxty);
}

void Pat2H2( const EdbPattern &p, EdbH2 &h2p)
{
  float xbin=50, ybin=50;
  float minx = p.Xmin();
  float maxx = p.Xmax();
  float miny = p.Ymin();
  float maxy = p.Ymax();
  int nx = (maxx-minx)/xbin;
  int ny = (maxy-miny)/ybin;
  
  h2p.InitH2(nx, minx, maxx, ny, miny, maxy);
  
  int n= p.N();
  for(int i=0; i<n; i++) 
  {
    EdbSegP *s = p.GetSegment(i);
    if(s->Flag()>=0) 
      if( abs(s->TX())<0.1 && abs(s->TY())<0.1 )
	h2p.Fill( s->X(),s->Y() );
  }
}

void FindPeaks(EdbH2 &h2p, EdbPattern &ptag, EdbMosaicIO &omio)
{
  int npmax=10;
  EdbPeak2 pf(h2p);
  gOH.spectrOriginal = pf.DrawSpectrumN( Form("spO%d_%d_%d",ptag.Plate(), ptag.Side(), ptag.ID()),"original spectrum");
  pf.Smooth();
  gOH.xy     = pf.DrawH2N( Form("xy%d_%d_%d",ptag.Plate(), ptag.Side(), ptag.ID()),Form("plate:%d side:%d id:%d microtracks xy plot",ptag.Plate(), ptag.Side(), ptag.ID()));
  gOH.spectrSmooth = pf.DrawSpectrumN( Form("spS%d_%d_%d",ptag.Plate(), ptag.Side(), ptag.ID()),"smoothed spectrum");

  int ir[2] = {2,2};
  pf.ProbPeaks(npmax);
  pf.Print();
  for(int i=0; i<npmax; i++) 
  {
    ptag.AddSegment( i, pf.eXpeak[i], pf.eYpeak[i],0,0, pf.ePeak[i]-pf.eMean[i] );
  }
  gOH.spectrProc = pf.DrawSpectrumN( Form("spP%d_%d_%d",ptag.Plate(), ptag.Side(), ptag.ID()),"processed spectrum");
}

void DrawEllipse(const EdbPattern &ptag, int col )
{
  int np = ptag.N();
  TText t(0,0,"a");
  for(int i=0; i<np; i++) {
    EdbSegP *s = ptag.GetSegment(i);
    TEllipse  *el = new TEllipse( s->X(), s->Y(), 300,300);
    el->SetFillStyle(0);
    el->SetLineColor(col);
    el->Draw();
    t.DrawText(s->X(), s->Y()+300, Form("%d",i) );
  }
}
