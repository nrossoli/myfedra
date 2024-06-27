//-- Author :  Valeri Tioukov   15/02/2024

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbMosaicIO                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TCanvas.h>
#include <TH2F.h>
#include "EdbMosaicIO.h"
#include "EdbLog.h"

ClassImp(EdbMosaicIO);

//-----------------------------------------------------------------------
void EdbMosaicIO::Init( const char *file, Option_t* option)
{
  eFile = new TFile( file, option );
}

//-----------------------------------------------------------------------
char *EdbMosaicIO::FileName(int brick, int plate, int major, int minor, const char *pref, const char *suff)
{
  return Form("p%3.3d/%s%d.%d.%d.%d%s",plate,pref,brick,plate,major,minor,suff);
}

//-----------------------------------------------------------------------
void EdbMosaicIO::SaveFragment(EdbPattern &p)
{
  if(eFile) 
  {
    eFile->cd();
    p.Write( Form("p%d_%d_%d", p.Plate(), p.Side(), p.ID() ) );
  }
}

//-----------------------------------------------------------------------
void EdbMosaicIO::SaveFragmentTag(TObject *ob, int plate, int side, int id, const char *pref)
{
  if(eFile) 
  {
    eFile->cd();
    ob->Write( Form("%s%d_%d_%d", pref,plate, side, id) );
  }
}

//-----------------------------------------------------------------------
EdbPattern *EdbMosaicIO::GetFragment(int plate, int side, int id, bool do_corr )
{
  EdbLayer *mapside = GetCorrMap( plate, side );
  EdbLayer *l = mapside->Map().GetLayer( id );
  char *name = Form("p%d_%d_%d", plate, side, id);
  Log(2,"EdbMosaicIO::GetFragment","%s",name);
  if(eFile)
  {
    EdbPattern *p = (EdbPattern *)(eFile->Get( name ));
    p->SetSide(side);
    p->SetID(id);
    if(do_corr) 
      if(l) 
      {
	p->Transform(l->GetAffineXY());
	p->TransformA(l->GetAffineTXTY());
	Log(2,"EdbMosaicIO::GetFragment","AffXY  : %s", l->GetAffineXY()->AsString() );
 	Log(2,"EdbMosaicIO::GetFragment","AffTXTY: %s", l->GetAffineTXTY()->AsString() );
      }
    return p;
  }
  else return 0;
}

//-----------------------------------------------------------------------
void EdbMosaicIO::SaveCorrMap(int plate, int side, EdbLayer &l)
{
  if(eFile)
  {
    eFile->cd();
    l.Write( Form("map_p%d_%d", plate,side ) );
  }
}

//-----------------------------------------------------------------------
void EdbMosaicIO::SaveCorrMap(int plate, int side, EdbLayer &l, const char *file)
{
  TFile f(file, "UPDATE");
  l.Write( Form("map_p%d_%d", plate,side ) );
  f.Close();
}

//-----------------------------------------------------------------------
EdbLayer *EdbMosaicIO::GetCorrMap(int plate, int side)
{
  if(eFile) 
    return (EdbLayer *)(eFile->Get( Form("map_p%d_%d", plate, side) ));
  else 
    return 0;
}

//-----------------------------------------------------------------------
void EdbMosaicIO::DrawFragment(EdbPattern &p)
{
  TH2F *h2xy = new TH2F("hxy","hxy",1000, -10000, 10000, 1000, -10000,10000);
  for(int i=0; i<p.N(); i++)
  {
    EdbSegP *s = p.GetSegment(i);
    h2xy->Fill( s->X(), s->Y() );
  }
  TCanvas *c = new TCanvas("cdf","cdf",800,800);
  //c->cd(1);
  h2xy->Draw("colz");
}

