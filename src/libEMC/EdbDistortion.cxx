//-- Author :  Valeri Tioukov   22.06.2011

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDistortionCorr - use overlapped frames as "cormtx"                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "EdbLog.h"
#include "EdbCluster.h"
#include "EdbSegment.h"
#include "EdbDistortion.h"

using namespace TMath;

ClassImp(EdbClDist)
ClassImp(EdbDistortionCorr)
ClassImp(EdbDistortionMap)

//____________________________________________________________________________________
EdbDistortionCorr::EdbDistortionCorr()
{
  eNClMin      = 20;
  eR2CenterMax = 15.;
  eNCenterMin  =  1;
  eRmax        =  1.;
  eOutputFile  =  0;
  
  eCl.SetClass("EdbClDist",10); eCl.ExpandCreateFast(100000000);
  eGr.SetClass("EdbSegment",10); eGr.Expand(10000000);
  
  eXpix  = 0; //0.30625;
  eYpix  = 0; //0.30714;
  eNXpix = 0; //1280;
  eNYpix = 0; //1024;
  
  eCorrectionMatrixStepX = 10.; //microns - will be tuned in setpar
  eCorrectionMatrixStepY = 10.; //microns - will be tuned in setpar
  
  eAreaMin  =0;  eAreaMax  =90000000;
  eVolumeMin=0;  eVolumeMax=90000000;
  
}

//____________________________________________________________________________________
void EdbDistortionCorr::SetPar(TEnv &env)
{
  eNClMin      = env.GetValue("viewdist.NClMin" , 20);
  eR2CenterMax = env.GetValue("viewdist.R2CenterMax" , 15.);
  eRmax        = env.GetValue("viewdist.Rmax" ,  1.);
  
  eXpix  = env.GetValue("viewdist.Xpix" , 0.); // 0.30625);
  eYpix  = env.GetValue("viewdist.Ypix" , 0.); //0.30714);
  eNXpix = env.GetValue("viewdist.NXpix", 0);  //1280);
  eNYpix = env.GetValue("viewdist.NYpix", 0);  //1024);
  eDirX  = env.GetValue("viewdist.DirX",  0);  // -1,0,1
  eDirY  = env.GetValue("viewdist.DirY",  0);  // -1,0,1
  
  eCorrectionMatrixStepX = env.GetValue("viewdist.MatrixStepX", 10);
  eCorrectionMatrixStepY = env.GetValue("viewdist.MatrixStepY", 10);
  
  eDumpGr = env.GetValue("viewdist.DumpGr", 0);
  
  sscanf( env.GetValue("viewdist.ClusterAreaLimits", "0 90000000"), "%f %f", &eAreaMin, &eAreaMax);
  sscanf( env.GetValue("viewdist.ClusterVolumeLimits", "0 90000000"), "%f %f", &eVolumeMin, &eVolumeMax);
  
  printf("\n----------------------- Processing Parameters ---------------------------\n");
  printf("eNClMin\t %d\n",eNClMin);
  printf("eR2CenterMax\t %6.2f\n",eR2CenterMax);
  printf("eRmax\t %f6.2\n",eRmax);
  printf("Pixel:  %9.7f x %9.7f \n", eXpix, eYpix );
  printf("Matrix: %d x %d pixels, \t  %15.7f x %15.7f microns", eNXpix, eNYpix, eXpix*eNXpix, eYpix*eNYpix);
  printf("Clusters Area limits: %7.0f  %7.0f \n",  eAreaMin, eAreaMax );
  printf("Clusters Volume limits: %7.0f  %7.0f \n",  eVolumeMin, eVolumeMax );
  printf("-------------------------------------------------------------------------\n\n");
  
}

//____________________________________________________________________________________
void EdbDistortionCorr::InitCorrMap()
{
  eCorrMap.InitMap(eNXpix, eNYpix, eXpix, eYpix, eCorrectionMatrixStepX, eCorrectionMatrixStepY);
  eCorrMap0.InitMap(eNXpix, eNYpix, eXpix, eYpix, eCorrectionMatrixStepX, eCorrectionMatrixStepY);
  eCorrMapTot.InitMap(eNXpix, eNYpix, eXpix, eYpix, eCorrectionMatrixStepX, eCorrectionMatrixStepY);
}

//____________________________________________________________________________________
void EdbDistortionCorr::InitGMap()
{
  float bin=5;
  float  mi[2] = { -eNXpix*Abs(eXpix)/2-bin, -eNYpix*Abs(eYpix)/2-bin };
  float  ma[2] = {  1.5*eNXpix*Abs(eXpix)+bin,  1.5*eNYpix*Abs(eYpix)+bin };
  int    n[2] = { int((ma[0]-mi[0])/bin), int((ma[1]-mi[1])/bin) };
  eGMap.InitCell(250, n, mi, ma);
}

//____________________________________________________________________________________
void EdbDistortionCorr::AddCluster( EdbClDist *c)
{
  float v[2]={ c->eX+c->eXv,  c->eY+c->eYv };
  int  ir[2]={1,1};
  TObjArray arr;
  EdbSegment *s0 = 0;
  int n = eGMap.SelectObjectsC( v, ir, arr);
  if(n>100) printf("n=%d\n",n);
  if(n) { 
    float r, rmin=2.*eRmax;
    for(int i=0; i<n; i++) {
      EdbSegment *s = (EdbSegment *)arr.At(i);
      r = Sqrt( (s->GetX0()-v[0])*(s->GetX0()-v[0]) + (s->GetY0()-v[1])*(s->GetY0()-v[1]) );
      if(r<eRmax) if(r<rmin) { rmin=r; s0 = s; }
    }
  }
  if(s0) { s0->AddElement(c); s0->SetX0(v[0]); s0->SetY0(v[1]); s0->SetZ0(c->eZ); }
  else  {
    int inds = eGr.GetEntriesFast();
    EdbSegment *s = new(eGr[inds]) EdbSegment(v[0],v[1],c->eZ,0,0);
    s->AddElement(c);
    eGMap.AddObject(v, s);
  }
  
}

//____________________________________________________________________________________
void EdbDistortionCorr::CalculateGrRefMean()
{
  // take as grain reference position the mean of grains closer then eR2CenterMax to the view center
  
  int ngr = eGr.GetEntries();
  for(int i=0; i<ngr; i++) {
    EdbSegment *s = ((EdbSegment*)eGr.UncheckedAt(i));
    s->SetSide(0);
    int nc = s->GetNelements();
    double x0=0, y0=0, r0=0;
    int n0=0;
    for( int ic=0; ic<nc; ic++ ) {
      EdbClDist *c = (EdbClDist *)(s->GetElements()->UncheckedAt(ic));
      float r = Sqrt( c->eX*c->eX + c->eY*c->eY );                                     // distance to the view center
      if(r<=eR2CenterMax) { x0 += (c->eXv+c->eX); y0+=(c->eYv+c->eY); r0+=r; n0++; }
    }
    if( n0 >= eNCenterMin)   // select as a reference
    { 
      x0 /= n0;
      y0 /= n0;
      r0 /= n0;
      s->SetX0( x0 ); s->SetY0( y0 ); s->SetDz(r0);
      s->SetSide( n0 );   // use as a counter
    }
  }
}

//____________________________________________________________________________________
void EdbDistortionCorr::CutGrRef()
{
  // remove marginal grains from calculation
  int ngr = eGr.GetEntries();
  for(int i=0; i<ngr; i++) {
    EdbSegment *s = ((EdbSegment*)eGr.UncheckedAt(i));
    int nc = s->GetNelements();
    if(nc<eNClMin)              s->GetElements()->Clear();
    if(s->GetSide()<1)          s->GetElements()->Clear();
  }
}

//____________________________________________________________________________________
TNtuple *EdbDistortionCorr::DumpGr(const char *name)
{
  TNtuple *nt = new TNtuple( name,"grains dump","ig:n:xg:yg:nc:ic:xc:yc:xcv:ycv:dx:dy:w:ncenter:iscenter");
  int ngr = eGr.GetEntries();
  for(int ig=0; ig<ngr; ig++) {
    EdbSegment *s = ((EdbSegment*)eGr.UncheckedAt(ig));
    int nc = s->GetNelements();                                  if(!nc) continue;
    float dx=0, dy=0;
    for( int ic=0; ic<nc; ic++ ) {
      EdbClDist *c = (EdbClDist *)(s->GetElements()->UncheckedAt(ic));
      int j = eCorrMap.Jcell(c->eX,c->eY);
      int w = eCorrMap.Bin(j);
      if(w) {
	dx = eCorrMap.DX(j);
	dy = eCorrMap.DY(j);
      } else {dx=dy=0;}
      nt->Fill( ig, s->GetNelements(), s->GetX0(), s->GetY0(), //s->GetDz(), 
		nc, ic, c->eX, c->eY, c->eXv, c->eYv, dx,dy,w, s->GetSide(), c->eIsCenter);
    }
  }
  return nt;
}

//____________________________________________________________________________________
void EdbDistortionCorr::CalculateCorr()
{
  int ngr = eGr.GetEntries();
  Log(2,"EdbDistortionCorr::CalculateCorr","%d grains used",ngr);
  for(int i=0; i<ngr; i++) {
    EdbSegment *s = ((EdbSegment*)eGr.UncheckedAt(i));
    int nc = s->GetNelements();                                  if(!nc) continue;
    float x0 = s->GetX0(), y0= s->GetY0();
    for( int ic=0; ic<nc; ic++ ) {
      EdbClDist *c = (EdbClDist *)(s->GetElements()->UncheckedAt(ic));
      float dx = c->eX+c->eXv - x0;
      float dy = c->eY+c->eYv - y0;
      eCorrMap.Fill(c->eX, c->eY, dx,dy);
    }
  }
  eCorrMap.Norm();
}

//____________________________________________________________________________________
void EdbDistortionCorr::SaveCorrMap()
{
  eCorrMap.Save(eOutputFile, "diff");
  eCorrMapTot.Save(eOutputFile, "tot");

  bool batch = gROOT->IsBatch();
  gROOT->SetBatch();
  eGMap.DrawH2("eGMap","entries/bin")->Write();
  gROOT->SetBatch(batch);
  
  eCorrMapTot.GenerateCorrectionMatrix("correction_matrix.txt");
}
//____________________________________________________________________________________
void EdbDistortionCorr::Print()
{
  printf("\n-----------------------------------------------------------------------------------------------------\n");
  int ncl = eCl.GetEntries();
  int ngr = eGr.GetEntries();
  int icgr=0,  iccl=0;
  for(int i=0; i<ngr; i++) {
    int nc = ((EdbSegment*)eGr.At(i))->GetNelements();
    if(nc>=eNClMin) { icgr++;  iccl+=nc; }
  }
  printf("\n%d grains found\n",ngr);
  printf("\n%d grains selected with ncl > %d and closer then %.2f to view center \n",icgr, eNClMin, eR2CenterMax );
  printf("\n%d clusters used in selected grains, mean: %d clusters/grain \n",iccl, (int)(iccl/icgr) );
  printf("Matrix definition: %d %d  %f %f view size: %.2f %.2f\n", eNXpix, eNYpix, eXpix, eYpix, eNXpix*eXpix, eNYpix*eYpix );
  
  eGMap.PrintStat();
  //eCorrMap.PrintStat();
  printf("-----------------------------------------------------------------------------------------------------\n");
}

//____________________________________________________________________________________
int EdbDistortionCorr::MakeViewDirectionList( EdbRun &run, int dirx, int diry, TArrayI &entries)
{
  // Select views in desired direction dirx, diry {-1,0,1}
  //    -1 negative, +1 positive, 0 any direction
  
  struct {int id; float x; float y; } sv1,sv2;
  
  float tolerance = 2; // tolerance in microns
  EdbView *v=run.GetView();
  int nentr = run.GetEntries();
  entries.Set(nentr);
  int count=0;
  for(int ie=0; ie<nentr-1; ie++ ) {
    run.GetEntry(ie  , 1,0,0,0,0);
    sv1.id = v->GetViewID();
    sv1.x  = v->GetXview();
    sv1.y  = v->GetYview();
    run.GetEntry(ie+1  , 1,0,0,0,0);
    sv2.id = v->GetViewID();
    sv2.x  = v->GetXview();
    sv2.y  = v->GetYview();
    
    float dx = sv2.x-sv1.x;
    float dy = sv2.y-sv1.y;
    if(dirx) {
      if( abs(dx)<tolerance) continue;
      else if( dirx*dx < 0 ) continue;
    }
    if(diry) {
      if( abs(dy)<tolerance) continue;
      else if( diry*dy < 0 ) continue;
    }
    //printf("%d %d   %f %f\n", sv1.id, sv2.id, sv2.x-sv1.x, sv2.y-sv1.y);
    entries[count++]=ie+1;
  }
  return count;
}

//____________________________________________________________________________________
void EdbDistortionCorr::ReadClusters( const char *fname, TClonesArray &arr )
{
  EdbRun run(fname);
  TArrayI entries;
  int nsel = MakeViewDirectionList( run, eDirX, eDirY, entries );
  printf("nsel = %d\n", nsel);
  
  int      n = run.GetEntries();
  EdbView *v = run.GetView();
  run.GetEntry(0,1);
  int count=0;
  for(int i=0; i<nsel; i++) {
    run.GetEntry( entries[i] , 1,1);
    int ncl   = v->Nclusters();
    float xv  = v->GetXview();
    float yv  = v->GetYview();
    int view  = v->GetViewID();
    printf("%6d ", ncl);
    for(int ic=0; ic<ncl; ic++) {
      EdbCluster *c = v->GetCluster(ic);
      if(c->GetArea()>eAreaMin&&c->GetArea()<eAreaMax) {
	EdbClDist *cm = (EdbClDist*)(arr[count++]);
	cm->eX=c->eX; cm->eY=c->eY; cm->eZ=c->eZ; 
	cm->eXv=xv; cm->eYv=yv;
	cm->eView=view; cm->eFrame=c->GetFrame();
	AddCluster(cm);
      }
    }
  }
  printf("\nread %d clusters from %s\n",count,fname);
}

//____________________________________________________________________________________
void EdbDistortionCorr::MakeDistortionMap( const char *fname, TEnv &cenv, const char *usefile, const char *addfile )
{
  SetPar(cenv);
  
  InitGMap();
  InitCorrMap();
  
  ReadClusters( fname, eCl );
  
  if(usefile) {
    eCorrMap0.ReadMatrix2Map(usefile);
    eCorrMapTot.Add(eCorrMap0);
    int ncl = eCl.GetEntries();
    for(int i=0; i<ncl; i++) {
      eCorrMap0.ApplyCorr( *((EdbClDist*)(eCl[i])) );
    }
  }
  
  CalculateGrRefMean();
  CutGrRef();
  CalculateCorr();          // calculate the differential map eCorrMap
  eCorrMapTot.Add(eCorrMap);
  
  eOutputFile = new TFile("map.root","RECREATE");
  SaveCorrMap();
  if(eDumpGr) { TNtuple *ntgr = DumpGr("grcut"); ntgr->Write(); }
  eOutputFile->Close();

  //Print();
}

//========================= EdbDistortionMap class implementation =====================================================
void EdbDistortionMap::ApplyCorr(EdbClDist &c)
{
  int j=Jcell(c.eX,c.eY);
  if(j>-1) {
    TArrayD *dxy = (TArrayD*)( eMap.GetObject(j, 0 ) );
    if(dxy) {
      c.eX -= (*dxy)[0];
      c.eY -= (*dxy)[1];
    }
  }
}

void EdbDistortionMap::Save(TFile *file, const char *suffix)
{
  if(file) {
    bool batch = gROOT->IsBatch();
    gROOT->SetBatch();
    DrawCorrMap( file, Form("corrmap_%s",suffix) );
    GetH2dX()->Write( Form("hdx_%s",suffix) );
    GetH2dY()->Write( Form("hdy_%s",suffix) );
    gROOT->SetBatch(batch);
  }
}

void EdbDistortionMap::DrawCorrMap(TFile *file, const char *name)
{
  bool batch = gROOT->IsBatch();
  if(file) {
    Log(2,"EdbDistortionMap::DrawCorrMap","Save to file %s", file->GetName());
    gROOT->SetBatch();
  }
  
  TCanvas *cc = new TCanvas(Form("c_%s",name),"Corrections map",1200,900);
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  
  float margin=10;
  double minXborder = eMap.Xmin() - margin;
  double maxXborder = eMap.Xmax() + margin;
  double minYborder = eMap.Ymin() - margin;
  double maxYborder = eMap.Ymax() + margin;
  
  TH2F *hh = new TH2F(Form("hh_%s",name),"Corrections map",100,minXborder,maxXborder,100,minYborder,maxYborder);
  hh->GetXaxis()->SetTitle("X (#mum)");
  hh->GetYaxis()->SetTitle("Y (#mum)");
  hh->Draw();
  
  TBox *plate = new TBox(eMap.Xmin(),eMap.Ymin(),eMap.Xmax(),eMap.Ymax());
  plate->SetFillColor(16);
  plate->SetFillStyle(3001);
  plate->Draw();
  
  Double_t meanx=0, meany=0, meanr=0, wtot=0;
  float scale = 15.;
  int nc=eMap.Ncell();
  for( int i=0; i<nc; i++ ) {
    int w = eMap.Bin(i);    if(!w) continue;
    TArrayD *dxy = (TArrayD*)( eMap.GetObject(i, 0 ) );
    float x = eMap.Xj(i);
    float y = eMap.Yj(i);
    float dx = (*dxy)[0];
    float dy = (*dxy)[1];
    
    TArrow *arrow = new TArrow(x,y,x+scale*dx,y+scale*dy,0.01);
    arrow->SetLineWidth(1);
    arrow->Draw();
    
    meanx += dx;
    meany += dy;
    meanr += sqrt(dx*dx+dy*dy);
    wtot  += w;
  }
  
  meanx /= wtot;
  meany /= wtot;
  meanr /= wtot;
  printf("\nmeanR = %g   meanx = %g   meany = %g  wtot = %f\n", meanr, meanx, meany, wtot);

  if(file) {
   const char *mapname = "corr_map";
   if(name) mapname = name;
   eMap.DrawH2(Form("eMap_%s",name),"entries/bin")->Write();
   cc->Write(name);
  }
  gROOT->SetBatch(batch);
}


//____________________________________________________________________________________
void EdbDistortionMap::ReadMatrix2Map(const char *file)
{
  Log(1,"EdbDistortionMap::ReadMatrix2Map","%s ",file);
  FILE *f = fopen(file,"r");
  if(!f) { Log(1,"ReadMatrix2Map","file %s not found!",file); return; }
  char str[256]; 
  fgets(str, 256, f);
  int nxpix, nypix;
  float xpix, ypix;
  sscanf(str,"%d %d  %f %f", &nxpix, &nypix, &xpix, &ypix);
  if( nxpix!=eNXpix||nypix!=eNYpix||abs(xpix-eXpix)>0.000001||abs(ypix-eYpix)>0.000001 ) 
  { Log(1,"ReadMatrix2Map","ERROR: matrix definition is different!",file); 
    printf("%d %d  %f %f\n", nxpix, nypix, xpix, ypix);
    return; 
  }
  
  int i0,j0;
  float x0,y0,dx0,dy0;
  for(int i=0; i<eNXpix; i++) {
    for(int j=0; j<eNYpix; j++) {
      if(fgets(str, 256, f)==NULL) Log(1,"GenerateCorrectionMatrix","ERROR: input file is not correct");
      if( sscanf(str, "%d %d %f %f %f %f", &i0, &j0, &x0, &y0, &dx0, &dy0) == 6 )  
      {
	TArrayD *dxy = (TArrayD *)eMap.GetObject(x0,y0,0);
	(*dxy)[0]+=dx0;
	(*dxy)[1]+=dy0;
	eMap.Fill(x0, y0);
      }
    }    
  }
  fclose(f);
  
  int nc=eMap.Ncell();
  for( int i=0; i<nc; i++ ) {
    TArrayD *dxy = (TArrayD*)( eMap.GetObject(i, 0 ) );    if(!dxy) continue;
    int w = eMap.Bin(i);                                   if(!w) continue;
    (*dxy)[0] /= w;
    (*dxy)[1] /= w;
  }
}

//____________________________________________________________________________________
void EdbDistortionMap::Fill(float x, float y, float dx, float dy)
{
  TArrayD *dxy = (TArrayD *)eMap.GetObject(x,y,0);
  if(dxy) {
    (*dxy)[0]+=dx;
    (*dxy)[1]+=dy;
    eMap.Fill(x, y);
  }
}

void EdbDistortionMap::Norm()
{
  int nc=eMap.Ncell();
  for( int i=0; i<nc; i++ ) {
    TArrayD *dxy = (TArrayD*)( eMap.GetObject(i, 0 ) );    if(!dxy) continue;
    int w = eMap.Bin(i);                                   if(!w) continue;
    (*dxy)[0] /= w;
    (*dxy)[1] /= w;
  }
}

void EdbDistortionMap::SetDX( int j, double dx )
{
  TArrayD *dxy = (TArrayD*)( eMap.GetObject(j, 0 ) );    
  if(dxy) (*dxy)[0] = dx;
}

void EdbDistortionMap::SetDY( int j, double dy )
{
  TArrayD *dxy = (TArrayD*)( eMap.GetObject(j, 0 ) );    
  if(dxy) (*dxy)[1] = dy;
}

void EdbDistortionMap::AddDX( int j, double dx )
{
  TArrayD *dxy = (TArrayD*)( eMap.GetObject(j, 0 ) );    
  if(dxy) (*dxy)[0] += dx;
}

void EdbDistortionMap::AddDY( int j, double dy )
{
  TArrayD *dxy = (TArrayD*)( eMap.GetObject(j, 0 ) );    
  if(dxy) (*dxy)[1] += dy;
}

double EdbDistortionMap::DX(int j) const
{
  TArrayD *dxy = (TArrayD*)( eMap.GetObject(j, 0 ) );
  if(dxy) return (*dxy)[0]; else return 0;
}

double EdbDistortionMap::DY(int j) const
{
  TArrayD *dxy = (TArrayD*)( eMap.GetObject(j, 0 ) );
  if(dxy) return (*dxy)[1]; else return 0;
}

//____________________________________________________________________________________
void EdbDistortionMap::Smooth(int n, Option_t *opt)
{
  TH2F *h2x = GetH2dX();
  TH2F *h2y = GetH2dY();
  for(int i=0; i<n; i++) {
    h2x->Smooth(1,opt);
    h2y->Smooth(1,opt);
  }
  PutH2dX(*h2x);
  PutH2dY(*h2y);
}

//____________________________________________________________________________________
void EdbDistortionMap::PutH2dX( const TH2F &h2 )
{
  // TODO: no any checks yet
  int nc = h2.GetNcells();
  int ix,iy,iz;
  for( int j=0; j<nc; j++ ) {
    double dx = h2.GetBinContent(j);
    h2.GetBinXYZ(j,ix,iy,iz);
    float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
    float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
    int jm = eMap.Jcell(x,y);
    SetDX( jm, dx );
  }
}

//____________________________________________________________________________________
void EdbDistortionMap::PutH2dY( const TH2F &h2 )
{
  // TODO: no any checks yet
  int nc = h2.GetNcells();
  int ix,iy,iz;
  for( int j=0; j<nc; j++ ) {
    double dy = h2.GetBinContent(j);
    h2.GetBinXYZ(j,ix,iy,iz);
    float x = ((TAxis*)h2.GetXaxis())->GetBinCenter(ix);
    float y = ((TAxis*)h2.GetYaxis())->GetBinCenter(iy);
    int jm = eMap.Jcell(x,y);
    SetDY( jm, dy );
  }
}

//____________________________________________________________________________________
TH2F *EdbDistortionMap::GetH2dX(const char *name)
{
  TH2F *hd = eMap.DrawH2(name,"CorrectionMatrix dX");
  hd->Reset();
  int nc=eMap.Ncell();
  for( int j=0; j<nc; j++ )
    hd->Fill( eMap.Xj(j), eMap.Yj(j), DX(j) );
  return hd;
}

//____________________________________________________________________________________
TH2F *EdbDistortionMap::GetH2dY(const char *name)
{
  TH2F *hd = eMap.DrawH2(name,"CorrectionMatrix dY");
  hd->Reset();
  int nc=eMap.Ncell();
  for( int j=0; j<nc; j++ )
    hd->Fill( eMap.Xj(j), eMap.Yj(j), DY(j) );
  return hd;
}

//____________________________________________________________________________________
void EdbDistortionMap::Add(const EdbDistortionMap &map, float k)
{
  int nc=eMap.Ncell();
  for( int j=0; j<nc; j++ ) {
    int jm = map.Jcell( eMap.Xj(j), eMap.Yj(j) );
    TArrayD *dxy = (TArrayD *)eMap.GetObject(j,0);
    (*dxy)[0] += k*map.DX(jm);
    (*dxy)[1] += k*map.DY(jm);
    eMap.Fill(j);
  }
}

//____________________________________________________________________________________
void EdbDistortionMap::Scale(const float k)
{
  int nc=eMap.Ncell();
  for( int j=0; j<nc; j++ ) {
    SetDX( j, DX(j)*k );
    SetDY( j, DY(j)*k );
  }
}
//____________________________________________________________________________________
void EdbDistortionMap::InitMap(int nxpix,int nypix, float xpix, float ypix, float stepX, float stepY)
{
  // matrix of the corrections in each entry is TArrayD with 2 entries:(dx,dy)
  eNXpix = nxpix;
  eNYpix = nypix;
  eXpix  = xpix;
  eYpix  = ypix;
  
  float  mi[2] = { -eNXpix*Abs(eXpix)/2., -eNYpix*Abs(eYpix)/2. };
  float  ma[2] = {  eNXpix*Abs(eXpix)/2.,  eNYpix*Abs(eYpix)/2. };
  int    n[2] = { int((ma[0]-mi[0])/stepX), int((ma[1]-mi[1])/stepY) };
  
  stepX = (ma[0]-mi[0])/n[0];
  stepY = (ma[1]-mi[1])/n[1];
  n[0] = int((ma[0]-mi[0]+1.)/stepX);
  n[1] = int((ma[1]-mi[1]+1.)/stepY);
  
  eMap.InitCell(2, n, mi, ma);
  //eMap.PrintStat();
  
  int nc=eMap.Ncell();
  for( int i=0; i<nc; i++ ) {
    TArrayD *a = new TArrayD(2);
    eMap.AddObject(i, (TObject*)a );
    eMap.SetBin(i,0);
  }
}

//____________________________________________________________________________________
void EdbDistortionMap::GenerateCorrectionMatrix(const char *file)
{
  int nentries = eMap.Ncell();
  Double_t *vx   = new Double_t[nentries];
  Double_t *vy   = new Double_t[nentries];
  Double_t *vdx  = new Double_t[nentries];
  Double_t *vdy  = new Double_t[nentries];
  
  for(int i=0; i<nentries; i++) {
    int w = eMap.Bin(i);    if(!w) continue;
    TArrayD *arr = (TArrayD*)eMap.GetObject(i,0);
    if(!arr) { Log(1,"EdbDistortionMap::GenerateCorrectionMatrix","ERROR: missed bin");  break; }
    vx[i]   = eMap.Xj(i);
    vy[i]   = eMap.Yj(i);
    vdx[i]  = arr->At(0);
    vdy[i]  = arr->At(1);
  }
  
  TGraph2D *gdx  = new TGraph2D("graphDX","graphDX",nentries,vx,vy,vdx);
  TGraph2D *gdy  = new TGraph2D("graphDY","graphDY",nentries,vx,vy,vdy);
  
  FILE *fmatr = fopen(file,"w");
  fprintf(fmatr,"%d %d  %f %f\n", eNXpix, eNYpix, eXpix, eYpix);
  
  for(int i=0; i<eNXpix; i++) {
    for(int j=0; j<eNYpix; j++) {
      float x = (i - 0.5*(eNXpix-1))*eXpix;
      float y = (j - 0.5*(eNYpix-1))*eYpix;
      //float dx = gdx->Interpolate(x,y);
      //float dy = gdy->Interpolate(x,y);
      float dx = ((TArrayD *)eMap.GetObject(x,y,0))->At(0);
      float dy = ((TArrayD *)eMap.GetObject(x,y,0))->At(1);
      fprintf(fmatr,"%5d %5d %12.6f %12.6f  %11.6f %11.6f\n", i, j, x,y,dx,dy);
    }
  }
  fclose(fmatr);
}

