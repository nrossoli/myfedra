//-- Author :  Valeri Tioukov   9.06.2003

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDataSet                                                           //
//                                                                      //
// scanned raw data set and correspondent parameters                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "EdbDataSet.h"
#include "EdbSegment.h"

ClassImp(EdbSegmentCut)
ClassImp(EdbLayer)
ClassImp(EdbDataPiece)
ClassImp(EdbDataSet)
ClassImp(EdbDataProc)

///______________________________________________________________________________
EdbSegmentCut::EdbSegmentCut(int xi, float var[10])
{
  eXI=xi; 
  for(int i=0;i<5;i++) {
    eMin[i]=var[i*2]; 
    eMax[i]=var[i*2+1];
  }
}

///______________________________________________________________________________
int EdbSegmentCut::PassCut(float var[5])
{
  if     (eXI==0)  return PassCutX(var);
  else if(eXI==1)  return PassCutI(var);
  return 0;
}

///______________________________________________________________________________
int EdbSegmentCut::PassCutX(float var[5])
{
  // exclusive cut: if var is inside cut volume - return 0

  //printf("passcut:\n");
  for(int i=0; i<5; i++) {
    //printf("%f \t%f %f \n",var[i],eMin[i],eMax[i]);
    if(var[i]<eMin[i])  return 1;
    if(var[i]>eMax[i])  return 1;
  }
  return 0;
}

///______________________________________________________________________________
int EdbSegmentCut::PassCutI(float var[5])
{
  // inclusive cut: if var is inside cut volume - return 1

  for(int i=0; i<5; i++) {
    if(var[i]<eMin[i])  return 0;
    if(var[i]>eMax[i])  return 0;
  }
  return 1;
}

///______________________________________________________________________________
void EdbSegmentCut::Print()
{
  printf("min: "); for(int i=0;i<5;i++) printf("\t %f",eMin[i]); printf("\n");
  printf("max: "); for(int i=0;i<5;i++) printf("\t %f",eMax[i]); printf("\n");
}

///==============================================================================
void EdbLayer::Print()
{
  printf("Layer:\t %d\n",eID);
  printf("ZLAYER\t %f %f %f\n",eZ,eZmin,eZmax);
  printf("SHR\t %f\n",eShr);
  printf("AFFXY\t");        eAffXY.Print();
  printf("AFFTXTY\t");      eAffTXTY.Print();
}

///==============================================================================
EdbDataPiece::EdbDataPiece()
{
  Set0();
}

///______________________________________________________________________________
EdbDataPiece::EdbDataPiece(int plate, int piece, char* file, int flag)
{
  Set0();
  ePlate=plate;
  ePiece=piece;
  eFileName=file;
  eFlag=flag;
}

///______________________________________________________________________________
EdbDataPiece::~EdbDataPiece()
{
  for(int i=0; i<3; i++)  if(eLayers[i])   delete eLayers[i];
  for(int i=0; i<3; i++)  if(eAreas[i])    delete eAreas[i];
  for(int i=0; i<3; i++)  if(eCond[i])     delete eCond[i];
}

///______________________________________________________________________________
void EdbDataPiece::Set0()
{
  ePlate=0;
  ePiece=0;
  eFileName="";
  eFlag=0;

  for(int i=0; i<3; i++) eLayers[i]=0;
  for(int i=0; i<3; i++) eAreas[i]=0;
  for(int i=0; i<3; i++) eCond[i]=0;
  eRun    = 0;
}

///______________________________________________________________________________
void EdbDataPiece::Print()
{
  printf("Piece: %s\n",GetName());
  printf("%d %d %s %d\n", ePlate,ePiece,eFileName.Data(),eFlag);
  for(int i=0; i<3; i++)  if(eLayers[i])  eLayers[i]->Print();
  for(int i=0; i<3; i++)  if(eCond[i])    eCond[i]->Print();
  for(int i=0; i<3; i++)  
    if(eCuts[i])
      for(int j=0; j<NCuts(i); j++)  GetCut(i,j)->Print();
}

///______________________________________________________________________________
const char *EdbDataPiece::MakeName()
{
  char name[8];
  sprintf(name,"%2.2d_%3.3d", ePlate,ePiece);
  SetName(name);
  return GetName();
}


///______________________________________________________________________________
EdbLayer *EdbDataPiece::GetMakeLayer(int id)
{
  if(id<0) return 0;
  if(id>2) return 0;
  if(!GetLayer(id))  eLayers[id] = new EdbLayer();
  return GetLayer(id);
}

///______________________________________________________________________________
void EdbDataPiece::AddSegmentCut(int layer, int xi, float var[10])
{
  if(!eCuts[layer])  eCuts[layer] = new TObjArray();
  eCuts[layer]->Add( new EdbSegmentCut(xi,var) );
}

///______________________________________________________________________________
int EdbDataPiece::NCuts(int layer)
{
  if(!eCuts[layer])  return 0;
  return eCuts[layer]->GetEntries();
}

///______________________________________________________________________________
EdbScanCond *EdbDataPiece::GetMakeCond(int id)
{
  if(id<0) return 0;
  if(id>2) return 0;
  if(!GetCond(id))  eCond[id] = new EdbScanCond();
  return GetCond(id);
}

///______________________________________________________________________________
int EdbDataPiece::TakePiecePar(const char *dir)
{
  TString file=dir;
  file += MakeName();
  file += ".par";
  return ReadPiecePar(file);
}

///______________________________________________________________________________
int EdbDataPiece::ReadPiecePar(const char *file)
{
  char            buf[256];
  char            key[256];
  char            name[256];

  FILE *fp=fopen(file,"r");
  if (fp==NULL)   {
    printf("ERROR open file: %s \n", file);
    return(-1);
  }else
    printf( "\nRead piece parameters from file: %s\n\n", file );

  int id;
  float z,zmin,zmax,shr;
  float a11,a12,a21,a22,b1,b2;
  float x1,x2,x3,x4;
  float var[10];

  while( fgets(buf,256,fp)!=NULL ) {

    for( int i=0; i<(int)strlen(buf); i++ ) 
      if( buf[i]=='#' )  {
	buf[i]='\0';                       // cut out comments starting from #
	break;
      }

   if( sscanf(buf,"%s",key)!=1 )                             continue;

   if      ( !strcmp(key,"INCLUDE")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	ReadPiecePar(name);
      }
   else if ( !strcmp(key,"ZLAYER")   )
      {
	sscanf(buf+strlen(key),"%d %f %f %f",&id,&z,&zmin,&zmax);
	GetMakeLayer(id)->SetZlayer(z,zmin,zmax);
      }
    else if ( !strcmp(key,"SHRINK")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&shr);
	GetMakeLayer(id)->SetShrinkage(shr);
      }
    else if ( !strcmp(key,"AFFXY")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f",&id,&a11,&a12,&a21,&a22,&b1,&b2);
	GetMakeLayer(id)->SetAffXY(a11,a12,a21,a22,b1,b2);
      }
    else if ( !strcmp(key,"AFFTXTY")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f",&id,&a11,&a12,&a21,&a22,&b1,&b2);
	GetMakeLayer(id)->SetAffTXTY(a11,a12,a21,a22,b1,b2);
      }
    else if ( !strcmp(key,"SIGMA0")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f",&id,&x1,&x2,&x3,&x4);
	GetMakeCond(id)->SetSigma0(x1,x2,x3,x4);
      }
    else if ( !strcmp(key,"DEGRAD")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&x1);
	GetMakeCond(id)->SetDegrad(x1);
      }
    else if ( !strcmp(key,"BINS")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f",&id,&x1,&x2,&x3,&x4);
	GetMakeCond(id)->SetBins(x1,x2,x3,x4);
      }
    else if ( !strcmp(key,"RAMP0")  )
      {
	sscanf(buf+strlen(key),"%d %f %f",&id,&x1,&x2);
	GetMakeCond(id)->SetPulsRamp0(x1,x2);
      }
    else if ( !strcmp(key,"RAMP04")  )
      {
	sscanf(buf+strlen(key),"%d %f %f",&id,&x1,&x2);
	GetMakeCond(id)->SetPulsRamp04(x1,x2);
      }
    else if ( !strcmp(key,"CHI2MAX")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&x1);
	GetMakeCond(id)->SetChi2Max(x1);
      }
    else if ( !strcmp(key,"CHI2PMAX")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&x1);
	GetMakeCond(id)->SetChi2PMax(x1);
      }
    else if ( !strcmp(key,"XCUT")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f %f %f %f %f",&id,
	       var,var+1,var+2,var+3,var+4,var+5,var+6,var+7,var+8,var+9);
	AddSegmentCut(id,0,var);
      }
    else if ( !strcmp(key,"ICUT")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f %f %f %f %f",&id,
	       var,var+1,var+2,var+3,var+4,var+5,var+6,var+7,var+8,var+9);
	AddSegmentCut(id,1,var);
      }

  }

  fclose(fp);

  return 0;
}

///______________________________________________________________________________
int EdbDataPiece::PassCuts(int id, float var[5])
{
  int nc = NCuts(id);
  for(int i=0; i<nc; i++)
    if( !(GetCut(id,i)->PassCut(var)) )  return 0;
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::TakeRawSegment(EdbView *view, int id, EdbSegP &segP, int side)
{
  EdbSegment *seg = view->GetSegment(id);

  float var[5];
  var[0] = seg->GetX0();
  var[1] = seg->GetY0();
  var[2] = seg->GetTx();
  var[3] = seg->GetTy();
  var[4] = seg->GetPuls();

  if( !PassCuts(side,var) )     return 0;

  EdbLayer  *layer=GetLayer(side);
  float x,y,tx,ty,puls;

  tx   = seg->GetTx()/layer->Shr();
  ty   = seg->GetTy()/layer->Shr();
  x    = (seg->GetX0()+view->GetXview()); // + dz*tx;
  y    = (seg->GetY0()+view->GetYview()); // + dz*ty;
  puls = seg->GetPuls();

  segP.Set( seg->GetID(),x,y,tx,ty);
  segP.SetZ( layer->Z() );
  segP.SetW( puls );

  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::GetAreaData(EdbPVRec *ali, int aid, int side)
{
  TIndexCell *elist   = eAreas[side]->At(aid);
  if(!elist)  return 0;
  EdbPattern *pat = new EdbPattern( 0.,0., GetLayer(side)->Z() );
  pat->SetID(side);


  EdbSegP    segP;
  EdbView    *view = eRun->GetView();
  int   nseg=0, nrej=0;
  int   entry;
  int   nsegV;

  int niu=elist->N();
  for(int iu=0; iu<niu; iu++) {
    entry = elist->At(iu)->Value();
    view = eRun->GetEntry(entry);

    nsegV = view->Nsegments();

    for(int j=0;j<nsegV;j++) {
      if(!TakeRawSegment(view,j,segP,side)) {
	nrej++;
	continue;
      }
      nseg++;
      segP.SetVid(entry,j);
      pat->AddSegment( segP);
    }
  }

  printf("Area: %d ( %d%%)  %d \t viwes: %d \t nseg: %d \t rejected: %d\n", 
	 aid,100*aid/eAreas[side]->N(1),side,niu,nseg, nrej );
  ali->AddPattern(pat);
  return nseg;
}

///______________________________________________________________________________
int EdbDataPiece::MakeLinkListArea()
{
  if(!eRun)  eRun =  new EdbRun( eFileName.Data(),"READ" );
  for(int i=0; i<3; i++) {
    if(eAreas[i])    delete eAreas[i];
    eAreas[i]= new TIndexCell();
  }

  int nentr = eRun->GetEntries();
  printf("Make views entry map,  nentr = %d\n",nentr);

  Long_t v[2];   // areaID,entry
  EdbViewHeader *head=0;

  for(int iv=0; iv<nentr; iv++ ) {
    head = eRun->GetEntryHeader(iv);
    v[0]=head->GetAreaID();
    v[1]=iv;
    if(head->GetNframesTop()==0) {         // fill down views
      eAreas[2]->Add(2,v);
    }
    else if(head->GetNframesBot()==0) {    // fill up views
       eAreas[1]->Add(2,v);
    }
  }

  eAreas[1]->Sort();
  eAreas[2]->Sort();
  eAreas[1]->PrintStat();
  eAreas[2]->PrintStat();
  return eAreas[1]->N(1);
}

///______________________________________________________________________________
int EdbDataPiece::MakeLinkListCoord()
{
  if(!eRun)  eRun =  new EdbRun( eFileName.Data(),"READ" );
  for(int i=0; i<3; i++) {
    if(eAreas[i])    delete eAreas[i];
    eAreas[i]= new TIndexCell();
  }

  int nentr = eRun->GetEntries();
  printf("Make views coordinate map,  nentr = %d\n",nentr);

  TIndexCell upc;
  TIndexCell downc;
  Long_t v[3];   // x,y,entry
  EdbViewHeader *head=0;

  Long_t xx=0, yy=0;
  float cx = 2000., cy = 2000.;  // 2x2 mm cells
  float dx = 400. , dy = 400.;   // 400 microns margins
  float xv,yv;
  int mx[9] = {0, 0, 0,-1,  1, -1,-1, 1, 1};
  int my[9] = {0,-1, 1, 0,  0, -1, 1,-1, 1};

  for(int iv=0; iv<nentr; iv++ ) {

    head = eRun->GetEntryHeader(iv);
    yv = head->GetXview();
    xv = head->GetYview();
    xx = (Long_t)(xv/cx);
    yy = (Long_t)(yv/cx);
    v[0] = xx;
    v[1] = yy;
    v[2] = iv;

    if(head->GetNframesBot()==0) {              // fill up views
      upc.Add(3,v);
    }
    else if(head->GetNframesTop()==0) {         // fill down views
      downc.Add(3,v);                           // add center
      for( int im=1; im<5; im++ ) {             // add sides margins
        v[0] = (Long_t)(( xv+dx*mx[im] ) / cx);
        v[1] = (Long_t)(( yv+dy*my[im] ) / cy);
        if( (v[0] != xx) || (v[1] != yy) )
            downc.Add(3,v);
      }
      for( int im=5; im<9; im++ ) {             // add angles margins
        v[0] = (Long_t)(( xv+dx*mx[im] ) / cx);
        v[1] = (Long_t)(( yv+dy*my[im] ) / cy);
        if( (v[0] != xx) && (v[1] != yy) )
          //      if(!downc.Find(3,v))
            downc.Add(3,v);
      }

    }

  }

  upc.Sort();
  downc.Sort();

  //upc.SetName("x:y:entry");
  //downc.SetName("x:y:entry");
  //upc.PrintStat();
  //downc.PrintStat();

  TIndexCell *clx=0;
  TIndexCell *cly=0;

  int areac=0;
  int nix,niy,nie;
  nix=upc.N(1);
  for(int ix=0; ix<nix; ix++) {
    clx = upc.At(ix);
    xx = clx->Value();
    niy=clx->N(1);
    for(int iy=0; iy<niy; iy++) {
      cly  = clx->At(iy);
      yy = cly->Value();
      areac++;

      nie = cly->N(1);
      for(int ie=0; ie<nie; ie++) {
        v[0]=areac;
        v[1]  = cly->At(ie)->Value();
        eAreas[1]->Add(2,v);
      }

      cly = downc.Find(xx)->Find(yy);
      nie=cly->N(1);
      for(int ie=0; ie<nie; ie++) {
        v[0] = areac;
        v[1] = cly->At(ie)->Value();
        eAreas[2]->Add(2,v);
      }
    }
  }

  eAreas[1]->Sort();
  eAreas[2]->Sort();
  eAreas[1]->PrintStat();
  eAreas[2]->PrintStat();

  return eAreas[1]->N(1);
}

///==============================================================================
EdbDataSet::EdbDataSet()
{
  Set0();
}

///-------------------------------------------------------------------------------
EdbDataSet::EdbDataSet(const char *file)
{
  Set0();
  ReadDataSetDef(file);
}

///______________________________________________________________________________
EdbDataSet::~EdbDataSet()
{
  if(ePieces.GetEntries()) ePieces.Delete();
}

///______________________________________________________________________________
void EdbDataSet::Set0()
{
  eInputList  = "runs.lst";
  eAnaDir     = "./";
  eParDir     = "./";
  eDBFileName = "pieces_DB.root";
  eDBFile     = 0;
  ePieces.SetOwner();
}

///______________________________________________________________________________
int EdbDataSet::ReadDataSetDef(const char *file)
{
  char            buf[256];
  char            key[256];
  char            name[256];

  FILE *fp=fopen(file,"r");
  if (fp==NULL)   {
    printf("ERROR open file: %s \n", file);
    return(-1);
  }else
    printf( "\nRead Data Set Definitions from: %s\n\n", file );


  while( fgets(buf,256,fp)!=NULL ) {
    for( int i=0; i<(int)strlen(buf); i++ ) 
      if( buf[i]=='#' )  {
	buf[i]='\0';                       // cut out comments starting from #
	break;
      }

   if( sscanf(buf,"%s",key)!=1 )                             continue;

   if      ( !strcmp(key,"OUTPUT_DATA_DIR")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eAnaDir = name;
      }
   else if ( !strcmp(key,"INPUT_RUNS_LIST")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eInputList=name;
	GetRunList(name);
      }
   else if ( !strcmp(key,"PARAMETERS_DIR")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eParDir=name;
      }
   else if ( !strcmp(key,"DBFILNAME")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eDBFileName=name;
      }
  }
  fclose(fp);
  return 0;
}

///______________________________________________________________________________
void EdbDataSet::PrintRunList()
{
  for(int i=0; i<ePieces.GetEntries(); i++){
    ((EdbDataPiece*)ePieces.At(i))->Print();
  }
}

///______________________________________________________________________________
void EdbDataSet::WriteRunList()
{
  if(!eDBFile) eDBFile = new TFile(eDBFileName.Data(),"UPDATE");
  //else eDBFile->Cd();

  for(int i=0; i<ePieces.GetEntries(); i++){
    ((EdbDataPiece*)ePieces.At(i))->MakeName();
  }

  ePieces.Write("pieces",1);
}

///______________________________________________________________________________
int EdbDataSet::GetRunList(const char *file)
{
  char            buf[256];
  char            filename[256];
  int             nrun=0;

  FILE *fp=fopen(file,"r");
  if (fp==NULL)   {
    printf("ERROR open file: %s \n", file);
    return(-1);
  }else
    printf( "\nRead runs list from file: %s\n\n", file );

  EdbDataPiece *piece=0;
  int plateID,pieceID,flag;

  int ntok=0;
  while( fgets(buf,256,fp)!=NULL ) {
    ntok = sscanf(buf,"%d %d %s %d",&plateID,&pieceID,filename,&flag);
    if(ntok!=4) break;
    nrun++;
    piece = new EdbDataPiece(plateID,pieceID,filename,flag);
    ePieces.Add(piece);
  }
  fclose(fp);

  return nrun;
}

///==============================================================================
EdbDataProc::EdbDataProc(const char *file)
{
  eDataSet = new EdbDataSet(file);
}

///------------------------------------------------------------------------------
EdbDataProc::~EdbDataProc()
{
  if(eDataSet) delete eDataSet;
}

///------------------------------------------------------------------------------
int EdbDataProc::Process()
{
  if(!eDataSet) return 0;
  EdbDataPiece *piece;
  int np=eDataSet->N();
  for(int i=0; i<np; i++) {
    piece = eDataSet->GetPiece(i);
    piece->TakePiecePar(eDataSet->GetParDir());
    piece->Print();
    if(piece->Flag()==1)    Link(*piece);
  }
  return np;
}

///==============================================================================
TTree *EdbDataProc::InitCouplesTree(const char *tree_name,
				    const char *file_name,
				    const char *mode)
{
  TTree *tree=0;

  //tree = (TTree*)gROOT->FindObject(tree_name);
   if (!tree) {
      TFile *f = new TFile(file_name,mode);
      //if (f)  tree = (TTree*)f->Get(tree_name);
      if(!tree) {

	f->cd();
        tree = new TTree(tree_name,tree_name);

        int pid1,pid2;
        EdbSegCouple *cp=0;
        EdbSegP      *s1=0;
        EdbSegP      *s2=0;
        EdbSegP      *s=0;

        tree->Branch("pid1",&pid1,"pid1/I");
        tree->Branch("pid2",&pid2,"pid2/I");
        tree->Branch("cp","EdbSegCouple",&cp,32000,99);
        tree->Branch("s1.","EdbSegP",&s1,32000,99);
        tree->Branch("s2.","EdbSegP",&s2,32000,99);
        tree->Branch("s." ,"EdbSegP",&s,32000,99);
      }
   }
   tree->Write();
   tree->SetAutoSave(2000000);
   return tree;
}

///______________________________________________________________________________
void EdbDataProc::FillCouplesTree( TTree *tree, EdbPVRec *al, int fillraw=0 )
{
  tree->GetDirectory()->cd();

  EdbPatCouple *patc=0;

  int pid1,pid2;
  EdbSegCouple *cp=0;
  EdbSegP      *s1=0;
  EdbSegP      *s2=0;
  EdbSegP      *s=0;

  tree->SetBranchAddress("pid1",&pid1);
  tree->SetBranchAddress("pid2",&pid2);
  tree->SetBranchAddress("cp"  ,&cp);
  tree->SetBranchAddress("s1." ,&s1);
  tree->SetBranchAddress("s2." ,&s2);
  tree->SetBranchAddress("s."  ,&s );

  if(fillraw) {
    // **** fill tree with raw segments ****
    EdbPattern *pat=0;
    int nic;
    int nip=al->Npatterns();
    for( int ip=0; ip<nip; ip++ ) {
      pat  = al->GetPattern(ip);
      pid1 = pat->ID();
      pid2 = -1;
      nic=pat->N();
      for( int ic=0; ic<nic; ic++ ) {
        s1 = pat->GetSegment(ic);
        tree->Fill();
      }
    }
  }

  // **** fill tree with found couples ****
  s = new EdbSegP();

  int nip=al->Ncouples();
  for( int ip=0; ip<nip; ip++ ) {
    patc = al->GetCouple(ip);
    pid1 = patc->Pat1()->ID();
    pid2 = patc->Pat2()->ID();

    int nic=patc->Ncouples();
    for( int ic=0; ic<nic; ic++ ) {
      //      printf("%d %d\n",ip,ic);
      cp = patc->GetSegCouple(ic);
      s1 = patc->Pat1()->GetSegment(cp->ID1());
      s2 = patc->Pat2()->GetSegment(cp->ID2());
      //s = new EdbSegP(*s1);
      //*s += *s2;
        s->Set( ic, s1->X(), s1->Y(),
              (s1->X()-s2->X())/(s1->Z()-s2->Z()),
              (s1->Y()-s2->Y())/(s1->Z()-s2->Z()),
              s1->W()+s2->W());
        s->SetZ( s1->Z() );

//        s->Set( ic, (s1->X()+s2->X())/2., (s1->Y()+s2->Y())/2.,
//            (s1->TX()+s2->TX())/2., (s1->TY()+s2->TY())/2.,
//                    s1->W()+s2->W());

//        s->SetZ( (s1->Z()+s2->Z())/2. );

      tree->Fill();
      //delete s;
    }
  }

}

///______________________________________________________________________________
void EdbDataProc::CloseCouplesTree(TTree *tree)
{
  tree->AutoSave();
  tree->GetCurrentFile()->Close();
  //if(tree) delete tree;
}

///______________________________________________________________________________
int EdbDataProc::Link(EdbDataPiece &piece)
{
  EdbPVRec  *ali;

  TString name=eDataSet->GetAnaDir();
  if(!piece.GetName()) name+="link";
  else name+=piece.GetName();
  name+=".root";
  TTree *cptree = InitCouplesTree("couples",
  				  name.Data(),
  				  "RECREATE");  // "RECREATE" "NEW" or "UPDATE" open modes

  EdbScanCond *cond = piece.GetCond(1);

  int nareas = piece.MakeLinkListCoord();
  for(int i=0; i<nareas; i++ ) {

    ali      = new EdbPVRec();
    ali->SetScanCond( cond );

    printf("\n");
    piece.GetAreaData(ali,i,1);
    piece.GetAreaData(ali,i,2);

    ali->SetSegmentsErrors();

    ali->SetCouplesAll();
    ali->SetChi2Max(cond->Chi2PMax());
    ali->Link();

    FillCouplesTree(cptree, ali,0);

    delete ali;
  }
  CloseCouplesTree(cptree);

  return nareas;
}
