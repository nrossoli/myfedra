#ifndef ROOT_EdbDisplay
#define ROOT_EdbDisplay

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDisplay                                                           //
//                                                                      //
// Class to display pattern volume in 3D                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TObject.h>

#include "TCanvas.h"
#include "TSlider.h"
#include "TView.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TClonesArray.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "EdbPVRec.h"
#include "TButton.h"
#include "EdbAffine.h"


#include <TGClient.h>
#include <TGButton.h>

#include <TTree.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TDiamond.h>
#include <TArc.h>
#include <TSliderBox.h>
#include <TGaxis.h>
#include <TVirtualX.h>
#include <TMath.h>

//_________________________________________________________________________
class EdbTrackG : public TPolyMarker3D {
 private:

  EdbTrackP *eTr;

 public:
  EdbTrackG() {eTr=0;}
  EdbTrackG(Int_t nhits):TPolyMarker3D(nhits) {eTr=0;}
  virtual ~EdbTrackG(){}

  void SetTrack(EdbTrackP *tr) {eTr=tr;}

  //virtual void          ExecuteEvent(Int_t event, Int_t px, Int_t py);
  //virtual void          InspectParticle(); // *MENU*
  virtual void          DumpTrack();    // *MENU*
  virtual void          InspectTrack(); // *MENU*

 ClassDef(EdbTrackG,1)  //EdbTrack graphics
};

//_________________________________________________________________________
class EdbSegG : public TPolyLine3D {
 private:

  EdbSegP *eSeg;

 public:
  EdbSegG() {eSeg=0;}
  EdbSegG(Int_t nhits):TPolyLine3D(nhits) {eSeg=0;}
  virtual ~EdbSegG(){}

  void SetSeg(EdbSegP *s) {eSeg=s;}

  //virtual void          ExecuteEvent(Int_t event, Int_t px, Int_t py);
  //virtual void          InspectParticle(); // *MENU*
  virtual void          DumpSegment(); // *MENU*
  virtual void          InspectSegment(); // *MENU*

 ClassDef(EdbSegG,1) //
};

//_________________________________________________________________________
class EdbDisplay: public TObject {
private:
  TCanvas* fCanvas;
  TButton* fPickButton;
  TPad* fCutpad;
  TPad* fImagepad;
  EdbAffine2D* fAff;
  TView*   fView;
  TPolyLine3D* pl;
  TSlider  *fXSlider,*fYSlider,*fTXSlider,*fTYSlider;

  Float_t cutX0,cutX1,cutY0,cutY1;
  Float_t cutTX0,cutTX1,cutTY0,cutTY1;
  Float_t cutW,cutCHI2low,cutCHI2high,cutN1,cutN2;

  EdbPVRec *ePVR;
  TObjArray *eArrSegP; //array of segments to be drawn
  TObjArray *eArrTr;   //array of tracks to be drawn

  TArrayI* colors;
  TArrayF* DZs;

  Float_t vx0,vy0,vz0,vx1,vy1,vz1;
  Int_t eNsegmax,eNsegmin,eNpieces;
  Bool_t eTr_Co; // Tracks/Couples switch
protected:
   Int_t             fZoomMode;             //=1 if in zoom mode
   Bool_t            fDrawAllViews;         //Flag True if AllViews selected
   Bool_t            fDrawParticles;        //Flag True if particles to be drawn
   Bool_t            fDrawHits;             //Flag True if Hits to be drawn
   Float_t           fPTcut;                //PT cut to display objects
   Float_t           fTheta;                //Viewing angle theta
   Float_t           fPhi;                  //Viewing angle phi
   Float_t           fPsi;                  //Viewving angle psi (rotation on display)
   Float_t           fRrange;               //Size of view in R
   Float_t           fZrange;               //Size of view along Z
   Float_t           fZoomX0[20];           //Low x range of zoom number i
   Float_t           fZoomY0[20];           //Low y range of zoom number i
   Float_t           fZoomX1[20];           //High x range of zoom number i
   Float_t           fZoomY1[20];           //High y range of zoom number i
   Int_t             fZooms;                //Number of zooms
   Int_t             fHitsCuts;             //Number of hits surviving cuts
 //  TCanvas          *fCanvas;               //Pointer to the display canvas
   TPad             *fTrigPad;              //Pointer to the trigger pad
   TPad             *fCutPad;               //Pointer to the momentum cut slider pad
   TPad             *fEtaPad;               //Pointer to the rapidity cut slider pad
   TPad             *fButtons;              //Pointer to the buttons pad
   TPad             *fPad;                  //Pointer to the event display main pad
   TSlider          *fCutSlider;            //Momentum cut slider
   TSlider          *fEtaSlider;            //Rapidity slider
   TSlider          *fRangeSlider;          //Range slider
//   TButton          *fPickButton;           //Button to activate Pick mode
   TButton          *fZoomButton;           //Button to activate Zoom mode
   TArc             *fArcButton;            //Gren/Red button to show Pick/Zoom mode
   TList            *fFruits;               //List for fruits

public:
Bool_t MaxChi;
         EdbDisplay();
         EdbDisplay(const char *title, Float_t x0, Float_t x1 , Float_t y0, Float_t y1, Float_t z0, Float_t z1);
    void Clears();
    void SetNsegmax(Int_t v){eNsegmax=v;};
    void SetNsegmin(Int_t v){eNsegmin=v;};
    void SetShowTracks(){eTr_Co=kTRUE;};
    void SetShowCouples(){eTr_Co=kFALSE;};
    void SetCuts(Float_t x0, Float_t x1 , Float_t y0, Float_t y1){cutX0=x0;cutX1=x1;cutY0=y0;cutY1=y1;}

    void Refresh();
    void SetPVR(EdbPVRec *pvr) {ePVR=pvr;}
    void SetArrSegP(TObjArray *arr) {eArrSegP=arr;}
    void SetArrTr(TObjArray *arr) {eArrTr=arr;}
    EdbPVRec *PVR() const {return ePVR;}
    void PatternDraw(EdbPattern &pat);
    void TrackDraw(EdbTrackP *tr);
    EdbSegG *SegLine(EdbSegP *seg);
    
    virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py);
	 ~EdbDisplay(){};
   virtual void      DisplayButtons();
   virtual void      SetView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      DrawAllViews();
   virtual void      DrawViewGL();
   virtual void      DrawViewX3D();
   virtual void      DrawTitle(Option_t *option="");
   virtual void      DrawView(Float_t theta, Float_t phi, Float_t psi=0);
   virtual void      SetPickMode();
   virtual void      SetRange(Float_t x0, Float_t x1 , Float_t y0, Float_t y1, Float_t z0, Float_t z1);
   virtual void      SetZoomMode();
   virtual void      UnZoom(); // *MENU*
   virtual void      Draw(Option_t *option="");
   virtual Int_t     DistancetoPrimitive(Int_t px, Int_t py);


  ClassDef(EdbDisplay,1) //class to display OPERA emulsion data
};



#endif /* ROOT_EdbDisplay */
