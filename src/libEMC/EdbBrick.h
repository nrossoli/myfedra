#ifndef ROOT_EdbBrick
#define ROOT_EdbBrick

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbBrick - OPERA Brick structure definition                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TObjArray.h"
#include "EdbLayer.h"

//______________________________________________________________________________
class EdbPlateP : public EdbLayer {

 private:
  TObjArray eLayers;      // 0-base, 1-up, 2-down

 public:
  EdbPlateP();
  ~EdbPlateP(){}

  void SetPlateLayout(float z0, float z1, float z2);
  EdbLayer *GetLayer(int i) {if(i<0) return 0; return (EdbLayer*)eLayers.At(i);}
  void SetDXDY(float dx, float dy);

  ClassDef(EdbPlateP,1)  // OPERA Plate structure definition
};

//______________________________________________________________________________
class EdbBrickP : public EdbLayer {

 private:
  TObjArray ePlates;
  TObjArray eSpacers;

 public:
  EdbBrickP();
  ~EdbBrickP(){}

  void SetPlatesLayout(float z0, float z1, float z2);
  void SetDXDY(float dx, float dy);
  int Npl() const {return ePlates.GetEntries();}
  void AddPlate(EdbPlateP *pl) { ePlates.Add(pl); }
  EdbPlateP *GetPlate(int i) {return (EdbPlateP*)ePlates.At(i);}
  void Print();

  ClassDef(EdbBrickP,1)  // OPERA Brick structure definition
};


#endif /* ROOT_EdbBrick */

