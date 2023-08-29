
//Marwa
#ifndef ROOTManager_h
#define ROOTManager_h 1

#include "globals.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "G4AnalysisManager.hh"
#include "G4RootAnalysisManager.hh"
#include <TFile.h>
#include <TTree.h>
#include <mutex>
namespace B1{

class ROOTManager
{

  static const int MaxNGenPart = 500;
  static const int MaxNScintHit = 500;
  
  struct ROOTTreeStruct_t {
//Generator
    Int_t Event;
    Int_t  NScintHit;  
    Int_t NGenPart;
    Int_t GenPartID[MaxNGenPart];
    Int_t GenPartPDG[MaxNGenPart];
      Float_t GenPartE;
    //  Float_t GenPartE[MaxNGenPart];
    // Float_t GenPartTheta[MaxNGenPart];
    Float_t GenPartTheta;
    Float_t GenPartThetaMinus[MaxNGenPart];
    Float_t GenPartPhi[MaxNGenPart];


//G4Particle
    Int_t ScintHitPrimaryID[MaxNScintHit];
    Float_t ScintHitE[MaxNScintHit];
    Float_t ScintHitPosX[MaxNScintHit];
    Float_t ScintHitPosY[MaxNScintHit];
    Float_t ScintHitPosZ[MaxNScintHit];
    Int_t ScintHitStation[MaxNScintHit];
    Int_t ScintHitModule[MaxNScintHit];
    Int_t HitPDG[MaxNScintHit];
    Int_t ScintHitBar[MaxNScintHit];
    Float_t ScintHitKe[MaxNScintHit];
    Float_t ScintHitPz[MaxNScintHit];
    Float_t ScintHitPhi[MaxNScintHit];
    Float_t ScintHitPosXWo[MaxNScintHit];
    Float_t ScintHitSmear[MaxNScintHit];
    Float_t ScintHitPosYWo[MaxNScintHit];
    Float_t ScintHitPosZWo[MaxNScintHit];
            
    Int_t ScintHitTrackID[MaxNScintHit];
    Float_t ScintHitTime[MaxNScintHit];

    Float_t ScintHitTrackLength[MaxNScintHit];
//my edits separation of station
   Int_t NScintHitX1; 
   Int_t NScintHitX2; 
   Int_t NScintHitX3; 
   Int_t NScintHitX4;   

   Int_t NScintHitY1; 
   Int_t NScintHitY2; 
   Int_t NScintHitY3; 
   Int_t NScintHitY4;   


    Float_t ScintHitPosXStationX4[MaxNScintHit];
    Float_t ScintHitPosYStationX4[MaxNScintHit];
    Float_t ScintHitPosZStationX4[MaxNScintHit];
    Int_t ScintHitBarStationX4[MaxNScintHit];
    Int_t ScintHitModuleStationX4[MaxNScintHit];
    Float_t ScintHitEStationX4[MaxNScintHit];
    Float_t ScintHitKeStationX4[MaxNScintHit];

    Float_t ScintHitPosXStationY4[MaxNScintHit];
    Float_t ScintHitPosYStationY4[MaxNScintHit];
    Float_t ScintHitPosZStationY4[MaxNScintHit];
    Int_t ScintHitBarStationY4[MaxNScintHit];
    Int_t ScintHitModuleStationY4[MaxNScintHit];
    Float_t ScintHitEStationY4[MaxNScintHit];
    Float_t ScintHitKeStationY4[MaxNScintHit];


    Float_t ScintHitPosXStationX1[MaxNScintHit];
    Float_t ScintHitPosYStationX1[MaxNScintHit];
    Float_t ScintHitPosZStationX1[MaxNScintHit];
    Int_t ScintHitBarStationX1[MaxNScintHit];
    Int_t ScintHitModuleStationX1[MaxNScintHit];
    Float_t ScintHitEStationX1[MaxNScintHit];
    Float_t ScintHitKeStationX1[MaxNScintHit];

    Float_t ScintHitPosXStationY1[MaxNScintHit];
    Float_t ScintHitPosYStationY1[MaxNScintHit];
    Float_t ScintHitPosZStationY1[MaxNScintHit];
    Int_t ScintHitBarStationY1[MaxNScintHit];
    Int_t ScintHitModuleStationY1[MaxNScintHit];
    Float_t ScintHitEStationY1[MaxNScintHit];
    Float_t ScintHitKeStationY1[MaxNScintHit];

    Float_t ScintHitPosXStationX2[MaxNScintHit];
    Float_t ScintHitPosYStationX2[MaxNScintHit];
    Float_t ScintHitPosZStationX2[MaxNScintHit];
    Int_t ScintHitBarStationX2[MaxNScintHit];
    Int_t ScintHitModuleStationX2[MaxNScintHit];
    Float_t ScintHitEStationX2[MaxNScintHit];
    Float_t ScintHitKeStationX2[MaxNScintHit];

    Float_t ScintHitPosXStationY2[MaxNScintHit];
    Float_t ScintHitPosYStationY2[MaxNScintHit];
    Float_t ScintHitPosZStationY2[MaxNScintHit];
    Int_t ScintHitBarStationY2[MaxNScintHit];
    Int_t ScintHitModuleStationY2[MaxNScintHit];
    Float_t ScintHitEStationY2[MaxNScintHit];
    Float_t ScintHitKeStationY2[MaxNScintHit];


    Float_t ScintHitPosXStationX3[MaxNScintHit];
    Float_t ScintHitPosYStationX3[MaxNScintHit];
    Float_t ScintHitPosZStationX3[MaxNScintHit];
    Int_t ScintHitBarStationX3[MaxNScintHit];
    Int_t ScintHitModuleStationX3[MaxNScintHit];
    Float_t ScintHitEStationX3[MaxNScintHit];
    Float_t ScintHitKeStationX3[MaxNScintHit];

    Float_t ScintHitPosXStationY3[MaxNScintHit];
    Float_t ScintHitPosYStationY3[MaxNScintHit];
    Float_t ScintHitPosZStationY3[MaxNScintHit];
    Int_t ScintHitBarStationY3[MaxNScintHit];
    Int_t ScintHitModuleStationY3[MaxNScintHit];
    Float_t ScintHitEStationY3[MaxNScintHit];
    Float_t ScintHitKeStationY3[MaxNScintHit];



 




   




}; 
 

public:
  
  ROOTManager();
  ~ROOTManager();
  
  static ROOTManager* Instance();   
    
  void Init();		
  void Save();
  void Fill();

  struct ROOTTreeStruct_t ROOTTreeStruct;
//edit
  struct ROOTTreeStruct_t ROOTTreeStruct_G4Particle;
  struct ROOTTreeStruct_t ROOTTreeStructStn_1;
  struct ROOTTreeStruct_t ROOTTreeStructStn_2;
  struct ROOTTreeStruct_t ROOTTreeStructStn_3;
  struct ROOTTreeStruct_t ROOTTreeStructStn_4;
private:

  static ROOTManager* fgInstance;
  //mtsafe
   static std::mutex fgMutex; // Mutex for thread safety
  //mtsafe
  TFile* ROOTFile;        
  TTree* ROOTTree;
//edit
  TTree* ROOTTreeG4Particle;
  TTree* ROOTTreeDefault;
  TTree* ROOTTreeStn_1;
  TTree* ROOTTreeStn_2;
  TTree* ROOTTreeStn_3;
  TTree* ROOTTreeStn_4;
};

  //mtsafe
  
  //mtsafe

}
#endif

