
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
  static const int MaxNRPCHit = 500;
  
  struct ROOTTreeStruct_t {
//Generator
    Int_t Event;
    Int_t  NRPCHit;  
    Int_t NGenPart;
    Int_t GenPartID[MaxNGenPart];
    //Int_t GenPartPDG[MaxNGenPart];
    Int_t GenPartPDG[MaxNGenPart];
      Float_t GenPartE[MaxNGenPart];
      Float_t GenPartM;
    //  Float_t GenPartE[MaxNGenPart];
    // Float_t GenPartTheta[MaxNGenPart];
    Float_t GenPartTheta[MaxNGenPart];
    Float_t GenPartPhi[MaxNGenPart];
    Float_t GenPartThetaMinus[MaxNGenPart];
   // Float_t GenPartPhi[MaxNGenPart];


//G4Particle
    Int_t RPCHitPrimaryID[MaxNRPCHit];
    Float_t RPCHitE[MaxNRPCHit];
    Float_t RPCHitPosX[MaxNRPCHit];
    Float_t RPCHitPosY[MaxNRPCHit];
    Float_t RPCHitPosZ[MaxNRPCHit];
    Int_t RPCHitStation[MaxNRPCHit];
    Int_t RPCHitModule[MaxNRPCHit];
    Int_t HitPDG[MaxNRPCHit];
    Int_t RPCHitBar[MaxNRPCHit];
    Float_t RPCHitKe[MaxNRPCHit];
    Float_t RPCHitPz[MaxNRPCHit];
    Float_t RPCHitPhi[MaxNRPCHit];
    Float_t RPCHitPosXWo[MaxNRPCHit];
    Float_t RPCHitSmear[MaxNRPCHit];
    Float_t RPCHitPosYWo[MaxNRPCHit];
    Float_t RPCHitPosZWo[MaxNRPCHit];
            
    Int_t RPCHitTrackID[MaxNRPCHit];
    Float_t RPCHitTime[MaxNRPCHit];

    Float_t RPCHitTrackLength[MaxNRPCHit];
//my edits separation of station
   Int_t NRPCHitX1; 
   Int_t NRPCHitX2; 
   Int_t NRPCHitX3; 
   Int_t NRPCHitX4;   

   Int_t NRPCHitY1; 
   Int_t NRPCHitY2; 
   Int_t NRPCHitY3; 
   Int_t NRPCHitY4;   


    Float_t RPCHitPosXStationX4[MaxNRPCHit];
    Float_t RPCHitPosYStationX4[MaxNRPCHit];
    Float_t RPCHitPosZStationX4[MaxNRPCHit];
    Int_t RPCHitBarStationX4[MaxNRPCHit];
    Int_t RPCHitModuleStationX4[MaxNRPCHit];
    Float_t RPCHitEStationX4[MaxNRPCHit];
    Float_t RPCHitKeStationX4[MaxNRPCHit];

    Float_t RPCHitPosXStationY4[MaxNRPCHit];
    Float_t RPCHitPosYStationY4[MaxNRPCHit];
    Float_t RPCHitPosZStationY4[MaxNRPCHit];
    Int_t RPCHitBarStationY4[MaxNRPCHit];
    Int_t RPCHitModuleStationY4[MaxNRPCHit];
    Float_t RPCHitEStationY4[MaxNRPCHit];
    Float_t RPCHitKeStationY4[MaxNRPCHit];


    Float_t RPCHitPosXStationX1[MaxNRPCHit];
    Float_t RPCHitPosYStationX1[MaxNRPCHit];
    Float_t RPCHitPosZStationX1[MaxNRPCHit];
    Int_t RPCHitBarStationX1[MaxNRPCHit];
    Int_t RPCHitModuleStationX1[MaxNRPCHit];
    Float_t RPCHitEStationX1[MaxNRPCHit];
    Float_t RPCHitKeStationX1[MaxNRPCHit];

    Float_t RPCHitPosXStationY1[MaxNRPCHit];
    Float_t RPCHitPosYStationY1[MaxNRPCHit];
    Float_t RPCHitPosZStationY1[MaxNRPCHit];
    Int_t RPCHitBarStationY1[MaxNRPCHit];
    Int_t RPCHitModuleStationY1[MaxNRPCHit];
    Float_t RPCHitEStationY1[MaxNRPCHit];
    Float_t RPCHitKeStationY1[MaxNRPCHit];

    Float_t RPCHitPosXStationX2[MaxNRPCHit];
    Float_t RPCHitPosYStationX2[MaxNRPCHit];
    Float_t RPCHitPosZStationX2[MaxNRPCHit];
    Int_t RPCHitBarStationX2[MaxNRPCHit];
    Int_t RPCHitModuleStationX2[MaxNRPCHit];
    Float_t RPCHitEStationX2[MaxNRPCHit];
    Float_t RPCHitKeStationX2[MaxNRPCHit];

    Float_t RPCHitPosXStationY2[MaxNRPCHit];
    Float_t RPCHitPosYStationY2[MaxNRPCHit];
    Float_t RPCHitPosZStationY2[MaxNRPCHit];
    Int_t RPCHitBarStationY2[MaxNRPCHit];
    Int_t RPCHitModuleStationY2[MaxNRPCHit];
    Float_t RPCHitEStationY2[MaxNRPCHit];
    Float_t RPCHitKeStationY2[MaxNRPCHit];


    Float_t RPCHitPosXStationX3[MaxNRPCHit];
    Float_t RPCHitPosYStationX3[MaxNRPCHit];
    Float_t RPCHitPosZStationX3[MaxNRPCHit];
    Int_t RPCHitBarStationX3[MaxNRPCHit];
    Int_t RPCHitModuleStationX3[MaxNRPCHit];
    Float_t RPCHitEStationX3[MaxNRPCHit];
    Float_t RPCHitKeStationX3[MaxNRPCHit];

    Float_t RPCHitPosXStationY3[MaxNRPCHit];
    Float_t RPCHitPosYStationY3[MaxNRPCHit];
    Float_t RPCHitPosZStationY3[MaxNRPCHit];
    Int_t RPCHitBarStationY3[MaxNRPCHit];
    Int_t RPCHitModuleStationY3[MaxNRPCHit];
    Float_t RPCHitEStationY3[MaxNRPCHit];
    Float_t RPCHitKeStationY3[MaxNRPCHit];



 




   




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

