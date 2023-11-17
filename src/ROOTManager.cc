
//Marwa
#include "ROOTManager.hh"

#include <TROOT.h>
#include <TFile.h>
#include "TBranch.h"
#include "TTree.h"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4UnitsTable.hh"
//#include"MuravesMessenger.hh"

namespace B1{
  // ROOTManager* ROOTManager::fgInstance = 0;
  ROOTManager* ROOTManager::fgInstance = nullptr;
  std::mutex ROOTManager::fgMutex;
  
  ROOTManager* ROOTManager::Instance()
{std::lock_guard<std::mutex> lock(fgMutex); // Lock the mutex

  if (fgInstance == nullptr) {
        fgInstance = new ROOTManager();
    }
  return fgInstance;
}      
  
ROOTManager::ROOTManager()
{ 
  fgInstance = this;
}

ROOTManager::~ROOTManager()
{
  // if ( ROOTFile ) delete ROOTFile;
  fgInstance = 0;  
}

void ROOTManager::Init()
{ 
 // Creating a tree container to handle histograms and ntuples.
 // This tree is associated to an output file.
 //
  // create ROOT file (dove alla fine salvare dati e grafici)
  //auto Label = MuravesMessenger::Instance()->GetRunLabel();
  //G4String input = "Muraves"+ Label + ".root";
  ROOTFile = new TFile("RPC.root", "RECREATE");
  if (!ROOTFile) {
    G4cout << " problem creating the ROOT TFile" << G4endl;
    return;
  }

  // creazione ttree p_incident_TOT
  ROOTTree = new TTree("Muraves_Generator", "Muraves_Generator");
   ROOTTreeG4Particle = new TTree("Muraves_G4Particle", "Muraves_G4Particle");
  ROOTTreeDefault = new TTree("Default", "Default");
  //tree for separetad station
  /*ROOTTreeStn_1 = new TTree("Station1", "Station1");
  ROOTTreeStn_2 = new TTree("Station2", "Station2");
  ROOTTreeStn_3 = new TTree("Station3", "Station3");
  ROOTTreeStn_4 = new TTree("Station4", "Station4");
  */




//Generator
  ROOTTree->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
   ROOTTree->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE[NGenPart]/F");
  ROOTTree->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta[NGenPart]/F");
 
  /* ROOTTree->Branch("Event", &ROOTTreeStruct.Event, "Event/I");
  ROOTTree->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
  ROOTTree->Branch("GenPartID", &ROOTTreeStruct.GenPartID, "GenPartID[NGenPart]/I");
 ROOTTree->Branch("GenPartPDG", &ROOTTreeStruct.GenPartPDG, "GenPartPDG[NGenPart]/I");
  ROOTTree->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE[NGenPart]/F");
  ROOTTree->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta[NGenPart]/F");
  ROOTTree->Branch("GenPartPhi", &ROOTTreeStruct.GenPartPhi, "GenPartPhi[NGenPart]/F");
   ROOTTree->Branch("NRPCHit", &ROOTTreeStruct.NRPCHit, "NRPCHit/I");
  */

//G4Particule
   ROOTTreeG4Particle->Branch("NRPCHit", &ROOTTreeStruct.NRPCHit, "NRPCHit/I");
   ROOTTreeG4Particle->Branch("RPCHitE", &ROOTTreeStruct.RPCHitE, "RPCHitE[NRPCHit]/F");
    ROOTTreeG4Particle->Branch("RPCHitPosX", &ROOTTreeStruct.RPCHitPosX, "RPCHitPosX[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitPosY", &ROOTTreeStruct.RPCHitPosY, "RPCHitPosY[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitPosZ", &ROOTTreeStruct.RPCHitPosZ, "RPCHitPosZ[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitE", &ROOTTreeStruct.RPCHitPz, "RPCHitPz[NRPCHit]/F");
   ROOTTreeG4Particle->Branch("RPCHitStation", &ROOTTreeStruct.RPCHitStation, "RPCHitStation[NRPCHit]/I");
    ROOTTreeG4Particle->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NRPCHit]/I");
   /*
  ROOTTreeG4Particle->Branch("RPCHitPrimaryID", &ROOTTreeStruct.RPCHitPrimaryID, "RPCHitPrimaryID[NRPCHit]/I");
  ROOTTreeG4Particle->Branch("RPCHitE", &ROOTTreeStruct.RPCHitE, "RPCHitE[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitPosX", &ROOTTreeStruct.RPCHitPosX, "RPCHitPosX[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitPosY", &ROOTTreeStruct.RPCHitPosY, "RPCHitPosY[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitPosZ", &ROOTTreeStruct.RPCHitPosZ, "RPCHitPosZ[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitStation", &ROOTTreeStruct.RPCHitStation, "RPCHitStation[NRPCHit]/I");
  ROOTTreeG4Particle->Branch("RPCHitModule", &ROOTTreeStruct.RPCHitModule, "RPCHitModule[NRPCHit]/I");
  ROOTTreeG4Particle->Branch("RPCHitBar", &ROOTTreeStruct.RPCHitBar, "RPCHitBar[NRPCHit]/I");
  ROOTTreeG4Particle->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NRPCHit]/I");
  ROOTTreeG4Particle->Branch("RPCHitTime", &ROOTTreeStruct.RPCHitTime, "RPCHitTime[NRPCHit]/F");

*/

  
 // ROOTTreeG4Particle->Branch("RPCHitKe", &ROOTTreeStruct.RPCHitKe, "RPCHitKe[NRPCHit]/F");  
 // ROOTTreeG4Particle->Branch("RPCHitPhi",&ROOTTreeStruct.RPCHitPhi,"RPCHitPhi[NRPCHit]/F");
  //ROOTTreeG4Particle->Branch("RPCHitPz",&ROOTTreeStruct.RPCHitPz,"RPCHitPz[NRPCHit]/F");


   /*
  ROOTTreeG4Particle->Branch("RPCHitPosXWo", &ROOTTreeStruct.RPCHitPosXWo, "RPCHitPosXWo[NRPCHit]/F");
  ROOTTreeG4Particle->Branch("RPCHitPosYWo", &ROOTTreeStruct.RPCHitPosYWo, "RPCHitPosYWo[NRPCHit]/F");
   */

  
 // ROOTTreeG4Particle->Branch("RPCHitSmear", &ROOTTreeStruct.RPCHitSmear, "RPCHitSmear[NRPCHit]/F");  
//Default
   // ROOTTreeDefault->Branch("Event", &ROOTTreeStruct.Event, "Event/I");
     ROOTTreeDefault->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
     ROOTTreeDefault->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE[NGenPart]/F");
     ROOTTreeDefault->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta[NGenPart]/F");
     ROOTTreeDefault->Branch("GenPartPhi", &ROOTTreeStruct.GenPartPhi, "GenPartPhi[NGenPart]/F");
     ROOTTreeDefault->Branch("GenPartPDG", &ROOTTreeStruct.GenPartPDG, "GenPartPDG[NGenPart]/I");
      ROOTTreeDefault->Branch("GenPartM", &ROOTTreeStruct.GenPartM, "GenPartM/F");


      ROOTTreeDefault->Branch("Event", &ROOTTreeStruct.Event, "Event/I");
      ROOTTreeDefault->Branch("NRPCHit", &ROOTTreeStruct.NRPCHit, "NRPCHit/I");
      ROOTTreeDefault->Branch("RPCHitPrimaryID", &ROOTTreeStruct.RPCHitPrimaryID, "RPCHitPrimaryID[NRPCHit]/I");
      ROOTTreeDefault->Branch("RPCHitTrackID", &ROOTTreeStruct.RPCHitTrackID, "RPCHitTrackID[NRPCHit]/I");
      ROOTTreeDefault->Branch("RPCHitTrackLength", &ROOTTreeStruct.RPCHitTrackLength, "RPCHitTrackLength[NRPCHit]/I");
      ROOTTreeDefault->Branch("RPCHitE", &ROOTTreeStruct.RPCHitE, "RPCHitE[NRPCHit]/F");
      ROOTTreeDefault->Branch("RPCHitPz", &ROOTTreeStruct.RPCHitPz, "RPCHitPz[NRPCHit]/F");
      ROOTTreeDefault->Branch("RPCHitPosX", &ROOTTreeStruct.RPCHitPosX, "RPCHitPosX[NRPCHit]/F");
      ROOTTreeDefault->Branch("RPCHitPosY", &ROOTTreeStruct.RPCHitPosY, "RPCHitPosY[NRPCHit]/F");
      ROOTTreeDefault->Branch("RPCHitPosZ", &ROOTTreeStruct.RPCHitPosZ, "RPCHitPosZ[NRPCHit]/F");
      ROOTTreeDefault->Branch("RPCHitStation", &ROOTTreeStruct.RPCHitStation, "RPCHitStation[NRPCHit]/I");
      ROOTTreeDefault->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NRPCHit]/I");
      
   /*
  ROOTTreeDefault->Branch("Event", &ROOTTreeStruct.Event, "Event/I");
  ROOTTreeDefault->Branch("NGenPart", &ROOTTreeStruct.NGenPart, "NGenPart/I");
  ROOTTreeDefault->Branch("GenPartID", &ROOTTreeStruct.GenPartID, "GenPartID[NGenPart]/I");
  ROOTTreeDefault->Branch("GenPartPDG", &ROOTTreeStruct.GenPartPDG, "GenPartPDG[NGenPart]/I");
  ROOTTreeDefault->Branch("GenPartE", &ROOTTreeStruct.GenPartE, "GenPartE[NGenPart]/F");
  ROOTTreeDefault->Branch("GenPartTheta", &ROOTTreeStruct.GenPartTheta, "GenPartTheta[NGenPart]/F");
  ROOTTreeDefault->Branch("GenPartPhi", &ROOTTreeStruct.GenPartPhi, "GenPartPhi[NGenPart]/F");


  ROOTTreeDefault->Branch("NRPCHit", &ROOTTreeStruct.NRPCHit, "NRPCHit/I");
  ROOTTreeDefault->Branch("RPCHitPrimaryID", &ROOTTreeStruct.RPCHitPrimaryID, "RPCHitPrimaryID[NRPCHit]/I");
  ROOTTreeDefault->Branch("RPCHitE", &ROOTTreeStruct.RPCHitE, "RPCHitE[NRPCHit]/F");
  ROOTTreeDefault->Branch("RPCHitPosX", &ROOTTreeStruct.RPCHitPosX, "RPCHitPosX[NRPCHit]/F");
  ROOTTreeDefault->Branch("RPCHitPosY", &ROOTTreeStruct.RPCHitPosY, "RPCHitPosY[NRPCHit]/F");
  ROOTTreeDefault->Branch("RPCHitPosZ", &ROOTTreeStruct.RPCHitPosZ, "RPCHitPosZ[NRPCHit]/F");
  ROOTTreeDefault->Branch("RPCHitStation", &ROOTTreeStruct.RPCHitStation, "RPCHitStation[NRPCHit]/I");
  ROOTTreeDefault->Branch("RPCHitModule", &ROOTTreeStruct.RPCHitModule, "RPCHitModule[NRPCHit]/I");
  ROOTTreeDefault->Branch("RPCHitBar", &ROOTTreeStruct.RPCHitBar, "RPCHitBar[NRPCHit]/I");

   

  
  ROOTTreeDefault->Branch("HitPDG", &ROOTTreeStruct.HitPDG, "HitPDG[NRPCHit]/I");
  ROOTTreeDefault->Branch("RPCHitPosXWo", &ROOTTreeStruct.RPCHitPosXWo, "RPCHitPosXWo[NRPCHit]/F");
  ROOTTreeDefault->Branch("RPCHitPosYWo", &ROOTTreeStruct.RPCHitPosYWo, "RPCHitPosYWo[NRPCHit]/F");  
  ROOTTreeDefault->Branch("RPCHitTime", &ROOTTreeStruct.RPCHitTime, "RPCHitTime[NRPCHit]/F");

*/
  
 // ROOTTreeDefault->Branch("RPCHitKe", &ROOTTreeStruct.RPCHitKe, "RPCHitKe[NRPCHit]/F");
 
//if you want to separet the station in the root output
// first station(1||9)
 /* ROOTTreeStn_1->Branch("NRPCHitX1", &ROOTTreeStruct.NRPCHitX1, "NRPCHitX1/I");
  ROOTTreeStn_1->Branch("RPCHitPosXStationX1", &ROOTTreeStruct.RPCHitPosXStationX1, "RPCHitPosXStationX1[NRPCHitX1]/F");
  ROOTTreeStn_1->Branch("RPCHitPosYStationX1", &ROOTTreeStruct.RPCHitPosYStationX1, "RPCHitPosYStationX1[NRPCHitX1]/F");
  ROOTTreeStn_1->Branch("RPCHitPosZStationX1", &ROOTTreeStruct.RPCHitPosZStationX1, "RPCHitPosZStationX1[NRPCHitX1]/F"); 
  ROOTTreeStn_1->Branch("RPCHitBarStationX1", &ROOTTreeStruct.RPCHitBarStationX1, "RPCHitBarStationX1[NRPCHitX1]/I");
  ROOTTreeStn_1->Branch("RPCHitModuleStationX1", &ROOTTreeStruct.RPCHitModuleStationX1, "RPCHitModuleStationX1[NRPCHitX1]/I");
  ROOTTreeStn_1->Branch("RPCHitEStationX1", &ROOTTreeStruct.RPCHitEStationX1, "RPCHitEStationX1[NRPCHitX1]/F");
  ROOTTreeStn_1->Branch("RPCHitKeStationX1", &ROOTTreeStruct.RPCHitKeStationX1, "RPCHitKeStationX1[NRPCHitX1]/F");

  ROOTTreeStn_1->Branch("NRPCHitY1", &ROOTTreeStruct.NRPCHitY1, "NRPCHitY1/I");
  ROOTTreeStn_1->Branch("RPCHitPosXStationY1", &ROOTTreeStruct.RPCHitPosXStationY1, "RPCHitPosXStationY1[NRPCHitY1]/F");
  ROOTTreeStn_1->Branch("RPCHitPosYStationY1", &ROOTTreeStruct.RPCHitPosYStationY1, "RPCHitPosYStationY1[NRPCHitY1]/F");
  ROOTTreeStn_1->Branch("RPCHitPosZStationY1", &ROOTTreeStruct.RPCHitPosZStationY1, "RPCHitPosZStationY1[NRPCHitY1]/F"); 
  ROOTTreeStn_1->Branch("RPCHitBarStationY1", &ROOTTreeStruct.RPCHitBarStationY1, "RPCHitBarStationY1[NRPCHitY1]/I");
  ROOTTreeStn_1->Branch("RPCHitModuleStationY1", &ROOTTreeStruct.RPCHitModuleStationY1, "RPCHitModuleStationY1[NRPCHitY1]/I");
  ROOTTreeStn_1->Branch("RPCHitEStationY1", &ROOTTreeStruct.RPCHitEStationY1, "RPCHitEStationY1[NRPCHitY1]/F");
  ROOTTreeStn_1->Branch("RPCHitKeStationY1", &ROOTTreeStruct.RPCHitKeStationY1, "RPCHitKeStationY1[NRPCHitY1]/F");

// station3
  ROOTTreeStn_3->Branch("NRPCHitX3", &ROOTTreeStruct.NRPCHitX3, "NRPCHitX3/I");
  ROOTTreeStn_3->Branch("RPCHitPosXStationX3", &ROOTTreeStruct.RPCHitPosXStationX3, "RPCHitPosXStationX3[NRPCHitX3]/F");
  ROOTTreeStn_3->Branch("RPCHitPosYStationX3", &ROOTTreeStruct.RPCHitPosYStationX3, "RPCHitPosYStationX3[NRPCHitX3]/F");
  ROOTTreeStn_3->Branch("RPCHitPosZStationX3", &ROOTTreeStruct.RPCHitPosZStationX3, "RPCHitPosZStationX3[NRPCHitX3]/F"); 
  ROOTTreeStn_3->Branch("RPCHitBarStationX3", &ROOTTreeStruct.RPCHitBarStationX3, "RPCHitBarStationX3[NRPCHitX3]/I");
  ROOTTreeStn_3->Branch("RPCHitModuleStationX3", &ROOTTreeStruct.RPCHitModuleStationX3, "RPCHitModuleStationX3[NRPCHitX3]/I");
  ROOTTreeStn_3->Branch("RPCHitEStationX3", &ROOTTreeStruct.RPCHitEStationX3, "RPCHitEStationX3[NRPCHitX3]/F");
  ROOTTreeStn_3->Branch("RPCHitKeStationX3", &ROOTTreeStruct.RPCHitKeStationX3, "RPCHitKeStationX3[NRPCHitX3]/F");

  ROOTTreeStn_3->Branch("NRPCHitY3", &ROOTTreeStruct.NRPCHitY3, "NRPCHitY3/I");
  ROOTTreeStn_3->Branch("RPCHitPosXStationY3", &ROOTTreeStruct.RPCHitPosXStationY3, "RPCHitPosXStationY3[NRPCHitY3]/F");
  ROOTTreeStn_3->Branch("RPCHitPosYStationY3", &ROOTTreeStruct.RPCHitPosYStationY3, "RPCHitPosYStationY3[NRPCHitY3]/F");
  ROOTTreeStn_3->Branch("RPCHitPosZStationY3", &ROOTTreeStruct.RPCHitPosZStationY3, "RPCHitPosZStationY3[NRPCHitY3]/F"); 
  ROOTTreeStn_3->Branch("RPCHitBarStationY3", &ROOTTreeStruct.RPCHitBarStationY3, "RPCHitBarStationY3[NRPCHitY3]/I");
  ROOTTreeStn_3->Branch("RPCHitModuleStationY3", &ROOTTreeStruct.RPCHitModuleStationY3, "RPCHitModuleStationY3[NRPCHitY3]/I");
  ROOTTreeStn_3->Branch("RPCHitEStationY3", &ROOTTreeStruct.RPCHitEStationY3, "RPCHitEStationY3[NRPCHitY3]/F");
  ROOTTreeStn_3->Branch("RPCHitKeStationY3", &ROOTTreeStruct.RPCHitKeStationY3, "RPCHitKeStationY3[NRPCHitY3]/F");

//station4
  ROOTTreeStn_4->Branch("NRPCHitX4", &ROOTTreeStruct.NRPCHitX4, "NRPCHitX4/I");
  ROOTTreeStn_4->Branch("RPCHitPosXStationX4", &ROOTTreeStruct.RPCHitPosXStationX4, "RPCHitPosXStationX4[NRPCHitX4]/F");
  ROOTTreeStn_4->Branch("RPCHitPosYStationX4", &ROOTTreeStruct.RPCHitPosYStationX4, "RPCHitPosYStationX4[NRPCHitX4]/F");
  ROOTTreeStn_4->Branch("RPCHitPosZStationX4", &ROOTTreeStruct.RPCHitPosZStationX4, "RPCHitPosZStationX4[NRPCHitX4]/F"); 
  ROOTTreeStn_4->Branch("RPCHitBarStationX4", &ROOTTreeStruct.RPCHitBarStationX4, "RPCHitBarStationX4[NRPCHitX4]/I");
  ROOTTreeStn_4->Branch("RPCHitModuleStationX4", &ROOTTreeStruct.RPCHitModuleStationX4, "RPCHitModuleStationX4[NRPCHitX4]/I");
  ROOTTreeStn_4->Branch("RPCHitEStationX4", &ROOTTreeStruct.RPCHitEStationX4, "RPCHitEStationX4[NRPCHitX4]/F");
  ROOTTreeStn_4->Branch("RPCHitKeStationX4", &ROOTTreeStruct.RPCHitKeStationX4, "RPCHitKeStationX4[NRPCHitX4]/F");

  ROOTTreeStn_4->Branch("NRPCHitY4", &ROOTTreeStruct.NRPCHitY4, "NRPCHitY4/I");
  ROOTTreeStn_4->Branch("RPCHitPosXStationY4", &ROOTTreeStruct.RPCHitPosXStationY4, "RPCHitPosXStationY4[NRPCHitY4]/F");
  ROOTTreeStn_4->Branch("RPCHitPosYStationY4", &ROOTTreeStruct.RPCHitPosYStationY4, "RPCHitPosYStationY4[NRPCHitY4]/F");
  ROOTTreeStn_4->Branch("RPCHitPosZStationY4", &ROOTTreeStruct.RPCHitPosZStationY4, "RPCHitPosZStationY4[NRPCHitY4]/F"); 
  ROOTTreeStn_4->Branch("RPCHitBarStationY4", &ROOTTreeStruct.RPCHitBarStationY4, "RPCHitBarStationY4[NRPCHitY4]/I");
  ROOTTreeStn_4->Branch("RPCHitModuleStationY4", &ROOTTreeStruct.RPCHitModuleStationY4, "RPCHitModuleStationY4[NRPCHitY4]/I");
  ROOTTreeStn_4->Branch("RPCHitEStationY4", &ROOTTreeStruct.RPCHitEStationY4, "RPCHitEStationY4[NRPCHitY4]/F");
  ROOTTreeStn_4->Branch("RPCHitKeStationY4", &ROOTTreeStruct.RPCHitKeStationY4, "RPCHitKeStationY4[NRPCHitY4]/F");

//station2(0||8)
  ROOTTreeStn_2->Branch("NRPCHitX2", &ROOTTreeStruct.NRPCHitX2, "NRPCHitX2/I");
  ROOTTreeStn_2->Branch("RPCHitPosXStationX2", &ROOTTreeStruct.RPCHitPosXStationX2, "RPCHitPosXStationX2[NRPCHitX2]/F");
  ROOTTreeStn_2->Branch("RPCHitPosYStationX2", &ROOTTreeStruct.RPCHitPosYStationX2, "RPCHitPosYStationX2[NRPCHitX2]/F");
  ROOTTreeStn_2->Branch("RPCHitPosZStationX2", &ROOTTreeStruct.RPCHitPosZStationX2, "RPCHitPosZStationX2[NRPCHitX2]/F");
  ROOTTreeStn_2->Branch("RPCHitBarStationX2", &ROOTTreeStruct.RPCHitBarStationX2, "RPCHitBarStationX2[NRPCHitX2]/I");
  ROOTTreeStn_2->Branch("RPCHitModuleStationX2", &ROOTTreeStruct.RPCHitModuleStationX2, "RPCHitModuleStationX2[NRPCHitX2]/I");
  ROOTTreeStn_2->Branch("RPCHitEStationX2", &ROOTTreeStruct.RPCHitEStationX2, "RPCHitEStationX2[NRPCHitX2]/F");
  ROOTTreeStn_2->Branch("RPCHitKeStationX2", &ROOTTreeStruct.RPCHitKeStationX2, "RPCHitKeStationX2[NRPCHitX2]/F");

  ROOTTreeStn_2->Branch("NRPCHitY2", &ROOTTreeStruct.NRPCHitY2, "NRPCHitY2/I");
  ROOTTreeStn_2->Branch("RPCHitPosXStationY2", &ROOTTreeStruct.RPCHitPosXStationY2, "RPCHitPosXStationY2[NRPCHitY2]/F");
  ROOTTreeStn_2->Branch("RPCHitPosYStationY2", &ROOTTreeStruct.RPCHitPosYStationY2, "RPCHitPosYStationY2[NRPCHitY2]/F");
  ROOTTreeStn_2->Branch("RPCHitPosZStationY2", &ROOTTreeStruct.RPCHitPosZStationY2, "RPCHitPosZStationY2[NRPCHitY2]/F");
  ROOTTreeStn_2->Branch("RPCHitBarStationY2", &ROOTTreeStruct.RPCHitBarStationY2, "RPCHitBarStationY2[NRPCHitY2]/I");
  ROOTTreeStn_2->Branch("RPCHitModuleStationY2", &ROOTTreeStruct.RPCHitModuleStationY2, "RPCHitModuleStationY2[NRPCHitY2]/I");
  ROOTTreeStn_2->Branch("RPCHitEStationY2", &ROOTTreeStruct.RPCHitEStationY2, "RPCHitEStationY2[NRPCHitY2]/F");
  ROOTTreeStn_2->Branch("RPCHitKeStationY2", &ROOTTreeStruct.RPCHitKeStationY2, "RPCHitKeStationY2[NRPCHitY2]/F");*/


}

void ROOTManager::Save()
{ 
  if (ROOTFile) {
    ROOTFile->Write();       
    ROOTFile->Close();       
    G4cout << "ROOT Tree closed" << G4endl;
  }
}

void ROOTManager::Fill()
{
  //ROOTTree->Fill();
  //ROOTTreeG4Particle->Fill();
   ROOTTreeDefault->Fill();
  //separeted station
 /* ROOTTreeStn_1->Fill();
  ROOTTreeStn_2->Fill();
  ROOTTreeStn_3->Fill();
  ROOTTreeStn_4->Fill();*/
}


}
